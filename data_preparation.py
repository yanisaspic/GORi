"""
Retrieves the Gene Ontology and Reactome data.
The GO OBO file is filtered according to the terms specific to a species found using the GAF file.
A Reactome OBO file is created using the hierarchy and the label Reactome files.
The Reactome OBO file is filtered according to the terms specific to a species found using the leaves Reactome file.

The files used are indicated on settings.py.

Gene Association File, GAF format (up to version 2.2):
http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/

OBO format (up to version 1.2):
https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_2.html

@ May 2021
@ Asloudj Yanis
"""

from goatools import obo_parser, associations
from Bio.UniProt.GOA import gafiterator
from settings import *
from scripts import ontology as onto
from scripts import common as cmn
from scripts import reactome as rc
from scripts import hpo


import time as tm
import pandas as pd


# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -

#_________________________________________ P R E P A R A T I O N __________________________________________

#_________________________________________ D O W N L O A D I N G

begin = tm.time()

associations.dnld_annotation(gaf_file)
cmn.replace_first_line(gaf_file, "!gaf-version: 2.1\n")
cmn.download_url(go_obo_url, go_obo_file)
cmn.download_url(reactome_hierarchy_url, reactome_hierarchy_file)
cmn.download_url(reactome_label_url, reactome_label_file)
cmn.download_url(reactome_annotation_url, reactome_annotation_file)

if species == "human":
    cmn.download_replace_first(hpo_annotation_url, hpo_annotation_file, "")
    cmn.download_url(hpo_obo_url, hpo_obo_file)

end_dl = tm.time()

#_________________________________________ L O A D I N G

## GO ANNOTATIONS
gene_go_annotation = {}
gene_symbol_id = set()
with open(gaf_file, 'rt') as gaf:
    for anno in gafiterator(gaf):
        gene_symbol_id.add((anno['DB_Object_Symbol'], anno['DB_Object_ID']))
        try:
            gene_go_annotation[anno['DB_Object_ID']].add(anno['GO_ID'])
        except KeyError:
            gene_go_annotation[anno['DB_Object_ID']] = set([anno['GO_ID']])
gene_symbol_id_dict = dict((x, y) for x, y in gene_symbol_id)

## GO ONTOLOGY
go_onto = obo_parser.GODag(go_obo_file, optional_attrs = "relationship")

## REACTOME ANNOTATIONS
gene_reactome_annotation = rc.load_reactome_annotation(reactome_annotation_file)

## REACTOME ONTOLOGY
reacterm = rc.load_reacterm_dict(reactome_hierarchy_file, reactome_label_file)
onto.save_as_obo(reacterm, reactome_obo_file, "ontology: reactome")
react_onto = obo_parser.GODag(reactome_obo_file)
species_genes = set(gene_go_annotation.keys())
species_genes.update(set(gene_reactome_annotation.keys()))

if species == "human":    
    ## HPO ANNOTATIONS
    gene_hpo_annotation = hpo.load_hpo_annotation(hpo_annotation_file)
    gene_uniprotkb_hpo_annotation = []
    # replace gene symbols by their UniProtKB IDs
    for symbol in gene_hpo_annotation.keys():
        try:
            uniprotkb_id = gene_symbol_id_dict[symbol]
            gene_uniprotkb_hpo_annotation.append((uniprotkb_id, gene_hpo_annotation[symbol]))
        except KeyError:
            pass
    gene_hpo_annotation = dict((x, y) for x, y in gene_uniprotkb_hpo_annotation)
    species_genes.update(set(gene_hpo_annotation.keys()))

## HPO ONTOLOGY
hpo_onto = obo_parser.GODag(hpo_obo_file)

# get the GO, Reactome (and if species == human, HPO) leaf nodes corresponding to each gene.
rdy2use_data = {}
for gene_id in species_genes:

    # count the annotation sources for this gene:
    n_onto = 0

    try:
        gene_leaf_goid = gene_go_annotation[gene_id]
        gene_aspect_leaf_goid = set()
        gene_aspect_goid = set()
        for goid in gene_leaf_goid:
            if go_onto[goid].namespace == aspect:
                gene_aspect_leaf_goid.add(goid)
                # gene_aspect_goid.update(go_onto[goid].get_all_upper())
        ancestry_goid = onto.get_ancestry_id(gene_id, gene_go_annotation, go_onto)
        ## If the entire tree is added to the JSON file
        # gene_goid = list(gene_aspect_goid)
        ## If only the leaves are added to the JSON file
        gene_goid = list(gene_aspect_leaf_goid)
        n_onto += 1

    except KeyError:
        gene_goid = set()

    try:
        # ancestry_reactid = onto.get_ancestry_id(gene_id, gene_reactome_annotation, react_onto)
        # gene_reactid = list(ancestry_reactid)
        gene_reactid = list(gene_reactome_annotation[gene_id])
        n_onto += 1
    except KeyError:
        gene_reactid = set()

    if species == "human":
        try:
            # ancestry_hpoid = onto.get_ancestry_id(gene_id, gene_hpo_annotation, hpo_onto)
            # gene_hpoid = list(ancestry_hpoid)
            gene_hpoid = list(gene_hpo_annotation[gene_id])
            n_onto += 1
        except KeyError:
            gene_hpoid = set()
    
    # if the gene can help find links between multiple ontologies, save the data
    if n_onto > 1:
        rdy2use_data[gene_id] = {'Gene Ontology': gene_goid, 'Reactome': gene_reactid}
        if species == "human":
            rdy2use_data[gene_id]['Human Phenotype Ontology'] = gene_hpoid

end_load = tm.time()

#_________________________________________ E X P O R T

## Export the generated data as rdy2use files.
cmn.save_as_json(rdy2use_data, "%s/%s_gene_annotation.json" % (data_path, species))
with open("%s/%s_gene_symbol.csv" % (data_path, species), 'wt') as csv:
    csv.write("symbol,id\n")
    for key in gene_symbol_id_dict.keys():
        csv.write("%s,%s\n" % (key, gene_symbol_id_dict[key]))

end_exp = tm.time()

print("\ndata preparation completed in %s seconds:\n \
    .downloading: %ss\n \
    .loading: %ss\n \
    .exporting: %ss\n"  % (
        (end_exp - begin), 
        (end_dl - begin), 
        (end_load - end_dl),
        (end_exp - end_load)))