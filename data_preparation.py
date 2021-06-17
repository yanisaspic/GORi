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
from collections import OrderedDict
from _scripts import ontology as onto, uniprot as up, common as cmn, reactome as rc, hpo


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
gaf_symbol_id_dict = dict((x, y) for x, y in gene_symbol_id)
species_genes = set(gene_go_annotation.keys())

## GO ONTOLOGY
go_onto = obo_parser.GODag(go_obo_file)

## REACTOME ANNOTATIONS
gene_reactome_annotation = rc.load_reactome_annotation(reactome_annotation_file)

## REACTOME ONTOLOGY
reacterm = rc.load_reacterm_dict(reactome_hierarchy_file, reactome_label_file)
onto.save_as_obo(reacterm, reactome_obo_file, "ontology: reactome")
react_onto = obo_parser.GODag(reactome_obo_file)
species_genes.update(set(gene_reactome_annotation.keys()))

## HPO ANNOTATIONS
if species == "human":

    gene_hpo_annotation = hpo.load_hpo_annotation(hpo_annotation_file)
    gene_hpo_annotation_with_uniprot = {}
    gene_hpo_annotation_with_symbol = {}
    
    # Try to replace the default gene symbols by their UniProtKB IDs according to the GAF file
    for symbol in gene_hpo_annotation.keys():
        try:
            uniprot_id = gaf_symbol_id_dict[symbol]
            gene_hpo_annotation_with_uniprot[uniprot_id] = gene_hpo_annotation[symbol]
        
        # if a corresponding id is not found from the GAF file, save the symbols and their values
        except KeyError:
            gene_hpo_annotation_with_symbol[symbol] = gene_hpo_annotation[symbol]

    # query UniProtKB to find the corresponding ids of the remaining symbols and add them
    mg_symbol_id_dict = up.get_symbol_dict(list(gene_hpo_annotation_with_symbol.keys()))
    for symbol in mg_symbol_id_dict.keys():
        uniprot_id = mg_symbol_id_dict[symbol]
        gene_hpo_annotation_with_uniprot[uniprot_id] = gene_hpo_annotation_with_symbol[symbol]

    species_genes.update(set(gene_hpo_annotation.keys()))

## HPO ONTOLOGY
hpo_onto = obo_parser.GODag(hpo_obo_file)

## save all the ontologies
ontologies = {
    'GO': go_onto, 
    'R-': react_onto
    }
if species == "human":
    ontologies['HP'] = hpo_onto

# get the GO, Reactome (and HPO if species == human) leaf nodes corresponding to each gene.
rdy2use_data = {}

annotations = [gene_go_annotation, gene_reactome_annotation]
if species == "human":
    annotations.append(gene_hpo_annotation)

for gene_id in species_genes:

    # count the annotation sources for this gene:
    n_onto = 0
    gene_terms_for_mult_onto = []

    for anno_file in annotations:
        gene_terms_for_one_onto = set()

        try:
            leaves_id = list(anno_file[gene_id])
            source = onto.get_term_ontology(leaves_id[0], ontologies)
            for leafid in leaves_id:
                term = source[leafid]

                # for the Gene Ontology, only keep the aspect of interests
                if id(source) == id(ontologies['GO']):
                    if term.namespace not in aspect:
                        continue
                
                # add the term and its parents
                gene_terms_for_one_onto.add(leafid)
                gene_terms_for_one_onto.update(term.get_all_parents())
            n_onto += 1

        except KeyError:
            pass
        gene_terms_for_mult_onto.append(gene_terms_for_one_onto)

    if n_onto > 1:
        rdy2use_data[gene_id] = {'Gene Ontology': gene_terms_for_mult_onto[0], 'Reactome': gene_terms_for_mult_onto[1]}
        if species == "human":
            rdy2use_data[gene_id]['Human Phenotype Ontology'] = gene_terms_for_mult_onto[2]

end_load = tm.time()

#_________________________________________ E X P O R T

## Export the generated data as rdy2use files.
cmn.save_as_json(rdy2use_data, "%s/%s_gene_annotation.json" % (data_path, species))
with open("%s/%s_gene_symbol.csv" % (data_path, species), 'wt') as csv:
    csv.write("symbol,id\n")
    for key in gaf_symbol_id_dict.keys():
        csv.write("%s,%s\n" % (key, gaf_symbol_id_dict[key]))
    for key in mg_symbol_id_dict.keys():
        csv.write("%s,%s\n" % (key, mg_symbol_id_dict[key]))

end_exp = tm.time()

print("\ndata preparation completed in %s seconds:\n \
    .downloading: %ss\n \
    .loading: %ss\n \
    .exporting: %ss\n"  % (
        (end_exp - begin), 
        (end_dl - begin), 
        (end_load - end_dl),
        (end_exp - end_load)))