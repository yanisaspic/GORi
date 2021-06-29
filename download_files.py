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

from goatools import obo_parser, associations, semantic
from Bio.UniProt.GOA import gafiterator
from settings import *
from collections import OrderedDict
from _scripts import ontology as onto, uniprot as up, common as cmn, reactome as rc, hpo

import time as tm
import pandas as pd
import math

def calculate_term_IC(terms_count, annotations_count):
    """
    # Description
    Returns a weight corresponding to the Information Content of the terms based on their frequencies.

    # Arguments
    ``terms_count`` (dict): the terms ID and their absolute frequency. \n
    ``annotations_count`` (dict): the prefix of an ontology terms as key and the number of genes
    annotated by the ontology as values.

    # Usage
    >>> terms_count = {'GO:0060974': 4, 'GO:0016477': 10, 'GO:2000566': 7, 'R-HSA-3656535': 4}
    >>> annotations_count = {'GO': 10, 'R-': 6}
    >>> print(calculate_term_IC(terms_count, annotations_count))
    ... {'GO:0060974': 1.0, 'GO:0016477': -0.0, 'GO:2000566': 0.3892595783536953, 'R-HSA-3656535': 1.0}
    """
    term_IC = {}
    max_IC = {'GO': 0, 'R-': 0, 'HP': 0}
    for term in terms_count.keys():
        IC_val = -math.log(terms_count[term] / annotations_count[term[:2]])
        term_IC[term] = IC_val
        max_IC[term[:2]] = max([term_IC[term], max_IC[term[:2]]])
    for term in term_IC.keys():
        term_IC[term] = term_IC[term] / max_IC[term[:2]]
    return term_IC

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
gene_go_annotation['n_genes_with_annotations'] = len(gene_go_annotation)

## GO ONTOLOGY
go_onto = obo_parser.GODag(go_obo_file)

## REACTOME ANNOTATIONS
gene_reactome_annotation = rc.load_reactome_annotation(reactome_annotation_file)
species_genes.update(set(gene_reactome_annotation.keys()))
gene_reactome_annotation['n_genes_with_annotations'] = len(gene_reactome_annotation)

## REACTOME ONTOLOGY
reacterm = rc.load_reacterm_dict(reactome_hierarchy_file, reactome_label_file)
onto.save_as_obo(reacterm, reactome_obo_file, "ontology: reactome")
react_onto = obo_parser.GODag(reactome_obo_file)

## save the ontologies
ontologies = {
    'GO': go_onto, 
    'R-': react_onto
    }

if species == "human":
    ## HPO ANNOTATIONS
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
    gene_hpo_annotation['n_genes_with_annotations'] = len(gene_hpo_annotation)

    ## HPO ONTOLOGY
    hpo_onto = obo_parser.GODag(hpo_obo_file)
    ontologies['HP'] = hpo_onto

# get the GO, Reactome (and HPO if species == human) leaf nodes corresponding to each gene.
rdy2use_annotations = {}

annotations = [gene_go_annotation, gene_reactome_annotation]
if 'HP' in ontologies:
    annotations.append(gene_hpo_annotation_with_uniprot)

for gene_id in species_genes:

    # count the annotation sources for this gene:
    n_onto = 0
    gene_terms_for_mult_onto = []

    for anno_dict in annotations:
        gene_terms_for_one_onto = set()

        try:
            leaves_id = list(anno_dict[gene_id])
            source = onto.get_term_ontology(leaves_id[0], ontologies)
            n_valid_leafid = 0

            for leafid in leaves_id:
                term = source[leafid]

                # for the Gene Ontology, only keep the terms from the aspects of interest
                if id(source) == id(ontologies['GO']) and term.namespace not in target_aspect:
                    continue
                n_valid_leafid += 1

                # add the term and its parents (RS excluded)
                ancestry = term.get_all_parents()
                gene_terms_for_one_onto.add(leafid)
                gene_terms_for_one_onto.update(ancestry)

            # if at least one term from the ontology is kept, update the counter of ontologies used
            if n_valid_leafid > 0:   
                n_onto += 1

        except KeyError:
            pass

        gene_terms_for_mult_onto.append(gene_terms_for_one_onto)

    if n_onto > 1:
        rdy2use_annotations[gene_id] = {'GO': gene_terms_for_mult_onto[0], 'R-': gene_terms_for_mult_onto[1]}
        if 'HP' in ontologies.keys():
            rdy2use_annotations[gene_id]['HP'] = gene_terms_for_mult_onto[2]

end_load = tm.time()

#_________________________________________ E X P O R T

## Export the generated data as rdy2use files.
cmn.save_as_json(rdy2use_annotations, "%s/%s_gene_annotation.json" % (data_path, species))
with open("%s/%s_gene_symbol.csv" % (data_path, species), 'wt') as csv:
    csv.write("symbol,id\n")
    for key in gaf_symbol_id_dict.keys():
        csv.write("%s,%s\n" % (key, gaf_symbol_id_dict[key]))
    if 'HP' in ontologies.keys():
        for key in mg_symbol_id_dict.keys():
            csv.write("%s,%s\n" % (key, mg_symbol_id_dict[key]))

end_exp = tm.time()

print("\ndata preparation completed in %s seconds:\n \
    .downloading: %ss\n \
    .loading: %ss\n \
    .exporting: %ss\n"  % (
        round(end_exp - begin, 1), 
        round(end_dl - begin, 1), 
        round(end_load - end_dl, 1),
        round(end_exp - end_load, 1)
        )
    )