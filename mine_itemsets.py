"""
For a given set of genes, identify the frequent itemsets involving multiple ontologies terms.
"""

from goatools import obo_parser
from settings import *

import csv
import json
import time as tm
import pandas as pd
from _scripts import common as cmn, ontology as onto, wofptree as wofp

def get_terms_frequencies(annotations):
    """
    # Description
    Returns a dict of terms as keys and their Information Content (IC) as values. The IC is calculated
    according to the frequencies of each term in the annotations.

    # Arguments
    ``annotations`` (dict of dict): a dict of genes as keys and their terms as values. The terms must 
    be grouped according to their source ontology.

    # Usage
    >>> annotations = {
        'GeneA': {
            'GO': ["GO:0009987", "GO:0008150"], 
            'R-': ["R-HSA-73857", "R-HSA-74160"],
            'HP': set()
            },
        'GeneB': {
            'GO': ["GO:0008150"], 
            'R-': ["R-HSA-74160"],
            'HP': ["HP:0000001"]
            },
        'GeneC': {
            'GO': ["GO:0009987", "GO:0008150", "GO:0030029"], 
            'R-': ["R-HSA-73857", "R-HSA-74160", "R-HSA-1643685"],
            'HP': ["HP:0000118", "HP:0000001"]
            }
        }
    >>> print(get_terms_frequencies(annotations))
    ... {'GO:0009987': 0.67, 'GO:0008150': 1.0, 'R-HSA-73857': 0.67, 'R-HSA-74160': 1.0, 
    'HP:0000001': 1.0, 'GO:0030029': 0.33, 'R-HSA-1643685': 0.33, 'HP:0000118': 0.5}
    """
    ## Reads the annotations to get the absolute frequencies of the terms and the number
    ## genes annotated by each ontology
    annotated_genes_count = {}
    terms_count = {}
    for gene in annotations.keys():
        one_gene_annotations_from_all_onto = annotations[gene]
        for ontology in one_gene_annotations_from_all_onto.keys():
            one_gene_annotations_from_one_onto = one_gene_annotations_from_all_onto[ontology]
            if len(one_gene_annotations_from_one_onto) > 0:
                annotated_genes_count[ontology] = annotated_genes_count.get(ontology, 0) + 1
                for term in one_gene_annotations_from_one_onto:
                    terms_count[term] = terms_count.get(term, 0) + 1

    # get the relative frequencies
    for term in terms_count.keys():
        term_prefix = term[:2]
        term_freq = round(terms_count[term] / annotated_genes_count[term_prefix], 2)
        terms_count[term] = term_freq

    return terms_count

def filter_association_rules_with_onto(asso_rules, ontologies, sep = " => ", mult_onto = True):
    """
    # Description
    Filters out the association rules based on two related ontology terms.

    # Arguments
    ``asso_rules`` (dict): separated terms as keys and dict of metrics (confidence, lift, coverage) 
    as values. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 
    characters as keys. \n
    ``sep`` (string): the string separating two terms in a key. \n
    ``mult_onto`` (boolean): should the associations between two terms of the same 
    ontology be filtered out as well ?

    # Usage
    >>> asso_rules = {
        'GO:0002429 => GO:0002757': {'confidence': 1.0, 'lift': 1.73, 'coverage': 0.58}, 
        'GO:0002757 => GO:0050865': {'confidence': 0.7, 'lift': 1.22, 'coverage': 0.58},
        'GO:0002757 => GO:0010646': {'confidence': 0.57, 'lift': 1.0, 'coverage': 0.58}}
    >>> ontologies = {
        'GO': obo_parser.GODag("go-basic.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> asso_rules = filter_association_rules_with_onto(asso_rules, ontologies)
    >>> print(asso_rules)
    ... {'GO:0002757 => GO:0050865': {'confidence': 0.7, 'lift': 1.22, 'coverage': 0.58}, 
        'GO:0002757 => GO:0010646': {'confidence': 0.57, 'lift': 1.0, 'coverage': 0.58}}
    """
    filtered_asso_rules = {}
    for rul in asso_rules.keys():
        body, head = rul.split(sep)
        if mult_onto and body[:2] == head[:2]:
            continue
        source = onto.get_term_ontology(body, ontologies)
        body = source[body]
        if head in body.get_all_upper():
            continue
        filtered_asso_rules[rul] = asso_rules[rul]
    return filtered_asso_rules

def rename_association_rules_with_onto(asso_rules, ontologies, sep = " => "):
    """
    # Description
    Renames the association rules with human readable labels.

    # Arguments
    ``asso_rules`` (dict): separated terms as keys and dict of metrics (confidence, lift, coverage) as values. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters as keys. \n
    ``sep`` (string): the string separating two terms in a key.

    # Usage
    >>> asso_rules = {
        'GO:0002429 => GO:0002757': {'confidence': 1.0, 'lift': 1.73, 'coverage': 0.58}, 
        'GO:0002757 => GO:0050865': {'confidence': 0.7, 'lift': 1.22, 'coverage': 0.58},
        'GO:0002757 => GO:0010646': {'confidence': 0.57, 'lift': 1.0, 'coverage': 0.58}}
    >>> ontologies = {
        'GO': obo_parser.GODag("go-basic.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> asso_rules = rename_association_rules_with_onto(asso_rules, ontologies)
    >>> print(asso_rules)
    ... {
        'immune response-activating cell surface receptor signaling pathway => immune response-activating signal transduction': {'confidence': 1.0, 'lift': 1.73, 'coverage': 0.58}, 
        'immune response-activating signal transduction => regulation of cell activation': {'confidence': 0.7, 'lift': 1.22, 'coverage': 0.58}, 
        'immune response-activating signal transduction => regulation of cell communication': {'confidence': 0.57, 'lift': 1.0, 'coverage': 0.58}
        }
    """
    renamed_asso_rules = {}
    for rul in asso_rules.keys():
        body, head = rul.split(sep)
        body_source = onto.get_term_ontology(body, ontologies)
        body_name = body_source[body].name
        head_source = onto.get_term_ontology(head, ontologies)
        head_name = head_source[head].name
        new_key = "%s (%s)%s%s (%s)" % (body_name, body, sep, head_name, head)
        renamed_asso_rules[new_key] = asso_rules[rul]
    return renamed_asso_rules

# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -

#_________________________________________ I T E M S E T S __________________________________________

#_________________________________________ L O A D I N G

begin = tm.time()

# Convert input gene ids into symbols
if symbol:
    gene_symbol_to_id_dict = {}
    reader = csv.reader(open("%s/%s_gene_symbol.csv" % (data_path, species), 'rt'))
    for row in reader:
        symbol, uniprot = row
        gene_symbol_to_id_dict[symbol] = uniprot
    target_genes = cmn.get_values_from_keys(keys=target_genes, dictio=gene_symbol_to_id_dict)

# load ontologies
ontologies = {
    'GO': obo_parser.GODag(go_obo_file, optional_attrs="relationship"), 
    'R-': obo_parser.GODag(reactome_obo_file, optional_attrs="relationship")
    }
if species == "human":
    ontologies['HP'] = obo_parser.GODag(hpo_obo_file, optional_attrs="relationship")

## extract the annotations corresponding to the selected genes:
with open("%s/%s_gene_annotation.json" % (data_path, species), 'rt') as anno:
    species_all_annotations = json.load(anno)

## extract the annotations of the target genes and convert them to transactions 
## regardless of the ontologies sources
transactions = []
known_target_genes = set(target_genes).intersection(set(species_all_annotations.keys()))

target_annotations = {}         # used to calculate the frequencies
ontology_terms_found = set()    # used to calculate the weights

for gene in known_target_genes:
    ti = []
    gene_annotations = species_all_annotations[gene]
    target_annotations[gene] = gene_annotations
    for source in ontologies.keys():
        ti.extend(gene_annotations[source])
        ontology_terms_found.update(gene_annotations[source])
    transactions.append(ti)
ontology_terms_found = list(ontology_terms_found)

end_load = tm.time()

#_________________________________________ W E I G H T S  &  F R E Q S 

#___________ Get the terms Information Content : 
# ontology feature-based (intrinsic)
intrinsic_IC = onto.get_all_intrinsic_IC(ontology_terms_found, ontologies)

#___________ Get the terms frequencies :
terms_freq = get_terms_frequencies(target_annotations)

end_weight = tm.time()

#_________________________________________ M I N I N G

##########################################
min_item_freq = 0.20                     #
min_item_weight = 0.95                   #
min_pattern_conf = 0.33                  #
##########################################

asso_rules = wofp.mine_association_rules(transactions, intrinsic_IC)
asso_rules = filter_association_rules_with_onto(asso_rules, ontologies)
asso_rules = rename_association_rules_with_onto(asso_rules, ontologies)
print(asso_rules)

end_min = tm.time()

print("\nitemset mining completed in %s seconds:\n \
    .loading: %ss\n \
    .weighting: %ss\n \
    .mining: %ss\n" % (
        (end_min - begin), 
        (end_load - begin),
        (end_weight - end_load),
        (end_min - end_weight)))