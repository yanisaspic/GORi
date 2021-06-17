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

def filter_association_rules_with_onto(asso_rules, ontologies, sep = " => "):
    """
    # Description
    Filters out the association rules based on two related ontology terms.

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
    >>> asso_rules = filter_association_rules_with_onto(asso_rules, ontologies)
    >>> print(asso_rules)
    ... {'GO:0002757 => GO:0050865': {'confidence': 0.7, 'lift': 1.22, 'coverage': 0.58}, 
        'GO:0002757 => GO:0010646': {'confidence': 0.57, 'lift': 1.0, 'coverage': 0.58}}
    """
    filtered_asso_rules = {}
    for rul in asso_rules.keys():
        body, head = rul.split(sep)
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
    reader = csv.reader(open("%s/%s_gene_symbol.csv" % (data_path, species), 'rt'))
    gene_symbol_id_dict = {}
    for row in reader:
        k, v = row
        gene_symbol_id_dict[k] = v
    genes = cmn.get_values_from_keys(genes, gene_symbol_id_dict)

# load ontologies
ontologies = {
    'GO': obo_parser.GODag("%s/go-basic.obo" % data_path, optional_attrs="relationship"), 
    'R-': obo_parser.GODag("%s/reactome.obo" % data_path, optional_attrs="relationship")
    }
if species == "human":
    ontologies['HP'] = obo_parser.GODag(hpo_obo_file, optional_attrs="relationship")

## extract the annotations corresponding to the selected genes:
with open("%s/%s_gene_annotation.json" % (data_path, species), 'rt') as anno:
    annotation = json.load(anno)

## extract the annotations corresponding to the selected ontologies:
transactions = []
known_genes = set(genes).intersection(set(annotation.keys()))
used_terms = set()
for gene in known_genes:
    ti = []
    gene_anno = annotation[gene]
    for source in target_ontology:
        ti.extend(gene_anno[source])
        used_terms.update(gene_anno[source])
    transactions.append(ti)

## calculate the items weights
item_weights = {}
for term in used_terms:
    item_weights[term] = onto.term_weight(term, ontologies)

end_load = tm.time()

#_________________________________________ M I N I N G

##########################################
min_item_freq = 0.25                     #
min_item_weight = 0.50                   #
min_pattern_conf = 0.75                  #
##########################################

asso_rules = wofp.mine_association_rules(transactions, item_weights, min_item_freq, min_item_weight, min_pattern_conf)
asso_rules = filter_association_rules_with_onto(asso_rules, ontologies)
asso_rules = rename_association_rules_with_onto(asso_rules, ontologies)
print(asso_rules)

end_min = tm.time()

print("\nitemset mining completed in %s seconds:\n \
    .loading: %ss\n \
    .mining: %ss\n" % (
        (end_min - begin), 
        (end_load - begin),
        (end_min - end_load)))