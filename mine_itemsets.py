"""
For a given set of genes, identify the frequent itemsets involving multiple ontologies terms.
"""

from goatools import obo_parser
from settings import *
from math import inf

import csv
import json
import time as tm
import pandas as pd
from _scripts import common as cmn, ontology as onto, fptree as fp

def prune_rules(rules, ontologies):
    """
    # Description
    Removes rules if a similar rule but with a deeper child term instead has been found.

    # Arguments
    ``rules`` (df): association rules with a body and a head. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters 
    as keys.

    # Usage
    >>> df = pd.DataFrame({'body': ['GO:0002755, 'GO:0002764, 'GO:0006955', 'GO:0000165'],
    'head': ['R-HSA-445095', 'R-HSA-445095', 'R-HSA-445095', 'R-HSA-445095']})
    >>> ontologies = ontologies = {
    'GO': obo_parser.GODag('go-basic.obo', optional_attrs="relationship"), 
    'R-': obo_parser.GODag('reactome.obo', optional_attrs="relationship")
    }
    >>> print(prune_rules(df, ontologies))
    ... body          head
    0  GO:0002755  R-HSA-445095
    3  GO:0000165  R-HSA-445095
    """
    rules_bodies = set(rules['body'])

    for body in rules_bodies:

        rules_with_body = rules.loc[rules['body'] == body]
        if len(rules_with_body) > 0:

            # get the ontology parent terms of the body
            source = onto.get_term_ontology(body, ontologies)
            term = source[body]
            parent_terms = term.get_all_upper()

            # get the heads associated with this body
            body_heads = list(rules_with_body['head'])

            # if a rule has a parent term as a body with a similar head, delete it
            rules = rules.drop(rules[ ( rules['body'].isin(parent_terms) ) &
                ( rules['head'].isin(body_heads) ) ].index)
            
            # for one body, if a rule has the parent term of an other one as a head, remove it
            heads_rules_to_ignore = set()

            for head in body_heads:

                # get the ontology parent terms of the head
                source = onto.get_term_ontology(head, ontologies)
                term = source[head]
                parent_terms = term.get_all_upper()

                found_parents = parent_terms.intersection(set(body_heads))
                heads_rules_to_ignore.update(found_parents)
            
            rules = rules.drop(rules[ ( rules['body'] == body ) &
                ( rules['head'].isin(heads_rules_to_ignore) ) ].index)

    return rules

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
    n_found_target_onto = 0
    ti = []
    gene_annotations = species_all_annotations[gene]
    target_annotations[gene] = gene_annotations

    # check if the annotation has terms corresponding to our target ontologies
    for source in target_onto:
        if len(gene_annotations[source]) > 0:
            ti.extend(gene_annotations[source])
            n_found_target_onto += 1

    # don't count uncomplete annotations in order to have valid frequencies during the FPTree
    if n_found_target_onto == len(target_onto):
        transactions.append(ti)
        ontology_terms_found.update(ti)
        
ontology_terms_found = list(ontology_terms_found)

end_load = tm.time()

#_________________________________________ W E I G H T S  &  F R E Q S 

#___________ Get the terms Information Content : 
# ontology feature-based (intrinsic)
intrinsic_IC = onto.get_all_intrinsic_IC(ontology_terms_found, ontologies)

end_weight = tm.time()

#_________________________________________ M I N I N G

rules_metrics = fp.mine_rules_metrics(transactions, intrinsic_IC)
print(len(rules_metrics))

# filter out the rules involving a body and a head sharing a relationship in an ontology
rules_metrics = rules_metrics[ 
    rules_metrics['body'].str[:2] != rules_metrics['head'].str[:2]
    ]
print(len(rules_metrics))

rules_metrics = prune_rules(rules_metrics, ontologies)
print(len(rules_metrics))

# plot the metrics together
bodies, heads = list(rules_metrics['body']), list(rules_metrics['head'])
rules_colors = fp.get_color_labels( [bodies, heads] )

rules_plot = fp.get_scatterplot(rules_metrics, 'conf', 'lift', 
"Confidence and lift values of the frequent patterns mined", rules_colors,
{'#ff00ff': 'GO + Reactome', '#00ffff': 'GO + HPO', '#ffff00': 'Reactome + HPO'})
rules_plot.show()
# set the threshold values
min_rule_conf = fp.get_numeric_input("minimum confidence", 0, 1)
min_rule_lift = fp.get_numeric_input("minimum lift", 1, inf)
# filter according to the thresholds
rules_metrics = fp.filter_confidence_and_lift(rules_metrics, min_rule_conf, min_rule_lift)
print(len(rules_metrics))

end_min = tm.time()

print("\nitemset mining completed in %s seconds:\n \
    .loading: %ss\n \
    .weighting: %ss\n \
    .mining: %ss\n" % (
        (end_min - begin), 
        (end_load - begin),
        (end_weight - end_load),
        (end_min - end_weight)))