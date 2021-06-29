"""
For a given set of genes, identify the frequent itemsets involving multiple ontologies terms.
"""

from _scripts import common as cmn, ontology as og, fptree as fp, analysis as ana
from goatools import obo_parser
from settings import *
from math import inf
from collections import OrderedDict

import csv
import json
import time as tm
import pandas as pd

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

#_________________________________________ I T E M I Z I N G

#___________ Get the terms Information Content : 
# ontology feature-based (intrinsic)
k = 0
intrinsic_IC, labels = og.get_all_intrinsic_IC_and_labels(ontology_terms_found, ontologies, index = 'Zhou', k)

# get a frequency, weights dataframe for items
items_metrics = fp.get_items_metrics(transactions, intrinsic_IC, labels)

# plot the metrics together
frequencies = OrderedDict()
weights = OrderedDict()
for index, row in items_metrics.iterrows():
    item_id = row['item'][:2]
    weights[item_id] = weights.get(item_id, []) + [ row['weight'] ]
    frequencies[item_id] = frequencies.get(item_id, []) + [ row['freq'] ]

hist_weights = list(weights.values())
hist_frequencies = list(frequencies.values())

#######################################################################################################
##### P L O T   T H E   W E I G H T S   I N   H I S T O G R A M S
#####   - incorpore ta fonction d'histogrammes ici ^ (l96)
#####   - utilise hist_weights et hist_frequenceis (2 listes de 3 sous-listes)
#####   - fait varier k dans le get_weight_zhou (de 0 à 1) [l80]
#####   - détermine si une valeur de k permet d'avoir des distributions similaires entre les ontologies
#######################################################################################################

scatter_items = fp.get_scatterplot(
    weights, frequencies, axes=['weight', 'freq'],
    title='Metrics of the ontology terms')
scatter_items.show()

# set the threshold values and filter
min_item_weight = fp.get_numeric_input("minimum weight", 0, 1)
min_item_freq = fp.get_numeric_input("minimum frequency", 0, 1)
items_metrics = fp.filter_frequency_and_weight(items_metrics, min_item_freq, min_item_weight)

end_item = tm.time()

#_________________________________________ M I N I N G

# construct the tree and extract rules metrics from it
tree, item_nodes = fp.construct_fptree(transactions, items_metrics)
rules_metrics = fp.get_rules_metrics(tree, item_nodes, items_metrics)

end_mine = tm.time()

#_________________________________________ P R U N I N G

# filter out the rules involving a body and a head sharing a relationship in an ontology
rules_metrics = rules_metrics[ 
    rules_metrics['body'].str[:2] != rules_metrics['head'].str[:2]
    ]
rules_metrics = og.prune_rules(rules_metrics, ontologies)

# plot the metrics together
confidences = OrderedDict()
lifts = OrderedDict()
coverages = OrderedDict()
for index, row in rules_metrics.iterrows():
    rule_elements = sorted([ row['body'][:2], row['head'][:2] ])
    rule_id = " & ".join(rule_elements)
    confidences[rule_id] = confidences.get(rule_id, []) + [ row['conf'] ]
    lifts[rule_id] = lifts.get(rule_id, []) + [ row['lift'] ]

scatter_rules = fp.get_scatterplot(
    confidences, lifts, coverages, axes=['conf', 'lift'], title="Metrics of the association rules")
scatter_rules.show()

# set the threshold values and filter
min_rule_conf = fp.get_numeric_input("minimum confidence", 0, 1)
min_rule_lift = fp.get_numeric_input("minimum lift", 0, inf)
rules_metrics = fp.filter_confidence_and_lift(rules_metrics, min_rule_conf, min_rule_lift)

end_prune = tm.time()

#_________________________________________ E X P O R T

# export the items and rules metrics used to draw a network

items_metrics = items_metrics.drop(columns=['weight'])
items_metrics.to_csv(items_file, index=False, header=True)
rules_metrics = rules_metrics.drop(columns=['cover'])
rules_metrics.to_csv(rules_file, index=False, header=True)

end_exp = tm.time()

print("\nitemset mining completed in %s seconds:\n \
    .loading: %ss\n \
    .itemizing: %ss\n \
    .mining: %ss\n \
    .pruning: %ss\n \
    .export: %ss\n" % (
        round(end_prune - begin, 1), 
        round(end_load - begin, 1),
        round(end_item - end_load, 1),
        round(end_mine - end_item, 1),
        round(end_prune - end_mine, 1),
        round(end_exp - end_prune, 1)
        )
    )
