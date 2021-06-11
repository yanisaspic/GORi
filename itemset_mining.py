"""
For a given set of genes, identify the frequent itemsets involving multiple ontologies terms.
"""

from goatools import obo_parser
from settings import *
from mlxtend.preprocessing.transactionencoder import TransactionEncoder

import csv
import json
import time as tm
import pandas as pd
from scripts import common as cmn
from scripts import ontology as onto
from scripts import wofptree as wofp

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
    'R-': obo_parser.GODag("%s/reactome.obo" % data_path)
    }
if species == "human":
    ontologies['HP'] = obo_parser.GODag(hpo_obo_file)

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
    for source in target_onto:
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
min_item_freq = 0.1                      #
min_trans_weight = 0.1                   #
##########################################

frequent_patterns = wofp.mine_frequent_itemsets(transactions, item_weights, min_item_freq, min_trans_weight)
print(frequent_patterns)

end_min = tm.time()

print("\nitemset mining completed in %s seconds:\n \
    .loading: %ss\n \
    .mining: %ss\n" % (
        (end_min - begin), 
        (end_load - begin),
        (end_min - end_load)))