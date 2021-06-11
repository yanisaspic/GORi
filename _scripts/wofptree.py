"""
Weighted Ordered FP-Tree Algorithm implementation.

WOFP:
Yuanyuan Li, Shaohong Yin (2020). 
Mining algorithm for weighted FP-tree frequent item sets based on two-dimensional table.
Journal of Physics: Conference Series.
DOI:10.1088/1742-6596/1453/1/012002

@ June 2021
@ Asloudj Yanis
"""

import pandas as pd
import time as tm
from treelib import Node, Tree
from collections import Counter, OrderedDict

def weigh_transaction(trans, weight_items):
    """
    # Description
    Calculates the average weight of a transaction according to the individual weights of the items.

    # Arguments
    ``trans`` (list): a transaction. Each element is an item. \n
    ``weight_items`` (dict): the items names as keys and their weights as values.

    # Usage
    >>> t = ["A", "B", "D", "C"]
    >>> weight = {'A': 2, 'B': 4, 'D': 8, 'C': 1}
    >>> print(weigh_transaction(t, weight))
    ... 3.75
    """
    total = 0
    for item in trans:
        total += weight_items[item]
    return total / len(trans)

def get_frequency_and_weight(trans_db, weight_items, min_item_freq = 0.1):
    """
    # Description
    Returns a data frame containing the frequencies of the frequent items only, in a descending order.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its items. \n
    ``min_item_freq`` (float): the minimum frequency necessary to keep an item. 0 < min_item_freq <= 1. 

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 2, 'B': 4, 'D': 8, 'C': 1, 'Y': 0.5, 'E': 23}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> print(freq_items)
    ... item  freq
    3    D   0.8
    0    A   0.6
    1    B   0.6
    2    C   0.4
    4    E   0.4
    >>> print(trans_weights)
    ... [3.75, 11.666666666666666, 9.666666666666666, 3.1666666666666665, 5.0]
    """
    item_freq = {}
    total_trans = len(trans_db)
    weight_trans = []
    for t in trans_db:
        weight_trans.append(weigh_transaction(t, weight_items))
        for item in t:
            try:
                item_freq[item] += 1/total_trans
            except KeyError:
                item_freq[item] = 1/total_trans
    freq_df = pd.DataFrame({'item': list(item_freq.keys()), 'freq': list(item_freq.values())})
    freq_df = freq_df[freq_df['freq'] > min_item_freq]
    freq_df = freq_df.sort_values(by=['freq', 'item'], ascending=[False, True])
    return freq_df, weight_trans

def order_transaction(trans, freq_items):
    """
    # Description
    Orders a transaction according to the frequency values of the items.

    # Arguments
    ``trans`` (list): a transaction. Each element is an item. \n
    ``freq_items`` (df): the items names and their frequencies, ordered by descending frequency and ascending name.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 2, 'B': 4, 'D': 8, 'C': 1, 'Y': 0.5, 'E': 23}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> print(order_transaction(trans[3], freq_items))
    ... ['D', 'C']
    """
    order = list(freq_items['item'])
    frequent_trans = [x for x in trans if x in order]
    return sorted(frequent_trans, key = order.index)

def construct_fptree(trans_db, freq_items, weight_trans):
    """
    # Description
    Returns an FP tree and a dict of items and their nids from the transactions stored in a database.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its items. \n
    ``freq_items`` (df): the items names and their frequencies, ordered by descending frequency and ascending name. \n
    ``weight_items`` (list): a list with the respective weights of each transactions.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 2, 'B': 4, 'D': 8, 'C': 1, 'Y': 0.5, 'E': 23}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, trans_weights)
    >>> tree.show(idhidden=False)
    ...
    [0]
    ├── A[7]
    │   └── B[8]
    │       └── E[9]
    └── D[1]
        ├── A[2]
        │   └── B[3]
        │       └── C[4]
        ├── B[5]
        │   └── E[6]
        └── C[10]
    >>> print(item_nodes)
    ... OrderedDict([('D', [1]), ('A', [2, 7]), ('B', [3, 5, 8]), ('C', [4, 10]), ('E', [6, 9])])
    """
    # ordered dict indicating the nodes ids corresponding to each item (one item can have multiple nids)
    item_nodes = OrderedDict()

    # init the tree
    tree = Tree()
    tree.create_node(tag="", identifier=0)
    nid = 1 # node id
    total_weight = sum(weight_trans)

    # transactions weights counter
    k = 0

    for t in trans_db:
        parent_nid = 0
        t = order_transaction(t, freq_items)

        for item in t:

            # get the children nids of the current parent
            children_nids = tree.is_branch(parent_nid)

            # if there is no node child corresponding to the item, add it on the tree and save the new nid in the dict
            if item not in item_nodes.keys() or not any(nid in item_nodes[item] for nid in children_nids):
                try:
                    item_nodes[item].append(nid)
                except KeyError:
                    item_nodes[item] = [nid]
                tree.create_node(tag=item, identifier=nid, parent=parent_nid, data = weight_trans[k] / total_weight)

                # generate a new unique nid and navigate to the new node
                parent_nid = nid
                nid += 1

            # if there is, update the weight of the node and navigate to it
            else:
                parent_nid = min(set(item_nodes[item]).intersection(set(children_nids)))
                tree[parent_nid].data += weight_trans[k] / total_weight
        k += 1

    return tree, item_nodes

def get_branch(fptree, nid):
    """
    # Description
    Returns a dict corresponding to the parents of a node and the data value of the child node

    # Arguments
    ``fptree`` (Tree Object): items nodes and their weights stored in a treelib Tree Object. \n
    ``nid`` (int): unique id corresponding to an existing node of the Tree Object.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 2, 'B': 4, 'D': 8, 'C': 1, 'Y': 0.5, 'E': 23}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, trans_weights)
    >>> print(get_branch(tree, 4))
    ... {'B': 0.11278195488721804, 'A': 0.11278195488721804, 'D': 0.11278195488721804}
    """
    branch = {}
    leaf_node = fptree.get_node(nid)
    value = leaf_node.data
    parent_node = fptree.get_node(leaf_node.predecessor(fptree._identifier))
    while parent_node.predecessor(fptree._identifier) != None:
        branch[parent_node.tag] = value
        parent_node = fptree.get_node(parent_node.predecessor(fptree._identifier))
    return branch

def get_condition_tree(item, fptree, item_nodes):
    """
    # Description
    Returns the conditional FP tree of a given item according to its nodes values on a tree.

    # Arguments
    ``item`` (string): the item you want the conditional FPtree of. \n
    ``fptree`` (Tree Object): items nodes and their weights stored in a treelib Tree Object. \n
    ``item_nodes`` (dict): items tags and their respective nodes ids.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 2, 'B': 4, 'D': 8, 'C': 1, 'Y': 0.5, 'E': 23}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, trans_weights)
    >>> print(get_condition_tree('E', tree, item_nodes))
    ... Counter({'B': 0.6416040100250626, 'D': 0.3508771929824561, 'A': 0.2907268170426065})
    """
    cond_tree = Counter()
    for nid in item_nodes[item]:
        cond_tree += Counter(get_branch(fptree, nid))
    return cond_tree

def get_association_rules(fptree, item_nodes, min_trans_weight = 0.20):
    """
    # Description
    Returns the items association rules by reading the FPtree.

    # Arguments
    ``fptree`` (Tree Object): items nodes and their weights stored in a treelib Tree Object. \n
    ``item_nodes`` (dict): items tags and their respective nodes ids. \n
    ``min_trans_weight`` (float): itemsets with a weight equal or higher to this value (0<v<1) are considered patterns.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 2, 'B': 4, 'D': 8, 'C': 1, 'Y': 0.5, 'E': 23}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, trans_weights)
    >>> print(get_association_rules(tree, item_nodes, min_trans_weight = 0.50))
    ... {'E;B': 0.64}
    """
    item_nodes.popitem(last=False)
    frequent_pattern = {}
    for item in item_nodes.keys():
        cond_tree = get_condition_tree(item, fptree, item_nodes)
        for parent_item in cond_tree.keys():
            if cond_tree[parent_item] > min_trans_weight:
                pattern = "%s;%s" % (item, parent_item)
                frequent_pattern[pattern] = round(cond_tree[parent_item], 2)
    return frequent_pattern

def filter_patterns(freq_pat, fix_length = 2):
    """
    # Description
    Returns a dict of frequent patterns after filtering according to items prefix or suffix.

    # Arguments
    ``freq_pat`` (dict): frequent patterns and their associated weight value. \n
    ``fix_length`` (int): the length of the suffix or prefix.

    # Usage
    >>> trans = [
            ["GO1", "GO3", "GO2", "R-1"],
            ["GO2", "HP1", "R-1"],
            ["HP1", "GO2", "GO1"],
            ["R-2", "GO3", "R-1"],
            ["R-1", "GO1"]
        ]
    >>> w_items = {'GO1': 1, 'GO2': 1, 'R-1': 1, 'GO3': 1, 'R-2': 1, 'HP1': 1}
    >>> freq_items, trans_weights = get_frequency_and_weight(trans_db = trans, weight_items = w_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, trans_weights)
    >>> asso_rules = get_association_rules(tree, item_nodes)
    >>> print(asso_rules)
    ... {'GO1;R-1': 0.4, 'GO2;GO1': 0.4, 'GO2;R-1': 0.4, 'GO3;R-1': 0.4, 'HP1;GO2': 0.4}
    >>> print(filter_patterns(asso_rules))
    ... {'GO1;R-1': 0.4, 'GO2;R-1': 0.4, 'GO3;R-1': 0.4, 'HP1;GO2': 0.4}
    """
    patterns_to_del = []
    for pat in freq_pat.keys():
        terms = pat.split(";")
        if terms[0][:fix_length] == terms[1][:fix_length]:
            patterns_to_del.append(pat)
    for k in patterns_to_del:
        freq_pat.pop(k)
    return freq_pat

def mine_frequent_itemsets(trans_db, weight_items, min_item_freq = 0.1, min_trans_weight = 0.2):
    """
    # Description
    Returns the frequent patterns in a trans database using a weighted ordered FP tree.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its items. \n
    ``weight_items`` (dict): the items names as keys and their weights as values. \n
    ``min_item_freq`` (float): the minimum frequency necessary to keep an item. 0 < min_item_freq <= 1. \n
    ``min_trans_weight`` (float): itemsets with a weight equal or higher to this value (0<v<1) are considered patterns.

    # Usage
    >>> trans = [
            ["GO1", "GO3", "GO2", "R-1"],
            ["GO2", "HP1", "R-1"],
            ["HP1", "GO2", "GO1"],
            ["R-2", "GO3", "R-1"],
            ["R-1", "GO1"]
        ]
    >>> w_items = {'GO1': 1, 'GO2': 1, 'R-1': 1, 'GO3': 1, 'R-2': 1, 'HP1': 1}
    >>> frequent_patterns = mine_frequent_itemsets(trans, w_items, 0.1, 0.2)
    >>> print(frequent_patterns)
    ... {'GO1;R-1': 0.4, 'GO2;R-1': 0.4, 'GO3;R-1': 0.4, 'HP1;GO2': 0.4}
    """
    freq_items, trans_weights = get_frequency_and_weight(trans_db, weight_items, min_item_freq)
    tree, item_nodes = construct_fptree(trans_db, freq_items, trans_weights)
    frequent_patterns = get_association_rules(tree, item_nodes, min_trans_weight)
    return filter_patterns(frequent_patterns)