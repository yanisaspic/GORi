"""
Ordered FP-Tree Algorithm implementation with weight filtering.

@ June 2021
@ Asloudj Yanis
"""

import pandas as pd
import time as tm
from treelib import Node, Tree
from collections import Counter, OrderedDict
from _scripts import ontology

class NodeData():
    """Data stored inside a node of the FP tree : node frequency and term weight."""
    def __init__(self, freq, weight):
        self.freq = freq
        self.weight = weight
    
    def __repr__(self):
        return 'freq: %s, weight: %s' % (self.freq, self.weight)

def get_frequency(trans_db):
    """
    # Description
    Returns a data frame containing the frequencies of every item.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its items. \n

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> freq_items = get_frequency(trans_db = trans)
    >>> print(freq_items)
    ...   item  freq
        0    A   0.6
        1    C   0.4
        2    B   0.6
        3    D   0.8
        4    E   0.4
        5    Y   0.2
    """
    item_freq = {}
    total_trans = len(trans_db)
    for t in trans_db:
        for item in t:
            try:
                item_freq[item] += 1/total_trans
            except KeyError:
                item_freq[item] = 1/total_trans
    freq_df = pd.DataFrame({'item': list(item_freq.keys()), 'freq': list(item_freq.values())})
    return freq_df

def filter_frequency(freq_items, min_item_freq):
    """
    # Description
    Filters out all the frequency below a threshold and orders the remaining ones descendingly.

    # Arguments
    ``freq_items`` (df): the items names and their frequencies. \n
    ``min_item_freq`` (float): the minimum frequency necessary to keep an item. 0 < min_item_freq <= 1. 

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> freq_items = get_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> print(freq_items)
    ...   item  freq
        3    D   0.8
        0    A   0.6
        2    B   0.6
        1    C   0.4
        4    E   0.4
    """
    freq_items = freq_items[freq_items['freq'] > min_item_freq]
    freq_items = freq_items.sort_values(by=['freq', 'item'], ascending=[False, True])
    return freq_items

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
    >>> freq_items = get_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> print(order_transaction(trans[3], freq_items))
    ... ['D', 'C']
    """
    order = list(freq_items['item'])
    frequent_trans = [x for x in trans if x in order]
    return sorted(frequent_trans, key = order.index)

def construct_fptree(trans_db, freq_items, weight_items):
    """
    # Description
    Returns an FP tree and a dict of items and their nids from the transactions stored in a database.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its items. \n
    ``freq_items`` (df): the items names and their frequencies, ordered by descending frequency and ascending name. \n
    ``weight_items`` (dict): the items names as keys and their weights as values.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 0.5, 'B': 0.1, 'D': 0.8, 'C': 0.1, 'Y': 0.5, 'E': 0.23}
    >>> freq_items = get_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, w_items)
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
    total_trans = len(trans_db)

    # init the tree
    tree = Tree()
    tree.create_node(tag="", identifier=0)
    nid = 1 # node id

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
                tree.create_node(
                    tag=item, identifier=nid, parent=parent_nid, 
                    data = NodeData(1/total_trans, weight_items[item])
                    )

                # generate a new unique nid and navigate to the new node
                parent_nid = nid
                nid += 1

            # if there is, update the frequency of the node and navigate to it
            else:
                parent_nid = min(set(item_nodes[item]).intersection(set(children_nids)))
                tree[parent_nid].data.freq += 1/total_trans
        k += 1

    return tree, item_nodes

def get_branch(fptree, nid, min_item_weight):
    """
    # Description
    Returns a dict of parent nodes with a weight equal or higher to the threshold. \n
    The tag acts as key and the frequency acts as value.

    # Arguments
    ``fptree`` (Tree Object): items nodes and their weights stored in a treelib Tree Object. \n
    ``nid`` (int): unique id corresponding to an existing node of the Tree Object. \n
    ``min_item_weight`` (float): items with a weight equal or higher to this value (0<v<1) are considered informative.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 0.5, 'B': 0.1, 'D': 0.8, 'C': 0.1, 'Y': 0.5, 'E': 0.23}
    >>> freq_items = get_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, w_items)
    >>> print(get_branch(tree, 4, 0.25))
    ... {'A': 0.2, 'D': 0.2}
    """
    branch = {}
    leaf_node = fptree.get_node(nid)
    pattern_frequency = leaf_node.data.freq
    parent_node = fptree.get_node(leaf_node.predecessor(fptree._identifier))
    while parent_node.predecessor(fptree._identifier) != None:
        if parent_node.data.weight >= min_item_weight:
            branch[parent_node.tag] = pattern_frequency
        parent_node = fptree.get_node(parent_node.predecessor(fptree._identifier))
    return branch

def get_condition_tree(item, fptree, item_nodes, min_item_weight):
    """
    # Description
    Returns the conditional FP tree of a given item according to its nodes values on a tree.

    # Arguments
    ``item`` (string): the item you want the conditional FPtree of. \n
    ``fptree`` (Tree Object): items nodes and their weights stored in a treelib Tree Object. \n
    ``item_nodes`` (dict): items tags and their respective nodes ids. \n
    ``min_item_weight`` (float): items with a weight equal or higher to this value (0<v<1) are considered informative.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 0.5, 'B': 0.1, 'D': 0.8, 'C': 0.1, 'Y': 0.5, 'E': 0.23}
    >>> freq_items = get_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, w_items)
    >>> print(get_condition_tree('C', tree, item_nodes, min_item_weight = 0.25))
    ... Counter({'D': 0.4, 'A': 0.2})
    """
    cond_tree = Counter()
    for nid in item_nodes[item]:
        cond_tree += Counter(get_branch(fptree, nid, min_item_weight))
    return cond_tree

def get_association_rules(fptree, item_nodes, freq_items, min_item_weight, min_pattern_conf):
    """
    # Description
    Returns the associtation rules by reading the FPtree.

    # Arguments
    ``fptree`` (Tree Object): items nodes and their weights stored in a treelib Tree Object. \n
    ``item_nodes`` (dict): items tags and their respective nodes ids. \n
    ``freq_items`` (df): the items names and their frequencies, ordered by descending frequency and ascending name. \n
    ``min_item_weight`` (float): items with a weight equal or higher to this value (0<v<1) are considered informative. \n
    ``min_pattern_conf`` (float): the confidence value used as a threshold to define an association rule.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 0.5, 'B': 0.1, 'D': 0.8, 'C': 0.1, 'Y': 0.5, 'E': 0.23}
    >>> freq_items = get_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, w_items)
    >>> print(get_association_rules(tree, item_nodes, freq_items, min_item_weight = 0.25, min_pattern_conf = 0.5))
    ... {'D => A': 0.5, 'A => B': 0.67, 'D => B': 0.5, 'D => C': 0.5}
    """
    item_nodes.popitem(last=False)
    
    frequent_pattern = {}
    for item in item_nodes.keys():
        cond_tree = get_condition_tree(item, fptree, item_nodes, min_item_weight)

        for parent_item in cond_tree.keys():

            # get the frequencies of the node and its parents
            parent_row = freq_items.loc[freq_items['item'] == parent_item]
            parent_freq = parent_row['freq'].iloc[0]
            item_row = freq_items.loc[freq_items['item'] == item]
            item_freq = item_row['freq'].iloc[0]

            # confidence(Y => X) = S(x ∩ y) / S(y)
            confidence = cond_tree[parent_item] / parent_freq
            if confidence >= min_pattern_conf:

                # lift > 1: rule body and rule head appear more often together than expected
                lift = confidence / item_freq
                if lift > 1:
                    # add the pattern and its metrics
                    pattern = "%s => %s" % (parent_item, item)
                    frequent_pattern[pattern] = {
                        'confidence': round(confidence, 2),
                        'lift': round(lift, 2),
                        'coverage': round(parent_freq, 2)}
    return frequent_pattern

def mine_association_rules(trans_db, weight_items, min_item_freq = 0.1, min_item_weight = 0.25, min_pattern_conf = 0.5):
    """
    # Description
    Returns the frequent patterns in a trans database using a weighted ordered FP tree.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its items. \n
    ``weight_items`` (dict): the items names as keys and their weights as values. \n
    ``min_item_freq`` (float): the minimum frequency necessary to keep an item. 0 < min_item_freq <= 1. \n
    ``min_item_weight`` (float): itemsets with a weight equal or higher to this value (0<v<1) are considered patterns.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> w_items = {'A': 0.5, 'B': 0.1, 'D': 0.8, 'C': 0.1, 'Y': 0.5, 'E': 0.23}
    >>> print(mine_association_rules(trans, w_items))
    ... {'D => A': 0.5, 'A => B': 0.67, 'D => B': 0.5, 'D => C': 0.5}
    """
    freq_items = get_frequency(trans_db)
    freq_items = filter_frequency(freq_items, min_item_freq)
    tree, item_nodes = construct_fptree(trans_db, freq_items, weight_items)
    return get_association_rules(tree, item_nodes, freq_items, min_item_weight, min_pattern_conf)

