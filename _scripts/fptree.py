"""
Ordered FP-Tree Algorithm implementation with weight filtering.

@ June 2021
@ Asloudj Yanis
"""

import pandas as pd
import time as tm
import matplotlib.pyplot as plt
import matplotlib.patches as mpa
from treelib import Node, Tree
from colorutils import Color
from collections import Counter, OrderedDict

class NodeData():
    """Data stored inside a node of the FP tree : node frequency and term weight."""
    def __init__(self, freq, weight):
        self.freq = freq
        self.weight = weight
    
    def __repr__(self):
        return 'freq: %s, weight: %s' % (self.freq, self.weight)

def get_item_metrics(trans_db, item_weights):
    """
    # Description
    Returns a data frame containing the frequencies and weights of every item.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its 
    items. \n
    ``item_weights`` (dict): the items names as keys and their weights as values. \n

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> item_metrics = get_item_metrics(trans, weights)
    >>> print(item_metrics)
    ...     item  freq  weight
        0    A   0.6     0.5
        1    C   0.4     0.9
        2    B   0.6     0.3
        3    D   0.8     0.7
        4    E   0.4     0.1
        5    Y   0.2     1.0
    """
    # ordered dicts are used to correctly align the frequencies and the weights values.
    ordered_item_freqs = OrderedDict()
    ordered_item_weights = OrderedDict()
    total_trans = len(trans_db)

    # get terms relative frequencies
    for t in trans_db:
        for item in t:
            ordered_item_freqs[item] = ordered_item_freqs.get(item, 0) + 1 / total_trans
            ordered_item_weights[item] = ordered_item_weights.get(item, item_weights[item])

    item_metrics = pd.DataFrame(
        {'item': list(ordered_item_freqs.keys()), 
        'freq': list(ordered_item_freqs.values()), 
        'weight': list(ordered_item_weights.values())})
    return item_metrics

def filter_frequency_and_weight(item_metrics, min_item_freq, min_item_weight):
    """
    # Description
    Filters out all the frequencies and weights below their respective threshold values and 
    orders the remaining ones descendingly (frequence) then ascendingly (alphabetically).

    # Arguments
    ``item_metrics`` (df): the items names, their weights and their frequencies. \n
    ``min_item_freq`` (float): the minimum frequency necessary to keep an item. 0 < min_item_freq <= 1. \n
    ``min_item_weight`` (float): the minimum weight necessary to keep an item. 0 <= min_item_freq <= 1. 

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> item_metrics = get_item_metrics(trans, weights)
    >>> item_metrics = filter_frequency_and_weight(item_metrics, 0.25, 0.5)
    >>> print(item_metrics)
    ...     item  freq  weight
        3    D   0.8     0.7
        0    A   0.6     0.5
        1    C   0.4     0.9
    """
    item_metrics = item_metrics[item_metrics['freq'] >= min_item_freq]
    item_metrics = item_metrics[item_metrics['weight'] >= min_item_weight]
    item_metrics = item_metrics.sort_values(by=['freq', 'item'], ascending=[False, True])
    return item_metrics

def order_transaction(ti, item_metrics):
    """
    # Description
    Orders a transaction according to the frequency values of the items.

    # Arguments
    ``ti`` (list): a transaction. Each element of the list is an item. \n
    ``item_metrics`` (df): the items names, their weights and their frequencies.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> item_metrics = get_item_metrics(trans, weights)
    >>> item_metrics = filter_frequency_and_weight(item_metrics, 0.25, 0.5)
    >>> print(order_transaction(trans[0], item_metrics))
    ... ['D', 'A', 'C']
    """
    order = list(item_metrics['item'])
    freq_ti = [x for x in ti if x in order]
    return sorted(freq_ti, key = order.index)

def construct_fptree(trans_db, item_metrics):
    """
    # Description
    Returns an FP tree with a dict of items names as keys and nodes ids (= nids) as values. 

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. 
    A sublist is a transaction with its items. \n
    ``item_metrics`` (df): the items names, their weights and their frequencies. \n

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> item_metrics = get_item_metrics(trans, weights)
    >>> item_metrics = filter_frequency_and_weight(item_metrics, 0.25, 0.5)
    >>> tree, item_nodes = construct_fptree(trans, item_metrics)
    >>> tree.show(idhidden=False)
    ... [0]
        ├── A[4]
        └── D[1]
            ├── A[2]
            │   └── C[3]
            └── C[5]
    >>> tree.show(idhidden=False, data_property="freq")
    ... [0]
        ├── 0.2[4]
        └── 0.8[1]
            ├── 0.4[2]
            │   └── 0.2[3]
            └── 0.2[5]
    >>> print(item_nodes)
    ... {'D': [1], 'A': [2, 4], 'C': [3, 5]}
    """
    # dict indicating the nodes ids corresponding to each item : one item can have multiple nids
    item_nodes = {}
    total_trans = len(trans_db)

    # init the tree
    tree = Tree()
    new_nid = 0
    tree.create_node(tag="", identifier=new_nid, data=NodeData("", ""))
    new_nid += 1  # increment the next node id value whenever a new node is created

    for ti in trans_db:

        parent_nid = 0  # root node
        ti = order_transaction(ti, item_metrics)

        for item in ti:

            # get the children nids of the current parent
            children_nids = tree.is_branch(parent_nid)

            # if there is a child nid corresponding to one of the known item's nid, 
            # navigate to the corresponding node and update its frequency.
            if item in item_nodes.keys() and any(nid in item_nodes[item] for nid in children_nids):
                parent_nid = min(set(item_nodes[item]) & (set(children_nids)))
                tree[parent_nid].data.freq += 1/total_trans
            
            # if there is none, create a new node, navigate to it and save its nid in a dict.
            else:
                item_nodes[item] = item_nodes.get(item, []) + [new_nid]
                tree.create_node(
                    tag=item, identifier=new_nid, parent=parent_nid, 
                    data = NodeData(
                        1/total_trans, 
                        item_metrics.loc[item_metrics['item'] == item, 'weight'].values[0])
                )

                # generate a new unique nid and navigate to the new node
                parent_nid = new_nid
                new_nid += 1

    return tree, item_nodes

def get_branch(fptree, nid):
    """
    # Description
    Returns a dict made of all the parent nodes of a nid. 
    The node tags act as keys and the frequencies act as values.

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
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> item_metrics = get_item_metrics(trans, weights)
    >>> item_metrics = filter_frequency_and_weight(item_metrics, 0.25, 0.5)
    >>> tree, item_nodes = construct_fptree(trans, item_metrics)
    >>> print(get_branch(tree, 3))
    ... {'A': 0.2, 'D': 0.2}
    """
    branch = {}
    leaf_node = fptree.get_node(nid)
    pattern_frequency = leaf_node.data.freq
    parent_node = fptree.get_node(leaf_node.predecessor(fptree._identifier))
    while parent_node.predecessor(fptree._identifier) != None:
        branch[parent_node.tag] = pattern_frequency
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
    >>> w_items = {'A': 0.5, 'B': 0.1, 'D': 0.8, 'C': 0.1, 'Y': 0.5, 'E': 0.23}
    >>> freq_items = get_items_frequency(trans_db = trans)
    >>> freq_items = filter_frequency(freq_items, min_item_freq = 0.25)
    >>> tree, item_nodes = construct_fptree(trans, freq_items, w_items)
    >>> print(get_condition_tree('C', tree, item_nodes))
    ... Counter({'D': 0.4, 'A': 0.2})
    """
    cond_tree = Counter()
    for nid in item_nodes[item]:
        cond_tree += Counter(get_branch(fptree, nid))
    return cond_tree

def get_rules_metrics(fptree, item_nodes, item_metrics):
    """
    # Description
    Returns the associtation rules and their metrics by reading the FPtree.

    # Arguments
    ``fptree`` (Tree Object): items nodes and their conditional frequencies stored in a treelib 
    Tree Object. \n
    ``item_nodes`` (dict): items tags and their respective nodes ids. \n
    ``item_metrics`` (df): the items names, their weights and their frequencies. \n

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> item_metrics = get_item_metrics(trans, weights)
    >>> item_metrics = filter_frequency_and_weight(item_metrics, 0.25, 0.5)
    >>> tree, item_nodes = construct_fptree(trans, item_metrics)
    >>> print(get_rules_metrics(tree, item_nodes, item_metrics))
    ... body head      conf      lift  cover
    0    D    A  0.500000  1.041667    0.8
    1    A    C  0.333333  1.388889    0.6
    2    D    C  0.500000  1.562500    0.8
    """
    rules_metrics = []

    for item in list(item_metrics['item'][1:]):
        cond_tree = get_condition_tree(item, fptree, item_nodes)
        for parent_item in cond_tree.keys():

            # get the frequencies of the node and its parents to calculate quality metrics
            parent_freq = item_metrics.loc[item_metrics['item'] == parent_item, 'freq'].values[0]
            item_freq = item_metrics.loc[item_metrics['item'] == item, 'freq'].values[0]

            # confidence(Y => X) = S(x ∩ y) / S(y)
            confidence = cond_tree[parent_item] / parent_freq
            # lift(Y => X) = conf(Y => X) / S(X) * S(Y)
            lift = confidence / (item_freq * parent_freq)

            rules_metrics.append(
                {'body': parent_item, 'head': item, 
                'conf': confidence, 'lift': lift, 'cover': parent_freq})

    rules_metrics = pd.DataFrame(rules_metrics)
    return rules_metrics

def filter_confidence_and_lift(rules_metrics, min_rule_conf, min_rule_lift):
    """
    # Description
    Filters out all the confidences and lifts below their respective threshold values and 
    orders the remaining ones descendingly (coverage) then ascendingly (alphabetically).

    # Arguments
    ``rules_metrics`` (df): association rules with a body, a head, confidence, lift and coverage values.
    ``min_rule_conf`` (float): the minimum confidence necessary to keep a rule. 0 < min_rule_conf <= 1. \n
    ``min_lift_conf`` (float): the minimum lift necessary to keep a rule. 1 <= min_rule_lift <= +Inf.

    # Usage
    >>> rules_metrics = pd.DataFrame({'body': ['D', 'D', 'A', 'C'], 'head': ['A', 'C', 'F', 'F'],
    'conf': [0.5, 0.5, 0.2, 0.4], 'lift': [1.0, 1.6, 0.8, 0.7], 'cover':[0.8, 0.8, 0.3, 0.2]})
    >>> print(filter_confidence_and_lift(rules_metrics, 0.5, 1))
    ... body head  conf  lift  cover
    0    D    A   0.5   1.0    0.8
    1    D    C   0.5   1.6    0.8
    """
    rules_metrics = rules_metrics[rules_metrics['conf'] >= min_rule_conf]
    rules_metrics = rules_metrics[rules_metrics['lift'] >= min_rule_lift]
    rules_metrics = rules_metrics.sort_values(by=['cover', 'body'], ascending=[True, True])
    return rules_metrics

def get_numeric_input(label, mini, maxi):
    """
    # Description
    Asks the user for an input value between min and max and returns basic messages if the value
    is invalid.

    # Arguments
    ``label`` (string): a label corresponding to the variable you expect.
    ``mini`` (float): the lowest threshold value.
    ``maxi`` (float): the highest threshold value.

    # Usage
    >>> score = get_numeric_input("score", 0, 20)
    ... Please set a valid score value (0 <= score <= 20): 30
    ... Wrong input: out of limits value.
    ... Please set a valid score value (0 <= score <= 20): A
    ... Wrong input: non numeric value.
    ... Please set a valid score value (0 <= score <= 20): 10
    >>> print(score)
    ... 10
    """
    var = input("Please set a valid %s value (%s <= %s <= %s): " % (label, mini, label, maxi))
    try:
        if float(var) < mini or float(var) > maxi:
            print("Wrong input: out of limits value.")
            var = get_numeric_input(label, mini, maxi)
    except ValueError:
        print("Wrong input: non numeric value.")
        var = get_numeric_input(label, mini, maxi)
    return float(var)

def get_color_labels(terms_lists, colors = {
    'R-': Color((255,0,0)), 
    'HP': Color((0,255,0)), 
    'GO': Color((0,0,255))
    }):
    """
    # Description
    Returns a list corresponding to the sum of the ontologies colors.
    e.g. for a list of 2 sublists, if the first element of the sublist is GO (blue) and HP (green),
    the first element of the resulting colors list is cyan).

    # Arguments
    ``terms_lists`` (list of sublists): a list of equally sized sublists.
    ``rules_metrics`` (df): a dataframe with rules bodies and heads. \n
    ``colors`` (dict of tuples): the color value corresponding to each prefix.

    # Usage
    >>> sub_list1 = ['GO:0001817', 'GO:0001817']
    >>> sub_list2 = ['GO:0010628', 'HP:0011024']
    >>> print(get_color_labels( [sub_list1, sub_list2] ))
    ... ['#0000ff', '#00ffff']
    """
    summed_colors = []
    len_subl = len(terms_lists[0])
    for i in range(len_subl):
        color_value = Color((0, 0, 0))
        for subl in terms_lists:
            term_prefix = subl[i][:2]
            color_value += colors[term_prefix]
        summed_colors.append(color_value.hex)
    return summed_colors

def get_scatterplot(data, xlab, ylab, title = None, colors = None, legend = None, mark = "."):
    """
    # Description
    Returns a scatterplot object using 2 columns of a dataframe as x and y axes.

    # Arguments
    ``data`` (df): a pandas dataframe with at least 2 columns of values. \n
    ``xlab`` (string): the name of the column you want to use as x values. \n
    ``ylab`` (string): the name of the column you want to use as y values. \n
    ``title`` (string): the title of your scatterplot. \n
    ``colors`` (list of strings): the ordered color labels corresponding to your points. \n
    ``legend`` (dict): colors as keys and their significance as values. \n
    ``marker`` (string): the appearance of your scatterplot points.

    # Usage
    >>> df = pd.DataFrame({'x': [0, 1, 2], 'y': [3, 4, 5]})
    >>> scatter = get_scatterplot(
    df, 'x', 'y', colors=['red', 'blue', 'green'], 
    legend={'red':'A', 'blue':'B', 'green':'C'}, mark='v')
    >>> scatter.show()
    """
    plt.scatter(data[xlab], data[ylab], c=colors, marker=mark)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    patches = []
    for k in legend.keys():
        patches.append(mpa.Patch(color=k, label=legend[k]))
    plt.legend(handles=patches, loc='lower left')
    return plt

def mine_rules_metrics(trans_db, item_weights):
    """
    # Description
    Returns the frequent patterns found in a transaction database and metrics using an ordered FP tree.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with 
    its items. \n
    ``item_weights`` (dict): the items names as keys and their weights as values.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> print(mines_rules_metrics(trans, weights))
    ... Please set a valid minimum frequency value (0 <= minimum frequency <= 1): 0.25
    ... Please set a valid minimum weight value (0 <= minimum weight <= 1): 0.5
    ... {
        'D => A': {'conf': 0.5, 'lift': 1.0416666666666665, 'cover': 0.8}, 
        'D => C': {'conf': 0.5, 'lift': 1.5624999999999998, 'cover': 0.8}
        }
    """
    # get a frequency, weights dataframe for items
    item_metrics = get_item_metrics(trans_db, item_weights)

    # plot the metrics together
    items = list(item_metrics['item'])
    items_colors = get_color_labels( [items] )

    items_plot = get_scatterplot(item_metrics, 'weight', 'freq', 
    "Relative weights and frequencies of the ontologies terms", items_colors, 
    {'#0000ff': 'Gene Ontology', '#ff0000': 'Reactome', '#00ff00': 'Human Phenotype Ontology'})
    items_plot.show()

    # set the threshold values
    min_item_weight = get_numeric_input("minimum weight", 0, 1)
    min_item_freq = get_numeric_input("minimum frequency", 0, 1)
    item_metrics = filter_frequency_and_weight(item_metrics, min_item_freq, min_item_weight)

    # construct the tree
    tree, item_nodes = construct_fptree(trans_db, item_metrics)

    # get a confidence, lift dataframe for rules with ascending coverage values
    rules_metrics = get_rules_metrics(tree, item_nodes, item_metrics)

    return rules_metrics