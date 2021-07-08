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
from collections import Counter, OrderedDict
from mpl_toolkits.mplot3d import Axes3D

class NodeData():
    """Data stored inside a node of the FP tree : node frequency and term weight."""
    def __init__(self, freq, weight):
        self.freq = freq
        self.weight = weight
    
    def __repr__(self):
        return 'freq: %s, weight: %s' % (self.freq, self.weight)

def get_items_metrics(trans_db, item_weights, item_labels):
    """
    # Description
    Returns a data frame containing the frequencies and weights of every item.

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. A sublist is a transaction with its 
    items. \n
    ``item_weights`` (dict): the items names as keys and their weights as values. \n
    ``item_labels`` (dict): the items names as keys and their labels as values.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> labels = {'A': 'Alpha', 'B': 'Beta', 'C': 'Gamma', 'D': 'Delta', 'E': 'Episode', 'Y': 'Yahtzee'}
    >>> items_metrics = get_items_metrics(trans, weights, labels)
    >>> print(items_metrics)
    ... item  freq  weight    label
    0    A   0.6     0.5    Alpha
    1    C   0.4     0.9    Gamma
    2    B   0.6     0.3     Beta
    3    D   0.8     0.7    Delta
    4    E   0.4     0.1  Episode
    5    Y   0.2     1.0  Yahtzee
    """
    # ordered dicts are used to correctly align the frequencies and the weights values.
    ordered_item_freqs = OrderedDict()
    ordered_item_weights = OrderedDict()
    ordered_item_labels = OrderedDict()
    total_trans = len(trans_db)
    test = []

    # get terms relative frequencies
    for t in trans_db:
        for item in t:
            ordered_item_freqs[item] = ordered_item_freqs.get(item, 0) + 1 / total_trans
            ordered_item_weights[item] = ordered_item_weights.get(item, item_weights[item])
            ordered_item_labels[item] = ordered_item_labels.get(item, item_labels[item])

    items_metrics = pd.DataFrame(
        {'item': list(ordered_item_freqs.keys()), 
        'freq': list(ordered_item_freqs.values()), 
        'weight': list(ordered_item_weights.values()),
        'label': list(ordered_item_labels.values())})
    return items_metrics

def filter_frequency_and_weight(items_metrics, min_item_freq, min_item_weight):
    """
    # Description
    Filters out all the frequencies and weights below their respective threshold values and 
    orders the remaining ones descendingly (frequence) then ascendingly (alphabetically).

    # Arguments
    ``items_metrics`` (df): the items names, their weights and their frequencies. \n
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
    >>> labels = {'A': 'Alpha', 'B': 'Beta', 'C': 'Gamma', 'D': 'Delta', 'E': 'Episode', 'Y': 'Yahtzee'}
    >>> items_metrics = get_items_metrics(trans, weights, labels)
    >>> items_metrics = filter_frequency_and_weight(items_metrics, 0.25, 0.5)
    >>> print(items_metrics)
    ...     item  freq  weight
        3    D   0.8     0.7
        0    A   0.6     0.5
        1    C   0.4     0.9
    """
    if min_item_freq < 0:
        items_metrics = items_metrics[items_metrics['freq'] <= -min_item_freq]
    else:
        items_metrics = items_metrics[items_metrics['freq'] >= min_item_freq]
    if min_item_weight < 0:
        items_metrics = items_metrics[items_metrics['weight'] <= -min_item_weight]
    else:
        items_metrics = items_metrics[items_metrics['weight'] >= min_item_weight]
    items_metrics = items_metrics.sort_values(by=['freq', 'item'], ascending=[False, True])
    return items_metrics

def order_transaction(ti, items_metrics):
    """
    # Description
    Orders a transaction according to the frequency values of the items.

    # Arguments
    ``ti`` (list): a transaction. Each element of the list is an item. \n
    ``items_metrics`` (df): the items names, their weights and their frequencies.

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> items_metrics = get_items_metrics(trans, weights)
    >>> items_metrics = filter_frequency_and_weight(items_metrics, 0.25, 0.5)
    >>> print(order_transaction(trans[0], items_metrics))
    ... ['D', 'A', 'C']
    """
    order = list(items_metrics['item'])
    freq_ti = [x for x in ti if x in order]
    return sorted(freq_ti, key = order.index)

def construct_fptree(trans_db, items_metrics):
    """
    # Description
    Returns an FP tree with a dict of items names as keys and nodes ids (= nids) as values. 

    # Arguments
    ``trans_db`` (list of sublists): your transaction database. 
    A sublist is a transaction with its items. \n
    ``items_metrics`` (df): the items names, their weights and their frequencies. \n

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> items_metrics = get_items_metrics(trans, weights)
    >>> items_metrics = filter_frequency_and_weight(items_metrics, 0.25, 0.5)
    >>> tree, item_nodes = construct_fptree(trans, items_metrics)
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
        ti = order_transaction(ti, items_metrics)

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
                        items_metrics.loc[items_metrics['item'] == item, 'weight'].values[0])
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
    >>> items_metrics = get_items_metrics(trans, weights)
    >>> items_metrics = filter_frequency_and_weight(items_metrics, 0.25, 0.5)
    >>> tree, item_nodes = construct_fptree(trans, items_metrics)
    >>> print(get_branch(tree, 3))
    ... {'A': 0.2, 'D': 0.2}
    """
    branch = {}
    leaf_node = fptree.get_node(nid)
    pattern_frequency = leaf_node.data.freq
    parent_node = fptree.get_node(leaf_node.predecessor(fptree._identifier))
    while parent_node.predecessor(fptree._identifier) != None:
        if leaf_node.tag[:2] != parent_node.tag[:2]:
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

def get_rule(pattern_freq, body_row, head_row):
    """
    # Description
    Returns a rule as a dictionary with a body, a head, a confidence, a lift and a coverage value.

    # Arguments
    ``pattern_freq`` (float): the frequency of a 2-items frequent pattern. \n
    ``body_row`` (df): the frequency dataframe row of the item considered as body. \n
    ``head_row`` (df): the frequency dataframe row of the item considered as head.

    # Usage
    >>> items_freq = pd.DataFrame({
        'item': ['A', 'B', 'C', 'D'], 'freq': [0.2, 0.4, 0.5, 0.6]
    })
    >>> B_row, C_row = items_freq[items_freq['item'] == 'B'], items_freq[items_freq['item'] == 'C']
    >>> freq_BC_pattern = 0.3
    >>> print( get_rule(freq_BC_pattern, B_row, C_row) )
    ... {'body': 'B', 'head': 'C', 'conf': 0.7499999999999999, 'lift': 3.749999999999999, 'cover': 0.4}
    >>> print( get_rule(freq_BC_pattern, C_row, B_row) )
    ... {'body': 'C', 'head': 'B', 'conf': 0.6, 'lift': 2.9999999999999996, 'cover': 0.5}
    """
    body_id, body_freq = body_row['item'].values[0], body_row['freq'].values[0]
    head_id, head_freq = head_row['item'].values[0], head_row['freq'].values[0]

    # confidence(Y => X) = S(X ∩ Y) / S(Y)
    confidence = pattern_freq / body_freq
    # lift(Y => X) = conf(Y => X) / S(X)
    lift = confidence / head_freq

    rule = {
        'body': body_id, 'head': head_id, 
        'conf': confidence, 'lift': lift, 'cover': body_freq}
    return rule

def get_rules_metrics(fptree, item_nodes, items_metrics):
    """
    # Description
    Returns the associtation rules and their metrics by reading the FPtree.

    # Arguments
    ``fptree`` (Tree Object): items nodes and their conditional frequencies stored in a treelib 
    Tree Object. \n
    ``item_nodes`` (dict): items tags and their respective nodes ids. \n
    ``items_metrics`` (df): the items names, their weights and their frequencies. \n

    # Usage
    >>> trans = [
            ["A", "C", "B", "D"],
            ["B", "E", "D"],
            ["E", "B", "A"],
            ["Y", "C", "D"],
            ["D", "A"]
        ]
    >>> weights = {'A': 0.5, 'B': 0.3, 'C':0.9, 'D': 0.7, 'E': 0.1, 'Y': 1}
    >>> labels = {'A': 'Alpha', 'B': 'Beta', 'C': 'Gamma', 'D': 'Delta', 'E': 'Episode', 'Y': 'Yahtzee'}
    >>> items_metrics = get_items_metrics(trans, weights, labels)
    >>> items_metrics = filter_frequency_and_weight(items_metrics, 0.25, 0.5)
    >>> tree, item_nodes = construct_fptree(trans, items_metrics)
    >>> print(get_rules_metrics(tree, item_nodes, items_metrics))
    ... body head      conf      lift  cover
    0    D    A  0.500000  1.041667    0.8
    1    A    D  0.666667  1.388889    0.6
    2    A    C  0.333333  1.388889    0.6
    3    C    A  0.500000  2.083333    0.4
    4    D    C  0.500000  1.562500    0.8
    5    C    D  1.000000  3.125000    0.4
    """
    rules_metrics = []

    for item in list(items_metrics['item'][1:]):

        cond_tree = get_condition_tree(item, fptree, item_nodes)

        for parent_item in cond_tree.keys():

            # get the frequencies of the node and its parents to calculate quality metrics
            parent_row = items_metrics[items_metrics['item'] == parent_item]
            item_row = items_metrics[items_metrics['item'] == item]
            parent_U_item_freq = cond_tree[parent_item]

            rules_metrics.extend([
                get_rule(parent_U_item_freq, parent_row, item_row),
                get_rule(parent_U_item_freq, item_row, parent_row)
            ])

    rules_metrics = pd.DataFrame(rules_metrics)
    return rules_metrics

def filter_confidence_and_lift(rules_metrics, min_rule_conf, min_rule_lift):
    """
    # Description
    Filters out all the rules with confidences and lifts below their respective threshold values and 
    orders the remaining ones descendingly (coverage) then ascendingly (alphabetically).

    # Arguments
    ``rules_metrics`` (df): association rules with a body, a head, confidence, lift and coverage values.
    ``min_rule_conf`` (float): the minimum confidence necessary to keep a rule. 0 < min_rule_conf <= 1. \n
    ``min_rule_lift`` (float): the minimum lift necessary to keep a rule. 1 <= min_rule_lift <= +Inf. \n

    # Usage
    >>> rules_metrics = pd.DataFrame({'body': ['D', 'D', 'A', 'C'], 'head': ['A', 'C', 'F', 'F'],
    'conf': [0.5, 0.5, 0.2, 0.4], 'lift': [1.0, 1.6, 0.8, 0.7], 'cover':[0.8, 0.8, 0.3, 0.2]})
    >>> print(filter_confidence_and_lift(rules_metrics, 0.5, 1))
    ... body head  conf  lift  cover
    0    D    A   0.5   1.0    0.8
    1    D    C   0.5   1.6    0.8
    """
    if min_rule_conf < 0:
        rules_metrics = rules_metrics[rules_metrics['conf'] <= -min_rule_conf]
    else:
        rules_metrics = rules_metrics[rules_metrics['conf'] >= min_rule_conf]
    if min_rule_lift < 0:
        rules_metrics = rules_metrics[rules_metrics['lift'] <= -min_rule_lift]
    else:
        rules_metrics = rules_metrics[rules_metrics['lift'] >= min_rule_lift]
    rules_metrics = rules_metrics.sort_values(by=['cover', 'body'], ascending=[True, True])
    return rules_metrics

def get_numeric_input(label, mini, maxi):
    """
    # Description
    Asks the user for an input value between min and max and returns basic messages if the value
    is invalid.

    # Arguments
    ``label`` (string): a label corresponding to the variable you expect. \n
    ``mini`` (float): the lowest threshold value. \n
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
    print("\n## You can add a minus (e.g. -0.1) to input a maximal %s value instead (%s <= 0.1)."
        % (label, label) )
    var = input("Please set a valid minimal %s value (%s <= %s <= %s): " % (
        label, mini, label, maxi))
    negative = False
    try:
        var = float(var)
        if var < 0:
            negative = True
            var = -var
        if var < mini or var > maxi:
            print("Wrong input: out of limits value.")
            var = get_numeric_input(label, mini, maxi)
    except ValueError:
        print("Wrong input: non numeric value.")
        var = get_numeric_input(label, mini, maxi)
    if negative:
        var = -var
    return var

def get_scatterplot(
    x, y, z = None, title = None, axes = ['x', 'y', 'z'], 
    colors = 'Set1', mark = ["*", ".", "p"]):
    """
    # Description
    Returns a scatterplot object using 2 columns of a darame as x and y axes.

    # Arguments
    ``x`` (dict): category identifier as key and list of numeric as values. \n
    ``y`` (dict): category identifier as key (same as x) and list of numeric as values. \n
    ``z`` (dict): category identifier as key (same as x and y) and list of numeric as values (3D). \n
    ``title`` (string): the title of your scatterplot. \n
    ``axes`` (list of strings): ordered names of the axes of the graph (x>y>z) \n
    ``colors`` (string): the name of the colormap used to represent each category. \n
    ``mark`` (list of strings): the markers used to represent each category.

    # Usage
    >>> x = {'A': [0, 1, 2], 'B': [5, 6, 7]}
    >>> y = {'A': [0.2, 0.3, 0.1], 'B': [0.9, 1.2, 0.1]}
    >>> z = {'A': [10, 11, 9], 'B': [3, 4, 8]}
    >>> scatter = get_scatterplot(x, y, title='Demo 2D scatterplot')
    >>> scatter.show()
    >>> scatter = get_scatterplot(x, y, z, title='Demo 3D scatterplot')
    >>> scatter.show()
    """
    categories = list(x.keys())
    plt.title(title)

    if z is None:
        layers = []
        labels = []

        for i in range(len(x)):
            cat = categories[i]
            labels.append(cat)
            lay = plt.scatter(x[cat], y[cat], cmap=colors, marker=mark[i], alpha=0.5)
            layers.append(lay)

        plt.xlabel(axes[0])
        plt.ylabel(axes[1])
        plt.legend(layers, labels, loc='lower left', fancybox=True, shadow=True)

    else:
        ax = plt.subplot(111, projection='3d')
        for i in range(len(x)):
            cat = categories[i]
            ax.plot(x[cat], y[cat], z[cat], mark[i], cmap=colors, label=cat)

        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])
        ax.set_zlabel(axes[2])
        plt.legend(loc='lower left', fancybox=True, shadow=True)

    return plt