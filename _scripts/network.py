"""
Data encoder and stylesheet for Dash Cytoscape.

@ ASLOUDJ Yanis
@ July 2021
"""

import pandas as pd
from math import pi, cos, sin

### CYTOSCAPE NETWORK STYLESHEET
style = [

    {
        'selector': 'node',
        'style': {
            # label-related
            'content': 'data(label)',
            "text-valign": "center",
            "text-halign": "center",
            'color': '#202020',
            'font-size': 'data(font)',
            'text-wrap': 'wrap',
            'text-max-width': 'data(line)',
            'font-weight': 'bolder',
            'text-background-color': 'data(color)',
            'text-background-opacity': 1,
            'text-background-padding': '1rem',
            'text-background-shape': 'round-rectangle',

            # node-related
            'width': 'data(size)',
            'height': 'data(size)',
            'background-color': 'data(color)'
        }
    },

    {
        'selector': 'edge',
        'style': {
            'width': 'data(size)',
            'line-color': 'data(color)'
        }
    }

]

def dictize_node(item_id, items_metrics_row, node_size, colors, links):
    """
    # Description
    Converts an items dataframe's row as a dict of parameters used by Cytoscape.

    # Arguments
    ``item_id`` (string): the unique identifier of the ontology term. \n
    ``items_metrics`` (df): the items id, weight, label and frequency. \n
    ``node_size`` (int): the size of a node. \n
    ``colors`` (dict): the keys are term prefixes and the values are the corresponding node colors. \n
    ``links`` (dict): the keys are term prefixes and the values are links to display the ontology.

    # Usage
    >>> data_items = pd.DataFrame(
        {'item': ['Alpha', 'Beta', 'Gamma', 'Delta'], 'freq': [0.8, 0.7, 0.3, 0.1],
        'label': ['A', 'B', 'C', 'D'], 'weight': [0.3, 0.3, 0.8, 0.95]})
    >>> target_item = 'Alpha'
    >>> items_colors = {
        'Al': [255, 0, 180], 'Be': [0, 0, 0], 'Ga': [75, 75, 75], 'De': [100, 200, 150] }
    >>> items_links = {'Al': 'https://', 'Be': 'www.', 'Ga': 'http://', 'De': 'www.placeholder/'}
    >>> data_items_rowA = data_items.loc[data_items['item'] == target_item]
    >>> print(dictize_node(target_item, data_items_rowA, 10, items_colors, items_links))
    ... {
        'id': 'Alpha', 'label': 'A', 'size': 8.0, 'font': 0.8, 'line': 7.2, 
        'color': [22.95, 0.0, 16.2], 'href': 'https://Alpha', 
        'weight': 0.3, 'freq': 0.8, 'connect': 0}
    """
    freq = items_metrics_row['freq'].values[0]
    weight = items_metrics_row['weight'].values[0]
    term_prefix = item_id[:2]

    node = {
            'id': item_id,
            'label': items_metrics_row['label'].values[0], 
            'size': round(node_size * freq, 1),
            'font': round(node_size * 0.1 * freq, 1),
            'line': round(node_size * 0.9 * freq, 1),
            'color': [ round(e*weight**2, 1) for e in colors[term_prefix] ],
            'href': links[term_prefix] + item_id,
            'weight': weight, 'freq': freq, 'connect': 0
    }
    return node

def dataframize_data(rules_metrics, items_metrics, node_size=200, colors = {
        'GO': [0, 255, 255],
        'R-': [255, 0, 255],
        'HP': [255, 255, 0]
    }, links = {
        'GO': 'ebi.ac.uk/QuickGO/term/',
        'R-': 'https://reactome.org/PathwayBrowser/#/',
        'HP': 'https://hpo.jax.org/app/browse/term/'
    }):
    """
    # Description
    Converts rules and items dataframes as nodes and edges dataframes with parameters used by Cytoscape.

    # Arguments
    ``rules_metrics`` (df): association rules with a body, a head, confidence and lift values. \n
    ``items_metrics`` (df): the items names and their weights, readable labels and frequencies. \n
    ``node_size`` (int): the default size of a node in your network. \n
    ``colors`` (dict): the keys are term prefixes and the values are the corresponding node colors. \n
    ``links`` (dict): the keys are term prefixes and the values are links to display the ontology.

    # Usage
    >>> data_rules = pd.DataFrame({
            'body': ['Alpha', 'Alpha', 'Beta', 'Gamma'],
            'head': ['Beta', 'Gamma', 'Gamma', 'Beta'],
            'conf': [0.8, 0.7, 0.6, 1], 'lift': [1.5, 0.3, 4, 9]})
    >>> data_items = pd.DataFrame({
            'item': ['Alpha', 'Beta', 'Gamma'], 
            'freq': [0.8, 0.7, 0.3],
            'weight': [0.3, 0.2, 0.8],
            'label': ['Al', 'Bl', 'Cl']})
    >>> nodes, edges = dataframize_data(
            data_rules, data_items, 
            colors={'Al': [0,0,0], 'Be': [200,200,200], 'Ga': [100,0,0]},
            links={'Al': 'https://whatever', 'Be': 'www.idk', 'Ga': 'http://ok'})
    >>> print(nodes)
    ... id label size font  line             color                   href weight freq connect
    0  Alpha    Al   16  1.6  14.4   [0.0, 0.0, 0.0]  https://whateverAlpha    0.3  0.8       2
    1   Beta    Bl   14  1.4  12.6   [8.0, 8.0, 8.0]            www.idkBeta    0.2  0.7       3
    2  Gamma    Cl    6  0.6   5.4  [64.0, 0.0, 0.0]         http://okGamma    0.8  0.3       3
    >>> print(edges)
    ...             color conf connect              id lift rev_conf size source target
    0  [-235, -235, -235]  0.8       1   Alpha => Beta  1.5     None  3.8  Alpha   Beta
    1  [-235, -235, -235]  0.7       1  Alpha => Gamma  0.3     None  2.9  Alpha  Gamma
    2   [0.0, 226.7, 0.0]  0.6       2   Beta => Gamma    4        1  0.6   Beta  Gamma
    """
    nodes, single_edges, two_way_edges, two_way_int_edges = {}, {}, {}, {}

    min_node_size = items_metrics['freq'].min() * node_size
    min_edge_size = min_node_size * 0.1
    
    two_way_lift = {
        'pos': [0, 255, 0],
        'ind': [0, 0, 0],
        'neg': [255, 0, 0]
    }

    for index, rule_row in rules_metrics.iterrows():

        body, head = rule_row['body'], rule_row['head'] 
        lift, conf = rule_row['lift'], rule_row['conf']

        # get the 2 corresponding nodes
        for term in [body, head]:
            term_metrics = items_metrics.loc[items_metrics['item'] == term]
            nodes[term] = nodes.get(term, dictize_node(term, term_metrics, node_size, colors, links))
            nodes[term]['connect'] += 1

        # get the corresponding edge
        edge_id = " => ".join( sorted([body, head]) )
        edge = {
                'id': edge_id, 'source': body, 'target': head, 
                'size': round(conf**2 * min_node_size, 1),
                'color': [ e - 235 for e in colors[body[:2]] ],
                'lift': lift, 'conf': conf, 'connect': 0, 'rev_conf': None, 'effect': 'hidden'
        }
        single_edges[edge_id] = single_edges.get(edge_id, edge)
        single_edges[edge_id]['connect'] += 1

        if single_edges[edge_id]['connect'] == 2:

            if lift < 0.5:
                mult_lift, two_way_sign = lift, 'neg'
                two_way_cat = two_way_int_edges
            elif lift > 2:
                mult_lift, two_way_sign = 1/lift, 'pos'
                two_way_cat = two_way_int_edges
            else:
                mult_lift, two_way_sign = 1, 'ind'
                two_way_cat = two_way_edges

            # remove the size calculated with confidence since there are 2 conf values for this edge
            single_edges[edge_id]['rev_conf'] = edge['conf']
            single_edges[edge_id]['effect'] = two_way_sign
            single_edges[edge_id]['size'] = min_edge_size
            two_way_cat[edge_id] = single_edges.pop(edge_id)

            # special edge color indicating if the item combination is likely (green) or unlikely (red)
            two_way_cat[edge_id]['color'] = [ 
                round(e - e * mult_lift, 1) for e in two_way_lift[two_way_sign] ]

    nodes_df = pd.DataFrame.from_records(nodes).transpose()
    nodes_df = nodes_df.reset_index(drop=True)
    edges_df = pd.concat([
        pd.DataFrame.from_records(single_edges).transpose(), 
        pd.DataFrame.from_records(two_way_edges).transpose(), 
        pd.DataFrame.from_records(two_way_int_edges).transpose() ])
    edges_df = edges_df.reset_index(drop=True)
    
    return nodes_df, edges_df

def encode_data(data):
    """
    # Description
    Returns a list of Cytoscape-readable data dict from a dataframe.

    # Arguments
    ``data`` (df): a pandas dataframe.

    # Usage
    >>> df = pd.DataFrame({'name': ['Anto', 'Beatrice'], 'age': [22, 34]})
    >>> print(encode_data(df))
    ... [{'data': {'name': 'Anto', 'age': 22}}, {'data': {'name': 'Beatrice', 'age': 34}}]
    """
    return [{'data': x} for x in data.to_dict('records')]

def subgraph(main_node_id, nodes, edges):
    """
    # Description
    Returns Dash Cytoscape readable data to draw a concentric graph centered around one node.

    # Arguments
    ``main_node_id`` (str): id of the existing central node. \n
    ``nodes`` (df): the nodes ids, labels, weights, frequencies, etc. \n
    ``edges`` (df): the edges ids, labels, lifts, confidences, colors, etc.
    """
    key_edges_df = edges.loc[ (edges['source'] == main_node_id) | (edges['target'] == main_node_id) ]
    key_edges_data = encode_data(key_edges_df)

    # main node manually added to display the node only in case there are no edges
    key_nodes = set(key_edges_df['source']) | set(key_edges_df['target']) | set([main_node_id])
    key_nodes_df = nodes.loc[ nodes['id'].isin(key_nodes) ]
    key_nodes_data = encode_data(key_nodes_df)

    return key_nodes_data + key_edges_data