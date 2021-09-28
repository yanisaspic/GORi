"""
Data encoder and stylesheet for Dash Cytoscape.

@ ASLOUDJ Yanis
@ July 2021
"""

import pandas as pd
import dash_html_components as html

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

def dictize_node(item_id, items_metrics_row, node_size, colors):
    """
    # Description
    Converts an items dataframe's row as a dict of parameters used by Cytoscape.

    # Arguments
    ``item_id`` (string): the unique identifier of the ontology term. \n
    ``items_metrics`` (df): the items id, weight, label and frequency. \n
    ``node_size`` (int): the size of a node. \n
    ``colors`` (dict): the keys are term prefixes and the values are the corresponding node colors. \n

    # Usage
    >>> data_items = pd.DataFrame(
        {'item': ['Alpha', 'Beta', 'Gamma', 'Delta'], 'freq': [0.8, 0.7, 0.3, 0.1],
        'label': ['A', 'B', 'C', 'D'], 'weight': [0.3, 0.3, 0.8, 0.95]})
    >>> target_item = 'Alpha'
    >>> items_colors = {
        'Al': [255, 0, 180], 'Be': [0, 0, 0], 'Ga': [75, 75, 75], 'De': [100, 200, 150] }
    >>> data_items_rowA = data_items.loc[data_items['item'] == target_item]
    >>> print(dictize_node(target_item, data_items_rowA, 10, items_colors))
    ... {
        'id': 'Alpha', 'label': 'A', 'size': 8.0, 'font': 0.8, 'line': 7.2, 
        'color': [22.95, 0.0, 16.2], 'weight': 0.3, 'freq': 0.8, 'connect': 0}
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
            'weight': weight, 'freq': freq, 'connect': 0
    }
    return node

def dataframize_data(rules_metrics, items_metrics, node_size=200, colors = {
        'GO': [0, 255, 255],
        'R-': [255, 0, 255],
        'HP': [255, 255, 0]
    }):
    """
    # Description
    Converts rules and items dataframes as nodes and edges dataframes with parameters used by Cytoscape.

    # Arguments
    ``rules_metrics`` (df): association rules with a body, a head, confidence and lift values. \n
    ``items_metrics`` (df): the items names and their weights, readable labels and frequencies. \n
    ``node_size`` (int): the default size of a node in your network. \n
    ``colors`` (dict): the keys are term prefixes and the values are the corresponding node colors. \n

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
            colors={'Al': [0,0,0], 'Be': [200,200,200], 'Ga': [100,0,0]})
    >>> print(nodes)
    ... id label size font line             color weight freq connect
    0  Alpha    Al  160   16  144   [0.0, 0.0, 0.0]    0.3  0.8       2
    1   Beta    Bl  140   14  126   [8.0, 8.0, 8.0]    0.2  0.7       3
    2  Gamma    Cl   60    6   54  [64.0, 0.0, 0.0]    0.8  0.3       3
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
            nodes[term] = nodes.get(term, dictize_node(term, term_metrics, node_size, colors))
            nodes[term]['connect'] += 1

        # get the corresponding edge
        edge_id = " => ".join( sorted([body, head]) )
        edge = {
                'id': edge_id, 'source': body, 'target': head, 
                'size': round(conf**2 * min_node_size, 1),
                'color': [ e - 150 for e in colors[body[:2]] ],
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

def zoom_node(main_node_id, nodes, edges):
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

def zoom_edge(main_edge_id, nodes, edges):
    """
    # Description
    Returns Dash Cytoscape readable data to draw a graph of two connected nodes.

    # Arguments
    ``main_edge_id`` (str): id of the existing edge. \n
    ``nodes`` (df): the nodes ids, labels, weights, frequencies, etc. \n
    ``edges`` (df): the edges ids, labels, lifts, confidences, colors, etc.
    """
    key_edge_row = edges[edges['id'] == main_edge_id]
    key_edge_data = encode_data(key_edge_row)

    # get the two connected nodes
    nodes_ids = main_edge_id.split(" => ")
    key_nodes_df = nodes[nodes['id'].isin(nodes_ids)]
    key_nodes_data = encode_data(key_nodes_df)

    return key_edge_data + key_nodes_data

def htmlize_elem(main_elem_id, nodes, edges):
    """
    # Description
    Returns the data of a network element as a list of HTML text elements.

    # Arguments:
    ``main_elem_id`` (str): id of the existing node or edge. \n
    ``nodes`` (df): the nodes ids, labels, weights, frequencies, etc. \n
    ``edges`` (df): the edges ids, labels, lifts, confidences, colors, etc.
    """
    links = {
        'GO': 'https://www.ebi.ac.uk/QuickGO/term/',
        'R-': 'https://reactome.org/PathwayBrowser/#/',
        'HP': 'https://hpo.jax.org/app/browse/term/'
    }

    if " => " not in main_elem_id:
        node_id = main_elem_id
        node_row = nodes.loc[ nodes['id'] == node_id ]
        onto_href = html.A(
            html.B("%s (%s)" % ( node_row['label'].values[0], node_row['id'].values[0]) ),
                href=links[node_id[:2]] + node_id, target="_blank", id="raw-data-elem")
        node_metrics = html.Div([
                html.P( [ html.B("Weight : "), round(node_row['weight'].values[0], 2) ]), 
                html.P( [ html.B("Frequency : "), round(node_row['freq'].values[0], 2) ]), 
                html.P( [ html.B("Connectivity : "), round(node_row['connect'].values[0] / 2, 2) ] ) 
            ], id="raw-data-metrics")

        data_children = [onto_href, node_metrics]
    
    return html.Div(data_children)