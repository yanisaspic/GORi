"""
Prepares the graph highlighting the rules terms, their confidences, lift and coverage values.

@ June 2021
@ ASLOUDJ Yanis
"""

from settings import *
import pandas as pd

# Dash interactive network
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import dash_cytoscape as cyto
import plotly.express as px
import dash

def encode_node(item_id, items_metrics_row, node_size):
    """
    # Description
    Converts an items dataframe's row as a readable dict used by Cytoscape to create a node.

    # Arguments
    ``item_id`` (string): the name 
    ``items_metrics`` (df): the items names and their frequencies. \n
    ``node_size`` (int): the size of a node.

    # Usage
    >>> data_items = pd.DataFrame(
        {'item': ['Alpha', 'Beta', 'Gamma', 'Delta'], 'freq': [0.8, 0.7, 0.3, 0.1],
        'label': ['A', 'B', 'C', 'D']})
    >>> target_item = 'Alpha'
    >>> data_items_rowA = data_items.loc[data_items['item'] == target_item]
    >>> print(encode_node(target_item, data_items_rowA, 10))
    ... {'data': {'id': 'Alpha'}, 'label': 'A (Alpha)', 'classes': 'Al', 'size': 8.0}
    """
    node = {
        'data': {'id': item_id},
        'label': "%s (%s)" % (items_metrics_row['label'].values[0], item_id),
        'classes': item_id[:2],
        'size': node_size * items_metrics_row['freq'].values[0]
    }
    return node

def encode_data(rules_metrics, items_metrics, node_size=50, edge_size=10):
    """
    # Description
    Converts a rules and an items dataframe as readable data dicts used by Cytoscape to create a network.

    # Arguments
    ``rules_metrics`` (df): association rules with a body, a head, confidence and lift values. \n
    ``items_metrics`` (df): the items names and their frequencies. \n
    ``node_size`` (int): the default size of a node in your network. \n
    ``edge_size`` (int): the default size of an edge in your network.

    # Usage
    >>> data_rules = pd.DataFrame({
        'body': ['Alpha', 'Alpha', 'Beta', 'Gamma', 'Delta'],
        'head': ['Beta', 'Gamma', 'Delta', 'Beta', 'Gamma'],
        'conf': [0.8, 0.7, 0.3, 0.6, 1], 'lift': [1.5, 0.3, 17, 4, 9]})
    >>> data_items = pd.DataFrame({
        'item': ['Alpha', 'Beta', 'Gamma', 'Delta'], 
        'freq': [0.8, 0.7, 0.3, 0.1],
        'label': ['A', 'B', 'C', 'D']})
    >>> print(encode_data(data_rules, data_items))
    ... [[{'data': {'id': 'Alpha'}, 'label': 'A (Alpha)', 'classes': 'Al', 'size': 40.0}, {'data': {'id': 'Beta'}, 'label': 'B (Beta)', 'classes': 'Be', 'size': 35.0}, {'data': {'id': 'Gamma'}, 'label': 'C (Gamma)', 'classes': 'Ga', 'size': 15.0}, {'data': {'id': 'Delta'}, 'label': 'D (Delta)', 'classes': 'De', 'size': 5.0}, {'data': {'source': 'Alpha', 'target': 'Beta'}, 'rel_lift': 0.08823529411764706, 'size': 8.0}, {'data': {'source': 'Alpha', 'target': 'Gamma'}, 'rel_lift': -1.0, 'size': 7.0}, {'data': {'source': 'Beta', 'target': 'Delta'}, 'rel_lift': 1.0, 'size': 3.0}, {'data': {'source': 'Gamma', 'target': 'Beta'}, 'rel_lift': 0.23529411764705882, 'size': 6.0}, {'data': {'source': 'Delta', 'target': 'Gamma'}, 'rel_lift': 0.5294117647058824, 'size': 10.0}]]
    """
    nodes = {}
    edges = []

    max_lift, min_lift = rules_metrics['lift'].max(), rules_metrics['lift'].min()
    print(max_lift, min_lift)

    for index, rule_row in rules_metrics.iterrows():
        body, head = rule_row['body'], rule_row['head']
        
        # get the 2 corresponding nodes
        body_metrics = items_metrics.loc[items_metrics['item'] == body]
        head_metrics = items_metrics.loc[items_metrics['item'] == head]
        nodes[body] = nodes.get(body, encode_node(body, body_metrics, node_size))
        nodes[head] = nodes.get(head, encode_node(head, head_metrics, node_size))

        # get the corresponding edge
        absolute_lift = rule_row['lift']
        if absolute_lift > 1:
            relative_lift = absolute_lift/max_lift
        else:
            relative_lift = -(1 - absolute_lift) / (1 - min_lift)
        
        edge = {
            'data': {'source': body, 'target': head}, 
            'rel_lift': relative_lift, 'size': edge_size * rule_row['conf']}
        edges.append(edge)
    
    return edges + list(nodes.values())

# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -

#_________________________________________ N E T W O R K __________________________________________

#_________________________________________ L O A D I N G

stylesheet = [

    # Group selectors
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)'
        }
    },

    # Class selectors
    {
        'selector': '.GO',
        'style': {
            'background-color': 'cyan'
        },

        'selector': '.R-',
        'style': {
            'background-color': 'magenta'
        },

        'selector': '.HP',
        'style': {
            'background-color': 'yellow'
        },
        'selector': 'edge',
        'style': {
            'size': 'data(size)'
        }
    }
]

items = pd.read_csv(items_file)
rules = pd.read_csv(rules_file)

network_data = encode_data(rules, items)

app = dash.Dash(__name__)

