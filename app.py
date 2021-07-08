"""
Prepares the graph highlighting the rules terms, their confidences, lift and coverage values.

@ June 2021
@ ASLOUDJ Yanis
"""

from settings import *
from math import inf
from _scripts.misc import truncate
import pandas as pd
import _scripts.network as nw
import json

# Dash interactive network
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import plotly.express as px
import dash

# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -
# - - - - - - - - -- - - - - - -- - - - - - -- - -- - - -- - -- - - -- - -- - -- - - -- - - -- - -  - - - -

#_________________________________________ N E T W O R K __________________________________________

#_________________________________________ L O A D I N G

cyto.load_extra_layouts()

data_items = pd.read_csv(items_file)
data_rules = pd.read_csv(rules_file)
nodes, edges = nw.dataframize_data(data_rules, data_items, 200)
interest_nodes, interest_edges = nodes.copy(), edges.copy()

# default values for the main network
data = nw.encode_data(interest_nodes) + nw.encode_data(interest_edges)    

# default node for the zoom network is the central one
nodes = nodes.sort_values(by='connect', ascending=False)
zoom_node_id = nodes.iloc[0]['id']
default_zoom_data = nw.subgraph(zoom_node_id, interest_nodes, interest_edges)

app = dash.Dash(__name__, title="GORi")

### WEBAPP SECTIONS ###

### 1. ONTOLOGY NETWORK 
network = cyto.Cytoscape(id='main_cytoscape', elements=data, layout={'name': 'concentric'}, style={}, stylesheet=nw.style)

### 2. WEBAPP MENU

### PARAMETERS

# Layout
layout_dropdown = dcc.Dropdown(
    id='dropdown-layout', className='long-dropdown', value='concentric', clearable=False,
    options=[ 
        {'label': name.capitalize(), 'value': name} 
        for name in ['cola', 'concentric', 'euler', 'klay', 'random'] ], style={}
)

# Ontology
elements_dropdown = dcc.Dropdown(
    id='dropdown-elements', className='long-dropdown', options=[
        {'label': 'Gene Ontology', 'value': 'id@GO'},
        {'label': 'Reactome', 'value': 'id@R-'},
        {'label': 'Human Phenotype Ontology', 'value': 'id@HP'},
        {'label': 'One-sided', 'value': 'effect@hidden'},
        {'label': 'Positive (two-sided)', 'value': 'effect@pos'},
        {'label': 'Negative (two-sided)', 'value': 'effect@neg'},
        {'label': 'Indefinite (two-sided)', 'value': 'effect@ind'}
    ], value = ['id@GO', 'id@R-', 'id@HP', 'effect@pos', 'effect@neg'],
    multi=True, style={}, clearable=False, placeholder = ""
)

sliders = []
df_metrics = {
    'nodes': [nodes, ['freq', 'weight'], ['Term\nFrequency', 'Term\nWeight']],
    'edges': [edges, ['conf', 'lift'], ['Association\nConfidence', 'Association\nLift']]
}
for elem in df_metrics.keys():
    df, categories, labels = df_metrics[elem][0], df_metrics[elem][1], df_metrics[elem][2]
    for i in range(len(categories)):
        cat = categories[i]
        mini, maxi = float( truncate(df[cat].min(), 2) ), float( truncate(df[cat].max(), 2) )

        slid = dcc.RangeSlider(
            id="slider-%s" % cat,
            className='cat-slider',
            min = mini,
            max = maxi,
            step = float( truncate( (maxi - mini)/10, 2) ),
            value = [mini, maxi],
            marks={
                mini: {'label': ""},
                maxi: {'label': ""}
            },
            allowCross=False
        )
        sliders.append(html.Div(id=cat, children=[ html.P( labels[i], className="slider-title" ), slid]) )
param_sliders = html.Div(id="sliders", children=sliders, style={})

parameters = html.Div(
    id="parameters", style={},
    children=[
        html.H2(html.Span('Parameters', style={}), style={}), 
        html.P('Select a layout to display the associations of interest :', style={}), layout_dropdown,
        html.P('Select adjectives describing the associations of interest :', style={}), elements_dropdown,
        html.P('Select metrics ranges including the associations of interest :', style={}), param_sliders])

### 3. WEBAPP ZOOM ###
zoom_network = cyto.Cytoscape(
    id='zoom_cytoscape', elements=default_zoom_data, 
    layout={'name': 'concentric', 'animate': False}, style={}, stylesheet=nw.style)
data = html.Div(
    id="data", style={},
    children=[
        html.H2(html.Span('Data', style={}), style={})
    ]
)

menu = html.Div(id='menu', children=[html.H1('GORi', style={}), parameters, data], style={})

### WEBPAGE INTERACTIONS

@app.callback(
    Output('main_cytoscape', 'elements'),
    Input('dropdown-elements', 'value'))
def update_elements(association_adjectives):
    global interest_nodes
    global interest_edges

    # separe the adjectives in categories
    adjectives = {
        'id': [],
        'effect': []
    }

    for adj in association_adjectives:
        split_adj = adj.split("@")
        adjectives[split_adj[0]].append(split_adj[1])
    valid_prefixes = "|".join(adjectives['id'])

    # filter the data for each category
    interest_nodes = nodes.loc[ nodes['id'].str.contains(valid_prefixes) ]
    interest_edges = edges.loc[ edges['source'].str.contains(valid_prefixes) & \
        edges['target'].str.contains(valid_prefixes)]
    interest_edges = interest_edges.loc[ interest_edges['effect'].isin(adjectives['effect']) ]

    # encode the data and use it as cytoscape elements
    data = nw.encode_data(interest_nodes) + nw.encode_data(interest_edges)
    return data

@app.callback(
    Output('zoom_cytoscape', 'elements'),
    Input('main_cytoscape', 'tapNodeData'),
    Input('zoom_cytoscape', 'tapNodeData'),
    Input('main_cytoscape', 'elements'))
def display_node(tapNodeData_main, tapNodeData_zoom, priority_holder):
    # last input added to make sure that the elements are filtered before running this callback
    global zoom_node_id

    context = dash.callback_context
    trigger = context.triggered[0]['prop_id']
    if trigger == 'main_cytoscape.tapNodeData':
        zoom_node_id = tapNodeData_main['id']
    elif trigger == 'zoom_cytoscape.tapNodeData':
        zoom_node_id = tapNodeData_zoom['id']

    return nw.subgraph(zoom_node_id, interest_nodes, interest_edges)

@app.callback(
    Output('main_cytoscape', 'layout'),
    Input('dropdown-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': True
    }

### APP LAUNCHER

app.layout = html.Div(id="webapp",
    children=[network, zoom_network, menu])

if __name__ == '__main__':
    app.run_server(debug=True)