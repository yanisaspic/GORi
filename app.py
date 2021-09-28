"""
Prepares the graph highlighting the rules terms, their confidences, lift and coverage values.

@ June 2021
@ ASLOUDJ Yanis
"""

from settings import *
from random import randint
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
last_tap_event = None

# default values for the main network
data = nw.encode_data(interest_nodes) + nw.encode_data(interest_edges)    

# default node for the zoom network is the central one
x = randint(0, nodes.shape[0] - 1)
zoom_key_id = nodes.iloc[x]['id']
default_zoom_data = nw.zoom_node(zoom_key_id, interest_nodes, interest_edges)

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
        for name in ['concentric', 'euler', 'klay', 'random'] ], style={}
)

# Ontology
elements_dropdown = dcc.Dropdown(
    id='dropdown-elements', className='long-dropdown', options=[
        {'label': 'Gene Ontology', 'value': 'id@GO'},
        {'label': 'Reactome', 'value': 'id@R-'},
        {'label': 'Human Phenotype Ontology', 'value': 'id@HP'},
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
    df, categories, labels = df_metrics[elem]
    for i in range(len(categories)):
        cat = categories[i]
        mini, maxi = df[cat].min(), df[cat].max()

        slid = dcc.RangeSlider(
            id="slider-%s" % cat,
            className='cat-slider',
            min = mini,
            max = maxi,
            step = (maxi - mini)/20,
            value = [mini, maxi],
            marks={
                mini: {'label': ""},
                maxi: {'label': ""}
            },
            allowCross=False,
        )
        sliders.append(html.Div(id=cat, children=[ 
            html.P( labels[i], className="slider-title" ), 
            slid,
            html.P(round(mini, 2), id="slider-min-%s" % cat, className="slider-value-min slider-value"),
            html.P(round(maxi, 2), id="slider-max-%s" % cat, className="slider-value-max slider-value") 
        ]))
param_sliders = html.Div(id="sliders", children=sliders, style={})

parameters = html.Div(
    id="parameters", style={},
    children=[
        html.H2(html.Span('Parameters', style={}), style={}), 
        html.P('Select a layout to display the associations of interest :', style={}), layout_dropdown,
        html.P('Select keywords describing the associations of interest :', style={}), elements_dropdown,
        html.P('Select metrics ranges including the associations of interest :', style={}), param_sliders])

data_raw = html.Div(id="data", children=[
    html.H2(html.Span('Data', style={}), style={}),
    html.Div(id="raw-data", children=nw.htmlize_elem(zoom_key_id, nodes, edges))]
    )

menu = html.Div(id='menu', children=[
    html.H1('GORi', style={}), parameters, data_raw], style={})

### 3. WEBAPP ZOOM ###
zoom_network = cyto.Cytoscape(
    id='zoom_cytoscape', elements=default_zoom_data, 
    layout={'name': 'concentric', 'animate': True}, style={}, stylesheet=nw.style)

### WEBPAGE INTERACTIONS

@app.callback(
    Output('main_cytoscape', 'elements'),
    Input('dropdown-elements', 'value'),
    Input('slider-freq', 'value'),
    Input('slider-weight', 'value'),
    Input('slider-conf', 'value'),
    Input('slider-lift', 'value'))
def update_main(association_adjectives, freq_lim, weight_lim, conf_lim, lift_lim):
    global interest_nodes, interest_edges

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
    interest_nodes = nodes[ nodes['id'].str.contains(valid_prefixes) ]
    for x, y in zip(['freq', 'weight'], [freq_lim, weight_lim]):
        interest_nodes = interest_nodes[ interest_nodes[x].between(y[0], y[1], inclusive=True) ]
    valid_nodes_ids = list(interest_nodes['id'])


    interest_edges = edges[ 
        edges['source'].isin(valid_nodes_ids) & edges['target'].isin(valid_nodes_ids)]
    interest_edges = interest_edges[
        interest_edges['lift'].between(lift_lim[0], lift_lim[1], inclusive=True) ]
    interest_edges = interest_edges[
        ( interest_edges['conf'].between(conf_lim[0], conf_lim[1], inclusive = True ) )
        | ( interest_edges['rev_conf'].between(conf_lim[0], conf_lim[1], inclusive = True) )
    ]
    interest_edges = interest_edges[ interest_edges['effect'].isin(adjectives['effect']) ]

    # encode the data and use it as cytoscape elements
    return nw.encode_data(interest_nodes) + nw.encode_data(interest_edges)

@app.callback(
    Output('zoom_cytoscape', 'elements'),
    Input('main_cytoscape', 'tapNodeData'),
    Input('zoom_cytoscape', 'tapNodeData'),
    Input('main_cytoscape', 'tapEdgeData'),
    Input('zoom_cytoscape', 'tapEdgeData'),
    Input('main_cytoscape', 'elements'))
def update_zoom(tapNodeData_main, tapNodeData_zoom, tapEdgeData_main, tapEdgeData_zoom, holder):
    
    # last input added to make sure that the elements are filtered before running this callback
    global zoom_key_id, last_tap_event
    context = dash.callback_context
    trigger = context.triggered[0]['prop_id']
    origin, event = trigger.split(".")

    if event == 'elements':
        if last_tap_event == 'tapEdgeData':
            return nw.zoom_edge(zoom_key_id, interest_nodes, interest_edges)
        else:
            return nw.zoom_node(zoom_key_id, interest_nodes, interest_edges)

    if event == 'tapNodeData':
        last_tap_event = 'tapNodeData'
        if origin == 'main_cytoscape':
            zoom_key_id = tapNodeData_main['id']
        else:
            zoom_key_id = tapNodeData_zoom['id']
        return nw.zoom_node(zoom_key_id, interest_nodes, interest_edges)

    if event == 'tapEdgeData':
        last_tap_event = 'tapEdgeData'
        if origin == 'main_cytoscape':
            zoom_key_id = tapEdgeData_main['id']
        else:
            zoom_key_id = tapEdgeData_zoom['id']
        return nw.zoom_edge(zoom_key_id, interest_nodes, interest_edges)

@app.callback(
    Output('main_cytoscape', 'layout'),
    Input('dropdown-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': True
    }

@app.callback(
    Output('raw-data', 'children'),
    Input('main_cytoscape', 'tapNodeData'),
    Input('zoom_cytoscape', 'tapNodeData')
)
def update_data(tapNodeData_main, tapNodeData_zoom):

    context = dash.callback_context
    trigger = context.triggered[0]['prop_id']
    origin, event = trigger.split(".")

    if origin == 'main_cytoscape':
        elem_id = tapNodeData_main['id']
    elif origin == 'zoom_cytoscape':
        elem_id = tapNodeData_zoom['id']
    else:
        elem_id = zoom_key_id
    return nw.htmlize_elem(elem_id, nodes, edges)

@app.callback(
    Output('slider-min-freq', 'children'),
    Output('slider-max-freq', 'children'),
    Input('slider-freq', 'value')
)
def update_slider(interval):
    return round(interval[0], 2), round(interval[1], 2)

@app.callback(
    Output('slider-min-weight', 'children'),
    Output('slider-max-weight', 'children'),
    Input('slider-weight', 'value')
)
def update_slider(interval):
    return round(interval[0], 2), round(interval[1], 2)

@app.callback(
    Output('slider-min-conf', 'children'),
    Output('slider-max-conf', 'children'),
    Input('slider-conf', 'value')
)
def update_slider(interval):
    return round(interval[0], 2), round(interval[1], 2)

@app.callback(
    Output('slider-min-lift', 'children'),
    Output('slider-max-lift', 'children'),
    Input('slider-lift', 'value')
)
def update_slider(interval):
    return round(interval[0], 2), round(interval[1], 2)



### APP LAUNCHER

app.layout = html.Div(id="webapp",
    children=[network, zoom_network, menu])

if __name__ == '__main__':
    app.run_server(debug=True)