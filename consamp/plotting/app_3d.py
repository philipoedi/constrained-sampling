import os
import glob


import numpy as np

import dash
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
import plotly.express as px
from dash.dependencies import Input, Output
import plotly.graph_objects as go

from app_utils import *

root_folder = os.path.abspath(os.path.join(os.path.dirname("__file__"),"..",".."))
results_folder = os.path.join(root_folder, "results")
global_selected_buttons = {"samples": False,
                    "seeds": False,
                    "surface": True,
                    "projections": False}
local_selected_buttons = global_selected_buttons.copy()
global_checklist = [{"label":k,"value":k} for k in global_selected_buttons.keys()]
local_checklist = [{"label":k,"value":k} for k in local_selected_buttons.keys()]

experiments_path = glob.glob(os.path.join(results_folder,"*/"))
experiments_name = [os.path.split(os.path.split(p)[0])[-1] for p in experiments_path] 
experiments_map = {experiments_name[i]: experiments_path[i] for i in range(len(experiments_name))}
experiments_dropdown_map = [{"label": e,"value": e} for e in experiments_name]


experiments_dropdown_map_right = experiments_dropdown_map.copy()
global_selected_buttons_right = global_selected_buttons.copy()
local_selected_buttons_right = local_selected_buttons.copy()
global_checklist_right = global_checklist.copy()
local_checklist_right = local_checklist.copy()


print(experiments_name)

app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP],
    meta_tags=[{"name":"viewport","content":"width=device-width, initial-scale=1"}])
"""
app.layout = html.Div(children=[
  html.Div(className="header", children=[
    html.P("experiment"),
    html.Div(className="dropdown", children=[
        dcc.Dropdown(id="experimentselect",
                     options=experiments_dropdown_map,
                     multi=False,
                     value="",
                     className="experimentselect"),
        dcc.Checklist(id="global_checklist",
                    options=global_checklist,
                    value=[k for k in global_selected_buttons.keys()]),
        dcc.Checklist(id="local_checklist",
                    options=local_checklist,
                    value=[k for k in local_selected_buttons.keys()])]),
  html.Div(className="graphs", children=[
         html.Div(className="graph1div", children=[dcc.Graph(id="graph")])]
         )
  ])])
"""

app.layout = dbc.Container(id="main-container", children=[
    dbc.Row(id="experiment-select-row", children=[
        dbc.Col([
            dcc.Dropdown(id="experiment-select-left", multi=False,
                options=experiments_dropdown_map,
                value="")], width={"size":6,"order":1}),
        dbc.Col([
            dcc.Dropdown(id="experiment-select-right", multi=False,
                options=experiments_dropdown_map_right,
                value="")], width={"size":6,"order":2}),
        ]),
    dbc.Row(id="show-options", children=[
        dbc.Col([
            dcc.Checklist(id="global_checklist_left",
                options=global_checklist,
                #value=[k for k in global_selected_buttons.keys()]),
                value=[]),
            dcc.Checklist(id="local_checklist_left",
                options=local_checklist,
                value=[])], width={"size":6}),
        dbc.Col([
            dcc.Checklist(id="global_checklist_right",
                options=global_checklist_right,
                value=[]),
            dcc.Checklist(id="local_checklist_right",
                options=local_checklist_right,
                value=[])], width={"size":6}),
        ]),
    dbc.Row(id="graph-left-row", children=[
        dbc.Col([
            dcc.Graph(id="graph-left")
            ], width={"size":6,"order":1}),
        dbc.Col([
            dcc.Graph(id="graph-right")
            ], width={"size":6,"order":2})
        ]),
    dbc.Row(id="hist-row", children=[
        dbc.Col([
            dcc.Graph(id="hist-left")
            ], width={"size":6,"order":1}),
        dbc.Col([
            dcc.Graph(id="hist-right")
            ], width={"size":6,"order":2})
        ]),
     dbc.Row(id="statistics-row", children=[
        dbc.Col([
            dcc.Graph(id="stats-left")
            ], width={"size":6,"order":1}),
        dbc.Col([
            dcc.Graph(id="stats-right")
            ], width={"size":6,"order":2})
        ])
 
    ],fluid=True)




@app.callback(Output("graph-left","figure"),
                [Input("experiment-select-left","value"),
                Input("global_checklist_left","value"),
                Input("local_checklist_left","value")])
def update_plot(experiment, global_checklist, local_checklist):
    experiment_name = os.path.join(results_folder,experiment)
    print(experiment_name) 
    samples, seeds, pdes = create_dataframe(experiment_name)
    print("finished loading data")
    
    plots = []
    try:
        if "samples" in global_checklist:
            plots.append(get_scatterplot(samples, local=False))
        
        if "seeds" in global_checklist:
            plots.append(get_scatterplot(seeds, local=False))
        
        if "samples" in local_checklist:
            plots.append(get_scatterplot(samples, local=True))
        
        if "seeds" in local_checklist:
            plots.append(get_scatterplot(seeds, local=True))
        if "projections" in global_checklist:
            plots.extend(get_projections(samples, seeds, local=False))

        fig = go.Figure(data=plots).update_layout(autosize=False,height=800,width=1600)
        if "surface" in global_checklist:
            fig.add_trace(get_surfaceplot(pdes)) 
    except:
        pdb.set_trace()
    fig.update_layout(width=1000)
    return fig


@app.callback(Output("graph-right","figure"),
                [Input("experiment-select-right","value"),
                Input("global_checklist_right","value"),
                Input("local_checklist_right","value")])
def update_plot2(experiment, global_checklist_right, local_checklist_right):
    experiment_name = os.path.join(results_folder,experiment)
    samples, seeds, pdes = create_dataframe(experiment_name)
    print("finished loading data")
    
    plots = []
    if "samples" in global_checklist_right:
        plots.append(get_scatterplot(samples, local=False))
    
    if "seeds" in global_checklist_right:
        plots.append(get_scatterplot(seeds, local=False))
    
    if "samples" in local_checklist_right:
        plots.append(get_scatterplot(samples, local=True))
    
    if "seeds" in local_checklist_right:
        plots.append(get_scatterplot(seeds, local=True))

    if "projections" in global_checklist_right:
        plots.extend(get_projections(samples, seeds, local=False))

    fig = go.Figure(data=plots).update_layout(autosize=False,height=800,width=1600)
    if "surface" in global_checklist_right:
        fig.add_trace(get_surfaceplot(pdes)) 
    
    fig.update_layout(width=1000)
    return fig


def load_prob_data(experiment):
    experiment_name = os.path.join(results_folder,experiment,"*probs.dat")
    files = glob.glob(experiment_name)
    data = np.loadtxt(files[0])
    return data 

def update_histogram(experiment):
    data = load_prob_data(experiment)
    plot = get_histogram(data)
    fig = go.Figure(data=[plot])
    return fig

def resub_entropy_estimate(data):
    return -np.log(data).sum() / len(data)


def get_iterations_data(experiment):
    experiment_name = os.path.join(results_folder,experiment,"*num_iterations.dat")
    files = glob.glob(experiment_name)
    data = []
    for f in files:
        sub_data = np.loadtxt(f)
        try:
            if len(sub_data) > 0:
                data.append(sub_data)
        except TypeError:
            data.append(np.array(int(sub_data)).reshape(1,))
    return np.concatenate(data) 

def avg_iterations(experiment):
    data = get_iterations_data(experiment)
    return np.mean(data)


def get_num_global_samples(experiment):
    experiment_name = os.path.join(results_folder,experiment,"*global_seeds.dat")
    files = glob.glob(experiment_name)
    data = np.loadtxt(files[0])
    return len(data)

def get_num_local_seeds(experiment):
    experiment_name = os.path.join(results_folder,experiment,"*num_iterations.dat")
    files = glob.glob(experiment_name)
    return len(files)

@app.callback(Output("hist-left","figure"),
     [Input("experiment-select-left","value")])
def update_histogram_left(experiment):
    return update_histogram(experiment)


@app.callback(Output("hist-right","figure"),
     [Input("experiment-select-right","value")])
def update_histogram_right(experiment):
    return update_histogram(experiment)

def update_statistics(experiment):      
    data = load_prob_data(experiment)
    variance = np.var(data)
    stdev = np.std(data)
    entropy = resub_entropy_estimate(data) 
    num_samples = len(data)
    avg_its = avg_iterations(experiment)
    global_samples = get_num_global_samples(experiment)
    local_seeds = get_num_local_seeds(experiment)
    filtered_seeds = global_samples - local_seeds
    headers = {"values" : ["Metric","Value"]}
    table = go.Table(header=headers,
        cells = {"values":
         [["variance","std","entropy","#samples","avg iterations","# global samples","# local seeds","# removed seeds"],
         [np.round(val,3) for val in [variance,stdev,entropy,num_samples,avg_its, global_samples, local_seeds, filtered_seeds]]]})
    fig = go.Figure(data=[table])
    return fig

@app.callback(Output("stats-left","figure"),
    [Input("experiment-select-left","value")])
def update_stats_left(experiment):
    return update_statistics(experiment)

@app.callback(Output("stats-right","figure"),
    [Input("experiment-select-right","value")])
def update_stats_right(experiment):
    return update_statistics(experiment)




if __name__ == "__main__":
    app.run_server(debug=True)
