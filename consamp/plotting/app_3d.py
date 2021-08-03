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
from plotly.subplots import make_subplots

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

app.layout = html.Div([
    dcc.Store(id="samples"),
    dcc.Store(id="samples2"),
    dcc.Store(id="seeds"),
    dcc.Store(id="seeds2"),
    dcc.Store(id="pdes"),
    dcc.Store(id="pdes2"),
    dbc.Container(id="main-container", children=[
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
    dbc.Row(id="local-number", children=[
        dbc.Col([
            dcc.Dropdown(id="local-left", multi=True,
                options=[],
                value="")], width={"size":6,"order":1}),
        dbc.Col([
            dcc.Dropdown(id="local-right", multi=True,
                options=[],
                value="")], width={"size":6,"order":2}),
        ]),
    dbc.Row(id="show-options", children=[
        dbc.Col([
            dcc.Checklist(id="global_checklist",
                options=global_checklist,
                #value=[k for k in global_selected_buttons.keys()]),
                value=[]),
            dcc.Checklist(id="local_checklist",
                options=local_checklist,
                value=[])], width={"size":6}),
         ]),
    dbc.Row(id="graph-left-row", children=[
        dbc.Col([
            dcc.Graph(id="graph-main", figure={})
            ], width={"size":12})
        ]),
    dbc.Row(id="hist-row", children=[
        dbc.Col([
            dcc.Graph(id="hist", figure={})
            ], width={"size":12}),
        ]),
     dbc.Row(id="statistics-row", children=[
        dbc.Col([
            dcc.Graph(id="stats-left")
            ], width={"size":6,"order":1}),
        dbc.Col([
            dcc.Graph(id="stats-right")
            ], width={"size":6,"order":2})
        ])
    ],fluid=True)])



def create_dropdown_options(data):
    opt = [{"label": str(l),"value": int(l)} for l in data["global"].unique()]
    opt.append({"label":"all","value": -1}) 
    return sorted(opt, key=lambda x: x["value"])


@app.callback(Output("local-left","options"),
            Output("local-right","options"),
            [Input("samples","data"),
            Input("samples2","data")])
def update_locals(samplesj, samples2j):
    samples = pd.read_json(samplesj, orient="split")
    samples2 = pd.read_json(samples2j, orient="split")
    opt1 = create_dropdown_options(samples)
    opt2 = create_dropdown_options(samples2)
    return opt1, opt2 

@app.callback(Output("samples","data"),
            Output("samples2","data"),
            Output("seeds","data"),
            Output("seeds2","data"),
            Output("pdes","data"),
            Output("pdes2","data"),
            [Input("experiment-select-left","value"),
            Input("experiment-select-right","value")])
def store_data(experiment,experiment2):
    experiment_name = os.path.join(results_folder,experiment)
    experiment_name2 = os.path.join(results_folder,experiment2)
    samples, seeds, pdes = create_dataframe(experiment_name)
    samples2, seeds2, pdes2 = create_dataframe(experiment_name2)
    samples = samples.to_json(orient="split") 
    samples2 = samples2.to_json(orient="split") 
    seeds = seeds.to_json(orient="split") 
    seeds2 = seeds2.to_json(orient="split") 
    pdes = pdes.to_json(orient="split") 
    pdes2 = pdes2.to_json(orient="split") 
    return samples, samples2, seeds, seeds2, pdes, pdes2


@app.callback(Output("graph-main","figure"),
                [Input("global_checklist","value"),
                Input("local_checklist","value"),
                Input("samples","data"),
                Input("samples2","data"),
                Input("seeds","data"),
                Input("seeds2","data"),
                Input("pdes","data"),
                Input("pdes2","data"),
                Input("local-left","value"),
                Input("local-right","value"),
                ])
def update_plot(global_checklist, local_checklist, samplesj, samples2j, seedsj, seeds2j, pdesj, pdes2j, local_left, local_right):
    samples = pd.read_json(samplesj, orient="split")
    samples2 = pd.read_json(samples2j, orient="split")
    seeds = pd.read_json(seedsj, orient="split")
    seeds2 = pd.read_json(seeds2j, orient="split")
    pdes = pd.read_json(pdesj, orient="split")
    pdes2 = pd.read_json(pdes2j, orient="split")

    plots = []
    plots2 = []
    max_col = max([pdes["pdes"].max(),pdes2["pdes"].max()])
    if -1 not in local_left:
        samples_plot = samples[np.isin(samples["global"].values, np.array(local_left))]
        seeds_plot = seeds[np.isin(seeds["global"].values, np.array(local_left))]
    else:
        samples_plot = samples
        seeds_plot = seeds

    if -1 not in local_right:
        samples2_plot = samples2[np.isin(samples2["global"].values,np.array(local_right))]
        seeds2_plot = seeds2[np.isin(seeds2["global"].values, np.array(local_right))]
    else:
        samples2_plot = samples2
        seeds2_plot = seeds2

    dim_l = pdes.shape[1] - 1
    dim_r = pdes2.shape[1] - 1
    if dim_l == 3:
        if "samples" in global_checklist:
            plots.append(get_scatterplot(samples_plot, local=False))
        
        if "seeds" in global_checklist:
            plots.append(get_scatterplot(seeds_plot, local=False))
        
        if "samples" in local_checklist:
            plots.append(get_scatterplot(samples_plot, local=True))
        
        if "seeds" in local_checklist:
            plots.append(get_scatterplot(seeds_plot, local=True))

        if "projections" in global_checklist:
            plots.extend(get_projections(samples_plot, seeds_plot, local=False))

        if "surface" in global_checklist:
            plots.append(get_surfaceplot(pdes,max_col)) 
    
    if dim_l  == 2:
        if "samples" in global_checklist:
            plots.append(get_scatterplot2(samples_plot, local=False))
        
        if "seeds" in global_checklist:
            plots.append(get_scatterplot2(seeds_plot, local=False))
        
        if "samples" in local_checklist:
            plots.append(get_scatterplot2(samples_plot, local=True))
        
        if "seeds" in local_checklist:
            plots.append(get_scatterplot2(seeds_plot, local=True))

        if "projections" in global_checklist:
            plots.extend(get_projections2(samples_plot, seeds_plot, local=False))

        if "surface" in global_checklist:
            plots.append(get_surfaceplot2(pdes,max_col)) 
  
    if dim_r == 3:
        if "samples" in global_checklist:
            plots2.append(get_scatterplot(samples2_plot, local=False))
        
        if "seeds" in global_checklist:
            plots2.append(get_scatterplot(seeds2_plot, local=False))
        
        if "samples" in local_checklist:
            plots2.append(get_scatterplot(samples2_plot, local=True))
        
        if "seeds" in local_checklist:
            plots2.append(get_scatterplot(seeds2_plot, local=True))

        if "projections" in global_checklist:
            plots2.extend(get_projections(samples2_plot, seeds2_plot, local=False))

        if "surface" in global_checklist:
            plots2.append(get_surfaceplot(pdes2,max_col)) 
    
    if dim_r == 2:
        if "samples" in global_checklist:
            plots2.append(get_scatterplot2(samples2_plot, local=False))
        
        if "seeds" in global_checklist:
            plots2.append(get_scatterplot2(seeds2_plot, local=False))
        
        if "samples" in local_checklist:
            plots2.append(get_scatterplot2(samples2_plot, local=True))
        
        if "seeds" in local_checklist:
            plots2.append(get_scatterplot2(seeds2_plot, local=True))

        if "projections" in global_checklist:
            plots2.extend(get_projections2(samples2_plot, seeds2_plot, local=False))

        if "surface" in global_checklist:
            plots2.append(get_surfaceplot2(pdes2,max_col)) 
  


    fig = make_subplots(rows=1,cols=2,specs=[[{"type":"surface"},{"type":"surface"}]])
    for c , plot_data in enumerate([plots, plots2]):
        for p in plot_data:
            fig.add_trace(p, row=1, col=c+1)    
    fig.update_layout(autosize=False,height=800,width=1600)
    return fig



def load_prob_data(experiment):
    experiment_name = os.path.join(results_folder,experiment,"*probs.dat")
    files = glob.glob(experiment_name)
    data = np.loadtxt(files[0])
    return data 

def get_histogram_plot(experiment):
    data = load_prob_data(experiment)
    plot = get_histogram(data)
    return plot 

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

@app.callback(Output("hist","figure"),
     [Input("experiment-select-left","value"),
     Input("experiment-select-right","value")])
def update_histogram(experiment,experiment2):
    p = get_histogram_plot(experiment)
    p2 = get_histogram_plot(experiment2)
    fig = go.Figure()
    fig.add_trace(p)
    fig.add_trace(p2)
    fig.update_layout(barmode="overlay")
    fig.update_traces(opacity=0.75)
    return fig 


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
         [np.round(val,5) for val in [variance,stdev,entropy,num_samples,avg_its, global_samples, local_seeds, filtered_seeds]]]})
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
