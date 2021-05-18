import os
import glob


import numpy as np

import dash
import dash_html_components as html
import dash_core_components as dcc
import plotly.express as px
from dash.dependencies import Input, Output
import plotly.graph_objects as go

from app_utils import *

root_folder = os.path.abspath(os.path.join(os.path.dirname("__file__"),"..",".."))
results_folder = os.path.join(root_folder, "results")
selected_buttons = {"samples": True,
                    "seeds": True,
                    "pdes": True}
checklist_map = [{"label":k,"value":k} for k in selected_buttons.keys()]

app = dash.Dash(__name__)

app.layout = html.Div(children=[
      html.Div(className="header", children=[
        html.P("experiment"),
        html.Div(className="dropdown", children=[
            dcc.Dropdown(id="experimentselect",
                         options=get_experiments(),
                         multi=False,
                         value="",
                         className="experimentselect"),
            dcc.Checklist(id="checklist",
                        options=checklist_map,
                        value=["samples", "seeds", "pdes"])])]),
      html.Div(className="graphs", children=[
        html.Div(className="graph1div", children=[dcc.Graph(id="graph",
                  #animate=True,
                  #config={"displayModeBar":False},
                  #figure=px.line([2,3,4,5,1,2,3],
                #  )]),
                )]),
      #html.Div(className="graph2", children=[
        html.Div(className="graph2div2", children=[dcc.Graph(id="graph2")])]),
      html.Div(id="metadata")])


@app.callback(Output("graph","figure"),
                Output("graph2","figure"),
                Output("metadata","children"),
                [Input("experimentselect","value"),
                Input("checklist","value")])
def update_plot(experiment, checklist):
    data_map = experiment_dict(experiment)
    data = np.loadtxt(data_map["samples"])
    plots = []
    plots2 = []
    if "samples" in checklist:
        samples_data = np.loadtxt(data_map["samples"])
        samples_fig = go.Scatter(x=samples_data[:,0],y=samples_data[:,1],mode="markers")
        plots2.append(samples_fig)
    if "seeds" in checklist:
        seeds_data = np.loadtxt(data_map["seeds"])
        seeds_fig = go.Scatter(x=seeds_data[:,0],y=seeds_data[:,1],mode="markers")
        plots2.append(seeds_fig)
    if "pdes" in checklist:
        pdes_data = np.loadtxt(data_map["pdes"])
        x = np.unique(pdes_data[:,0])
        y = np.unique(pdes_data[:,1])
        z = pdes_data[:,2].reshape(len(x),len(y))
        pdes_fig = go.Surface(x=x,y=y,z=z)
        plots.append(pdes_fig)
    fig = go.Figure(data=plots).update_layout(autosize=False,height=800,width=1600)
    fig2 = go.Figure(data=plots2).update_layout(autosize=False,height=800,width=800)
    meta= [] 
    with open(data_map["meta"],"r") as f:
        for line in f:
            meta.append(html.P(line))
    return fig, fig2, meta


if __name__ == "__main__":
    app.run_server(debug=True)
