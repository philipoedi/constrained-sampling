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
global_selected_buttons = {"samples": True,
                    "seeds": True,
                    "surface": True,
                    "projections": True}
local_selected_buttons = global_selected_buttons.copy()
global_checklist = [{"label":k,"value":k} for k in global_selected_buttons.keys()]
local_checklist = [{"label":k,"value":k} for k in local_selected_buttons.keys()]

experiments_path = glob.glob(os.path.join(results_folder,"*/"))
experiments_name = [os.path.split(os.path.split(p)[0])[-1] for p in experiments_path] 
experiments_map = {experiments_name[i]: experiments_path[i] for i in range(len(experiments_name))}
experiments_dropdown_map = [{"label": e,"value": e} for e in experiments_name]


print(experiments_name)

app = dash.Dash(__name__)

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


@app.callback(Output("graph","figure"),
                [Input("experimentselect","value"),
                Input("global_checklist","value"),
                Input("local_checklist","value")])
def update_plot(experiment, global_checklist, local_checklist):
    experiment_name = os.path.join(results_folder,experiment)
    samples, seeds, pdes = create_dataframe(experiment_name)
    print("finished loading data")
    
    plots = []
    if "samples" in global_checklist:
        plots.append(get_scatterplot(samples, local=False))
    
    if "seeds" in global_checklist:
        plots.append(get_scatterplot(seeds, local=False))
    
    if "samples" in local_checklist:
        plots.append(get_scatterplot(samples, local=True))
    
    if "seeds" in local_checklist:
        #local_seeds_data = seeds[seeds["local"] != 0] 
        #local_seeds_fig = go.Scatter3d(x=local_seeds_data["x"],y=local_seeds_data["y"],z=local_seeds_data["z"])
        plots.append(get_scatterplot(seeds, local=True))
    if "projections" in global_checklist:
        plots.extend(get_projections(samples, seeds, local=False))

    fig = go.Figure(data=plots).update_layout(autosize=False,height=800,width=1600)
    if "surface" in global_checklist:
       # plots.append(get_surfaceplot(pdes)) 
        fig.add_trace(get_surfaceplot(pdes)) 
    

    #if "seeds" in global_checklist:
    #    seeds_data = np.loadtxt(data_map["seeds"])
    #    seeds_fig = go.Scatter(x=seeds_data[:,0],y=seeds_data[:,1],mode="markers")
    #    plots2.append(seeds_fig)
    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
