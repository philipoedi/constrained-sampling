import pdb
import os
import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

root_folder = os.path.abspath(os.path.join(os.path.dirname("__file__"), "..", ".."))
results_folder = os.path.join(root_folder, "results")


def experiment_name(path):
    """
    takes in path to .meta file and returns experiment name
    """
    return os.path.split(path)[1].split(".")[0]

def experiment_dict(name):
    """
    takes in experiment name and returns dictionary with paths to data
    """
    exp_dict = {
        "meta": os.path.join(results_folder, name+".meta"),
        "samples": os.path.join(results_folder, name+"_results.dat"),
            "seeds": os.path.join(results_folder, name+"_samples.dat"),
        "pdes": os.path.join(results_folder, name+"_pdes.dat")}
    return exp_dict

def get_experiments():
    files = glob.glob(os.path.join(results_folder,"*.meta"))
    #experiments = [{"label": experiment_name(f), "value": experiment_dict(experiment_name(f))} for f in files if "meta" in f]
    experiments = [{"label": experiment_name(f), "value": experiment_name(f)} for f in files if "meta" in f]
    return experiments

def create_dataframe(folder_path):
    files = glob.glob(os.path.join(results_folder, folder_path,"*"))
    sample_files = [f for f in files if "samples" in f]
    seed_files = [f for f in files if "seeds" in f]
    pdes_file = [f for f in files if "pdes" in f]
    samples = load_data(sample_files)
    seeds = load_data(seed_files)
    pdes = load_pdes(pdes_file[0])
    return samples, seeds, pdes

def load_pdes(filename):
    data = np.loadtxt(filename)
    df = pd.DataFrame(data)
    df.columns = ["x","y","z","pdes"]
    return df


def read_file(filename):
    """
    reads a .dat file containing coordinates in 3d space. Will 
    transform the data to a dataframe containint x,y,z columns
    
    Parameter
    ---------
    filename: str

    Return
    ---------
    data: pd.DataFrame
    """
    data = []
    with open(filename) as f:
        for r in f:
            row = r.split(" ")
            data.append([float(num) for num in row[:-1]])
    data = pd.DataFrame(data, columns=["x","y","z"])
    return data


def load_data(files):
    """
    consecutively reads in .dat files with 3d coordinates and creates dataframes 
    from these. Given the .dat file name a mapping of local points to their 
    global roots is added as columns 'global' and 'local' (local being the local id)

    Parameter
    ---------
    files: list

    Return:
    ---------
    data: pd.DataFrame
    """
    data = []
    for f in files:
        f_data = read_file(f)
        metadata = get_metadata(f)
        if metadata["level"] == "global":
            f_data["global"] = np.arange(len(f_data))
            f_data["local"] = -1 
        else:
            f_data["global"] = metadata["root"]
            f_data["local"] = np.arange(len(f_data))
        f_data.columns = ["x","y","z","global","local"]
        data.append(f_data)
    data = pd.concat(data,0)
    return data

def get_metadata(filename):
    """
    Takes a filename and transforms it to a dictionary holding
    all specifications for the contents of the file, such as "date", "global_sampler" etc.

    Parameter
    ---------
    filename: str

    return:
    metadata: dict
    """
    details = os.path.basename(filename).split("_")
    metadata = {"date": details[0],
                "global_sampler": details[1],
                "global_optimizer": details[2],
                "local_sampler": details[3],
                "local_optimizer": details[4],
                "level": details[5],
                "root": None if len(details) == 7 else int(details[-2]),
                "type": details[-1].split(".")[0]}
    return metadata

def get_plotdata(data, local):
    if local:
        return data[data["local"] >= 0]
    else:
        return data[data["local"] == -1]

def get_scatterplot(data, local):
    plotdata = get_plotdata(data, local)
    return go.Scatter3d(x=plotdata["x"], y=plotdata["y"], z=plotdata["z"],mode="markers",marker={"size":2})

def get_projections(samples, seeds, local):
    sample_data = get_plotdata(samples, local)
    seeds_data = get_plotdata(seeds, local)
    figs = []
    count = 0
    for index, row in seeds_data.iterrows():
        start = row[["x","y","z"]]
        end = sample_data.iloc[index][["x","y","z"]]
        line_data = pd.concat([start,end],axis=1).T 
        figs.append(go.Scatter3d(x=line_data["x"], y=line_data["y"], z=line_data["z"],mode="lines", line={"color":"#ffe476"}))
    return figs



def get_surfaceplot(data):
    #data = data.head(200)
    #data = data.sort_values("x")
    z = data[["z"]]
    z["ww"] = data["z"]
    #return go.Surface(x=data["x"],y=data["y"],z=z.values)#z=data[["z","pdes"]])
    #return go.Scatter3d(x=data["x"],y=data["y"],z=data["z"],mode="markers")#z=data[["z","pdes"]])
    return go.Mesh3d(x=data["x"],y=data["y"],z=data["z"], alphahull=0, intensity=data["pdes"])



def seeds_plot(path):
    pass


def samples_plot():
    pass


def pde_plot():
    pass

