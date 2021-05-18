import os
import glob

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

def seeds_plot(path):
    pass

def samples_plot():
    pass

def get_metadata():
    pass

def pde_plot():
    pass

