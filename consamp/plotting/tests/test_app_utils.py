import context

from app_utils import *

import os
from os.path import dirname

import glob


ROOT_FOLDER = os.path.abspath(os.path.join(dirname(__file__),"..","..",".."))
TEST_FOLDER = os.path.join(ROOT_FOLDER,"results","test")  
data = create_dataframe(TEST_FOLDER)
def test_get_metadata():
    meta = get_metadata("222_uniform_biased_rrt_biased_global_samples.dat")
    assert meta["date"] == "222"
    assert meta["global_sampler"] == "uniform"
    assert meta["global_optimizer"] == "biased"
    assert meta["local_sampler"] == "rrt"
    assert meta["local_optimizer"] == "biased"
    assert meta["level"] == "global"
    assert meta["root"] == None
    assert meta["type"] == "samples"

    meta = get_metadata("222_uniform_biased_rrt_biased_global_1_samples.dat")
    assert meta["date"] == "222"
    assert meta["global_sampler"] == "uniform"
    assert meta["global_optimizer"] == "biased"
    assert meta["local_sampler"] == "rrt"
    assert meta["local_optimizer"] == "biased"
    assert meta["level"] == "global"
    assert meta["root"] == 1
    assert meta["type"] == "samples"

def test_read_file():
    filename = os.path.join(TEST_FOLDER, "202106291012_uniform_biased_RRT_biased_local_29_seeds.dat")
    data = read_file(filename)
    assert all(data.columns == ["x","y","z"])


def test_load_data():
    files = glob.glob(os.path.join(TEST_FOLDER,"*"))
    data = load_data(files[:2])
    assert all(data.columns == ["x","y","z","global","local"])
    assert all(data.dtypes == [float, float, float, int, int])

def test_create_dataframe():
    samples, seeds = create_dataframe(TEST_FOLDER)
    assert len(samples) <= len(seeds)
    assert len(seeds) >= 2000

    





