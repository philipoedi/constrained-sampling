import context

from app_utils import *

import os
from os.path import dirname
import pdb

import glob


ROOT_FOLDER = os.path.abspath(os.path.join(dirname(__file__),"..","..",".."))
TEST_FOLDER = os.path.join(ROOT_FOLDER,"consamp","plotting","tests")  
TEST_FOLDER_GLOBAL = os.path.join(TEST_FOLDER,"data","global_data")

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

    meta = get_metadata("222_uniform_biased___global_samples.dat")
    assert meta["date"] == "222"
    assert meta["global_sampler"] == "uniform"
    assert meta["global_optimizer"] == "biased"
    assert meta["local_sampler"] == ""
    assert meta["local_optimizer"] == ""
    assert meta["level"] == "global"
    assert meta["type"] == "samples"


def test_read_file():
    filename = glob.glob(os.path.join(TEST_FOLDER_GLOBAL,"*samples.dat"))[0]
    data = read_file(filename)
    assert all(data.columns == ["x","y"])


def test_load_data():
    files = glob.glob(os.path.join(TEST_FOLDER_GLOBAL,"*samples.dat"))
    data = load_data(files)
    assert all(data.columns == ["x","y","global","local"])
    assert all(data.dtypes == [float, float, int, int])

def test_create_dataframe():
    samples, seeds, pdes = create_dataframe(TEST_FOLDER_GLOBAL)
    assert len(samples)/2 <= len(seeds)
    assert len(seeds) >= 100
    assert len(pdes) >= 100

    





