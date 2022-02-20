# Constrained Sampling

by Philip Oedi  
philip.oedi@gmail.com  
340896  
Computational Engineering Science  
TU Berlin - Learning and Intelligent Systems  
Prof. Marc Toussaint  
2022/02/22  


## About

This repository holds all code and data for my master thesis on Constrained Sampling at TU Berlins research group Learning and Intelligent Systems under Prof. Marc Toussaint. Jung-Su Ha and Joaquim Ortiz de Haro have been supervisors to this project.  

This README will give an overview of the data collected in this repository and how the code can be run to conduct the experiments described in the thesis, as  well as their evaluation.

## Overview

The structure of the 

### Code

The code to generate the experiment data, evaluate and visualize is split into three modules

1. Experiment creation

All code necessary is written in c++. The `*.hpp` and `*.cpp` files are located in the src folder. The test folder holds various test and more importantly `run_experiment.cpp` and example_experiment.cpp. How to specify, compile and run these will be explained later.  

Only the unit line experiments have been implemented in `python 3.8.5`

2. Experiment evaluation

The evaluation of the experiment data for entropy, variance, AIPS (Average Iterations per Sample) and ADTS (Average Distance to Sample Set) is written in `python 3.8.5`.
Specific jupter-notebooks (`evaluation.ipnb`) can be found in the folder `experiments_3d_sphere` and `robotics_evaluation`.


3. Experiment visualization

A dash-plotly webapp has been programmed to visualize the experiment results. The webapp can display for two experiments, samples, seeds and the probability density surface (only sphere experiments).
It is in `consamp`. A seperate readme can be found in consamp that holds all instructions on how to run the webapp.  

The static figures that are displayed in the thesis have been generated using python and plotly. All data for these can be found in figures folder.


#### Dependencies

The experiment creation program is written in `c++17` and depends on `Eigen` and `nlopt`.  

For the dependecies need to run all python code for the visualization and evaluation can be taken from `requirements.txt`.  

Alternatively, using anaconda the environment can be set up the following using the terminal:  

1. Create a new environment

```
conda create -n constrained_sampling
```

2. Activate new environment

```
conda activate constrained_sampling
```

3. Install necessary packages

```
conda install --file requirements.txt
``` 


### Data

* consamp: Visualization webapp
* experiments_1d_line: unit line experiment implementation and results
* experiments_2d_line: experiment results
* experiments_3d_sphere: experiment results. `Experiment Evaluation.ipynb` to evaluate metrics on results. `*.csv` with the evaluated metrics.
* figures: `figures.ipynb` to create all figures for thesis. Other data is auxilliary to the figure creation. 
* results: Experiment folders that will be used in the webapp
* robotics_evaluation: `Robotics Evaluation.ipynb` to evaluate the metrics-
* src: c++ header files 
* tests: tests and c++ main `run_experiment.cpp` and `example_experiment.cpp` 

#### Experiment Naming

An experiment folder holds all data for a single experiment run

This is an example of an experiment folder:

```
202201112158_uniform_biased_grid-walk_biased_center_connected_single_0
```

It can be interpredeted as the following algorithm and scenario choices:

```
date_globalsampler_globaloptimizer_localsampler_localoptimizer_scenario_description_#samplingrun
```

If `reference` is contained in the naming, then the data in this folder corresponds to the rejection sampled referenec dataset. If a folder has naming similar like `202201112159_uniform_biased___center_connected_36`, there are no local optimzer and sampler, then this is data for a purely i.i.d-Biased-Sampler (not a composite sampler).


The files in a folder is either local or global and has any of the following suffixes:

* samples: space seperated coordinates of the feasible samples
* seeds: space seperated coordinates of the corresponding seeds
* probs: KDE estimate at sample location
* pdes: KDE estimated over a grid for visualization
* num_iterations: #iterations for solving a single biased optimization problem


### Running the experiments

The file `example_experiment.cpp` contains code for a single experiment run. The comments in the code explain how the scenariom, parameters and algorithms can be specified.

The file `run_experiments.cpp` contains the code to run all sphere experiments. `experiment_line_2d.cpp` runs the experiment on the line in 2d with i.i.d-Biased-Sampler.



Compile these as follows:

```
g++ -O3 -std=c++17 -I ../src ../src/utils.cpp ../src/filter.cpp example_experiment.cpp -o example_experiment -lnlopt
```

and run:

```
./example_experiment
```

The unit-line experiments can be run using python and ipython-notebook from the `experiments_1d_folder` folder.

```
conda activate constrained_sampling
jupyter-lab
```

Then select `cellular sampler.ipynb` or `line experiment.ipynb`.

### Experiment Evaluation

The evaluation has been done using python and ipython notebooks. Within each folder there is an evaluation `*.ipynb`.  
These can be run using jupyter-lab as above shown.






