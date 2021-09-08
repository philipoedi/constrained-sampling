import numpy as np
import pandas as pd
import glob
import os

from numpy.linalg import norm


from app_utils import create_dataframe, results_folder


print(glob.glob(os.path.join(results_folder,"*/")))

cols_keep = ["x","y","z"]

def rename_cols(cols, pre):
    return {c:pre+"_"+c for c in cols}

def drop_and_rename(data,pre):
    new_data = data[cols_keep]    
    new_data.rename(columns=rename_cols(cols_keep,pre), inplace=True)
    return new_data


def create_summary_table():
    
    for simulation, reference in  ....:
        # variance

        # entropy

        # avg iterations

        # coverage

        # coverage / avg ref distance






folder = "202108171458_uniform_biased__" 

samples, seeds, pdes = create_dataframe(folder)
num_iterations = 
#samples = drop_and_rename(samples,"samples")

dists = np.zeros(len(samples))
for i in range(len(samples)):
    row_samp = samples.iloc[i,:3].values
    row_seed = seeds.iloc[i,:3].values
    dists[i] = norm(row_samp - row_seed)



