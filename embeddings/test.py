import numpy as np
import matplotlib.pyplot as plt
from map import *


data = np.loadtxt("gridwalk.dat")
ref = np.loadtxt("reference.dat")
npole = np.array([0,0,1])
spole = np.array([0,0,-1])


ref_pdes = np.loadtxt("rename_pdes.dat")
gw_pdes = np.loadtxt("pdes_grid.dat")
upper_cart , lower_cart = sphere_embeddings(gw_pdes)
