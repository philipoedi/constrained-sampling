import numpy as np
import pdb

from scipy.interpolate import griddata
import plotly.graph_objects as go


def get_spherical(x):
    theta = np.arccos(x[2])
    phi = np.arctan2(x[1],x[0])
    return np.array([theta,phi])


north_pole = get_spherical(np.array([0,0,1]))
south_pole = get_spherical(np.array([0,0,-1])) 

def have_dist(x,y):
    # (theta, phi)
    # (lat, long)
    y_lat = y[0]
    x_lat = x[0]
    y_long = y[1]
    x_long = y[1]
    return 2*np.arcsin(np.sqrt(np.sin((y_lat - x_lat) / 2)**2 + np.cos(x_lat) * np.cos(y_lat) * np.sin((y_long - x_long)/2)**2))

def dist_from_spherical(x,y):
    return haversine_distances(x.reshape(1,2), y.reshape(1,2))[0][0]
    
def dist_north_pole(x):
    return have_dist(x, north_pole)

def dist_south_pole(x):
    return have_dist(x, south_pole)

def polar_to_cartesian(p):
    x = p[0]*np.cos(p[1])
    y = p[0]*np.sin(p[1])
    return np.array([x,y])

def sphere_embeddings(data):
    data_sp = data[:,:3]
    data_pdes = data[:,3]
    le_0 = data_sp[:,2] <= 0
    ge_0 = data_sp[:,2] >= 0
    lower = data_sp[le_0]
    upper = data_sp[ge_0]
    if len(lower) > 0:
        lower_sp = np.asarray(list(map(get_spherical, lower)))
        lower_sp[:,0] = np.asarray(list(map(dist_south_pole, lower_sp))) 
        lower_cart = np.asarray(list(map(polar_to_cartesian, lower_sp)))
        lower_vals = np.zeros(len(lower_cart)*3).reshape(len(lower_cart),3)
        lower_vals[:,:2] = lower_cart
        lower_vals[:,2] = data_pdes[le_0]
    else:
        lower_vals = None
    if len(upper) > 0:
        upper_sp = np.asarray(list(map(get_spherical, upper)))
        upper_sp[:,0] = np.asarray(list(map(dist_north_pole, upper_sp))) 
        upper_cart = np.asarray(list(map(polar_to_cartesian, upper_sp)))
        upper_vals = np.zeros(len(upper_cart)*3).reshape(len(upper_cart),3)
        upper_vals[:,:2] = upper_cart
        upper_vals[:,2] = data_pdes[ge_0]
    else: 
        upper_vals = None
    return upper_vals, lower_vals


def projection_plot(data, pdes):
    b = np.pi/2*1.1
    x,y = np.mgrid[-b:b:1000j, -b:b:1000j]
    z = griddata(data,pdes, (x,y),method="cubic")
    x = np.unique(x)
    y = np.unique(y)
    return go.Heatmap(x=x,y=y,z=z)
    #use griddata for interpolation
    #use imshow




