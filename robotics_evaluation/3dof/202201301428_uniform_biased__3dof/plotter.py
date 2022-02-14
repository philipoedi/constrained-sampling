import numpy as np
import matplotlib.pyplot as plt
import glob 


fig = plt.figure()
ax = fig.add_subplot(projection="3d")

files = glob.glob("*local*samples*")
for f in files:
	data = np.loadtxt(f)
	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
	ax.scatter(x,y,z, marker="o")
plt.show()
