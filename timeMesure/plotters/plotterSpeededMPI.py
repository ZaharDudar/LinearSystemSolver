import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

import matplotlib.pylab as plt
import numpy as np

names = []
data = []
with open(dir_path+"/TimeSpeededMPI.csv", "r") as f:

    names = f.readline().strip().split(",")
    for line in f.readlines():
        data.append([float(i.strip()) for i in line.strip().split(",")])

data = np.array(data)
fig, axs = plt.subplots(1, 2)
print(names)
for i in range(1,6,2):
    axs[0].plot(data[:,0], data[:,i],"--", label=names[i])
    axs[1].plot(data[:,i+1], data[:,i],"--", label=names[i])
axs[0].set_yscale("log")
axs[1].set_yscale("log")
axs[0].grid()
axs[1].grid()
axs[0].legend()
axs[1].legend()
axs[0].set_ylabel("Error")
axs[1].set_ylabel("Error")
axs[0].set_xlabel("N Iter")
axs[1].set_xlabel("Time (μs)")
plt.show()
