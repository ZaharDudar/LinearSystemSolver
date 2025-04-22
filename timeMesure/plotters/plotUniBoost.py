import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

import matplotlib.pylab as plt
import numpy as np

names = []
data = []
with open(dir_path+"/UniversalBoost.csv", "r") as f:
    names = f.readline().strip().split(",")
    for line in f.readlines():
        data.append([float(i.strip()) for i in line.strip().split(",")])

data = np.array(data)
fig, axs = plt.subplots(1, 2)
print(names)
flag=True
for i in range(1,len(names)-1,2):
    if flag:
        axs[0].plot(data[np.where(data[:,i]>1e-7)][:,0], data[np.where(data[:,i]>1e-7)][:,i],"--", label=" ".join(names[i].split("_")[:-1]))
        axs[1].plot(data[np.where(data[:,i]>1e-7)][:,i+1], data[np.where(data[:,i]>1e-7)][:,i],"--", label=" ".join(names[i].split("_")[:-1]))
    else:
        axs[0].plot(data[np.where(data[:,i]>1e-7)][:,0], data[np.where(data[:,i]>1e-7)][:,i],"-.", label=" ".join(names[i].split("_")[:-1]))
        axs[1].plot(data[np.where(data[:,i]>1e-7)][:,i+1], data[np.where(data[:,i]>1e-7)][:,i],"-.", label=" ".join(names[i].split("_")[:-1]))
    flag=not flag
    
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
