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
for i in range(1,8,2):
    y = data[:,i]
    xt = data[:,i+1]
    x = data[:,0]
    for dot in range(1,len(y)):
        print(abs(y[dot]-y[dot-1]), names[i])
        if abs(y[dot]-y[dot-1]) <= 1e-9:
            y = data[:dot,i]
            xt = data[:dot,i+1]
            x = data[:dot,0]
            break

    axs[0].plot(x, y,"--", label=names[i])
    axs[1].plot(xt, y,"--", label=names[i])

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
