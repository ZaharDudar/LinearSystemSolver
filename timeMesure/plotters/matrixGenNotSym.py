import scipy as sp
import numpy as np
import math
N = 256
mul = 10

matrix = sp.sparse.rand(N,N,density=0.5).toarray()
matrix += np.eye(N) * np.max(np.sum(matrix,axis=1))
matrix *= mul
matrix = matrix.astype(int)
w, v = np.linalg.eig(matrix)
print(np.max(w), np.min(w))
print(matrix)
with open(f"./{N}x{N}matrixDOKNotSym.txt", "w") as f:
    for i in range(N):
        for j in range(N):
            if matrix[i,j]!=0:
                f.write(f"{i},{j},{matrix[i,j]}\n")

# np.concatenate(matrix_lines).tofile(f"./{N}x{N}matrix.txt", sep=",")