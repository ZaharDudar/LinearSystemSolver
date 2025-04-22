import scipy as sp
import numpy as np
import math
N, M = 16,16 #сетка, размер матрицы = N*M
# if(int(math.sqrt(N)) - math.sqrt(N) != 0): raise Exception("N is not a square")
# n = int(math.sqrt(N))
block_size = M

central_block = np.eye(block_size,block_size)*4
central_block[1:, :-1] -= np.eye(block_size-1,block_size-1)
central_block[:-1, 1:] -= np.eye(block_size-1,block_size-1)

sub_block = -np.eye(block_size,block_size)
zero_block = np.zeros_like(sub_block)

matrix_lines = []
for i in range(N):
    line = [zero_block]*N
    line[i] = central_block
    try:
        if(i-1>=0):
            line[i-1] = sub_block
    except:
        pass
    try:
        line[i+1] = sub_block
    except:
        pass
    matrix_lines.append(np.concatenate(line, axis=1))

print(np.concatenate(matrix_lines))
matrix = np.concatenate(matrix_lines)
w, v = np.linalg.eig(matrix)
print(np.max(w), np.min(w))
print(w)
with open(f"./{N*M}x{N*M}matrixDOK.txt", "w") as f:
    for i in range(N*M):
        for j in range(N*M):
            if matrix[i,j]!=0:
                f.write(f"{i},{j},{matrix[i,j]}\n")

# np.concatenate(matrix_lines).tofile(f"./{N}x{N}matrix.txt", sep=",")