import scipy as sp
import numpy as np
rng = np.random.default_rng()
N=100
S = sp.sparse.random(N, N, density=0.5) * 5
# S = sp.sparse.eye(N) * 200
S1 = S.toarray()
S1 = S1.astype(int)
S2 = np.triu(S1)+ np.triu(S1,k=1).T
for i in range(N):
    S2[i,i]+=np.sum(S2[i,:])+1
print(S2)
# file = open("matrix.txt", "w", encoding="utf-8")
(S2).tofile("./matrix.txt",sep=",")
# file.close()