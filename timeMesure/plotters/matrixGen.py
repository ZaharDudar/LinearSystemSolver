import scipy as sp
import numpy as np
rng = np.random.default_rng()
N=10
S = sp.sparse.random(N, N, density=0.65) * 5
S += sp.sparse.eye(N) * 20
S1 = S.toarray()
S1 = S1.astype(int)
print(np.triu(S1)+ np.triu(S1,k=1).T)

# file = open("matrix.txt", "w", encoding="utf-8")
(np.triu(S1)+ np.triu(S1,k=1).T).tofile("./matrix.txt",sep=",")
# file.close()