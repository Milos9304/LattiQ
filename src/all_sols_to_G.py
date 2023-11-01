import numpy as np

G=np.array([[49,0,1],[0,1,0],[0,0,1]])
for x1 in [-1,0,1,2]:
    for x2 in [-1,0,1,2]:
        for x3 in [-1,0,1,2]:
            x=np.array([x1,x2,x3])
            r=np.matmul(np.matmul(x.transpose(), G), x)
            if r == 1:
                print(x)