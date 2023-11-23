import glob
import numpy as np
from numpy import linalg as lg

import matplotlib.pyplot as plt

from qiskit.quantum_info import SparsePauliOp

result=[]

files=glob.glob("*.ising")
for file_name in files:
    #print(file_name)
    f=open(file_name)
    fg=open(file_name[0:-5]+"gram")
    G=fg.readlines()
    G=np.array(list(map(lambda y:list(map(int,y.split())),map(lambda x: x[:-1] if x[-1] == '\n' else x,G))))
    
    hamil_list=[]
    for line in f.readlines():
        c, s = line.split(" ")
        if s[-1] == '\n':
            s=s[0:-1][::-1]
        hamil_list.append((s,c))

hamiltonian = SparsePauliOp.from_list(hamil_list)
evals, evects = lg.eig(hamiltonian)
minimum_indices = np.where(evals == min(evals))

files=glob.glob("*.many_params")
for file_name in files:
    #print(file_name)
    f=open(file_name)
    G=f.readlines()
    for line in G:
        row=list((map(float,line.split())))
        result.append(row)

bars=plt.bar(range(len(result[0])), np.mean(result, axis=0), yerr=np.std(result, axis=0), ecolor='black', capsize=10)#, fmt='-o')
for sol in list(minimum_indices[0]):
    bars[sol].set_color('red')
for i in range(len(bars)):
    yval = bars[i].get_height()
    plt.text(bars[i].get_x(), yval + .005, int(np.real(evals[i])),fontsize='x-large')
plt.title(file_name)
plt.show()