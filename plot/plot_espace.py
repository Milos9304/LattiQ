import glob
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import dtype

ax_rows=8;
ax_cols=8;

sorted_by = "natural"#"gray" #"natural"

def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)
files=glob.glob("*_espace")
#fig, axs = plt.subplots(ax_rows, ax_cols)#len(files))

def grayCode(n):
    res=[]
    for i in range(0, 1<<n):
        res.append(i^(i>>1))
    return res

grayCodeCalculated=[]

i=0
pxs=[]
pys=[]

ys_averaged=[]

counter=0
for file_name in files:
    
    print(float(counter)/len(files)*100)
    counter+=1
    
    ys = []
    xs = []
    f=open(file_name)
    for line in f.readlines():
        index, energy = line.split(" ")
        xs.append(int(index))
        ys.append(int(energy))
    
    if sorted_by == "gray":    
        if len(xs) > len(grayCodeCalculated):
            grayCodeCalculated = grayCode(math.floor(math.log2(len(xs))))
            s=grayCodeCalculated[0:len(xs)]
        old_xs=xs
        old_ys=ys
        xs=list(np.linspace(0,len(xs)-1,len(xs),dtype=int))

        ys=[ys[xs.index(i)] for i in s]
    
    elif sorted_by == "natural":
        s=argsort(xs)
        xs=[xs[i] for i in s]
        ys=[ys[i] for i in s]
        ys_averaged.append(ys)

    else:
        print("Invalid sorted_by value!")
        throw          
    
    sv_vect_min=min(ys)
   
    if xs == pxs and ys == pys:
        print("SAAME")
  
    
    """axs[int(i/ax_cols)][i%ax_cols].bar(xs,ys)#,color=color)
    axs[int(i/ax_cols)][i%ax_cols].get_xaxis().set_visible(False)
    opt_xs=[]
    for j in range(0,len(xs)):
        if sorted_by == "gray": 
            if ys[j] == sv_vect_min:
                #print(j)
                #k
                axs[int(i/ax_cols)][i%ax_cols].axvline(x=xs[j],color="red")   
                opt_xs.append(s.index(j)) 
        elif sorted_by == "natural":    
            if ys[j] == sv_vect_min:
                axs[int(i/ax_cols)][i%ax_cols].axvline(x=xs[j],color="red")   
                opt_xs.append(xs[j]) 
        else:
            print("Invalid sorted_by value!")
            throw  
    axs[int(i/ax_cols)][i%ax_cols].title.set_text("SV="+str(sv_vect_min)+" "+str(opt_xs))
    #axs[i].ylim(0,max(ys)+10)
    #plt.show()"""
    pxs=xs
    pys=ys
    i+=1
print("Let's see")
ys_averaged=np.array(ys_averaged)
np.save("ys_averaged", ys_averaged)
#plt.show()

#x = np.array([1, 2, 3, 4, 5])
#y = np.power(x, 2) # Effectively y = x**2
#e = np.array([1.5, 2.6, 3.7, 4.6, 5.5])

plt.errorbar(xs, np.mean(ys_averaged, axis=0), np.std(ys_averaged, axis=0), linestyle='None', fmt='-o')
#plt.bar(xs, np.mean(ys_averaged, axis=0))
plt.show()
