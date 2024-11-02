import glob
import matplotlib.pyplot as plt
import matplotlib as mpl

iis={}
mins={}

legends=[]
ps=[]

fig, ax = plt.subplots(2,3)

class Instance:
    energies = []
    nbQubits = 0

    def __init__(self, energies, nbQubits):
        self.energies = energies
        self.nbQubits = nbQubits

for file in glob.glob("../experiments/performance_experiment/test_interim_energies_*"):
    f = open(file)
    
    idn=file.split("_q")[-1]
    legends.append(idn)
    num_qs=int(idn.split("_")[0])
    p=int(idn.split("_p")[-1])
    s=f.read().split("\n")
    if len(s) == 1:
        continue
    ps.append(p-1)
    
    if p not in mins:
        iis[p]=[]
        mins[p]=[]
    
    iis[p].append(Instance(list(map(float, s[1].split(" ")[:-1])), num_qs))
    mins[p].append(iis[p][-1].energies[-1])

mmins={}
cmap = mpl.colormaps['viridis']
max_qs=0
for p in ps:
    mmins[p+1]=min(mins[p+1])
    max_qs=max(max_qs, max(list(map(lambda x:x.nbQubits,iis[p+1]))))
for p in ps:
    for i in range(len(iis[p+1])):
        delta=mmins[p+1]-(iis[p+1][i].energies[-1])
        iis[p+1][i].energies = list( map(lambda x: x+delta, iis[p+1][i].energies))
        
        if p > 2:
            row = 1
        else:
            row = 0
        ax[row][p % 3].plot(iis[p+1][i].energies[:],c=cmap(1-(iis[p+1][i].nbQubits-1)/(max_qs-1)))
        ax[row][p % 3].set_title("p="+str(p+1))
        ax[row][p % 3].set_xlim([0,50*(p+1)+50])
    #plt.legend(legends)
plt.show()
