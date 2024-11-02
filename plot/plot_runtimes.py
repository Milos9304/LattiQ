import glob
import matplotlib.pyplot as plt
import matplotlib as mpl

iis={}

legends=[]
ps=[]

fig, ax = plt.subplots(2,3)

class Instance:
    duration = []
    nbQubits = 0

    def __init__(self, duration, nbQubits):
        self.duration = duration
        self.nbQubits = nbQubits

for file in glob.glob("../experiments/performance_experiment/test_interim_energies_*"):
    f = open(file)
    
    idn=file.split("_q")[-1]
    legends.append(idn)
    num_qs=int(idn.split("_")[0])
    p=int(idn.split("_p")[-1])
    s=f.read().split("\n")[0].split("=")
    if len(s) == 1:
        continue
    
    ps.append(p-1)
    
    if p not in iis:
        iis[p]=[]

    iis[p].append(Instance(int(s[1])/60.0, num_qs))

cmap = mpl.colormaps['viridis']
max_qs=0
for p in ps:
    max_qs=max(max_qs, max(list(map(lambda x:x.nbQubits,iis[p+1]))))
for p in ps:
    for i in range(len(iis[p+1])):
        if p > 2:
            row = 1
        else:
            row = 0
        #print(iis[p+1][i].nbQubits," ", iis[p+1][i].duration)
        ax[row][p % 3].plot(iis[p+1][i].nbQubits,iis[p+1][i].duration,'o',c=cmap((iis[p+1][i].nbQubits-1)/(max_qs)))
        ax[row][p % 3].set_title("p="+str(p+1))
        ax[row][p % 3].annotate(str(iis[p+1][i].nbQubits), (iis[p+1][i].nbQubits, iis[p+1][i].duration))
        #ax[row][p % 3].set_xlim([0,50*(p+1)])
    #plt.legend(legends)
plt.show()
