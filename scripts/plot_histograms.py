import numpy as np
import matplotlib.pyplot as plt
import sys

#with open('histogram_cmqaoa'+app, 'r') as file:
#    cmqaoa = file.readline().split()
    
#with open('histogram_qaoa'+app, 'r') as file:
#    qaoa = file.readline().split()

count=50 #number of graphs top show
rows=10
cols=5
assert(rows*cols == count)

fig, ax = plt.subplots(nrows=rows,ncols=cols)

jj=0

#for typ in ["cmqaoa","qaoa"]:
    
for i in range(count):
    
    average_data=[]
    average_data2=[]


    with open('histograms/hist_cmqaoa_'+str(i), 'r') as file:
        cmqaoa = file.readline().split()

    with open('histograms/hist_qaoa_'+str(i), 'r') as file:
        qaoa = file.readline().split()

    cmqaoa = list(map(float, cmqaoa))
    qaoa = list(map(float, qaoa))

    rnd=2**(-14)

    data=cmqaoa
    window = 100
    ssum=0
    for ind in range(len(data)):
        #if typ == "h":
        #    average_data.append(np.mean(data[ind]))
        #else:
        ssum+=data[ind]
        average_data.append(ssum)

    data=qaoa
    ssum=0
    for ind in range(len(data)):
        #if typ == "h":
        #    average_data2[app].append(np.mean(data[ind:ind+window]))
        #else:
        ssum+=data[ind]
        average_data2.append(ssum)

    #fig, ax = plt.subplots(nrows=2, ncols=1)

    """print(sum(cmqaoa))
    print(sum(qaoa))

    print(len(cmqaoa))
    print(len(qaoa))
    """
    """if typ == "h":

        average_data=average_data[app]
        average_data2=average_data2[app] 

        x = np.linspace(0, len(cmqaoa)-1, len(cmqaoa))
        plt.scatter(x,cmqaoa, color='blue', s=1, alpha=0.5)
        plt.scatter(x,qaoa, color='yellow', s=1) #alpha=0.4)
        plt.plot(average_data, color='black')
        plt.plot(average_data2, color='brown')

        plt.axhline(y=rnd, color='r', linestyle='-')

    elif typ == "i":
    """
    print(jj,int(jj/cols), jj-cols*int(jj/cols))
    ax[int(jj/cols)][jj-cols*int(jj/cols)].plot(average_data, color="blue", label="CM")
    ax[int(jj/cols)][jj-cols*int(jj/cols)].plot(average_data2, color='red', label="QAOA")
    ax[int(jj/cols)][jj-cols*int(jj/cols)].plot([0,2**14],[0,1], color='green', linewidth=1, label="Random guess")


    jj+=1

plt.legend()
plt.show()
