import numpy as np
import matplotlib.pyplot as plt
import sys

if sys.argv[1] == "9":
    app = "_p9"
else:
    app = ""    

if sys.argv[2] == "h":
    typ = "h" #histogram
    window=100
elif sys.argv[2] == "i":
    typ = "i" #integral
    window=1
else:
    raise("Error in input")

average_data={}
average_data2={}

Range = [""]#,"_p9"]

for app in Range:
    
    average_data[app]=[]
    average_data2[app]=[]

    #with open('histogram_cmqaoa'+app, 'r') as file:
    #    cmqaoa = file.readline().split()
        
    #with open('histogram_qaoa'+app, 'r') as file:
    #    qaoa = file.readline().split()

    with open('histogram_cmqaoa_p3_1inst', 'r') as file:
        cmqaoa = file.readline().split()

    with open('histogram_qaoa_p3_1inst', 'r') as file:
        qaoa = file.readline().split()


    cmqaoa = list(map(float, cmqaoa))
    qaoa = list(map(float, qaoa))

    rnd=2**(-14)

    data=cmqaoa
    window = 100
    ssum=0
    for ind in range(len(data) - window + 1):
        if typ == "h":
            average_data[app].append(np.mean(data[ind:ind+window]))
        else:
            ssum+=data[ind]
            average_data[app].append(ssum)

    data=qaoa
    ssum=0
    for ind in range(len(data) - window + 1):
        if typ == "h":
            average_data2[app].append(np.mean(data[ind:ind+window]))
        else:
            ssum+=data[ind]
            average_data2[app].append(ssum)

    #fig, ax = plt.subplots(nrows=2, ncols=1)

    print(sum(cmqaoa))
    print(sum(qaoa))

    print(len(cmqaoa))
    print(len(qaoa))

if typ == "h":

    average_data=average_data[app]
    average_data2=average_data2[app] 

    x = np.linspace(0, len(cmqaoa)-1, len(cmqaoa))
    plt.scatter(x,cmqaoa, color='blue', s=1, alpha=0.5)
    plt.scatter(x,qaoa, color='yellow', s=1) #alpha=0.4)
    plt.plot(average_data, color='black')
    plt.plot(average_data2, color='brown')

    plt.axhline(y=rnd, color='r', linestyle='-')

elif typ == "i":

    for a in Range:
        
        if a == "_p9":
            aa="9"
            alpha=1
        else:
            aa="3"
            alpha=0.5

        plt.plot(average_data[a], color="blue", alpha=alpha, label=aa+"CM")
        plt.plot(average_data2[a], color='red', alpha=alpha, label=aa+"QAOA")

plt.legend()
plt.show()
