import numpy as np
import matplotlib.pyplot as plt

with open('histogram_cmqaoa', 'r') as file:
    cmqaoa = file.readline().split()
    
with open('histogram_qaoa', 'r') as file:
    qaoa = file.readline().split()

cmqaoa = list(map(float, cmqaoa))
qaoa = list(map(float, qaoa))

rnd=2**(-14)

data=cmqaoa
window = 100
average_data = []
for ind in range(len(data) - window + 1):
    average_data.append(np.mean(data[ind:ind+window]))
    
data=qaoa
average_data2 = []
for ind in range(len(data) - window + 1):
    average_data2.append(np.mean(data[ind:ind+window]))

#fig, ax = plt.subplots(nrows=2, ncols=1)

print(sum(cmqaoa))
print(sum(qaoa))

print(len(cmqaoa))
print(len(qaoa))


x = np.linspace(0, len(cmqaoa)-1, len(cmqaoa))
plt.scatter(x,cmqaoa, color='blue', s=1, alpha=0.5)
plt.scatter(x,qaoa, color='yellow', s=1) #alpha=0.4)
plt.plot(average_data, color='black')
plt.plot(average_data2, color='brown')

plt.axhline(y=rnd, color='r', linestyle='-')

plt.show()
