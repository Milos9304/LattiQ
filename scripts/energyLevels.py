def mean(data):
  n = len(data)
  mean = sum(data) / n
  return mean

def variance(data):
   
  n = len(data)
   
  mean = sum(data) / n
   
  deviations = [(x - mean) ** 2 for x in data]
    
  variance = sum(deviations) / n
  return variance

def stdev(data):
  import math
  var = variance(data)
  std_dev = math.sqrt(var)
  return std_dev

import matplotlib.pyplot as plt
import numpy as np

f=open("../build/a4_6", "r")
s=f.read()
l=list(map(float, s.split()))

print(mean(l))
print(stdev(l))

plt.hist(l, bins=100)#, bins=np.arange(l.min(), l.max()+1))
plt.show()
