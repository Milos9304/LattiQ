import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f(x, a, b):
    return a*(x**b)

ms=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

#CM_QAOA
a=[1.00001,1.15984,2.28085,1.60549,2.1921,5.93133,9.88554,13.8006,27.4429,56.3166,90.6087,159.448,202.084,243.47,293.127,325.474,402.876,441.893,492.594,]

#QAOA
a2=[1.14711,1.18816,1.13386,1.68751,2.77557,2.99467,4.79571,7.40275,93.451,179.464,256.138,355.215,320.873,358.389,389.017,451.603,606.609,597.354,679.494,]

popt, pcov = curve_fit(f, ms, a, maxfev=30000)
print("popt: ", popt)
print("pcov: ", pcov)

popt2, pcov2 = curve_fit(f, ms, a2, maxfev=30000)
print("popt2: ", popt2)
print("pcov2: ", pcov2)


print("x cm qaoa line_cm line_qaoa")

for i in range(len(ms)):
  
    print(ms[i], a[i], a2[i], f(ms[i], popt[0], popt[1]), f(ms[i], popt2[0], popt2[1]))

# Now I try to fix the common constant

len_a=len(a)
len_a2=len(a2)

ms.extend(ms)

print(ms)

def f_com(x, k, b, c):
    """cost=0.
    for ii in x:
        if ii > 100:
            i = x - 100
        else:
            i = x
        cost += a*(i**b)+a*(i**c)
    return cost"""
    res = k*(x[:int(len(x)/2)]**b)
    return np.append(res, k*(x[int(len(x)/2):]**c))


a.extend(a2)

popt_com, pcov_com = curve_fit(f_com, ms, a, maxfev=30000)

print("popt_com: ", popt_com)
print("pcov_com: ", pcov_com)
print()
print("x cm qaoa line_cm line_qaoa")

for i in range(int(len(ms)/2)):

    print(ms[i], a[i], a2[i], popt_com[0]*(ms[i]**popt_com[1]), popt_com[0]*(ms[i]**popt_com[2]))
#

"""
plt.plot(ms, a, label="CM_QAOA approx factor")
plt.plot(ms, a2, label="QAOA approx factor")

plt.plot(ms, f(ms, *popt), label="CM_QAOA fit "+str(round(popt[0], 2))+"*x^"+str(round(popt[1], 2)))
plt.plot(ms, f(ms, *popt2), label="QAOA fit "+str(round(popt2[0], 2))+"*x^"+str(round(popt2[1], 2)))



plt.xticks(ms)
plt.title("approx_factors_plot.py")
plt.xlabel("lattice dimension")
plt.ylabel("approximation factor")
plt.legend()
plt.show()
"""
