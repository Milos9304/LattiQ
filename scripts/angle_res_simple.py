import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

#m: [mean, stdv]
#data=[ 0.0684507       ,    0.0324452   ,        0.0162824    ,      0.00841373    ,      0.00680977    ,      0.00600239     ,     0.00476806] #0.97
#xs=np.array([4,5,6,7,8,9,10])

#data=[0.106731   ,        0.0515169     ,      0.0207424     ,      0.0126368   ,         0.005107   ,       0.00519504     ,      0.0044405] #0.83
#xs=np.array([4,5,6,7,8,9,10])

#data=[0.102938    ,       0.0451122   ,        0.0091706    ,      0.00585934   ,       0.00539514    ,      0.00442476] #0.72
#xs=np.array([5,6,7,8,9,10])

#data=[0.0859569    ,       0.0323786    ,      0.00566087 ,         0.00542296        ,  0.00443021] #0.661
#xs=np.array([6,7,8,9,10])

#data=[0.0508255  ,         0.0153818    ,      0.00288475 ,       0.00273762] #0.678
#xs=np.array([7,8,9,10])

#data=[0.0364387   ,       0.00482137    ,      0.00144237] #0.667
#xs=np.array([8,9,10])

#data=[0.0213782      ,    0.00144802] #0.667
#xs=np.array([9,10])

#data=[0.0122817] #0.63
#xs=np.array([10])


data=[0.108526   ,        0.0586565  ,        0.00488591      ,    0.00575196    ,      0.00542205    ,      0.00443147   ,       0.00323217      ,    0.00218827    ,      0.00137737     ,    0.000826238        , 0.000485569     ,    0.000277501    ,     0.000156021     ,    8.67612e-05  ,       4.78607e-05    ,     2.61246e-05]
data=[0.120244,           0.0608056        ,  0.00739476       ,   0.00375513 ,         0.00383703       ,   0.00289825     ,     0.00193899     ,     0.00121443     ,    0.000729683      ,   0.000426075    ,     0.000243654   ,      0.000137135      ,   7.62219e-05   ,       4.1938e-05     ,    2.28826e-05   ,      1.23981e-05]
xs=np.array([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]) #0.7

#data=[2.61246e-05]
#xs=np.array([20]) #0.7

def f(n, alpha):
    return 1/2**(alpha*n)

ys=np.array(data)
popt, pcov = curve_fit(f, xs, ys, maxfev=30000)
#zd=np.polyfit(list(map(lambda x: int(x), headers)), data["depth"], 1)
print("alpha: ", popt)
print(pcov)

#a=popt[0]

#xs=np.array([3,4,5,6,7])
#ys=np.array([0.136877, 0.0267844, 0.036026, 0.00151015, 0.00646481])

#popt, pcov = curve_fit(f, xs, ys, maxfev=30000)
#zd=np.polyfit(list(map(lambda x: int(x), headers)), data["depth"], 1)
#print("alpha: ", popt)
#print(pcov)

#a=popt[0]

print(xs, xs*xs, np.sum(xs*xs))
print(np.log(ys))
b=-np.sum(np.log(ys)*xs)/np.sum(xs*xs)
print("b=",b)
print("nom=",-np.sum(np.log(ys)*xs))
print("den=",np.sum(xs*xs))


print(-np.sum(np.log(ys)*xs)/np.sum(xs*xs))

#print(a,b)
print(np.sum(np.square(ys-np.power(np.e,-a*xs))))
print(np.sum(np.square(ys-np.power(np.e,-b*xs))))
print()
print(np.sum(np.square(ys-np.power(np.e,-(b-0.01)*xs))))
print(np.sum(np.square(np.log(ys)-(b-0.01)*xs)))
print()
print(np.square(ys-np.power(np.e,-b*xs)))

"""
fig, ax = plt.subplots()
fig.canvas.draw()

ax.plot(m_range, plot_rnd_guess_notPnl)
ax.plot(m_range, plot_rnd_guess_pnl)
ax.plot(m_range, plot_means) 
ax.legend(["rnd_guess_notPnl","rnd_guess_pnl","qaoa_mean_guess"])

labels = [item.get_text() for item in ax.get_xticklabels()]
labels = list(map(lambda x: 'dim '+x+'='+str(2*int(x)-2)+' qubits', labels))

ax.set_xticklabels(labels)

plt.show()
"""
