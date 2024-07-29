import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

#m: [mean, stdv]
data_p2_i3000 = {4: [4.20837,3.75774],5:[7.96436,7.81426],6:[14.2568,15.8808 ],7:[23.9784,29.1775 ], 8:[45.509,62.2681 ] ,9:[75.5909,111.304], 10: [134.586,224.513]}
data_p3_i3000 = {4:[4.35712,4.13457], 5:[8.12003,8.3758], 6:[14.016,16.0804],  7:[23.6221,29.0158],  8:[44.7842,59.5323], 9:[77.6903,109.271], 10:[133.81,205.992]}

m_range = range(4, 10+1)
plot_rnd_guess_notPnl=[]
plot_rnd_guess_pnl=[]
plot_means=[]

means = []

for m in m_range:
    dim=m
    nqs=m*2-2
    rnd_guess_notPnl=1./(2**dim)
    rnd_guess_pnl=1./(2**nqs)

    plot_rnd_guess_notPnl.append(rnd_guess_notPnl)
    plot_rnd_guess_pnl.append(rnd_guess_pnl)
    plot_means.append(data_p3_i3000[m][0]*rnd_guess_pnl)

    #means.append(data[m][0])

def f(n, alpha, c):
    return 1/2**(alpha*2*(n-1)+c)

xs=np.array(list(m_range))
print(plot_means)
ys=np.array(plot_means)
popt, pcov = curve_fit(f, xs, ys, maxfev=30000)
#zd=np.polyfit(list(map(lambda x: int(x), headers)), data["depth"], 1)
print("alpha, c: ", popt)
print(pcov)

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
