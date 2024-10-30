import matplotlib.pyplot as plt
import glob
from statistics import median, mean

cm_vals={}
qaoa_vals={}

cm_average={}
qaoa_average={}

cm_median={}
qaoa_median={}

for file in glob.glob("appr_cm_*"):
    f=open(file)
    dim=int(file.split("_")[-1])
    vals=list(map(float, f.readlines()[0].split()))
    cm_vals[dim]=vals
    cm_average[dim]=mean(vals)
    cm_median[dim]=median(vals)

for file in glob.glob("appr_qaoa_*"):
    f=open(file)
    dim=int(file.split("_")[-1])
    vals=list(map(float, f.readlines()[0].split()))
    qaoa_vals[dim]=vals
    qaoa_average[dim]=mean(vals)
    qaoa_median[dim]=median(vals)

for i in range(4,22+1):
    for val in cm_vals[i]:
        plt.scatter(i-0.1, val, color='blue', s=1)
    for val in qaoa_vals[i]:
        plt.scatter(i+0.1, val, color='red', s=1)


dims=list(range(4,22+1))
cm_avg=list(map(lambda x: x[1], sorted(cm_average.items())))
qaoa_avg=list(map(lambda x: x[1], sorted(qaoa_average.items())))

cm_med=list(map(lambda x: x[1], sorted(cm_median.items())))
qaoa_med=list(map(lambda x: x[1], sorted(qaoa_median.items())))


plt.plot(dims,cm_avg,color='blue')
plt.plot(dims,qaoa_avg,color='red')

plt.plot(dims,cm_med,color='blue',linestyle='dashed')
plt.plot(dims,qaoa_med,color='red',linestyle='dashed')



plt.show()

