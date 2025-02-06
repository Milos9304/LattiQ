import matplotlib.pyplot as plt
import glob
from statistics import median, mean
import numpy as np

# Use PGF backend for LaTeX-style rendering
plt.rcParams.update({
    "pgf.texsystem": "pdflatex",  # Use LaTeX for rendering
    "text.usetex": True,  # Enable LaTeX
    "font.family": "serif",  # Use LaTeX font
    "pgf.rcfonts": False,  # Prevent Matplotlib from overriding fonts
    "axes.grid": True,  # Add a pgfplots-style grid
    "figure.figsize": (6, 4.5),  # ✅ Half of typical LaTeX text width (~7 inches)
    "figure.autolayout": True  # ✅ Adjust layout to prevent cutoff
    #"pgf.preamble": r"\renewcommand{\mathdefault}[1]{#1}"  # Fix for mathdefault error
    })

cm_vals={}
qaoa_vals={}

cm_average={}
qaoa_average={}

cm_median={}
qaoa_median={}

def ff(a):
    return np.quantile(a, 0.25)

for file in glob.glob("appr_cm_*"):
    f=open(file)
    dim=int(file.split("_")[-1])
    vals=list(map(float, f.readlines()[0].split()))
    cm_vals[dim]=vals[:100]
    cm_average[dim]=mean(vals)
    cm_median[dim]=ff(vals)

for file in glob.glob("appr_qaoa_*"):
    f=open(file)
    dim=int(file.split("_")[-1])
    vals=list(map(float, f.readlines()[0].split()))
    qaoa_vals[dim]=vals[:100]
    qaoa_average[dim]=mean(vals)
    qaoa_median[dim]=ff(vals)

"""
for i in range(4,22+1):
    for val in cm_vals[i]:
        plt.scatter(i-0.1, val, color='blue', s=1)
    for val in qaoa_vals[i]:
        plt.scatter(i+0.1, val, color='red', s=1)
"""

for i in range(4,22+1):
    plt.hlines(cm_vals[i], i-0.4, i+0.3, color="blue", linewidth=1, alpha=0.3)
    plt.hlines(qaoa_vals[i], i-0.3, i+0.4, color="red", linewidth=1, alpha=0.3)

dims=list(range(4,22+1))
cm_avg=list(map(lambda x: x[1], sorted(cm_average.items())))
qaoa_avg=list(map(lambda x: x[1], sorted(qaoa_average.items())))

cm_med=list(map(lambda x: x[1], sorted(cm_median.items())))
qaoa_med=list(map(lambda x: x[1], sorted(qaoa_median.items())))

plt.plot(dims,cm_avg,color='blue', label='Fixed-angle CM-QAOA mean a.f.', alpha=1)
plt.plot(dims,qaoa_avg,color='red', label='Fixed-angle QAOA mean a.f.', alpha=1)

plt.plot(dims,cm_med,color='blue',linestyle='dashed', label='cm median')
plt.plot(dims,qaoa_med,color='red',linestyle='dashed', label='qaoa median')

plt.xlabel(r"Lattice basis dimension", fontsize=12)
plt.ylabel(r"Approximation factor", fontsize=12)

ax=plt.gca()
ax.set_xticks(np.arange(4, 23, 1))  
ax.xaxis.set_major_formatter('{:.0f}'.format)  # Remove decimal points

plt.legend(    
    loc="lower center",           # Position the anchor point at the bottom
    bbox_to_anchor=(0.475, -0.3),   # Adjust the y-offset further down if needed
    ncol=2,
    fontsize=12,                # Adjust font size
    columnspacing=0.0,            # Reduce space between columns
    handlelength=1.0             # Reduce legend line length
)
#plt.savefig("approx_factors.pgf")  # Exports to a TikZ-compatible file
print("Not saving")
plt.show()
