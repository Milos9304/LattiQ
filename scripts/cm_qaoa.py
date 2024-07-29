import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

class colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

data="""1
2 0.728891 0.25 0.00267001 0.402741 
2 0.698801 0.25 0.0326716 0.40603 
2 0.519592 0.25 0.198941 0.0260714 
3 0.0312376 0.015625 0.0558158 0.0729225 
3 0.015625 0.015625 0.03125 0.0505795 
3 0.0446095 0.015625 0.0315811 0.0313254 
4 0.00179637 0.00390625 0.00553131 0.00657018 
4 0.000493654 0.00390625 0.00935024 0.0078125 
4 0.00379503 0.00390625 0.00533572 0.0208294 
2
2 0.936577 0.25 0.0451636 0.464303 
2 0.96532 0.25 0.019896 0.400004 
2 0.869821 0.25 0.0436346 0.298971 
3 0.0118131 0.015625 0.0164662 0.0696785 
3 0.0169913 0.015625 0.0155935 0.00751973 
3 0.0100635 0.015625 0.0347947 0.0910574 
4 0.000263423 0.00390625 0.00426993 0.00726747 
4 0.00184613 0.00390625 0.0165165 0.0161763 
4 0.0258258 0.00390625 0.00847047 0.025208 
3
2 0.860382 0.25 0.130765 0.691189 
2 0.98973 0.25 0.00682257 0.743829 
2 0.814786 0.25 0.00839911 0.68186 
3 0.0456591 0.015625 0.0452498 0.0572917 
3 0.0567083 0.015625 0.0369633 0.0916146 
3 0.0153883 0.015625 0.0997666 0.0495503 
4 0.00183617 0.00390625 0.0148235 0.00277628 
4 0.0104327 0.00390625 0.0157574 0.00306509 
4 0.0136405 0.00390625 0.0185777 0.00998913 
4
2 0.772069 0.25 0.201537 0.749976 
2 0.959696 0.25 0.0220731 0.749965 
2 0.967756 0.25 0.000243993 0.748451 
3 0.0668921 0.015625 0.0472338 0.0278294 
3 0.00184744 0.015625 0.0770256 0.0796467 
3 0.00385664 0.015625 0.0561148 0.015042 
4 0.0163985 0.00390625 0.0141291 0.0197953 
4 0.0012032 0.00390625 0.00415686 0.00129618 
4 0.00141544 0.00390625 0.0226238 0.00352091 
5
2 0.896832 0.25 0.0948165 0.749987 
2 0.979264 0.25 0.0170463 0.749994 
2 0.997081 0.25 0.000268137 0.738859 
3 0.0160565 0.015625 0.0170091 0.0411613 
3 0.0563038 0.015625 0.0664793 0.0189654 
3 0.0124285 0.015625 0.0756633 0.0777794 
4 0.00262334 0.00390625 0.0134484 0.00556395 
4 0.000475909 0.00390625 0.00152425 0.00701392 
4 0.00277265 0.00390625 0.00434552 0.0033788 
6
2 0.996728 0.25 0.00028327 0.749661 
2 0.970172 0.25 0.0175306 0.75 
2 0.990857 0.25 0.004548 0.749969 
3 0.123537 0.015625 0.0425462 0.081076 
3 0.0223283 0.015625 0.111188 0.0405938 
3 0.0422441 0.015625 0.0618582 0.0524234 
4 0.00150277 0.00390625 0.00849918 0.00539635 
4 0.00914574 0.00390625 0.00700732 0.011399 
4 0.000606327 0.00390625 0.00751151 0.00304072"""

from pathlib import Path
data = Path('m_2_5').read_text()
data=data.split('\n')[:-1]

def meann(data):
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
  var = variance(data)
  std_dev = math.sqrt(var)
  return std_dev

def mean(data, rnd_guess):
    return (meann(data)/rnd_guess, stdev(data)/rnd_guess)

plot_zero_qaoa=[[],[],[],[],[],[]]
plot_zero_cm=[[],[],[],[],[],[]]
plot_sv_qaoa=[[],[],[],[],[],[]]
plot_sv_cm=[[],[],[],[],[],[]]

zero_svs=[]
zero_cms=[]
zero_sv_qaoas=[]
sv_cm_qaoas=[]

prev_m=int(data[2][0])
iterator=0

set_ms=set()
for d in data:
    iterator+=1
    if len(d) == 1:
       p = int(d)
       ms=[]
       continue

    s=d.split(' ')
    print(d, s)
    m=int(s[0])

    zero_sv=float(s[1])
    zero_cm=float(s[2])
    sv_qaoa=float(s[3])
    sv_cm_qaoa=float(s[4])

    if prev_m != m or iterator==len(data):
        if len(ms) == 0:
            shift = -1
        else:
            shift = 0

        rnd_guess=(zero_cms[0],0)

        if p == 6:
            print("P:",p-1+shift,"   ",m,"    ",rnd_guess[0],"    ",zero_svs)

        plot_zero_qaoa[p-1+shift].append(mean(zero_svs,rnd_guess[0]))
        plot_zero_cm[p-1+shift].append(rnd_guess)
        plot_sv_qaoa[p-1+shift].append(mean(zero_sv_qaoas,rnd_guess[0]))
        plot_sv_cm[p-1+shift].append(mean(sv_cm_qaoas,rnd_guess[0]))

        ms = []
        zero_svs=[]
        zero_cms=[]
        zero_sv_qaoas=[]
        sv_cm_qaoas=[]

    ms.append(m)
    set_ms.add(m)
    zero_svs.append(zero_sv)
    zero_cms.append(zero_cm) 
    zero_sv_qaoas.append(sv_qaoa)
    sv_cm_qaoas.append(sv_cm_qaoa)

    prev_m=m



fig, ax = plt.subplots(2,3)
print(set_ms)
for p in range(1,6+1):
    print(p, plot_zero_cm[p-1])
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].errorbar(np.array(list(set_ms))-0.15, list(map(lambda x: x[0], plot_zero_qaoa[p-1])), list(map(lambda x: x[1], plot_zero_qaoa[p-1])), linestyle='None', marker='o', color='blue', capsize=5, label="zero_QAOA")
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].errorbar(np.array(list(set_ms))-0.05, list(map(lambda x: x[0], plot_zero_cm[p-1])), list(map(lambda x: x[1], plot_zero_cm[p-1])), linestyle='None', marker='o', color='red', capsize=5, label="zero_QAOA2")
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].errorbar(np.array(list(set_ms))+0.05, list(map(lambda x: x[0], plot_sv_qaoa[p-1])), list(map(lambda x: x[1], plot_sv_qaoa[p-1])), linestyle='None', marker='x', color='blue', capsize=5, label="sv_QAOA")
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].errorbar(np.array(list(set_ms))+0.15, list(map(lambda x: x[0], plot_sv_cm[p-1])), list(map(lambda x: x[1], plot_sv_cm[p-1])), linestyle='None', marker='x', color='red', capsize=5, label="sv_QAOA2")
    #fig.legend(["zero_qaoa","zero_cm","sv_qaoa","sv_cm"])
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].set_title("p="+str(p))
#plt.show()
plt.legend()
plt.show()

def f(n, alpha, c):
    return 1/2**(n*alpha*np.ceil(np.log2(n))+c)

def f_simpl(n, alpha):
    return 1/2**(n*alpha*np.ceil(np.log2(n)))

for p in range(1,6+1):
    
    #print(colors.ENDC+"p="+str(p))

    #print("zero_qaoa")
    popt, pcov = curve_fit(f, [2,3,4,5], list(map(lambda x: x[0], plot_zero_qaoa[p-1])), maxfev=30000)
    alpha_zero=popt[0]
    c_zero=popt[1]
    #print("alpha, c: ", popt, "\nCOV:",pcov)
    #print("zero_cm")
    popt, pcov = curve_fit(f, [2,3,4,5], list(map(lambda x: x[0], plot_zero_cm[p-1])), maxfev=30000)
    alpha_zero_cm=popt[0]
    c_zero_cm=popt[1]
    #print("alpha, c: ", popt, "\nCOV:",pcov)

    #print("sv_qaoa")
    popt, pcov = curve_fit(f, [2,3,4,5], list(map(lambda x: x[0], plot_sv_qaoa[p-1])), maxfev=30000)
    alpha_sv=popt[0]
    c_sv=popt[1]
    #print("alpha, c: ", popt, "\nCOV:",pcov)
    #print("sv_cm")
    popt, pcov = curve_fit(f, [2,3,4,5], list(map(lambda x: x[0], plot_sv_cm[p-1])), maxfev=30000)
    alpha_sv_cm=popt[0]
    c_sv_cm=popt[1]

    popt, pcov = curve_fit(f_simpl, [2,3,4,5], list(map(lambda x: x[0], plot_zero_qaoa[p-1])), maxfev=30000)
    s_alpha_zero=popt[0]
    #print("alpha, c: ", popt, "\nCOV:",pcov)
    #print("zero_cm")
    popt, pcov = curve_fit(f_simpl, [2,3,4,5], list(map(lambda x: x[0], plot_zero_cm[p-1])), maxfev=30000)
    s_alpha_zero_cm=popt[0]
    #print("alpha, c: ", popt, "\nCOV:",pcov)

    #print("sv_qaoa")
    popt, pcov = curve_fit(f_simpl, [2,3,4,5], list(map(lambda x: x[0], plot_sv_qaoa[p-1])), maxfev=30000)
    s_alpha_sv=popt[0]
    #print("alpha, c: ", popt, "\nCOV:",pcov)
    #print("sv_cm")
    popt, pcov = curve_fit(f_simpl, [2,3,4,5], list(map(lambda x: x[0], plot_sv_cm[p-1])), maxfev=30000)
    s_alpha_sv_cm=popt[0]
    #

    x_min = 2  
    x_max = 8                                #min/max values for x axis
    x_fit = np.linspace(x_min, x_max, 100) 
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].plot(x_fit, f(x_fit, alpha_sv, c_sv),label="qaoa_fit")
    ax[int((p-1)/3)][(p-1)-int((p-1)/3)*3].plot(x_fit, f(x_fit, alpha_sv_cm, c_sv_cm),label="cm_qaoa_fit")
    #fig.legend(["zero_qaoa","zero_cm","sv_qaoa","sv_cm", "qaoa_fit", "cm_qaoa_fit"])
    print(colors.ENDC + "{: >30} {: >30} {: >30} {: >30} {: >30}".format(*["p="+str(p), "QAOA", "CM_QAOA", "QAOA_simple_fit", "CM_QAOA_simple_fit"]))
    print(colors.RED + "{: >30} {: >30} {: >30} {: >30} {: >30}".format(*["ZERO",
            "2^-("+format(alpha_zero, '.2f')+"nlog(n)"+('+' if c_zero > 0 else '')+format(c_zero, '.2f')+")",
            "2^-("+format(alpha_zero_cm, '.2f')+"nlog(n)"+('+' if c_zero_cm > 0 else '')+format(c_zero_cm, '.2f')+")",
            "2^-("+format(s_alpha_zero, '.2f')+"nlog(n))",
            "2^-("+format(s_alpha_zero_cm, '.2f')+"nlog(n))"]))
    print(colors.BLUE + "{: >30} {: >30} {: >30} {: >30} {: >30}".format(*["GROUND_STATE",
            "2^-("+format(alpha_sv, '.2f')+"nlog(n)"+('+' if c_sv > 0 else '')+format(c_sv, '.2f')+")",
            "2^-("+format(alpha_sv_cm, '.2f')+"nlog(n)"+('+' if c_sv_cm > 0 else '')+format(c_sv_cm, '.2f')+")",
            "2^-("+format(s_alpha_sv, '.2f')+"nlog(n))",
            "2^-("+format(s_alpha_sv_cm, '.2f')+"nlog(n))"]))
    print()
print("I RESCALED EVERYTHING WITH RAND_GUESS!")
plt.legend(loc="upper right")
plt.show()
