import math

training_until=10

opt = False
simple_alpha=False
p_to_print = 3

c = 7

def get_c_alpha(ll):
    
    if simple_alpha:
    
        nom=0.
        den=0.
        for (l, dim) in ll:
            nom += math.log2(l)*dim
            den += dim*dim

        return (0, nom/den)       

    else:

        sum_yi=0.
        sum_xi=0.
        sum_xi2=0.
        sum_xi_yi=0.
        
        for (l, dim) in ll:
            sum_yi += math.log2(l)
            sum_xi+=dim
            sum_xi2+=dim*dim
            sum_xi_yi+=math.log2(l)*dim

        a = (sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(len(ll)*sum_xi2-sum_xi*sum_xi)
        b = (len(ll)*sum_xi_yi-sum_xi*sum_yi)/(len(ll)*sum_xi2-sum_xi*sum_xi)
    
    return (a,b)


if opt:
    pass
   
else:
 
  ms=[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

  if p_to_print == 3:
    
    if c == 1.7:
        ys={"CMQAOA": [0.0693292716715608,0.0478684307863023,0.02503258715718,0.013010958675684,0.00756397051018621,0.00404356014344025,0.00291991858919842,0.00192547400383744,0.00101099447042996,0.000467167756089037,0.000309612904156142,0.000158459439530568,6.70735204985571e-05,4.89165529224329e-05,2.64208971469158e-05,1.7675640133377e-05,6.75079000732792e-06,2.97459837880133e-06,1.81140015937729e-06],
            "QAOA non_pen": [0.165984611064018,0.105151678820791,0.0649284904747558,0.0531402319921266,0.0300615103406024,0.0180973354668209,0.0193490097476256,0.00782804389346197,0.00465994676806367,0.00511659640271405,0.00197546403001666,0.0010626163639798,0.000953802779592921,0.000458528408198807,0.000233537762325203,0.000186498827705746,8.10705197331645e-05,5.73023784526785e-05,3.94192989995113e-05]}
        #ys={"CMQAOA": [(0.069329271671561, 0.052979587718531), (0.047868430786302, 0.037381545336635), (0.025032587157180, 0.023652483637020), (0.012698733015311, 0.010940024752456), (0.006938764219226, 0.007202205463881), (0.003528632980062, 0.004269744970614), (0.002098228447024, 0.002432489298470), (0.00
    
    elif c == 3:
        
        
    
    elif c == 7:      
        ys={"CMQAOA": [0.449801453789462,0.460661287544699,0.826437516254708,0.984554977103171,0.99609375,0.998046875,0.9990234375,0.99951171875,0.999755859375,0.9998779296875,0.99993896484375,0.999969482421875,0.999984741210936,0.999992370605468,0.999996185302734,0.999998092651368,0.999999046325687,0.999999523162845,0.999999761581416],
        "QAOA non_pen": [0.583215092105652,0.620732658674798,0.908621169280694,0.972284842752436,0.972676344049525,0.987181970335781,0.98999260051623,0.994022368958332,0.997250448603223,0.998530373787053,0.998851934147088,0.99904653866116,0.999431644524317,0.999821725582917,0.999847436572696,0.999847666625252,0.999944598374456,0.999950459518479,0.999961140381018]}
           

cm_overlaps = []
qaoa_overlaps = []

print("x cm_overlap_sv qaoa_nonpen_overlap_sv alpha_line_cm alpha_line_qaoa fixed_used_in_training fixed_used_in_training2 alpha_cm alpha_qaoa")
for i in range(len(ms)):
    cm_overlap=ys["CMQAOA"][i]
    qaoa_overlap=ys["QAOA non_pen"][i]
    
    cm_overlaps.append((cm_overlap, ms[i]))
    qaoa_overlaps.append((qaoa_overlap, ms[i]))

line_cm=get_c_alpha(cm_overlaps)
line_qaoa=get_c_alpha(qaoa_overlaps)

for i in range(len(ms)):
    cm_overlap=ys["CMQAOA"][i]
    qaoa_overlap=ys["QAOA non_pen"][i]
    
    training2= 1
    training = 1
    if ms[i] > training_until:
        training = 0
    if ms[i] > training_until-1:
        training2 = 0
    
        
    print(ms[i], cm_overlap, qaoa_overlap, 2**(line_cm[0]+ms[i]*line_cm[1]), 2**(line_qaoa[0]+ms[i]*line_qaoa[1]), training, training2, round(line_cm[1],3), round(line_qaoa[1],3))
print(31, 0, 0, 2**(line_cm[0]+31*line_cm[1]), 2**(line_qaoa[0]+31*line_qaoa[1]), 0, 0, round(line_cm[1],3), round(line_qaoa[1],3))




    
