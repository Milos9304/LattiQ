import numpy as np

"""
    m=3, qs=4
    const int loglevel = 1;
    const int round_decimals = 5; //-1 undefined
    const int opt_strategy = 0;      //0=trivially, 1=rank_reduction
    const int num_steps = 15;
    const int ansatz_depth = 2;
    const double xtol = 10e-5;
    const double catol = 0.0002;
    const bool classical_esolver_compare = false;
    const bool outputLogToFile = false;
    const bool checkHessian = true;
    const bool printGroundStateOverlap = false;
    const bool print_eps = false;
    const int eval_limit_step = 600; //max iterations per step

    acceleratorOptions.log_level = loglevel;
    acceleratorOptions.logFileName = "aqc_pqc_log.txt";
    acceleratorOptions.roundDecimalPlaces = round_decimals;
    acceleratorOptions.optStrategy = opt_strategy;
    acceleratorOptions.accelerator_type = "quest";
    acceleratorOptions.nbSteps = num_steps;
    acceleratorOptions.ansatz_name = "Ry_Cz_nn_Ry";//"Ry_CNOT_nn_Rz_CNOT_Rz"
    acceleratorOptions.ansatz_depth = ansatz_depth;
    acceleratorOptions.xtol = xtol;
    acceleratorOptions.catol = catol;
    acceleratorOptions.compareWithClassicalEigenSolver = classical_esolver_compare;
    acceleratorOptions.outputLogToFile = outputLogToFile;
    acceleratorOptions.checkHessian = checkHessian;
    acceleratorOptions.printGroundStateOverlap = printGroundStateOverlap;
    acceleratorOptions.printEpsilons = print_eps;
    acceleratorOptions.eval_limit_step = eval_limit_step;
    acceleratorOptions.initialGroundState = FastVQA::InitialGroundState::PlusState;
    
"""
dat1=np.array([1,
0.938198,
1.63846e-05,
0.0463168,
0.00860362,
1.30969e-09,
1.51838e-06,
0.00170384,
4.08472e-11,
1,
1.28865e-06,
0.999614,
1.28865e-06,
4.35045e-06,
1.28865e-06,
0.999999,
0.0414371,
0.823495,
0.993975,
6.20475e-07,
0.0463168,
0.992856,
1.32755e-07,
0.312363,
9.06273e-07,
1.28865e-06,
1.30969e-09,
0.485674,
0.0105432,
1.15657e-05,
1.15657e-05,
5.7876e-07,
0.0269374,
1,
1.28865e-06,
0.83285,
0.0205302,
9.06273e-07,
0.823495,
0.992686,
0.823495,
0.938198,
0.936581,
0.701474,
0.999614,
0.312363,
1.63846e-05,
0.998869,
5.59816e-10,
0.701474,
0.83285,
0.0376059,
1,
1.30969e-09,
3.21677e-05,
1.32755e-07,
0.312363,
1,
2.39454e-05,
0.544519,
1.15657e-05,
0.0376059,
1,
1,
0.774869,
4.35045e-06,
4.35045e-06,
9.06273e-07,
0.938198,
9.16085e-09,
0.83285,
1,
5.59816e-10,
0.0105432,
5.59816e-10,
0.00860362,
3.21677e-05,
1.15657e-05,
1.30969e-09,
0.992686])

"""
    m=3, qs=4
    const int loglevel = 1;
    const int round_decimals = 5; //-1 undefined
    const int opt_strategy = 0;      //0=trivially, 1=rank_reduction
    const int num_steps = 20;
    const int ansatz_depth = 2;
    const double xtol = 10e-5;
    const double catol = 0.0002;
    const bool classical_esolver_compare = false;
    const bool outputLogToFile = false;
    const bool checkHessian = true;
    const bool printGroundStateOverlap = false;
    const bool print_eps = false;
    const int eval_limit_step = 600; //max iterations per step

    acceleratorOptions.log_level = loglevel;
    acceleratorOptions.logFileName = "aqc_pqc_log.txt";
    acceleratorOptions.roundDecimalPlaces = round_decimals;
    acceleratorOptions.optStrategy = opt_strategy;
    acceleratorOptions.accelerator_type = "quest";
    acceleratorOptions.nbSteps = num_steps;
    acceleratorOptions.ansatz_name = "Ry_Cz_nn_Ry";//"Ry_CNOT_nn_Rz_CNOT_Rz"
    acceleratorOptions.ansatz_depth = ansatz_depth;
    acceleratorOptions.xtol = xtol;
    acceleratorOptions.catol = catol;
    acceleratorOptions.compareWithClassicalEigenSolver = classical_esolver_compare;
    acceleratorOptions.outputLogToFile = outputLogToFile;
    acceleratorOptions.checkHessian = checkHessian;
    acceleratorOptions.printGroundStateOverlap = printGroundStateOverlap;
    acceleratorOptions.printEpsilons = print_eps;
    acceleratorOptions.eval_limit_step = eval_limit_step;
    acceleratorOptions.initialGroundState = FastVQA::InitialGroundState::PlusState;
    
"""
dat2=np.array([0.999696,
0.0031397,
7.22685e-05,
0.999999,
0.999988,
0.257978,
5.03812e-08,
1.80923e-05,
4.21564e-07,
1,
6.66419e-08,
0.0130603,
6.66419e-08,
1.53729e-07,
6.66419e-08,
0.0208859,
0.00572148,
0.999999,
1,
1.41095e-07,
0.999999,
0.999528,
8.95894e-08,
0.999975,
7.32035e-07,
6.66419e-08,
0.257978,
0.990299,
3.90414e-07,
2.91557e-09,
2.91557e-09,
7.89318e-09,
0.808332,
0.999984,
6.66419e-08,
0.998987,
2.76602e-05,
7.32035e-07,
0.999999,
0.342713,
0.999999,
0.0031397,
0.999999,
0.99454,
0.0130603,
0.999975,
7.22685e-05,
0.329955,
2.44199e-08,
0.99454,
0.998987,
0.999994,
0.999696,
0.257978,
6.32113e-08,
8.95894e-08,
0.999975,
1,
0.000225836,
0.999998,
2.91557e-09,
0.999994,
0.464712,
0.999984,
0.998964,
1.53729e-07,
1.53729e-07,
7.32035e-07,
0.0031397,
1.82726e-06,
0.998987,
1,
2.44199e-08,
3.90414e-07,
2.44199e-08,
0.999988,
6.32113e-08,
2.91557e-09,
0.257978,
0.342713])

"""
    m=4, qs=6
    const int loglevel = 1;
    const int round_decimals = 5; //-1 undefined
    const int opt_strategy = 0;      //0=trivially, 1=rank_reduction
    const int num_steps = 20;
    const int ansatz_depth = 2;
    const double xtol = 10e-5;
    const double catol = 0.0002;
    const bool classical_esolver_compare = false;
    const bool outputLogToFile = false;
    const bool checkHessian = true;
    const bool printGroundStateOverlap = false;
    const bool print_eps = false;
    const int eval_limit_step = 600; //max iterations per step

    acceleratorOptions.log_level = loglevel;
    acceleratorOptions.logFileName = "aqc_pqc_log.txt";
    acceleratorOptions.roundDecimalPlaces = round_decimals;
    acceleratorOptions.optStrategy = opt_strategy;
    acceleratorOptions.accelerator_type = "quest";
    acceleratorOptions.nbSteps = num_steps;
    acceleratorOptions.ansatz_name = "Ry_Cz_nn_Ry";//"Ry_CNOT_nn_Rz_CNOT_Rz"
    acceleratorOptions.ansatz_depth = ansatz_depth;
    acceleratorOptions.xtol = xtol;
    acceleratorOptions.catol = catol;
    acceleratorOptions.compareWithClassicalEigenSolver = classical_esolver_compare;
    acceleratorOptions.outputLogToFile = outputLogToFile;
    acceleratorOptions.checkHessian = checkHessian;
    acceleratorOptions.printGroundStateOverlap = printGroundStateOverlap;
    acceleratorOptions.printEpsilons = print_eps;
    acceleratorOptions.eval_limit_step = eval_limit_step;
    acceleratorOptions.initialGroundState = FastVQA::InitialGroundState::PlusState;
    
"""
dat3=np.array([7.16326e-05,
5.7081e-08,
0.000749623,
0.404922,
0.390136,
0.28661,
0.0254302,
3.13195e-06,
0.000130265,
0.000329007,
1.02331e-06,
2.42596e-10,
0.450907,
1.09312e-06,
0.280528,
4.33044e-07,
0.799416,
0.000329007,
7.60121e-07,
0.00553196,
0.000130265,
0.914344,
5.7081e-08,
8.33697e-12,
0.000599617,
0.914344,
0.00353751,
4.05546e-07,
2.12967e-14,
3.13136e-06,
1.09312e-06,
0.28661,
1,
6.23933e-06,
0.000963878,
7.60121e-07,
0.000312627,
1.60587e-19,
0.00353751,
0.904409,
0.000963878,
0.409777,
4.33452e-08,
0.00232136,
0.404922,
5.10063e-13,
0.000688935,
3.35383e-12,
2.47775e-09,
1.98781e-10,
7.54061e-06,
3.00448e-21,
0.308321,
2.12967e-14,
3.00448e-21,
0.0018398,
6.23933e-06,
0.00230632,
0.680814,
2.7758e-13,
0.132948,
0.404922,
2.7758e-13,
0.0013429,
0.964421,
0.0254302,
1.97603e-08,
3.22061e-09,
0.0013429,
0.280528,
3.22061e-09,
8.97994e-12,
1.02331e-06,
1.98781e-10,
1.33298e-17,
7.16326e-05,
1.31699e-06,
1.2471e-10,
0.000599617,
0.138355])

threshold=0.1
print("m=3")
print()
print("dat1")
print(" mean: ", np.mean(dat1))
print(" stdev: ", np.std(dat1))
print(" ", len(np.where(dat1 > threshold)[0]),"/",len(dat1)," solved with threshold=",threshold)

print("dat2")
print(" mean: ", np.mean(dat2))
print(" stdev: ", np.std(dat2))
print(" ", len(np.where(dat2 > threshold)[0]),"/",len(dat2)," solved with threshold=",threshold)

print("m=4")
print()
print("dat3")
print(" mean: ", np.mean(dat3))
print(" stdev: ", np.std(dat3))
print(" ", len(np.where(dat3 > threshold)[0]),"/",len(dat3)," solved with threshold=",threshold)

