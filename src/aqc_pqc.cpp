#include "experiment_runner.h"

void AqcPqcExperiment::run(){

	const int loglevel = 1;
	const int round_decimals = 5; //-1 undefined
	const int opt_strategy = 0;	  //0=trivially, 1=rank_reduction
	const int num_steps = 20;
	const int ansatz_depth = 1;
	const double xtol = 10e-5;
	const double catol = 0.0002;
	const bool classical_esolver_compare = false;
	const bool outputLogToFile = true;
	const bool checkHessian = true;
	const bool printGroundStateOverlap = true;
	const bool print_eps = false;
	const int eval_limit_step = 600; //max iterations per step


	FastVQA::AqcPqcAcceleratorOptions acceleratorOptions;

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

	this->mapOptions->penalty = 0;

	for(int m = m_start; m <= m_end; ++m){

		std::vector<Instance> dataset = _generate_dataset(3, m, false); //3 is arbitrary

		int i = 0;
		for(auto &instance: dataset){

			if(i > 0){
				loge("aqcpqc.cpp breaking too soon");
				break;
			}

			std::vector<long long int> solutions;
			for(auto &sol: instance.zero_solutions){
				solutions.push_back(sol.index);
			}
			acceleratorOptions.solutions = solutions;

			FastVQA::AqcPqcAccelerator accelerator(acceleratorOptions);

			loge(std::to_string(acceleratorOptions.solutions.size()));

			FastVQA::PauliHamiltonian h0(instance.h.nbQubits);
			h0.initializeSumMinusSigmaXHamiltonian();

			accelerator.initialize(&h0, &instance.h);
			accelerator.run();

			i++;
		}


	}

}
