#include "experiment_runner.h"
#include <algorithm>

void AqcPqcExperiment::run(){

	/*
	 * DEFAULTS
	 *
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
	const int eval_limit_step = 600; //max iterations per step*/

	const int loglevel = 1;
	const int round_decimals = 5; //-1 undefined
	const int opt_strategy = 0;	  //0=trivially, 1=rank_reduction
	const int num_steps = 20;
	const int ansatz_depth = 2;
	const double xtol = 10e-5;
	const double catol = 0.0002;
	const bool classical_esolver_compare = false;
	const bool outputLogToFile = false;
	const bool checkHessian = true;
	const bool printGroundStateOverlap = false;
	const bool print_eps = false;
	const int eval_limit_step = 6000; //max iterations per step


	FastVQA::AqcPqcAcceleratorOptions acceleratorOptions;

	/*
	 * DEFAULTS
	 *
	 * acceleratorOptions.log_level = loglevel;
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
	acceleratorOptions.initialGroundState = FastVQA::InitialGroundState::PlusState;*/

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

	//this->mapOptions->penalty = 0;

	std::vector<int> num_iters;
	std::vector<double> final_overlap;

	for(int m = m_start; m <= m_end; ++m){

		//std::vector<Instance> dataset = _generate_dataset(3, m, true); //3 is arbitrary, true is for penalise
		std::vector<Instance> dataset = _generate_dataset(1, 3, true);
		loge("very small instance");

		int i = 0;
		for(auto &instance: dataset){

			if(i > 0){
				logw("aqcpqc.cpp breaking too soon");
				break;
			}

			std::vector<long long int> solutions;
			for(auto &sol: instance.sv_solutions){
				solutions.push_back(sol.index);
			}
			acceleratorOptions.solutions = solutions;

			FastVQA::AqcPqcAccelerator accelerator(acceleratorOptions);

			/*logw("Num sols: " + std::to_string(acceleratorOptions.solutions.size()));
			for(const auto &sol: solutions){
				std::cerr<<"solution index="<<sol<<" with value="<<instance.sv1Squared<<std::endl;
			}
			logw("SV1Squared: " + std::to_string(instance.sv1Squared));
			for(const auto &sol: instance.zero_solutions){
				std::cerr<<sol.value<<" "<<sol.index<<std::endl;
			}*/

			FastVQA::PauliHamiltonian h0(instance.h.nbQubits);


			Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> matrix = instance.h.getMatrixRepresentation2(true);
			std::vector<double> energies;
			/*for(int j = 0; j < matrix.cols(); ++j){
				std::cerr<<j<<":  "<<matrix(j,j)<<std::endl;
				energies.push_back(matrix(j,j));
			}

			std::sort(energies.begin(), energies.end());
			for(auto &a:energies)
				std::cerr<<a<<" ";*/

			h0.initializeSumMinusSigmaXHamiltonian();

			accelerator.initialize(&h0, &instance.h);
			accelerator.run();

			i++;
		}

		logw("Breaking after first m");
		break;
	}

}
