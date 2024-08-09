#include "experiment_runner.h"

void AqcPqcExperiment::run(){

	const int round_decimals = 5; //-1 undefined
	const int opt_strategy = 0;	  //0=trivially, 1=rank_reduction

	FastVQA::AqcPqcAcceleratorOptions acceleratorOptions;
	FastVQA::AqcPqcAccelerator accelerator(acceleratorOptions);

	acceleratorOptions.log_level = qaoaOptions->accelerator->log_level;
	acceleratorOptions.roundDecimalPlaces = round_decimals;
	acceleratorOptions.optStrategy = opt_strategy;
	acceleratorOptions.accelerator_type = "quest";
/*	acceleratorOptions.nbSteps = num_steps->value();
	acceleratorOptions.ansatz_name = "Ry_Cz_nn_Ry";//"Ry_CNOT_nn_Rz_CNOT_Rz"
	acceleratorOptions.ansatz_depth = ansatz_depth->value();
	acceleratorOptions.xtol = xtol->value();
	acceleratorOptions.catol = catol->value();
	acceleratorOptions.compareWithClassicalEigenSolver = classical_esolver->is_set();
	acceleratorOptions.outputLogToFile = true;
	acceleratorOptions.checkHessian = true;
	acceleratorOptions.printGroundStateOverlap = true;
	acceleratorOptions.printEpsilons = eps_print->is_set();
	acceleratorOptions.eval_limit_step = eval_limit_step->value();
	acceleratorOptions.initialGroundState = FastVQA::InitialGroundState::PlusState;*/

}
