/*
 * experiment_runner.h
 *
 *  Created on: Nov 20, 2023
 *      Author: Milos Prokop
 */

#ifndef SRC_EXPERIMENTRUNNER_H_
#define SRC_EXPERIMENTRUNNER_H_

#include "FastVQA/fastVQA.h"
#include <string>

class ExperimentSetup{
public:

	std::string experiment_type;
	FastVQA::PauliHamiltonian* hamiltonian;
	FastVQA::QAOAOptions* qaoaOptions;
	qreal minimum_energy;

	//case specific
	int num_rand_params;
};

class AngleSearchExperiment{

public:
	int q = 97;
	int n = 1;
	int m = 4;

	int max_num_instances = 100;



	AngleSearchExperiment();

	void run();

private:
	const double test_ratio = 0.2;
	int num_instances;
	int num_train_instances;
	int num_test_instances;

	void _generate_dataset();
};

void experiment_runner(ExperimentSetup*, std::string experiment_name);

#endif /* SRC_EXPERIMENTRUNNER_H_ */
