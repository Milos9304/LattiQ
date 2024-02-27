/*
 * experiment_runner.h
 *
 *  Created on: Nov 20, 2023
 *      Author: Milos Prokop
 */

#ifndef SRC_EXPERIMENTRUNNER_H_
#define SRC_EXPERIMENTRUNNER_H_

#include "FastVQA/fastVQA.h"
#include "lattice/lattice.h"
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

	struct Cost{
		double mean;
		double stdev;
		double mean_num_of_sols;

		Cost(){}

		Cost(double mean, double stdev, double mean_num_of_sols){
			this->mean = mean;
			this->stdev = stdev;
			this->mean_num_of_sols=mean_num_of_sols;
		}
	};

	int loglevel = 1;

	int q = 97;
	int n = 1;
	int m = 4;

	int max_num_instances = 100;

	FastVQA::QAOAOptions* qaoaOptions;

	AngleSearchExperiment(int loglevel, FastVQA::QAOAOptions*, MapOptions*);

	void run();

private:

	struct Instance{
		FastVQA::PauliHamiltonian h;
		FastVQA::RefEnergies solutions;
		qreal min_energy;
	};

	const double test_ratio = 0.2;
	int num_instances;
	int num_train_instances;
	int num_test_instances;

	int nbQubits;
	int num_params;

	std::vector<Instance> train_set;
	std::vector<Instance> test_set;

	void _generate_dataset(MapOptions*);
	Cost _cost_fn(std::vector<Instance>*, double *angles);
	void run_p1();
	void run_p2();
	void run_p2_test();

};

void experiment_runner(ExperimentSetup*, std::string experiment_name, int loglevel);

#endif /* SRC_EXPERIMENTRUNNER_H_ */
