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
#include "io/sql_io.h"
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

class CmQaoaExperiment{
public:

	const int q = 97;
	const int n = 1;

	const int m_start = 2;
	const int m_end   = 5;//8;

	const int p_start = 1;
	const int p_end   = 6;

	const int num_instances = /*9*/6;

	FastVQA::QAOAOptions* qaoaOptions;
	MapOptions* mapOptions;

	CmQaoaExperiment(FastVQA::QAOAOptions*, MapOptions*, Database*, int loglevel=1);

	void run();


private:

	int loglevel=1;
	Database* database;

	struct Instance{
		FastVQA::PauliHamiltonian h;
		FastVQA::RefEnergies zero_solutions;
		FastVQA::RefEnergies sv_solutions;
		//FastVQA::RefEnergies eigenspace; //for debug
		//qreal min_energy;
		//qreal random_guess;

		double volume;
		//int sv1Squared;

		int num_qubits_per_dim;
		int q,m,n;
	};

	struct Cost{
		double prob_zero;
		double prob_first_excited_state;

		long long int first_excited_state;
		int degeneracy;

		//Cost(){}
	};

	Cost _cost_fn(Instance*, bool use_database=false);


};

class AngleExperimentBase{
public:

	FastVQA::QAOAOptions* qaoaOptions;

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
protected:

	Database* database;

	struct Instance{
			FastVQA::PauliHamiltonian h;
			FastVQA::RefEnergies solutions;
			FastVQA::RefEnergies eigenspace; //for debug
			qreal min_energy;
			qreal random_guess;

			double volume;
			int sv1Squared;

			int q,m,n;
	};

	Cost _cost_fn(std::vector<Instance>*, const double *angles, bool use_database=false);
};



class AngleResultsExperiment : AngleExperimentBase{

public:

	int q = 97;

	int m_start = 4;
	int m_end = 10;

	int max_num_instances = 3000;

	//const std::vector<double> angles{0.4,0.48,5.56,0.28};
	const std::vector<double> angles{0.4,0.5,0.36,0.2,0.04,0.5};

	int loglevel = 1;
	AngleResultsExperiment(int loglevel, FastVQA::QAOAOptions*, MapOptions*, Database*);

	MapOptions* mapOptions;

	void run();

private:
	std::vector<Instance> _generate_dataset(int n, int m);

};

class AngleSearchExperiment : AngleExperimentBase{

public:

	int loglevel = 1;

	int q = 97;
	int n = 1;
	int m = 5;//7;

	int max_num_instances = 100;

	AngleSearchExperiment(int loglevel, FastVQA::QAOAOptions*, MapOptions*);

	void run();

private:
	long long int num_instances;
	int nbQubits;
	const double test_ratio = 0.2;
	int num_train_instances;
	int num_test_instances;

	int num_params;

	std::vector<Instance> train_set;
	std::vector<Instance> test_set;

	void _generate_dataset(MapOptions*);
	void run_p1();
	void run_p2();
	void run_p2_full_bruteforce();
	void run_p3_full_bruteforce();
	void run_p2_test();
	void run_cobyla_p2();

};

void experiment_runner(ExperimentSetup*, std::string experiment_name, int loglevel);

#endif /* SRC_EXPERIMENTRUNNER_H_ */
