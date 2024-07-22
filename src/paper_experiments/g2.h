/*
 * g1.h
 *
 *  Created on: May 20, 2023
 *      Author: Milos Prokop
 */

#ifndef SRC_G2_H_
#define SRC_G2_H_

#include "FastVQA/fastVQA.h"
#include "lattice/lattice.h"
//#include "io/sql_io.h"
#include <string>

class G2{

	public:
		int loglevel = 1;

	G2(int loglevel, FastVQA::QAOAOptions*, MapOptions*);

	void run();

	private:

	int p = 6;
	int q = 97;
	int n = 2;
	int m_start = 3;
	int m_end = 20;
	int max_num_instances = 100;//1000;

	int num_params = p*2;

	//it all starts with 3
	//std::pair<int, int>(4, 10) means training is done on 3 and 4

	struct AlphaMinimizationExperimentInstance{
		FastVQA::PauliHamiltonian h;
		FastVQA::RefEnergies zero_solutions;
		FastVQA::RefEnergies sv_solutions;

		double volume;

		int num_qubits_per_dim;
		int q,m,n;
	};

	double _cost_fn(std::vector<AlphaMinimizationExperimentInstance>, const double *angles, bool use_database=false);

	FastVQA::QAOAOptions* qaoaOptions;
	MapOptions* mapOptions;
};

#endif /*SRC_G2_H_*/
