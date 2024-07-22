/*
 * g1.h
 *
 *  Created on: May 20, 2023
 *      Author: Milos Prokop
 */

#ifndef SRC_G1_H_
#define SRC_G1_H_

#include "FastVQA/fastVQA.h"
#include "lattice/lattice.h"
//#include "io/sql_io.h"
#include <string>

class G1{

	public:
		int loglevel = 1;

	G1(int loglevel, FastVQA::QAOAOptions*, MapOptions*);

	void run();

	private:

	int p = 6;
	int q = 97;
	int n = 2;
	int m_start = 3;
	int m_end = 15;
	int max_num_instances = 1000;//1000;
	double test_ratio = 0.2;

	int num_params = p*2;

	//it all starts with 3
	//std::pair<int, int>(4, 10) means training is done on 3 and 4
	std::vector<std::pair<int, int>> train_dim_ratios { std::pair<int, int>(4, m_end), std::pair<int, int>(5, m_end), std::pair<int, int>(6, m_end), std::pair<int, int>(7, m_end), std::pair<int, int>(8, m_end), std::pair<int, int>(9, m_end) };

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

#endif /*SRC_G1_H_*/
