/*
 * averaged_svp.h
 *
 *  Created on: 17 Oct, 2023
 *      Author: Milos Prokop
 */

#ifndef AVERAGED_SVP_H_
#define AVERAGED_SVP_H_

#include "FastVQA/fastVQA.h"
#include <Eigen/Dense>

//typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> HamiltonianD;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Hamiltonian;
typedef Eigen::Vector<int, Eigen::Dynamic> DiagonalHamiltonian;

class HamiltonianWrapper{
public:
	std::string name;
	Hamiltonian hamiltonian;
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> K;

	HamiltonianWrapper(Hamiltonian hamiltonian, std::string name){
		this->hamiltonian=hamiltonian;
		this->name=name;
	}
};

/*class HamiltonianWrapperD{
public:
	std::string name;
	HamiltonianD hamiltonian;
	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> K;

	HamiltonianWrapperD(HamiltonianD hamiltonianD, std::string name){
		this->hamiltonian=hamiltonianD;
		this->name=name;
	}
};*/

#include "io/logger.h"

struct GeneratorParam{

	//all these are lambdas squared
	int lambda1 = 10;
	int lambda2 = 12;
	int lambda3_lb = 15;
	int lambda_ub = 17;

	int num_instances = 100;

	int n, m;
	int q;

	bool shuffle;
	int seed=0;

	int sol_elem_bound = 1;

	int cutoff = -1;

	/*GeneratorParam(int q, int n, bool shuffle=false, int shuffle_seed=0){ //diagonal

		this->__diagonal = true;
		this->q = q;

		if(n < 3)
			throw_runtime_error("n should be larger than 2");

		if(lambda3_lb > lambda_ub || lambda2 >= lambda3_lb || lambda1 >= lambda2 || lambda1 < 1)
			throw_runtime_error("Invalid lambda values");

		this->n = n;

		this->shuffle = shuffle;
		this->shuffle_seed = shuffle_seed;
	}*/

	GeneratorParam(int q, int n, int m, bool shuffle=false, int shuffle_seed=0, int cutoff=-1){ //nondiagonal

		this->__diagonal = false;
		this->q = q;
		this->n = n;
		this->m = m;

		this->shuffle = shuffle;
		this->seed = shuffle_seed;

		this->cutoff = cutoff;

		if(this->num_instances > cutoff && cutoff > 0){
			this->num_instances = cutoff;
		}

	}

		bool __diagonal;
};
typedef std::function<std::vector<DiagonalHamiltonian>(GeneratorParam)> DiagonalInstanceGenerator;
typedef std::function<std::vector<HamiltonianWrapper>(GeneratorParam)> InstanceGenerator;

extern DiagonalInstanceGenerator generateDiagonalUniform;
extern DiagonalInstanceGenerator generateDiagonalExtensive;
extern InstanceGenerator generateQaryUniform;
extern InstanceGenerator generateQaryUniformFPLLLWay;
extern InstanceGenerator generateFromEvalDecomposition;

void calculateAverage(int n, DiagonalHamiltonian*);

void saveEigenspaceToFile(std::string filename, FastVQA::RefEnergies eigenspace);

#endif /* AVERAGED_SVP_H_ */
