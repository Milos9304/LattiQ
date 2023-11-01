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

#include "io/logger.h"

struct GeneratorParam{

	//all these are lambdas squared
	int lambda1 = 10;
	int lambda2 = 12;
	int lambda3_lb = 15;
	int lambda_ub = 17;

	int num_instances = 100;

	int n, m;
	int q = 7;

	GeneratorParam(int n){ //diagonal

		this->__diagonal = true;

		if(n < 3)
			throw_runtime_error("n should be larger than 2");

		if(lambda3_lb > lambda_ub || lambda2 >= lambda3_lb || lambda1 >= lambda2 || lambda1 < 1)
			throw_runtime_error("Invalid lambda values");

		this->n = n;
	}

	GeneratorParam(int n, int m){ //nondiagonal

		this->__diagonal = false;
		this->n = n;
		this->m = m;
	}

		bool __diagonal;
};
typedef std::function<std::vector<DiagonalHamiltonian>(GeneratorParam)> DiagonalInstanceGenerator;
typedef std::function<std::vector<HamiltonianWrapper>(GeneratorParam)> InstanceGenerator;

extern DiagonalInstanceGenerator generateDiagonalUniform;
extern DiagonalInstanceGenerator generateDiagonalExtensive;
extern InstanceGenerator generateQaryUniform;
extern InstanceGenerator generateQaryUniformFPLLLWay;

void calculateAverage(int n, DiagonalHamiltonian*);

void saveEigenspaceToFile(std::string filename, FastVQA::RefEnergies eigenspace);

#endif /* AVERAGED_SVP_H_ */
