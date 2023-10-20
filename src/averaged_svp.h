/*
 * averaged_svp.h
 *
 *  Created on: 17 Oct, 2023
 *      Author: Milos Prokop
 */

#ifndef AVERAGED_SVP_H_
#define AVERAGED_SVP_H_

#include <Eigen/Dense>
typedef Eigen::Vector<int, Eigen::Dynamic> DiagonalHamiltonian;

#include "io/logger.h"

struct GeneratorParam{

	//all these are lambdas squared
	int lambda0 = 0;
	int lambda1 = 10;
	int lambda2 = 12;
	int lambda3_lb = 15;
	int lambda_ub = 17;

	int num_instances = 100;

	int n;

	GeneratorParam(int n){
		if(n < 4)
			throw_runtime_error("n should be larger than 3");

		if(lambda3_lb > lambda_ub || lambda2 >= lambda3_lb || lambda1 >= lambda2 || lambda1 < 1)
			throw_runtime_error("Invalid lambda values");

		this->n = n;
	}
};
typedef std::function<std::vector<DiagonalHamiltonian>(GeneratorParam)> InstanceGenerator;
extern InstanceGenerator generateDiagonalUniform;
extern InstanceGenerator generateDiagonalExtensive;

void calculateAverage(int n, DiagonalHamiltonian*);

#endif /* AVERAGED_SVP_H_ */
