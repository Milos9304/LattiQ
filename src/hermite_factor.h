/*
 * qaoa.h
 *
 *  Created on: Apr 10, 2024
 *      Author: Milos Prokop
 */

#ifndef SRC_HERMITE_FACTOR_H_
#define SRC_HERMITE_FACTOR_H_

#include "FastVQA/fastVQA.h"

double calculate_hermite_factor(int n, double volume_1n, std::shared_ptr<Qureg>, FastVQA::RefEnergies* refEnergies);

#endif /* SRC_HERMITE_FACTOR_H_ */
