/*
 * test_runner.h
 *
 *  Created on: Nov 29, 2023
 *      Author: Milos Prokop
 */

#ifndef SRC_TEST_RUNNER_H_
#define SRC_TEST_RUNNER_H_

#include "io/sql_io.h"

void test_variable_substitution(MapOptions* mapOptions);
void test_execution_time(FastVQA::QAOAOptions* qaoaOptions, Database* database);

#endif /* SRC_TEST_RUNNER_H_ */
