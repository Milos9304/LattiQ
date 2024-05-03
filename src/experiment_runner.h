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
			FastVQA::RefEnergies zero_solutions;
			FastVQA::RefEnergies sv_solutions;
			//FastVQA::RefEnergies eigenspace; //for debug
			qreal min_energy;
			qreal random_guess;

			double volume;
			int sv1Squared;

			int q,m,n;
	};

	Cost _cost_fn(std::vector<Instance>*, const double *angles, bool use_database=false);
};

/*
(0.4, 0.7, 2.55, 0.6) m=1.20085 std=0.735449       ] [01m:23s<08m:21s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.4, 0.7, 2.55, 0.7) m=1.20233 std=0.763791       ] [01m:23s<08m:21s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.6, 2.1, 0.5) m=1.20024 std=0.714068       ] [01m:31s<08m:14s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.6, 2.1, 0.6) m=1.2135 std=0.741774        ] [01m:31s<08m:14s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.6, 2.1, 0.7) m=1.22119 std=0.774281       ] [01m:31s<08m:14s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.6, 2.6, 0.7) m=1.20686 std=0.690428       ] [01m:31s<08m:14s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 0.2, 0.6) m=1.21386 std=0.772585       ] [01m:31s<08m:13s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 0.2, 0.7) m=1.22436 std=0.792591       ] [01m:31s<08m:13s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 2.1, 0.4) m=1.20024 std=0.751449       ] [01m:32s<08m:13s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 2.1, 0.5) m=1.22297 std=0.776773       ] [01m:32s<08m:13s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 2.1, 0.6) m=1.2407 std=0.806287        ] [01m:32s<08m:13s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 2.1, 0.7) m=1.25273 std=0.838836       ] [01m:32s<08m:13s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 2.6, 0.6) m=1.2053 std=0.728198        ] [01m:32s<08m:12s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.45, 0.7, 2.6, 0.7) m=1.2153 std=0.746282        ] [01m:32s<08m:12s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.4, 1.8, 0.7) m=1.20513 std=0.719734        ] [02m:16s<07m:33s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.5, 1.1, 0.7) m=1.20276 std=0.772552        ] [02m:17s<07m:32s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.5, 1.8, 0.6) m=1.20255 std=0.704427        ] [02m:17s<07m:32s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.5, 1.8, 0.7) m=1.20999 std=0.726975        ] [02m:17s<07m:32s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.6, 1.1, 0.6) m=1.20976 std=0.794958        ] [02m:18s<07m:31s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.6, 1.1, 0.7) m=1.21748 std=0.816357        ] [02m:18s<07m:31s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.6, 1.45, 0.6) m=1.20263 std=0.835759       ] [02m:18s<07m:31s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.6, 1.45, 0.7) m=1.20994 std=0.861156       ] [02m:18s<07m:31s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.6, 1.8, 0.7) m=1.20376 std=0.745728        ] [02m:18s<07m:31s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.7, 1.1, 0.5) m=1.20279 std=0.825044        ] [02m:19s<07m:30s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.7, 1.1, 0.6) m=1.21699 std=0.851496        ] [02m:19s<07m:30s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.7, 1.1, 0.7) m=1.2251 std=0.872055         ] [02m:19s<07m:30s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.7, 1.45, 0.6) m=1.2099 std=0.864439        ] [02m:19s<07m:30s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(0.7, 0.7, 1.45, 0.7) m=1.21761 std=0.887912       ] [02m:19s<07m:30s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.2, 0.7, 1.6, 0.7) m=1.20571 std=0.72866         ] [03m:55s<05m:58s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.4, 1.4, 0.6) m=1.20742 std=0.657616       ] [05m:16s<04m:37s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.4, 1.4, 0.7) m=1.21249 std=0.693342       ] [05m:16s<04m:37s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.5, 1.4, 0.5) m=1.21913 std=0.650394       ] [05m:17s<04m:36s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.5, 1.4, 0.6) m=1.23374 std=0.679792       ] [05m:17s<04m:36s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.5, 1.4, 0.7) m=1.2402 std=0.709574        ] [05m:17s<04m:36s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.6, 0.8, 0.7) m=1.2044 std=0.719789        ] [05m:18s<04m:35s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.6, 1.4, 0.4) m=1.21706 std=0.670801       ] [05m:19s<04m:35s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.6, 1.4, 0.5) m=1.24227 std=0.689541       ] [05m:19s<04m:35s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.6, 1.4, 0.6) m=1.25945 std=0.711656       ] [05m:19s<04m:35s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.6, 1.4, 0.7) m=1.26792 std=0.736485       ] [05m:19s<04m:35s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.6, 2.05, 0.7) m=1.20424 std=0.736236      ] [05m:19s<04m:35s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.5, 0.6) m=1.20901 std=0.7865         ] [05m:19s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.5, 0.7) m=1.21552 std=0.804463       ] [05m:19s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.8, 0.5) m=1.20372 std=0.73826        ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.8, 0.6) m=1.21731 std=0.756439       ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.8, 0.7) m=1.22429 std=0.775797       ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.85, 0.6) m=1.20646 std=0.727297      ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 0.85, 0.7) m=1.21282 std=0.752763      ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 1.4, 0.4) m=1.23375 std=0.711285       ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 1.4, 0.5) m=1.262 std=0.728491         ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 1.4, 0.6) m=1.28186 std=0.748646       ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 1.4, 0.7) m=1.29254 std=0.771486       ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 2.05, 0.5) m=1.20549 std=0.754843      ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 2.05, 0.6) m=1.21927 std=0.773254      ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.65, 0.7, 2.05, 0.7) m=1.22636 std=0.793157      ] [05m:20s<04m:34s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.8, 0.7, 0.35, 0.5) m=1.20028 std=0.720313       ] [05m:48s<04m:06s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(1.8, 0.7, 0.35, 0.6) m=1.20082 std=0.715097       ] [05m:48s<04m:06s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.25, 0.5, 0.25, 0.7) m=1.20666 std=0.657876      ] [07m:11s<02m:44s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.25, 0.6, 0.25, 0.6) m=1.20515 std=0.63879       ] [07m:12s<02m:43s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.25, 0.6, 0.25, 0.7) m=1.21892 std=0.659179      ] [07m:12s<02m:43s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.25, 0.7, 0.25, 0.6) m=1.20862 std=0.658873      ] [07m:14s<02m:41s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.25, 0.7, 0.25, 0.7) m=1.22345 std=0.673932      ] [07m:14s<02m:41s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.5, 0.3, 0.6) m=1.20176 std=0.732425        ] [08m:00s<01m:57s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.5, 0.3, 0.7) m=1.20282 std=0.743403        ] [08m:00s<01m:57s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 0.3, 0.4) m=1.20712 std=0.728205        ] [08m:02s<01m:56s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 0.3, 0.5) m=1.22285 std=0.740456        ] [08m:02s<01m:56s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 0.3, 0.6) m=1.23223 std=0.751652        ] [08m:02s<01m:56s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 0.3, 0.7) m=1.23491 std=0.76248         ] [08m:02s<01m:56s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 1.05, 0.5) m=1.21361 std=0.715643       ] [08m:02s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 1.05, 0.6) m=1.222 std=0.732047         ] [08m:02s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 1.05, 0.7) m=1.22409 std=0.750495       ] [08m:02s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 1.55, 0.5) m=1.20143 std=0.678285       ] [08m:02s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 1.55, 0.6) m=1.20852 std=0.684767       ] [08m:02s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.6, 1.55, 0.7) m=1.20983 std=0.69564        ] [08m:02s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 0.3, 0.3) m=1.20187 std=0.738367        ] [08m:03s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 0.3, 0.4) m=1.22631 std=0.748232        ] [08m:03s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 0.3, 0.5) m=1.24455 std=0.758426        ] [08m:03s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 0.3, 0.6) m=1.25586 std=0.770143        ] [08m:03s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 0.3, 0.7) m=1.25978 std=0.784244        ] [08m:03s<01m:55s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.05, 0.4) m=1.2128 std=0.742194        ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.05, 0.5) m=1.2287 std=0.75507         ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.05, 0.6) m=1.23831 std=0.770027       ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.05, 0.7) m=1.24123 std=0.787226       ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.55, 0.4) m=1.21326 std=0.721561       ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.55, 0.5) m=1.22924 std=0.725548       ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.55, 0.6) m=1.2389 std=0.732482        ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.5, 0.7, 1.55, 0.7) m=1.24186 std=0.743778       ] [08m:03s<01m:54s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.6, 0.7, 0.3, 0.4) m=1.20565 std=0.654758        ] [08m:23s<01m:36s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.6, 0.7, 0.3, 0.5) m=1.21052 std=0.675124        ] [08m:23s<01m:36s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.6, 0.7, 0.3, 0.6) m=1.20857 std=0.700061        ] [08m:23s<01m:36s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.75, 0.6, 1.9, 0.7) m=1.20289 std=0.684776=>     ] [08m:53s<01m:08s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.75, 0.7, 1.45, 0.7) m=1.20884 std=0.803466>     ] [08m:54s<01m:07s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.75, 0.7, 1.9, 0.6) m=1.22173 std=0.702966=>     ] [08m:54s<01m:07s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.75, 0.7, 1.9, 0.7) m=1.24354 std=0.719761=>     ] [08m:54s<01m:07s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.9, 0.4, 0.15, 0.7) m=1.20056 std=0.727066===>   ] [09m:20s<00m:43s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.9, 0.7, 0.95, 0.7) m=1.20829 std=0.767989===>   ] [09m:24s<00m:39s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.95, 0.7, 2.6, 0.6) m=1.20516 std=0.783387====>  ] [09m:34s<00m:28s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(2.95, 0.7, 2.6, 0.7) m=1.21541 std=0.80946=====>  ] [09m:34s<00m:28s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3, 0.6, 2.2, 0.6) m=1.21207 std=0.745026========> ] [09m:43s<00m:20s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3, 0.6, 2.2, 0.7) m=1.23031 std=0.773396========> ] [09m:43s<00m:20s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3, 0.7, 2.2, 0.5) m=1.21139 std=0.765012========> ] [09m:44s<00m:19s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3, 0.7, 2.2, 0.6) m=1.24102 std=0.79716=========> ] [09m:44s<00m:19s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3, 0.7, 2.2, 0.7) m=1.26175 std=0.833805========> ] [09m:44s<00m:19s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.6, 2.8, 0.4) m=1.20894 std=0.756895======>] [09m:54s<00m:10s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.6, 2.8, 0.5) m=1.2131 std=0.77812========>] [09m:54s<00m:10s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.6, 2.8, 0.6) m=1.20997 std=0.796467======>] [09m:54s<00m:10s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.7, 2.8, 0.2) m=1.20596 std=0.745482======>] [09m:55s<00m:09s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.7, 2.8, 0.3) m=1.23324 std=0.772544======>] [09m:55s<00m:09s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.7, 2.8, 0.4) m=1.25328 std=0.796601======>] [09m:55s<00m:09s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.7, 2.8, 0.5) m=1.26528 std=0.817729======>] [09m:55s<00m:09s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.7, 2.8, 0.6) m=1.26877 std=0.83619=======>] [09m:55s<00m:09s] Running Angle Search Experiment FULL BRUTEFORCE p=2
(3.05, 0.7, 2.8, 0.7) m=1.2636 std=0.85229========>] [09m:55s<00m:09s] Running Angle Search Experiment FULL BRUTEFORCE p=2
[==================================================] [10m:05s<00m:00s] Running Angle Search Experiment FULL BRUTEFORCE p=2
*/
class AngleResultsExperiment : AngleExperimentBase{

public:

	int q = 97;

	int m_start = 4;
	int m_end = 20; //10;

	int max_num_instances = 10;//3000;

	//const std::vector<double> angles{0.4,0.48,5.56,0.28}; //work interesting

	//WORKS GREAT
	const std::vector<double> angles{2.2287, 2.13607 ,2.23057 ,2.17114, 0.736168, -0.775112, -1.17933, -2.84297, -1.32729 ,-0.37258 ,1.90813 ,-0.247229};
/*
 *  q=97
   Averages:
 n \ m                   4                   5                   6                   7                   8                   9                  10
  1              0.0641671           0.0333538           0.0157872          0.00804493          0.00401294          0.00211669          0.00100611
  2              0.0684507           0.0324452           0.0162824          0.00841373          0.00680977          0.00600239          0.00476806
  3               0.106731           0.0515169           0.0207424           0.0126368            0.005107          0.00519504           0.0044405
  4                      x            0.102938           0.0451122           0.0091706          0.00585934          0.00539514          0.00442476
  5                      x                   x           0.0859569           0.0323786          0.00566087          0.00542296          0.00443021
  6                      x                   x                   x           0.0508255           0.0153818          0.00288475          0.00273762
  7                      x                   x                   x                   x           0.0364387          0.00482137          0.00144237
  8                      x                   x                   x                   x                   x           0.0213782          0.00144802
  9                      x                   x                   x                   x                   x                   x           0.0122817
 *
 */


	//const std::vector<double> angles{0.583336, 2.16313 ,2.24903 ,2.18185};

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

	int max_num_instances = 1000;

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
	void run_p6_test();
	void run_cobyla();

};

class AlphaMinimizationExperiment{

	public:
		int loglevel = 1;

	AlphaMinimizationExperiment(int loglevel, FastVQA::QAOAOptions*, MapOptions*);

	void run();

	private:

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

void experiment_runner(ExperimentSetup*, std::string experiment_name, int loglevel);

#endif /* SRC_EXPERIMENTRUNNER_H_ */
