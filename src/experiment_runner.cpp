#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include <fstream>
#include <numeric>
#include <boost/tuple/tuple.hpp>
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"
#include <iomanip>
#include <cmath>

#define s(X) std::to_string(X)

const double pi = 3.141592654;

std::string nlopt_res_to_str(int result){
  switch(result)
  {
    case -1: return "FAILURE";
    case -2: return "INVALID_ARGS";
    case -3: return "OUT_OF_MEMORY";
    case -4: return "ROUNDOFF_LIMITED";
    case -5: return "FORCED_STOP";
    case 1: return "SUCCESS";
    case 2: return "STOPVAL_REACHED";
    case 3: return "FTOL_REACHED";
    case 4: return "XTOL_REACHED";
    case 5: return "MAXEVAL_REACHED";
    case 6: return "MAXTIME_REACHED";
    default: return NULL;
  }
}

CmQaoaExperiment::CmQaoaExperiment(FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database, int loglevel){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->database = database;
}

void CmQaoaExperiment::run(){

	int nbQubits_acc = -1;
	this->mapOptions->penalty = 0;

	for(int m = this->m_start; m <= this->m_end; ++m){

		int qs = ceil(log2(m));
		mapOptions->num_qbits_per_x = qs;

		GeneratorParam param(this->q, this->n, m, true, 97, this->num_instances); //q, n, m, shuffle, seed, cutoff
		std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

		for(int i = 0; i < this->num_instances; ++i){
			CmQaoaExperiment::Instance instance;
			Lattice l(gramian_wrappers[i].hamiltonian, gramian_wrappers[i].name);
			instance.q = this->q;
			instance.m = m;
			instance.n = this-> n;
			instance.h = l.getHamiltonian(mapOptions);

			if(nbQubits_acc < 0)
				nbQubits_acc = instance.h.nbQubits;
			else if(nbQubits_acc != instance.h.nbQubits){
				logw("Instances with different number of qubits found");
				//Need to destroy qureg, because next time experiments with different number of qubits will be run
				qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
				qaoaOptions->accelerator->finalize();
			}

			qaoaOptions->accelerator->initialize(&instance.h);
			qaoaOptions->accelerator->options.createQuregAtEachInilization = false;

			instance.solutions = qaoaOptions->accelerator->getSolutions();
			if(instance.solutions.size() != 1)
				throw_runtime_error("CmQaoaExperiment: Unimplemented, more than 1 solution marked");

			for(auto &sol: instance.solutions){
				if(sol.value != 0)
					throw_runtime_error("CmQaoaExperiment: Something else than 0 state marked as a solution");
			}
			long long int zero_index = instance.solutions[0].index;

			//instance.solutions = qaoaOptions->accelerator->getSolutions();

			for(int p = this->p_start; p <= this->p_end; ++p){
				this->qaoaOptions->p = p;

				logi("Running m="+s(m)+" p="+s(p)+" qs="+s(qs)+" index="+s(i)+"     Total Qubits: "+s(m*qs), this->loglevel);

			}

		}

	}

}

CmQaoaExperiment::Cost CmQaoaExperiment::_cost_fn(CmQaoaExperiment::Instance*, bool use_database){
	CmQaoaExperiment::Cost cost;



	return cost;
}


AngleResultsExperiment::AngleResultsExperiment(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->database = database;

	if(this->qaoaOptions->p != 3)
		throw_runtime_error("Angleres is only for p=5!");

}

std::vector<AngleResultsExperiment::Instance> AngleResultsExperiment::_generate_dataset(int n, int m){

	std::vector<AngleResultsExperiment::Instance> dataset;

	GeneratorParam param(q, n, m, true, 97, this->max_num_instances); //q, n, m, shuffle, seed, cutoff
	std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

	int nbQubits = -1;

	//logw("Saving eigensace which is not needed and very costly");
	//logw("!!!Random guess is now after the penalization!!!");

	long long int num_instances = pow(q,(m-n));
	if(max_num_instances < num_instances)
		num_instances = max_num_instances;

	for(int i = 0; i < num_instances; ++i){

		AngleResultsExperiment:Instance instance;

		Lattice l(gramian_wrappers[i].hamiltonian, gramian_wrappers[i].name);
		mapOptions->penalty = l.getSquaredLengthOfFirstBasisVector(); //penalty set to length of first vector squared

		instance.h = l.getHamiltonian(mapOptions);
		if(nbQubits < 0)
			nbQubits = instance.h.nbQubits;
		else if(nbQubits != instance.h.nbQubits)
			loge("Instances with different number of qubits found");

		qaoaOptions->accelerator->initialize(&instance.h);
		qaoaOptions->accelerator->options.createQuregAtEachInilization = false;
		instance.solutions = qaoaOptions->accelerator->getSolutions();

		//instance.eigenspace = qaoaOptions->accelerator->getEigenspace();//delete

		qreal min_energy = instance.solutions[0].value;
		int num_sols_with_min_energy = 0;
		for(const auto &sol: instance.solutions){
			if(sol.value < min_energy){
				min_energy = sol.value;
				num_sols_with_min_energy = 1;
			}else if(sol.value == min_energy)
				num_sols_with_min_energy++;
		}
		instance.min_energy = min_energy;

		//THIS CHOOSES WHICH RANDOM GUESS IS BEING USED
		instance.random_guess = (qreal)(1./pow(2, nbQubits)) * num_sols_with_min_energy;
		//instance.random_guess = l.get_random_guess_one_vect() * instance.h.custom_solutions.size();

		//std::cerr<<nbQubits<<" "<<instance.solutions.size()<<" "<<instance.random_guess<<std::endl;

		instance.volume	= l.getVolume();
		instance.sv1Squared= l.getSquaredLengthOfFirstBasisVector();
		instance.q = q;
		instance.m = m;
		instance.n = n;

		dataset.push_back(instance);
	}

	//Need to destroy qureg, because next time experiments with different number of qubits will be run
	qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
	qaoaOptions->accelerator->finalize();

	//logi("Experiment dataset generated", this->loglevel);

	return dataset;
}


void AngleResultsExperiment::run(){

	std::map<std::pair<int, int>, double> stdev_map;

	const int colWidth = 20;
	std::cout<<" q="<<this->q<<std::endl;
	std::cout<<"   Averages:"<<std::endl;
	std::cout << " n \\ m";
	for(int m = this->m_start; m <= this->m_end; ++m)
		std::cout << std::setw(colWidth) << std::internal << m;
	std::cout<<std::endl;

	for(int n = 1; n < m_end; ++n){
		std::cout << std::setw(3) << std::internal << n << "   ";
		for(int m = this->m_start; m <= m_end; ++m){
			if(n >= m){
				std::cout << std::setw(colWidth) << std::internal << "x";
				continue;
			}
			std::vector<AngleResultsExperiment::Instance> dataset = this->_generate_dataset(n, m);

			Cost cost = this->_cost_fn(&dataset, &this->angles[0], true);

			double mean = cost.mean;
			double stdev = cost.stdev;

			stdev_map.emplace(std::pair<int, int>(n,m), stdev);
			std::cout << std::setw(colWidth) << std::internal << mean << "/" << stdev << std::flush;
		}
		std::cout<<std::endl;
	}

	std::cout<<std::endl<<std::endl;
	std::cout<<"   Standard deviations:"<<std::endl;
	std::cout << " n \\ m";
		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout<<std::endl;

		for(int n = 1; n < m_end; ++n){
			std::cout << std::setw(3) << std::internal << n << "   ";
			for(int m = this->m_start; m <= m_end; ++m){
				if(n >= m){
					std::cout << std::setw(colWidth) << std::internal << "x";
					continue;
				}
				std::cout << std::setw(colWidth) << std::internal << stdev_map[std::pair<int, int>(n,m)];
			}
			std::cout<<std::endl;
		}

}

AngleSearchExperiment::AngleSearchExperiment(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;

	this->num_instances = pow(q,(m-n));
	if(max_num_instances < this->num_instances)
		this->num_instances = max_num_instances;
	this->num_test_instances = this->num_instances * this->test_ratio;
	this->num_train_instances = this->num_instances - this->num_test_instances;
	logi("q="+std::to_string(q)+" "+"n="+std::to_string(n)+" "+"m="+std::to_string(m)+" : "+std::to_string(this->num_instances) + " instances", this->loglevel);
	logi("Test ratio="+std::to_string(this->test_ratio)+" "+"train_size="+std::to_string(this->num_train_instances)+" "+"test_size="+std::to_string(this->num_test_instances), this->loglevel);

	//std::mt19937 gen(0); //rd() instead of 0 - seed
	//std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);

	logw("All qaoa instances started at the same initial point", loglevel);
	this->num_params = this->qaoaOptions->p*2;

	/*if(this->qaoaOptions->initial_params.size() > 0)
		this->qaoaOptions->initial_params.clear();
	for(int j = 0; j < this->num_params; ++j){
		this->qaoaOptions->initial_params.push_back(dis(gen));
	}*/

	this->_generate_dataset(mapOptions);
};

void AngleSearchExperiment::_generate_dataset(MapOptions* mapOptions){

	logi("Generating experiment dataset", this->loglevel);

	GeneratorParam param(q, n, m, true, 97, this->max_num_instances); //q, n, m, shuffle, seed, cutoff
	std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

	nbQubits = -1;

	logw("Saving eigensace which is not needed and very costly");
	logw("!!!Random guess is now after the penalization!!!");

	for(int i = 0; i < num_instances; ++i){

		Instance instance;

		Lattice l(gramian_wrappers[i].hamiltonian, gramian_wrappers[i].name);
		mapOptions->penalty = l.getSquaredLengthOfFirstBasisVector(); //penalty set to length of first vector squared
		instance.h = l.getHamiltonian(mapOptions);
		if(nbQubits < 0)
			nbQubits = instance.h.nbQubits;
		else if(nbQubits != instance.h.nbQubits)
			loge("Instances with different number of qubits found");
		qaoaOptions->accelerator->initialize(&instance.h);
		qaoaOptions->accelerator->options.createQuregAtEachInilization = false;
		instance.solutions = qaoaOptions->accelerator->getSolutions();

		instance.eigenspace = qaoaOptions->accelerator->getEigenspace();//delete

		qreal min_energy = instance.solutions[0].value;
		int num_sols_with_min_energy = 0;
		for(const auto &sol: instance.solutions){
			if(sol.value < min_energy){
				min_energy = sol.value;
				num_sols_with_min_energy = 1;
			}else if(sol.value == min_energy)
				num_sols_with_min_energy++;
		}
		instance.min_energy = min_energy;


		//THIS CHOOSES WHICH RANDOM GUESS IS BEING USED
		instance.random_guess = (qreal)(1./pow(2, nbQubits)) * num_sols_with_min_energy;
		//instance.random_guess = l.get_random_guess_one_vect() * instance.h.custom_solutions.size();

		//std::cerr<<nbQubits<<" "<<instance.solutions.size()<<" "<<instance.random_guess<<std::endl;

		if(i < num_train_instances)
			this->train_set.push_back(instance);
		else
			this->test_set.push_back(instance);
	}

	logi("Experiment dataset generated", this->loglevel);

}
// Function for calculating median
/*double median(std::vector<double> v, int n)
{
    // Sort the vector
    std::sort(v.begin(), v.end());

    // Check if the number of elements is odd
    if (n % 2 != 0)
        return (double)v[n / 2];

    // If the number of elements is even, return the average
    // of the two middle elements
    return (double)(v[(n - 1) / 2] + v[n / 2]) / 2.0;
}*/
AngleExperimentBase::Cost AngleExperimentBase::_cost_fn(std::vector<Instance>* dataset, const double *angles, bool use_database){

	std::vector<int> num_sols;

	int i = 0;

	double mean;
	double stdev;
	FastVQA::Qaoa qaoa_instance;

	std::vector<double> gs_overlaps;

	for(auto &instance: (*dataset)){

		if(instance.h.nbQubits != this->qaoaOptions->accelerator->getNumQubitsInQureg()){
			this->qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
			this->qaoaOptions->accelerator->finalize();
		}

		FastVQA::ExperimentBuffer buffer;
		buffer.storeQuregPtr = true;

		double ground_state_overlap = 0;
		if(use_database){
			Database::DatasetRow output_row;
			this->database->getOrCalculate_qary_with_fixed_angles(&buffer, angles, 6, &instance.h, &output_row, this->qaoaOptions, &qaoa_instance);
			for(auto &sol: instance.solutions){
				long long int index = sol.index;
				ground_state_overlap+=output_row.finalStateVectorMap[index].second;
			}
		}else{
			qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
			for(auto &sol: instance.solutions){
				long long int index = sol.index;
				ground_state_overlap+=buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
			}
		}

		num_sols.push_back(instance.h.custom_solutions.size());
		qreal improvement_ratio = ground_state_overlap / instance.random_guess;
		/*if(improvement_ratio < 1){
			std::cerr<<ground_state_overlap <<" "<< instance.random_guess << std::endl;
			logw("Improvement ratio < 1", loglevel);
		}*/
		gs_overlaps.push_back(improvement_ratio);
		//std::cerr<<ground_state_overlap<<" "<<instance.random_guess<<"\n";
		i++;

		this->qaoaOptions->accelerator->options.createQuregAtEachInilization = false;

	}

	double sum = std::accumulate(gs_overlaps.begin(), gs_overlaps.end(), 0.0);
	mean = sum / gs_overlaps.size();

	double sq_sum = std::inner_product(gs_overlaps.begin(), gs_overlaps.end(), gs_overlaps.begin(), 0.0);
	stdev = std::sqrt(sq_sum / gs_overlaps.size() - mean * mean);

	//mean = median(gs_overlaps, gs_overlaps.size());

	this->qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
	this->qaoaOptions->accelerator->finalize();


	double sum_mean_num_of_sols = std::accumulate(num_sols.begin(), num_sols.end(), 0.0);
	double mean_num_of_sols = sum_mean_num_of_sols / num_sols.size();

	return Cost(mean, stdev, mean_num_of_sols);

}

void AngleSearchExperiment::run(){

	//run_cobyla_p2();
	//run_p2_test();
	//run_p2();
	//return;

	if(this->qaoaOptions->p == 1)
		run_p1();
	else if(this->qaoaOptions->p == 2)
		run_p2_test();//run_p2_full_bruteforce();
		//run_p2();
	else if(this->qaoaOptions->p == 3)
		run_p3_full_bruteforce();
	else
		throw_runtime_error("Not implemented");

}

void AngleSearchExperiment::run_p3_full_bruteforce(){

	if(this->num_params != 6)
		throw_runtime_error("Unimplemented for other depth than 3");

	double *angles = (double*) malloc(this->num_params * sizeof(double));
		for(int j = 0; j < this->num_params; ++j){
			angles[j]=0;//dis(gen));
		}

	std::vector<std::tuple<double, double>> first_round_angles;

	double mean_threshold = 10;
	double stdev_threshold = 20;//0.019;

	double incr_beta = 0.1;///0.12;
	double incr_gamma = 0.04;//0.04;///0.01;

	double beta_min = incr_beta;//pi/16;
	double beta_max = pi/4;//pi;//;pi/8;

	double gamma_min = incr_gamma;//0;
	double gamma_max = /*2**/pi/4;//2*pi;

	double range_beta = beta_max - beta_min;
	double range_gamma = gamma_max - gamma_min;


	int axis_range_beta = ceil(range_beta / incr_beta);
	int axis_range_gamma = ceil(range_gamma / incr_gamma);
	int num_iters = axis_range_beta *  axis_range_beta * axis_range_beta * axis_range_gamma * axis_range_gamma * axis_range_gamma;

	int x = 0;
	int y = 0;

	std::vector<std::vector<double>> final_plot_means(axis_range_gamma);
	std::vector<std::vector<double>> final_plot_stdevs(axis_range_gamma);

	ProgressBar bar{
		option::BarWidth{50},
		option::MaxProgress{num_iters},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Running Angle Search Experiment FULL BRUTEFORCE p=3"},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::ForegroundColor{Color::yellow},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	double debug_highest_mean=6.8;
	double lgamma1, lbeta1, lgamma2, lbeta2, lgamma3, lbeta3;

	double i = 0;
	AngleSearchExperiment::Cost cost;
	for(double gamma1 = gamma_min; gamma1 < gamma_max; gamma1+=incr_gamma){
		for(double beta1 = beta_min; beta1 < beta_max; beta1+=incr_beta){
			for(double gamma2 = gamma_min; gamma2 < gamma_max; gamma2+=incr_gamma){
				for(double beta2 = beta_min; beta2 < beta_max; beta2+=incr_beta){
					for(double gamma3 = gamma_min; gamma3 < gamma_max; gamma3+=incr_gamma){
						for(double beta3 = beta_min; beta3 < beta_max; beta3+=incr_beta){
							bar.tick();
							angles[0] = gamma1;
							angles[1] = beta1;
							angles[2] = gamma2;
							angles[3] = beta2;
							angles[4] = gamma3;
							angles[5] = beta3;
							cost = this->_cost_fn(&this->train_set, angles);
							//final_plot_means[x].push_back(cost.mean);
							//final_plot_stdevs[x].push_back(cost.stdev);

							if(cost.mean > debug_highest_mean){
								debug_highest_mean = cost.mean;
								lgamma1=gamma1;
								lbeta1=beta1;
								lgamma2=gamma2;
								lbeta2=beta2;
								lgamma3=gamma3;
								lbeta3=beta3;

								std::cerr<<"("<<gamma1<<", "<<beta1<<", "<<gamma2<<", "<<beta2<<", "<<gamma3<<", "<<beta3<<") m="<<cost.mean<<" std="<<cost.stdev<<"\n";

							}

						}
					}
				}
			}
		}
		//x++;
	}

	std::cerr<<debug_highest_mean <<" "<<lgamma1<<" "<<lbeta1<<" "<<lgamma2<<" "<<lbeta2<<" "<<lgamma3<<" "<<lbeta3<<std::endl;

	//logi("First round size: " + std::to_string(first_round_angles.size()), loglevel);

}

void AngleSearchExperiment::run_p2_full_bruteforce(){

	if(this->num_params != 4)
		throw_runtime_error("Unimplemented for other depth than 2");

	double *angles = (double*) malloc(this->num_params * sizeof(double));
		for(int j = 0; j < this->num_params; ++j){
			angles[j]=0;//dis(gen));
		}

	std::vector<std::tuple<double, double>> first_round_angles;

	double mean_threshold = 1.34;
	double stdev_threshold = 20;//0.019;

	double beta_min = 0;//pi/16;
	double beta_max = pi/4;//pi;//;pi/8;

	double gamma_min = 0;//0;
	double gamma_max = /*2**/pi;//2*pi;

	double range_beta = beta_max - beta_min;
	double range_gamma = gamma_max - gamma_min;
	double incr_beta = 0.1;///0.12;
	double incr_gamma = 0.05;//0.04;///0.01;

	int axis_range_beta = ceil(range_beta / incr_beta);
	int axis_range_gamma = ceil(range_gamma / incr_gamma);
	int num_iters = axis_range_beta *  axis_range_beta * axis_range_gamma * axis_range_gamma;

	int x = 0;
	int y = 0;

	std::vector<std::vector<double>> final_plot_means(axis_range_gamma);
	std::vector<std::vector<double>> final_plot_stdevs(axis_range_gamma);

	ProgressBar bar{
		option::BarWidth{50},
		option::MaxProgress{num_iters},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Running Angle Search Experiment FULL BRUTEFORCE p=2"},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::ForegroundColor{Color::yellow},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	double debug_lowest_mean=100000;
	double lgamma1, lbeta1, lgamma2, lbeta2;

	double i = 0;
	AngleSearchExperiment::Cost cost;
	for(double gamma1 = gamma_min; gamma1 < gamma_max; gamma1+=incr_gamma){
		for(double beta1 = beta_min; beta1 < beta_max; beta1+=incr_beta){
			for(double gamma2 = gamma_min; gamma2 < gamma_max; gamma2+=incr_gamma){
				for(double beta2 = beta_min; beta2 < beta_max; beta2+=incr_beta){
					bar.tick();
					angles[0] = gamma1;
					angles[1] = beta1;
					angles[2] = gamma2;
					angles[3] = beta2;
					cost = this->_cost_fn(&this->train_set, angles);
					//final_plot_means[x].push_back(cost.mean);
					//final_plot_stdevs[x].push_back(cost.stdev);

					if(cost.mean < debug_lowest_mean){
						debug_lowest_mean = cost.mean;
						lgamma1=gamma1;
						lbeta1=beta1;
						lgamma2=gamma2;
						lbeta2=beta2;
					}

					if(cost.mean >= mean_threshold && cost.stdev <= stdev_threshold){
						//first_round_angles.push_back(std::tuple<double, double>(gamma, beta));

						std::cerr<<"("<<gamma1<<", "<<beta1<<", "<<gamma2<<", "<<beta2<<") m="<<cost.mean<<" std="<<cost.stdev<<"\n";

					}
				}
			}
		}
		//x++;
	}

	std::cerr<<debug_lowest_mean <<" "<<lgamma1<<" "<<lbeta1<<" "<<lgamma2<<" "<<lbeta2<<std::endl;

	//logi("First round size: " + std::to_string(first_round_angles.size()), loglevel);

}


void AngleSearchExperiment::run_cobyla_p2(){

	if(this->num_params != 4)
		throw_runtime_error("Unimplemented for other depth than 2");

	this->qaoaOptions->ftol = 1e-16;
	this->qaoaOptions->max_iters = 10000;

	ProgressBar bar{
			option::BarWidth{50},
			option::MaxProgress{this->qaoaOptions->max_iters},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{"Running Angle Search Experiment p=2 with COBYLA"},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::ForegroundColor{Color::yellow},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		};

	unsigned int iteration_i = 0;

	FastVQA::OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
		iteration_i++;
		bar.tick();
		std::vector<double> angles(x);
		//std::cerr<<this->_cost_fn(&this->train_set, &angles[0]).mean<<std::endl;
		return -this->_cost_fn(&this->train_set, &angles[0]).mean;
	}, this->num_params);

	std::vector<double> initial_params;
	std::mt19937 gen(0); //rd() instead of 0 - seed
	std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
	for(int i = 0; i < this->num_params/2; ++i){
		double param1 = dis(gen);
		double param2 = dis(gen);
		std::cerr<<param1<<" "<<param2<<std::endl;
		initial_params.push_back(param1);
		initial_params.push_back(param2);
	}

	std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
	std::vector<double> upperBounds(initial_params.size(), 3.141592654);

	logd("QAOA starting optimization", this->loglevel);
	FastVQA::OptResult result = this->qaoaOptions->optimizer->optimize(f, initial_params, this->qaoaOptions->ftol, this->qaoaOptions->max_iters, lowerBounds, upperBounds);
	logd("QAOA finishing optimization", this->loglevel);

	std::cerr<<"e: "<< -result.first.first <<"\n";
	std::cerr<<"num_iters: "<<iteration_i<<std::endl;
	for(auto &a: result.first.second){
		std::cerr<<a<<" ";
	}
	std::cerr<<"\n"<<nlopt_res_to_str(result.second)<<std::endl;

}

void AngleSearchExperiment::run_p1(){

	if(this->num_params != 2)
		throw_runtime_error("Depth must be 1");

	double *angles = (double*) malloc(this->num_params * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

	/*if(this->num_params != 2)
		throw_runtime_error("unimplemented");
	*/

	std::vector<std::tuple<double, double>> first_round_angles;

	double mean_threshold = 9.3;//1.2;
	double stdev_threshold = 2000000;//0.019;

	double beta_min = 0;//pi/16;
	double beta_max = pi/2;//pi;//;pi/8;

	double gamma_min = 0;//0;
	double gamma_max = pi;//2*pi;

	double range_beta = beta_max - beta_min;
	double range_gamma = gamma_max - gamma_min;
	double incr_beta = 0.04;///0.12;
	double incr_gamma = 0.04;//0.04;///0.01;

	int axis_range_beta = ceil(range_beta / incr_beta);
	int axis_range_gamma = ceil(range_gamma / incr_gamma);
	int num_iters = axis_range_beta * axis_range_gamma;

	int x = 0;
	int y = 0;

	std::vector<std::vector<double>> final_plot_means(axis_range_gamma);
	std::vector<std::vector<double>> final_plot_stdevs(axis_range_gamma);

	ProgressBar bar{
		option::BarWidth{50},
		option::MaxProgress{num_iters},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Running Angle Search Experiment p=1"},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::ForegroundColor{Color::yellow},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	double debug_lowest_mean=100000;
	double lgamma, lbeta;

	double i = 0;
	AngleSearchExperiment::Cost cost;
	for(double gamma = gamma_min; gamma < gamma_max; gamma+=incr_gamma){
		for(double beta = beta_min; beta < beta_max; beta+=incr_beta){
			bar.tick();
			angles[0] = gamma;
			angles[1] = beta;
			cost = this->_cost_fn(&this->train_set, angles);
			final_plot_means[x].push_back(cost.mean);
			final_plot_stdevs[x].push_back(cost.stdev);

			//loge("ahoj");

			if(cost.mean < debug_lowest_mean){
				debug_lowest_mean = cost.mean;
				lgamma=gamma;
				lbeta=beta;
			}

			if(cost.mean >= mean_threshold && cost.stdev <= stdev_threshold){
				first_round_angles.push_back(std::tuple<double, double>(gamma, beta));

				std::cerr<<"("<<gamma<<" "<<beta<<") m="<<cost.mean<<" std="<<cost.stdev<<"\n";

			}

		}
		x++;
	}

	std::cerr<<debug_lowest_mean <<" "<<lgamma<<" "<<lbeta<<std::endl;

	logi("First round size: " + std::to_string(first_round_angles.size()), loglevel);

	if(!false){
		Gnuplot gp;
		gp << "set multiplot layout 1,2 rowsfirst\n";
		gp << "unset key\n";
		gp << "set pm3d\n";
		gp << "set hidden3d\n";
		gp << "set palette rgb 33,13,10\n";
		gp << "set view map\n";
		gp << "set xrange [ 0 : " << axis_range_beta-1 <<" ] \n";
		gp << "set yrange [ 0 : " << axis_range_gamma-1 <<" ] \n";
		gp << "set xlabel 'Beta [0,pi/4]' \n";
		gp << "set ylabel 'Gamma [-pi,pi]' \n";
		//gp << "set title 'range="<< axis_range <<" nbQs="<< this->nbQubits <<" rg="<< 1./(pow(2,this->nbQubits)) <<" mean_nbSols="<<cost.mean_num_of_sols<<" rdmOverlap="<<1./(pow(2,this->nbQubits))*cost.mean_num_of_sols<<"'\n";
		gp << "set title 'LHS: means RHS: stdev'\n";

		gp << "splot '-'\n";
		gp.send2d(final_plot_means);

		//gp << "set palette model HSV\n";
		gp << "set palette rgb 13,10,33\n";

		gp << "splot '-'\n";
		gp.send2d(final_plot_stdevs);

		gp.flush();
	}

	/*double mx, my;int  mb;
	 * while(true){
		gp.getMouse(mx, my, mb, "");
		printf("You pressed mouse button %d at x=%f y=%f\n", mb, mx, my);
		//printf("%f", final_plot[(int)mx][(int)my]);
	}*/

	free(angles);

};

void AngleSearchExperiment::run_p2(){

	if(this->num_params != 4)
		throw_runtime_error("Depth must be 2");

/*
(0.4 0.4) m=9.30907 std=8.70639                    ] [00m:07s<00m:50s] Running Angle Search Experiment p=1
(0.4 0.44) m=9.40437 std=8.96823                   ] [00m:07s<00m:50s] Running Angle Search Experiment p=1
(0.4 0.48) m=9.44947 std=9.1068                    ] [00m:07s<00m:50s] Running Angle Search Experiment p=1
(0.4 0.52) m=9.44697 std=9.14503                   ] [00m:07s<00m:50s] Running Angle Search Experiment p=1
(0.4 0.56) m=9.4025 std=9.11298                    ] [00m:07s<00m:50s] Running Angle Search Experiment p=1
(0.4 0.6) m=9.32416 std=9.04324                    ] [00m:07s<00m:50s] Running Angle Search Experiment p=1
[==================================================] [01m:06s<00m:00s] Running Angle Search Experiment
 *
 */

	std::vector<std::tuple<double, double>> first_round_angles;
	first_round_angles.push_back(std::tuple<double, double>(0.4, 0.4));
	first_round_angles.push_back(std::tuple<double, double>(0.4, 0.48));
	first_round_angles.push_back(std::tuple<double, double>(0.4, 0.6));


	double *angles = (double*) malloc(4 * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

		double mean_threshold = 9.5;
		double stdev_threshold = 1000;
		double beta_min = 0;
		double beta_max = pi/2;

		double gamma_min = 0;
		double gamma_max = 2*pi;

		double range_beta = beta_max - beta_min;
		double range_gamma = gamma_max - gamma_min;
		double incr_beta = 0.04;//0.12;
		double incr_gamma = 0.04;//0.01;

		int axis_range_beta = ceil(range_beta / incr_beta);
		int axis_range_gamma = ceil(range_gamma / incr_gamma);
		int num_iters = axis_range_beta * axis_range_gamma;

		int x = 0;
		int y = 0;

		std::vector<std::vector<double>> final_plot_means(axis_range_gamma);
		std::vector<std::vector<double>> final_plot_stdevs(axis_range_gamma);

		ProgressBar bar{
			option::BarWidth{50},
			option::MaxProgress{num_iters},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{"Running Angle Search Experiment p=2"},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::ForegroundColor{Color::yellow},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		};

		double i = 0;
		AngleSearchExperiment::Cost cost;
		for(double gamma2 = gamma_min; gamma2 < gamma_max; gamma2+=incr_gamma){
			for(double beta2 = beta_min; beta2 < beta_max; beta2+=incr_beta){
				bar.tick();
				double max_mean = -1;

				double max_gamma1=-1;
				double max_beta1=-1;

				for(auto &first_round : first_round_angles){

					/*if(max_mean != -1)
						break;*/
					double gamma1 = std::get<0>(first_round);
					double beta1 = std::get<1>(first_round);
					angles[0] = gamma1;
					angles[1] = beta1;
					angles[2] = gamma2;
					angles[3] = beta2;
					cost = this->_cost_fn(&this->train_set, angles);
					if(max_mean == -1){
						final_plot_means[x].push_back(cost.mean);
						final_plot_stdevs[x].push_back(cost.stdev);
					}
					if(cost.mean > max_mean){
						final_plot_means[x][final_plot_means[x].size()-1]=cost.mean;
						final_plot_stdevs[x][final_plot_stdevs[x].size()-1]=cost.stdev;
						max_gamma1=gamma1;
						max_beta1=beta1;
						max_mean = cost.mean;
					}
				}
				if(max_mean >= mean_threshold && final_plot_stdevs[x][final_plot_stdevs[x].size()-1] <= stdev_threshold){
					//first_round_angles.push_back(std::tuple<double, double>(gamma, beta));
					std::cerr<<"("<<max_gamma1<<","<<max_beta1<<","<<gamma2<<","<<beta2<<") m="<<max_mean<<" std="<<final_plot_stdevs[x][final_plot_stdevs[x].size()-1]<<"\n";
				}

			}
			x++;
		}

		Gnuplot gp;
			gp << "set multiplot layout 1,2 rowsfirst\n";
			gp << "unset key\n";
			gp << "set pm3d\n";
			gp << "set hidden3d\n";
			gp << "set palette rgb 33,13,10\n";
			gp << "set view map\n";
			gp << "set xrange [ 0 : " << axis_range_beta-1 <<" ] \n";
			gp << "set yrange [ 0 : " << axis_range_gamma-1 <<" ] \n";
			gp << "set xlabel 'Beta [0,pi/4]' \n";
			gp << "set ylabel 'Gamma [-pi,pi]' \n";
			//gp << "set title 'range="<< axis_range <<" nbQs="<< this->nbQubits <<" rg="<< 1./(pow(2,this->nbQubits)) <<" mean_nbSols="<<cost.mean_num_of_sols<<" rdmOverlap="<<1./(pow(2,this->nbQubits))*cost.mean_num_of_sols<<"'\n";
			gp << "set title 'LHS: means RHS: stdev'\n";

			gp << "splot '-'\n";
			gp.send2d(final_plot_means);

			//gp << "set palette model HSV\n";
			gp << "set palette rgb 13,10,33\n";

			gp << "splot '-'\n";
			gp.send2d(final_plot_stdevs);

			gp.flush();

			/*double mx, my;int  mb;
			 * while(true){
				gp.getMouse(mx, my, mb, "");
				printf("You pressed mouse button %d at x=%f y=%f\n", mb, mx, my);
				//printf("%f", final_plot[(int)mx][(int)my]);
			}*/


	free(angles);

}

void AngleSearchExperiment::run_p2_test(){

	struct Angl{
		double g1,b1,g2,b2;
		double train_mean;
		double train_std;
		Angl(double g1, double b1, double g2, double b2, double train_mean, double train_std){
			this->g1=g1;this->b1=b1;this->g2=g2;this->b2=b2;this->train_mean=train_mean;this->train_std=train_std;
		}
	};

/*	 */

	std::vector<Angl> second_round_angles;

	/*
	 * (0.4,0.48,0.32,0.04) m=9.50037 std=9.18755         ] [00m:23s<07m:20s] Running Angle Search Experiment p=2
(0.4,0.48,0.32,0.08) m=9.50689 std=9.17675         ] [00m:23s<07m:20s] Running Angle Search Experiment p=2
(0.4,0.48,0.36,0.04) m=9.5999 std=9.31231          ] [00m:26s<07m:17s] Running Angle Search Experiment p=2
(0.4,0.48,0.36,0.08) m=9.68412 std=9.43483         ] [00m:26s<07m:17s] Running Angle Search Experiment p=2
(0.4,0.48,0.36,0.12) m=9.69658 std=9.46659         ] [00m:26s<07m:17s] Running Angle Search Experiment p=2
(0.4,0.48,0.36,0.16) m=9.63625 std=9.40428         ] [00m:26s<07m:17s] Running Angle Search Experiment p=2
(0.4,0.48,0.36,0.2) m=9.50667 std=9.24941          ] [00m:26s<07m:17s] Running Angle Search Experiment p=2
(0.4,0.48,1.84,0.04) m=9.5038 std=9.10692          ] [02m:17s<05m:33s] Running Angle Search Experiment p=2
(0.4,0.48,1.84,0.08) m=9.51616 std=8.96549         ] [02m:17s<05m:33s] Running Angle Search Experiment p=2
(0.4,0.48,1.92,0.08) m=9.51976 std=8.81561         ] [02m:23s<05m:27s] Running Angle Search Experiment p=2
(0.4,0.48,1.92,0.12) m=9.51127 std=8.56281         ] [02m:23s<05m:27s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.08) m=9.54653 std=9.23614          ] [02m:32s<05m:19s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.12) m=9.61485 std=9.25598          ] [02m:32s<05m:19s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.16) m=9.65381 std=9.21634          ] [02m:32s<05m:19s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.2) m=9.66276 std=9.11257           ] [02m:32s<05m:19s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.24) m=9.64084 std=8.94321          ] [02m:32s<05m:19s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.28) m=9.58708 std=8.71132          ] [02m:32s<05m:18s] Running Angle Search Experiment p=2
(0.4,0.6,2.04,0.32) m=9.50059 std=8.42482          ] [02m:32s<05m:18s] Running Angle Search Experiment p=2
(0.4,0.48,2.32,0.04) m=9.55894 std=9.20883         ] [02m:53s<04m:58s] Running Angle Search Experiment p=2
(0.4,0.48,2.32,0.08) m=9.61263 std=9.19993         ] [02m:53s<04m:58s] Running Angle Search Experiment p=2
(0.4,0.48,2.32,0.12) m=9.60944 std=9.08295         ] [02m:53s<04m:58s] Running Angle Search Experiment p=2
(0.4,0.48,2.32,0.16) m=9.55174 std=8.86525         ] [02m:53s<04m:58s] Running Angle Search Experiment p=2
(0.4,0.48,2.72,0.04) m=9.51443 std=9.23339         ] [03m:23s<04m:28s] Running Angle Search Experiment p=2
(0.4,0.48,2.72,0.08) m=9.52894 std=9.26695         ] [03m:23s<04m:28s] Running Angle Search Experiment p=2
(0.4,0.48,2.84,0.04) m=9.51777 std=9.13546         ] [03m:32s<04m:19s] Running Angle Search Experiment p=2
(0.4,0.48,2.84,0.08) m=9.51957 std=9.05822         ] [03m:32s<04m:19s] Running Angle Search Experiment p=2
(0.4,0.48,3.6,0.04) m=9.64408 std=9.20377          ] [04m:29s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.08) m=9.88224 std=9.57883           ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.12) m=10.0794 std=9.80599           ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.16) m=10.2137 std=9.99014           ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.2) m=10.2839 std=10.12              ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.24) m=10.2931 std=10.1862           ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.28) m=10.2478 std=10.1827           ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.32) m=10.1574 std=10.1085           ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.36) m=10.033 std=9.96919            ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.4) m=9.88659 std=9.7783             ] [04m:30s<03m:23s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.44) m=9.72935 std=9.55716           ] [04m:30s<03m:22s] Running Angle Search Experiment p=2
(0.4,0.6,3.6,0.48) m=9.57116 std=9.33294           ] [04m:30s<03m:22s] Running Angle Search Experiment p=2
(0.4,0.48,5.08,0.04) m=9.51092 std=9.0222>         ] [06m:27s<01m:34s] Running Angle Search Experiment p=2
(0.4,0.48,5.08,0.08) m=9.51696 std=8.82539         ] [06m:27s<01m:34s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.04) m=9.7166 std=9.19306====>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.08) m=9.94835 std=9.27002===>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.12) m=10.1402 std=9.35094===>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.16) m=10.2898 std=9.44933===>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.2) m=10.3968 std=9.57607====>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.24) m=10.4627 std=9.73631===>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.28) m=10.4903 std=9.92731===>     ] [07m:09s<00m:58s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.32) m=10.4835 std=10.138====>     ] [07m:09s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.36) m=10.4465 std=10.3505===>     ] [07m:09s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.4) m=10.3833 std=10.5436====>     ] [07m:09s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.44) m=10.2976 std=10.6967===>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.48) m=10.1923 std=10.793====>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.48,5.56,0.52) m=10.0696 std=10.8221===>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.4,5.56,0.56) m=9.93796 std=10.5708====>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.4,5.56,0.6) m=9.80455 std=10.4657=====>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.4,5.56,0.64) m=9.65913 std=10.3145====>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
(0.4,0.4,5.56,0.68) m=9.50213 std=10.1239====>     ] [07m:10s<00m:57s] Running Angle Search Experiment p=2
[==================================================] [08m:10s<00m:00s] Running Angle Search Experiment p=2  */
	second_round_angles.push_back(Angl(0.4,0.48,5.56,0.28, 10.4903, 9.92731));

	for(auto &setting : second_round_angles){
		double angles[4];
		angles[0] = setting.g1;
		angles[1] = setting.b1;
		angles[2] = setting.g2;
		angles[3] = setting.b2;
		AngleSearchExperiment::Cost cost = this->_cost_fn(&this->test_set, angles);

		std::cerr<<"Train/Test mean="<<setting.train_mean<<"/"<<cost.mean<<"   std="<<setting.train_std<<"/"<<cost.stdev<<"\n";

		//loge("Break from second_roung_angles");
		//break;

	}


}
void experiment_runner(ExperimentSetup* setup, std::string experiment_name, int loglevel){

	FastVQA::Qaoa qaoa_instance;

	if(setup->experiment_type == "angleSearch"){
		//this is handled by another function
	}else if(setup->experiment_type == "manyParams"){

		std::vector<std::vector<double>> results;

		std::mt19937 gen(0); //rd() instead of 0 - seed
		std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);

		int num_params = setup->num_rand_params;
		for(int i = 0; i < setup->num_rand_params; ++i){
			std::cerr<<setup->qaoaOptions->p<<std::endl;
			setup->qaoaOptions->initial_params.clear();

			for(int j = 0; j < setup->qaoaOptions->p*2; ++j){
				setup->qaoaOptions->initial_params.push_back(dis(gen));
			}

			FastVQA::ExperimentBuffer buffer;
			buffer.storeQuregPtr = true;
			qaoa_instance.run_qaoa(&buffer, setup->hamiltonian, setup->qaoaOptions);
			if(buffer.opt_val != setup->minimum_energy){
				logw("Solution returned has incorrect energy", loglevel);
			}
			/*for(auto &param: buffer.initial_params){
				std::cerr<<"INIT: "<<std::get<0>(param)<<" "<<std::get<1>(param)<<std::endl;
			}*/

			setup->qaoaOptions->accelerator->options.createQuregAtEachInilization = false; //we do not need to re-initialize qureg fir further runs

			std::vector<double> result(buffer.stateVector->numAmpsTotal);
			for(long long int j = 0; j < buffer.stateVector->numAmpsTotal; ++j){
				result[j] = buffer.stateVector->stateVec.real[j]*buffer.stateVector->stateVec.real[j]+buffer.stateVector->stateVec.imag[j]*buffer.stateVector->stateVec.imag[j];
			}
			results.push_back(result);
		}
		std::ofstream f;
		f.open(experiment_name+".many_params");

		for(auto &row : results){
			for(auto &x : row){
				f << x << " ";
			}f << "\n";
		}
		f.close();

	}else if(setup->experiment_type == "singleInstance"){

		FastVQA::ExperimentBuffer buffer;
		qaoa_instance.run_qaoa(&buffer, setup->hamiltonian, setup->qaoaOptions);
		if(buffer.opt_val != setup->minimum_energy){
			logw("Solution returned has incorrect energy", loglevel);
		}
		for(auto &param: buffer.initial_params){
			std::cerr<<std::get<0>(param)<<" "<<std::get<1>(param)<<"\n";
		}

		//hit_rates.push_back(buffer.getTotalHitRate());

	}else{
		loge("Wrong Experiment Name provided");
	}

}
