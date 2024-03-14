#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include <fstream>
#include <numeric>
#include <boost/tuple/tuple.hpp>
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"

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

	loge("Saving eigenspaaace which is not needed and very costly");

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
		for(const auto &sol: instance.solutions){
			if(sol.value < min_energy)
				min_energy = sol.value;
		}
		instance.min_energy = min_energy;
		instance.random_guess = l.get_random_guess_one_vect() * instance.h.custom_solutions.size();

		if(i < num_train_instances)
			this->train_set.push_back(instance);
		else
			this->test_set.push_back(instance);
	}

	logi("Experiment dataset generated", this->loglevel);

}

AngleSearchExperiment::Cost AngleSearchExperiment::_cost_fn(std::vector<Instance>* dataset, double *angles){

	FastVQA::Qaoa qaoa_instance;
	std::vector<double> gs_overlaps;
	std::vector<int> num_sols;

	int i = 0;
	for(auto &instance: (*dataset)){

		FastVQA::ExperimentBuffer buffer;
		buffer.storeQuregPtr = true;
		qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
		double ground_state_overlap = 0;
		for(auto &sol: instance.solutions){
			long long int index = sol.index;

			/*if(instance.h.getMatrixRepresentation2(true)(index, index) != instance.min_energy+1)
				logw("An outcome with different energy marked as a solution!");*/

			/*if(i == 73){
				std::cerr << index<<" "<<std::bitset<5>(index)<<"\n";
			}*/

			ground_state_overlap+=buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
		}

		//for(auto &e: instance.eigenspace){
		//	std::cerr<<e.first<<" "<<e.second<<"\n";
		//}

		//std::cerr<<"degeneracy = " << instance.h.custom_solutions.size() << "\n";
		//std::cerr<<"GS_overlap = " << ground_state_overlap << "\n";
		//std::cerr<<"random guess = " << instance.random_guess << "\n";

		/*if(i == 73){
			std::cerr<<ground_state_overlap / instance.random_guess<<std::endl;
			std::cerr<<ground_state_overlap<<std::endl; //0.00119142
			std::cerr<<instance.random_guess<<std::endl;

			for(int j = 0; j < buffer.stateVector->numAmpsTotal; ++j){
				double c = buffer.stateVector->stateVec.real[j]*buffer.stateVector->stateVec.real[j]+buffer.stateVector->stateVec.imag[j]*buffer.stateVector->stateVec.imag[j];
				if(c > 0.01)
					std::cerr<<j<<"   "<<std::bitset<8>(j)<<":   "<<c<<std::endl;
			}
		}*/

		num_sols.push_back(instance.h.custom_solutions.size());
		gs_overlaps.push_back(ground_state_overlap / instance.random_guess);
		i++;

		//loge("Tried only one instance");
		//break;

	}

	double debug_min=10, debug_i;

	/*i=0;
	for(const auto &o: gs_overlaps){
		if(o < debug_min){
			debug_min = o;
			debug_i = i;
		}i++;
	}
	std::cerr<<debug_min<<" "<<debug_i<<"\n";*/

	double sum = std::accumulate(gs_overlaps.begin(), gs_overlaps.end(), 0.0);
	double mean = sum / gs_overlaps.size();

	double sq_sum = std::inner_product(gs_overlaps.begin(), gs_overlaps.end(), gs_overlaps.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / gs_overlaps.size() - mean * mean);

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
		run_p2_full_bruteforce();//run_p2();
	else
		throw_runtime_error("Not implemented");

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

	double mean_threshold = 1.2;
	double stdev_threshold = 2;//0.019;

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
(0.4 0.56) m=1.35427 std=1.17801                   ] [00m:08s<00m:59s] Running Angle Search Experiment
(0.4 1.2) m=1.2138 std=1.06883                     ] [00m:09s<00m:58s] Running Angle Search Experiment
(0.84 0.56) m=1.20434 std=1.22635                  ] [00m:18s<00m:50s] Running Angle Search Experiment
(1.24 0.64) m=1.21043 std=0.977747                 ] [00m:27s<00m:41s] Running Angle Search Experiment
(1.24 0.68) m=1.21611 std=0.97851                  ] [00m:27s<00m:41s] Running Angle Search Experiment
(1.24 0.72) m=1.21713 std=0.974095                 ] [00m:27s<00m:41s] Running Angle Search Experiment
[==================================================] [01m:06s<00m:00s] Running Angle Search Experiment
 *
 */

	std::vector<std::tuple<double, double>> first_round_angles;
	first_round_angles.push_back(std::tuple<double, double>(0.4, 0.56));
	first_round_angles.push_back(std::tuple<double, double>(0.84, 0.56));
	first_round_angles.push_back(std::tuple<double, double>(1.24, 0.72));


	double *angles = (double*) malloc(4 * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

		double mean_threshold = 1.4;
		double stdev_threshold = 100;
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

/*	(0.4,0.56,3.6,0.08) m=1.42399 std=1.23366          ] [04m:37s<03m:29s] Running Angle Search Experiment
(0.4,0.56,3.6,0.12) m=1.4455 std=1.25496           ] [04m:37s<03m:29s] Running Angle Search Experiment
(0.4,0.56,3.6,0.16) m=1.45727 std=1.27018          ] [04m:38s<03m:29s] Running Angle Search Experiment
(0.4,0.56,3.6,0.2) m=1.4595 std=1.27838            ] [04m:38s<03m:29s] Running Angle Search Experiment

(0.4,0.56,3.6,0.24) m=1.45298 std=1.27902          ] [04m:38s<03m:29s] Running Angle Search Experiment
(0.4,0.56,3.6,0.28) m=1.43904 std=1.272            ] [04m:38s<03m:29s] Running Angle Search Experiment
(0.4,0.56,3.6,0.32) m=1.41933 std=1.2577           ] [04m:38s<03m:29s] Running Angle Search Experiment
(0.4,0.56,5.56,0.08) m=1.41008 std=1.18321===>     ] [07m:10s<00m:58s] Running Angle Search Experiment
(0.4,0.56,5.56,0.12) m=1.42966 std=1.18758===>     ] [07m:10s<00m:58s] Running Angle Search Experiment

(0.4,0.56,5.56,0.16) m=1.44393 std=1.19591===>     ] [07m:10s<00m:58s] Running Angle Search Experiment
(0.4,0.56,5.56,0.2) m=1.45345 std=1.20976====>     ] [07m:10s<00m:58s] Running Angle Search Experiment
(0.4,0.56,5.56,0.24) m=1.45897 std=1.22983===>     ] [07m:10s<00m:58s] Running Angle Search Experiment

(0.4,0.56,5.56,0.28) m=1.4613 std=1.25567====>     ] [07m:10s<00m:58s] Running Angle Search Experiment
(0.4,0.56,5.56,0.32) m=1.4612 std=1.28581====>     ] [07m:10s<00m:58s] Running Angle Search Experiment
(0.4,0.56,5.56,0.36) m=1.45931 std=1.31806===>     ] [07m:10s<00m:57s] Running Angle Search Experiment
(0.4,0.56,5.56,0.4) m=1.45608 std=1.35009====>     ] [07m:10s<00m:57s] Running Angle Search Experiment
(0.4,0.56,5.56,0.44) m=1.45168 std=1.37988===>     ] [07m:10s<00m:57s] Running Angle Search Experiment

(0.4,0.56,5.56,0.48) m=1.44609 std=1.40612===>     ] [07m:10s<00m:57s] Running Angle Search Experiment
(0.4,0.56,5.56,0.52) m=1.43906 std=1.42829===>     ] [07m:10s<00m:57s] Running Angle Search Experiment
(0.4,0.56,5.56,0.56) m=1.43023 std=1.44644===>     ] [07m:11s<00m:57s] Running Angle Search Experiment
(0.4,0.56,5.56,0.6) m=1.41918 std=1.46087====>     ] [07m:11s<00m:57s] Running Angle Search Experiment
(0.4,0.56,5.56,0.64) m=1.40553 std=1.47168===>     ] [07m:11s<00m:57s] Running Angle Search Experiment
[==================================================] [08m:08s<00m:00s] Running Angle Search Experiment */

	std::vector<Angl> second_round_angles;
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.08, 1.42399, 1.23366));
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.12, 1.4455 ,1.25496  ));
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.16,1.45727 ,1.27018 ));
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.2,1.4595, 1.27838));
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.24, 1.45298, 1.27902));
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.28,1.43904 ,1.272 ));
	second_round_angles.push_back(Angl(0.4,0.56,3.6,0.32, 1.41933, 1.2577 ));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.08, 1.41008, 1.18321));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.12, 1.42966 ,1.18758));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.16, 1.44393 ,1.19591));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.2, 1.45345 ,1.20976));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.24, 1.45897 ,1.22983));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.28, 1.4613, 1.25567));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.32, 1.4612, 1.28581));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.36, 1.45931, 1.31806));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.4,1.45608, 1.35009));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.44,1.45168, 1.37988));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.48, 1.44609, 1.40612));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.52,1.43906, 1.42829));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.56,1.43023, 1.44644));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.6, 1.41918, 1.46087));
	second_round_angles.push_back(Angl(0.4,0.56,5.56,0.64,1.40553, 1.4716));

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
