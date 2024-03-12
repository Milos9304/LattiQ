#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include <fstream>
#include <numeric>
#include <boost/tuple/tuple.hpp>
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"

const double pi = 3.141592654;

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

			ground_state_overlap+=buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
		}

		//for(auto &e: instance.eigenspace){
		//	std::cerr<<e.first<<" "<<e.second<<"\n";
		//}

		//std::cerr<<"degeneracy = " << instance.h.custom_solutions.size() << "\n";
		//std::cerr<<"GS_overlap = " << ground_state_overlap << "\n";
		//std::cerr<<"random guess = " << instance.random_guess << "\n";

		num_sols.push_back(instance.h.custom_solutions.size());
		gs_overlaps.push_back(ground_state_overlap / instance.random_guess);
		i++;

		//loge("Tried only one instance");
		//break;

	}

	double sum = std::accumulate(gs_overlaps.begin(), gs_overlaps.end(), 0.0);
	double mean = sum / gs_overlaps.size();

	double sq_sum = std::inner_product(gs_overlaps.begin(), gs_overlaps.end(), gs_overlaps.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / gs_overlaps.size() - mean * mean);

	double sum_mean_num_of_sols = std::accumulate(num_sols.begin(), num_sols.end(), 0.0);
	double mean_num_of_sols = sum_mean_num_of_sols / num_sols.size();

	return Cost(mean, stdev, mean_num_of_sols);

}

void AngleSearchExperiment::run(){

	//run_p2_test();
	run_p1();

	return;

	if(this->qaoaOptions->p == 1)
		run_p1();
	else if(this->qaoaOptions->p == 2)
		run_p2();
	else
		throw_runtime_error("Not implemented");

}

void AngleSearchExperiment::run_p1(){

	double *angles = (double*) malloc(this->num_params * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

	/*if(this->num_params != 2)
		throw_runtime_error("unimplemented");
	*/

	std::vector<std::tuple<double, double>> first_round_angles;

	double mean_threshold = 0.032;
	double stdev_threshold = 0.019;

	double beta_min = 0;//pi/16;
	double beta_max = pi;//;pi/8;

	double gamma_min = 0;
	double gamma_max = 2*pi;

	double range_beta = beta_max - beta_min;
	double range_gamma = gamma_max - gamma_min;
	double incr_beta = 0.04;///0.12;
	double incr_gamma = 0.04;///0.01;

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
		option::PostfixText{"Running Angle Search Experiment"},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::ForegroundColor{Color::yellow},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

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

			if(cost.mean >= mean_threshold && cost.stdev <= stdev_threshold){
				first_round_angles.push_back(std::tuple<double, double>(gamma, beta));

				//std::cerr<<"("<<gamma<<" "<<beta<<") m="<<cost.mean<<" std="<<cost.stdev<<"\n";

			}

		}
		x++;
	}

	logi("First round size: " + std::to_string(first_round_angles.size()), loglevel);


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

};

void AngleSearchExperiment::run_p2(){

	std::vector<std::tuple<double, double>> first_round_angles;
	first_round_angles.push_back(std::tuple<double, double>(0.12, 0.31635));
	first_round_angles.push_back(std::tuple<double, double>(0.13, 0.31635));
	first_round_angles.push_back(std::tuple<double, double>(0.52, 0.31635));
	first_round_angles.push_back(std::tuple<double, double>(0.93, 0.31635));
	first_round_angles.push_back(std::tuple<double, double>(1.68, 0.31635));
	first_round_angles.push_back(std::tuple<double, double>(2.09, 0.31635));
	first_round_angles.push_back(std::tuple<double, double>(2.89, 0.31635));

	double *angles = (double*) malloc(4 * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

	double mean_threshold = 0.04;
		double stdev_threshold = 0.019;
		double beta_min = 0;
		double beta_max = pi/4;

		double gamma_min = -pi;
		double gamma_max = pi;

		double range_beta = beta_max - beta_min;
		double range_gamma = gamma_max - gamma_min;
		double incr_beta = /*0.04*/0.12;
		double incr_gamma = /*0.04*/0.01;

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
			option::PostfixText{"Running Angle Search Experiment"},
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

/*	(0.12,0.31635,-0.121593,0.12) m=0.0430018 std=0.01738521m:22s<01m:28s] Running Angle Search Experiment
	(0.12,0.31635,-0.121593,0.24) m=0.0472991 std=0.01866071m:22s<01m:28s] Running Angle Search Experiment
	(1.68,0.31635,1.41841,0.12) m=0.0403283 std=0.0175522[02m:04s<00m:47s] Running Angle Search Experiment
	(1.68,0.31635,1.42841,0.12) m=0.0417132 std=0.0188931[02m:04s<00m:46s] Running Angle Search Experiment
	(0.93,0.31635,1.79841,0.12) m=0.0401401 std=0.017432 [02m:14s<00m:36s] Running Angle Search Experiment
	(0.13,0.31635,1.79841,0.24) m=0.041746 std=0.0188806 [02m:14s<00m:36s] Running Angle Search Experiment
	(0.93,0.31635,2.16841,0.12) m=0.0423418 std=0.0184396[02m:25s<00m:26s] Running Angle Search Experiment
	(0.52,0.31635,2.58841,0.12) m=0.0402696 std=0.0185866[02m:36s<00m:15s] Running Angle Search Experiment
	(0.52,0.31635,2.59841,0.12) m=0.0426105 std=0.0182783[02m:36s<00m:14s] Running Angle Search Experiment
	(0.52,0.31635,2.60841,0.12) m=0.0413229 std=0.0182011[02m:36s<00m:14s] Running Angle Search Experiment
	(0.13,0.31635,2.96841,0.12) m=0.0406699 std=0.0188591[02m:46s<00m:04s] Running Angle Search Experiment
	(0.13,0.31635,2.97841,0.12) m=0.0401275 std=0.017309 [02m:46s<00m:04s] Running Angle Search Experiment
	(0.13,0.31635,2.98841,0.12) m=0.0421467 std=0.0175874[02m:47s<00m:04s] Running Angle Search Experiment
	(0.13,0.31635,2.98841,0.24) m=0.0473511 std=0.0178743[02m:47s<00m:04s] Running Angle Search Experiment
	(0.13,0.31635,2.98841,0.36) m=0.0455804 std=0.0187591[02m:47s<00m:04s] Running Angle Search Experiment
	(0.12,0.31635,3.00841,0.12) m=0.0425464 std=0.0188244[02m:47s<00m:03s] Running Angle Search Experiment
	(0.12,0.31635,3.01841,0.12) m=0.0409526 std=0.0176841[02m:47s<00m:03s] Running Angle Search Experiment
	[==================================================] [02m:51s<00m:00s] Running Angle Search Experiment */

	std::vector<Angl> second_round_angles;

	second_round_angles.push_back(Angl(0.12,0.31635,-0.121593,0.12,0.0430018,0.01738521));
	second_round_angles.push_back(Angl(0.12,0.31635,-0.121593,0.24,0.0472991,0.01866071));
	second_round_angles.push_back(Angl(1.68,0.31635,1.41841,0.12,0.0403283,0.0175522));
	second_round_angles.push_back(Angl(1.68,0.31635,1.42841,0.12,0.0417132,0.0188931));
	second_round_angles.push_back(Angl(0.93,0.31635,1.79841,0.12,0.0401401,0.017432 ));
	second_round_angles.push_back(Angl(0.13,0.31635,1.79841,0.24,0.041746,0.0188806 ));
	second_round_angles.push_back(Angl(0.93,0.31635,2.16841,0.12,0.0423418,0.0184396));
	second_round_angles.push_back(Angl(0.52,0.31635,2.58841,0.12,0.0402696,0.0185866));
	second_round_angles.push_back(Angl(0.52,0.31635,2.59841,0.12,0.0426105,0.0182783));
	second_round_angles.push_back(Angl(0.52,0.31635,2.60841,0.12,0.0413229,0.0182011));
	second_round_angles.push_back(Angl(0.13,0.31635,2.96841,0.12,0.0406699,0.0188591));
	second_round_angles.push_back(Angl(0.13,0.31635,2.97841,0.12,0.0401275,0.017309 ));
	second_round_angles.push_back(Angl(0.13,0.31635,2.98841,0.12,0.0421467,0.0175874));
	second_round_angles.push_back(Angl(0.13,0.31635,2.98841,0.24,0.0473511,0.0178743));
	second_round_angles.push_back(Angl(0.13,0.31635,2.98841,0.36,0.0455804,0.0187591));
	second_round_angles.push_back(Angl(0.12,0.31635,3.00841,0.12,0.0425464,0.0188244));
	second_round_angles.push_back(Angl(0.12,0.31635,3.01841,0.12,0.0409526,0.0176841));

	for(auto &setting : second_round_angles){
		double angles[4];
		angles[0] = setting.g1;
		angles[1] = setting.b1;
		angles[2] = setting.g2;
		angles[3] = setting.b2;
		AngleSearchExperiment::Cost cost = this->_cost_fn(&this->train_set, angles);

		std::cerr<<"Train/Test mean="<<setting.train_mean<<"/"<<cost.mean<<"   std="<<setting.train_std<<"/"<<cost.stdev<<"\n";

		loge("Break from second_roung_angles");
		break;

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
