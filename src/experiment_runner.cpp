#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include <fstream>
#include <numeric>
#include <boost/tuple/tuple.hpp>
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"

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

	GeneratorParam param(q, n, m, true, 97, this->max_num_instances);
	std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

	nbQubits = -1;

	for(int i = 0; i < num_instances; ++i){

		Instance instance;

		Lattice l(gramian_wrappers[i].hamiltonian, gramian_wrappers[i].name);
		instance.h = l.getHamiltonian(mapOptions);
		if(nbQubits < 0)
			nbQubits = instance.h.nbQubits;
		else if(nbQubits != instance.h.nbQubits)
			loge("Instances with different number of qubits found");
		qaoaOptions->accelerator->initialize(&instance.h);
		qaoaOptions->accelerator->options.createQuregAtEachInilization = false;
		instance.solutions = qaoaOptions->accelerator->getSolutions();
		instance.min_energy = std::get<0>(instance.solutions[0]);


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
		//qaoaOptions->accelerator->initialize(&instance.h);
		qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
		double ground_state_overlap = 0;
		for(auto &sol: instance.solutions){
			if(instance.min_energy != std::get<0>(sol))
				logw("An outcome with different energy marked as a solution!");
			long long int index = std::get<1>(sol);

			/*if(instance.h.getMatrixRepresentation2(true)(index, index) != instance.min_energy+1)
				logw("An outcome with different energy marked as a solution!");*/

			ground_state_overlap+=buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
		}
		num_sols.push_back(instance.solutions.size());
		gs_overlaps.push_back(ground_state_overlap);
		i++;
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

	double *angles = (double*) malloc(this->num_params * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

	if(this->num_params != 2)
		throw_runtime_error("unimplemented");

	double range = 2*3.141592654;
	double incr = 0.04;

	int axis_range = ceil(range / incr);
	int num_iters = axis_range * axis_range;

	int x = 0;
	int y = 0;

	std::vector<std::vector<double>> final_plot_means(axis_range);
	std::vector<std::vector<double>> final_plot_stdevs(axis_range);

	double i = 0;
	AngleSearchExperiment::Cost cost;
	for(double a = -3.141592654; a < 3.141592654; a+=incr){
		for(double b = -3.141592654; b < 3.141592654; b+=incr){
			std::cerr<<(i++) / num_iters * 100 << "%" << std::endl;
			angles[0] = a;angles[1] = b;
			cost = this->_cost_fn(&this->train_set, angles);
			final_plot_means[x].push_back(cost.mean);
			final_plot_stdevs[x].push_back(cost.stdev);
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
	gp << "set xrange [ 0 : " << axis_range-1 <<" ] \n";
	gp << "set yrange [ 0 : " << axis_range-1 <<" ] \n";
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
