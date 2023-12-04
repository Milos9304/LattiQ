#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include <fstream>

AngleSearchExperiment::AngleSearchExperiment(){
	this->num_instances = pow(q,(m-n));
	if(max_num_instances < this->num_instances)
		this->num_instances = max_num_instances;
	this->num_test_instances = this->num_instances * this->test_ratio;
	this->num_train_instances = this->num_instances - this->num_test_instances;
	logi("q="+std::to_string(q)+" "+"n="+std::to_string(n)+" "+"m="+std::to_string(m)+" : "+std::to_string(this->num_instances) + " instances");
	logi("Test ratio="+std::to_string(this->test_ratio)+" "+"train_size="+std::to_string(this->num_train_instances)+" "+"test_size="+std::to_string(this->num_test_instances));

	this->_generate_dataset();
};

void AngleSearchExperiment::_generate_dataset(){

}

void AngleSearchExperiment::run(){

	int cutoff = 100;
	GeneratorParam param(q, n, m, true, 97, cutoff);
	std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

	std::cerr<<gramian_wrappers.size()<<"\n";

	for(int i = 0; i < cutoff; ++i)
		std::cerr<<gramian_wrappers[i].hamiltonian<<std::endl;

};

void experiment_runner(ExperimentSetup* setup, std::string experiment_name){

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
				logw("Solution returned has incorrect energy");
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
			logw("Solution returned has incorrect energy");
		}
		for(auto &param: buffer.initial_params){
			std::cerr<<std::get<0>(param)<<" "<<std::get<1>(param)<<"\n";
		}

		//hit_rates.push_back(buffer.getTotalHitRate());

	}else{
		loge("Wrong Experiment Name provided");
	}

}
