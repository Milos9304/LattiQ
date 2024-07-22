#include "g2.h"
#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include "hermite_factor.h"

#include <fstream>
#include <numeric>
#include <boost/tuple/tuple.hpp>
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"
#include <iomanip>
#include <cmath>

# if QuEST_PREC==1
	#define QREAL_MAX FLOAT_MAX
# elif QuEST_PREC==2
    #define QREAL_MAX DOUBLE_MAX
# elif QuEST_PREC==4
   #define QREAL_MAX 1e300
# endif

#define s(X) std::to_string(X)

const double pi = 3.141592654;

G2::G2(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){
	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;

	this->mapOptions->penalty = 0;
}

void G2::run(){

	std::vector<double> final_angles{-0.0232639, -0.916444, 0.735415 ,0.550755 ,0.754262, 0.282567, 0.865909 ,0.360061 ,0.875211 ,0.284794, 1.1415, 0.0951161};

	std::vector<double> xs, ys_test;
	std::vector<std::vector<AlphaMinimizationExperimentInstance>> test_dataset;

	for(int m = m_start; m <= m_end; ++m){

		std::vector<AlphaMinimizationExperimentInstance> m_test_instances;

		long long int num_instances = q;
		for(int i = 0; i < (m-n); ++i){ //pow
			num_instances *= q;
			if(max_num_instances < num_instances){
				num_instances = max_num_instances;
				break;
			}
		}

		if(max_num_instances < num_instances)
			num_instances = max_num_instances;

		int num_test_instances = num_instances;
		logi("Num_instances="+std::to_string(num_instances), this->loglevel);

		GeneratorParam param(q, n, m, true, 97, num_instances); //q, n, m, shuffle, seed, cutoff
		std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);
		int nbQubits_acc = -1;
		int counter = 0;
		for(auto &gw : gramian_wrappers){
			Lattice l(gw.hamiltonian, gw.name);
			AlphaMinimizationExperimentInstance instance;
			instance.h = l.getHamiltonian(mapOptions);
			instance.q = q;
			instance.n = n;
			instance.m = m;
			if(nbQubits_acc < 0)
				nbQubits_acc = instance.h.nbQubits;
			else if(nbQubits_acc != instance.h.nbQubits){
				//Need to destroy qureg, because next time experiments with different number of qubits will be run
				qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
				qaoaOptions->accelerator->finalize();
			}

			qaoaOptions->accelerator->initialize(&instance.h);
			qaoaOptions->accelerator->options.createQuregAtEachInilization = false;

			instance.zero_solutions = qaoaOptions->accelerator->getSolutions();
			if(instance.zero_solutions.size() != 1)
				throw_runtime_error("CmQaoaExperiment: Unimplemented, more than 1 solution marked");

			for(auto &sol: instance.zero_solutions){
				if(sol.value != 0)
					throw_runtime_error("CmQaoaExperiment: Something else than 0 state marked as a solution");
			}
			long long int zero_index = instance.zero_solutions[0].index;

			long long int numAmpsTotal = qaoaOptions->accelerator->getQuregPtr()->numAmpsTotal;

			FastVQA::RefEnergies refEnergies = qaoaOptions->accelerator->getEigenspace();

			qreal min = QREAL_MAX;//refEnergies[0].value;
			for(long long int j = 0; j < numAmpsTotal; ++j){

				if(refEnergies[j].value == min)
					instance.sv_solutions.push_back(FastVQA::RefEnergy(min, refEnergies[j].index, false));
				else if(refEnergies[j].value > 0 && refEnergies[j].value < min){
					instance.sv_solutions.clear();
					min = refEnergies[j].value;
					instance.sv_solutions.push_back(FastVQA::RefEnergy(min, refEnergies[j].index, false));
				}
			}

			m_test_instances.push_back(instance);
		}
		qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
		qaoaOptions->accelerator->finalize();

		test_dataset.push_back(m_test_instances);
		if(test_dataset.size() == 1)
			continue;
		double sum_yi=0;
		double sum_xi=0;
		double sum_xi2=0;
		double sum_xi_yi=0;
		for(auto &dim: test_dataset){

			double overlap = this->_cost_fn(dim, &final_angles[0]);
			//overlaps.push_back(log(overlap));
			sum_yi+=log(overlap);
			sum_xi+=dim[0].m;
			sum_xi2+=dim[0].m * dim[0].m;
			sum_xi_yi+=log(overlap)*dim[0].m;

			std::cerr<<dim[0].m<<" ";
		}
		double a=0,b=0;
		a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(test_dataset.size()*sum_xi2-sum_xi*sum_xi);
		b=(test_dataset.size()*sum_xi_yi-sum_xi*sum_yi)/(test_dataset.size()*sum_xi2-sum_xi*sum_xi);
		//std::cerr<<pow(2.71828, a)<<"e^n*"<<b<<std::endl;
		//std::cerr<<"e^"<<a<<"+n*"<<b<<std::endl;
		double alpha = -b;
		std::cerr<<alpha<<"\n";
		xs.push_back(m);
		ys_test.push_back(alpha);
	}


	std::cerr<<"xs=[";
	int i=0;
	for(auto &x: xs){
		std::cerr<<x;
		if(i++<xs.size())
			std::cerr<<",";
	}
	std::cerr<<"]"<<std::endl;

	std::cerr<<"ys_test=[";
	i=0;
	for(auto &x: ys_test){
		std::cerr<<x;
		if(i++<xs.size())
			std::cerr<<",";
	}
	std::cerr<<"]"<<std::endl;

}

double G2::_cost_fn(std::vector<AlphaMinimizationExperimentInstance> dataset, const double *angles, bool use_database){
	std::vector<int> num_sols;

		int i = 0;

		double mean;
		double stdev;
		FastVQA::Qaoa qaoa_instance;

		std::vector<double> gs_overlaps;

		/*loge("Overriding angles");
		bool print=false;//std::cerr<<angles[0] <<" " <<angles[1]<<"\n";
		if(angles[0] == 0.01 && angles[1]>0.3899 && angles[1]<0.3901){
			angles[0] = 2;angles[1] = 0;
			print=true;
		}*/

		for(auto &instance: dataset){

			if(instance.h.nbQubits != this->qaoaOptions->accelerator->getNumQubitsInQureg()){
				this->qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
				this->qaoaOptions->accelerator->finalize();
			}

			FastVQA::ExperimentBuffer buffer;
			buffer.storeQuregPtr = true;

			double ground_state_overlap = 0;
			if(use_database){throw;
				/*Database::DatasetRow output_row;
				this->database->getOrCalculate_qary_with_fixed_angles(&buffer, angles, 6, &instance.h, &output_row, this->qaoaOptions, &qaoa_instance);
				for(auto &sol: instance.sv_solutions){
					long long int index = sol.index;
					ground_state_overlap+=output_row.finalStateVectorMap[index].second;
				}*/
			}else{
				//qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
				qaoa_instance.run_cm_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles, instance.zero_solutions[0].index);
				//qaoa_instance.run_cm_qaoa(&buffer, &instance.h, this->qaoaOptions, instance.zero_solutions[0].index);

				/*for(auto &f: buffer.finalParams){
					std::cerr<<f<<" ";
				}std::cerr<<std::endl;*/

				for(auto &sol: instance.sv_solutions){
					long long int index   = sol.index;
					ground_state_overlap += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
				}
			}

			num_sols.push_back(instance.h.custom_solutions.size());
			qreal improvement_ratio = ground_state_overlap;// / instance.random_guess;
			//std::cerr<<improvement_ratio<<std::endl;throw;
			/*if(print){

				std::cerr<<angles[0]<<"   "<<angles[1]<<std::endl;

				std::cerr<<improvement_ratio<<" "<<ground_state_overlap <<" "<< instance.random_guess << std::endl;
				std::cerr<<"Zero_index: "<< instance.zero_solutions[0].index<<std::endl;
				double zero_prob=0;
				for(auto &z: instance.zero_solutions){
					long long int index   = z.index;
					zero_prob += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
				}
				std::cerr<<"Zero_probability: "<< zero_prob<<std::endl;

				ground_state_overlap = 0;
				for(auto &sol: instance.sv_solutions){
					double tmp;
					long long int index   = sol.index;
					tmp = buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
					std::cerr<<"Sol index: "<<index<<" with prob "<<tmp<<std::endl;
					ground_state_overlap += tmp;
				}

				std::cerr<<std::endl<<std::endl;

				FastVQA::RefEnergies refEnergies = qaoaOptions->accelerator->getEigenspace();//delete
				for(long long int j = 0; j < buffer.stateVector->numAmpsTotal; ++j){

					long long int index = refEnergies[j].index;
					double tmp = buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];

					std::cerr<<refEnergies[j].index<<" "<<refEnergies[j].value<<" "<<tmp<<std::endl;

				}

				throw;
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

		//std::cerr<<mean<<"  "<<stdev<<std::endl;

		return mean;

}
