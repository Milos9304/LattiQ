#include "g1.h"
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

std::string nlopt_res_to_str2(int result){
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

G1::G1(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){
	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;

	this->mapOptions->penalty = 0;
}

void G1::run(){

	std::vector<double> xs, ys_test, ys_train;

	std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset;
	std::vector<std::vector<AlphaMinimizationExperimentInstance>> test_dataset;

	for(int m = m_start; m <= m_end; ++m){

		std::vector<AlphaMinimizationExperimentInstance> m_train_instances;
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

		int num_test_instances = num_instances * test_ratio;
		int num_train_instances = num_instances - num_test_instances;
		logi("Num_instances="+std::to_string(num_instances)+" Test ratio="+std::to_string(test_ratio)+" "+"train_size="+std::to_string(num_train_instances)+" "+"test_size="+std::to_string(num_test_instances), this->loglevel);

		if(num_test_instances == 0)
			loge("Zero test instances!");

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
			if(counter++ < num_train_instances)
				m_train_instances.push_back(instance);
			else
				m_test_instances.push_back(instance);
		}
		qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
		qaoaOptions->accelerator->finalize();

		train_dataset.push_back(m_train_instances);
		test_dataset.push_back(m_test_instances);
	}

	logi("Dataset generated");

	for(auto &graph_line: train_dim_ratios){

		int m_train_end = graph_line.first;
		int m_test_end = graph_line.second;
		assert(m_test_end == m_end);

		this->qaoaOptions->ftol = 1e-10;
		this->qaoaOptions->max_iters = 500;

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

		std::pair<double, double> final_ab;
		FastVQA::OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
			iteration_i++;
			bar.tick();
			std::vector<double> angles(x);
			//std::cerr<<this->_cost_fn(&this->train_set, &angles[0]).mean<<std::endl;


			double nom=0, den=0;

			/*for(auto &dim: train_dataset){
				double overlap = this->_cost_fn(dim, &angles[0]);
				//std::cerr<<overlap<<std::endl;
				nom += log(overlap) * dim[0].m;
				den += dim[0].m * dim[0].m;
				//overlaps.push_back(this->_cost_fn(dim, &angles[0]));
			}
			double alpha = -nom/den;
			*/

			//std::vector<double> overlaps;
			double sum_yi=0;
			double sum_xi=0;
			double sum_xi2=0;
			double sum_xi_yi=0;
			for(auto &dim: train_dataset){

				if(dim[0].m > m_train_end)
					continue;

				double overlap = this->_cost_fn(dim, &angles[0]);
				//overlaps.push_back(log(overlap));
				sum_yi+=log(overlap);
				sum_xi+=dim[0].m;
				sum_xi2+=dim[0].m * dim[0].m;
				sum_xi_yi+=log(overlap)*dim[0].m;
			}
			double a=0,b=0;
			a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(/*train_dataset.size()*/(m_train_end-3+1)*sum_xi2-sum_xi*sum_xi);
			b=((m_train_end-3+1)*sum_xi_yi-sum_xi*sum_yi)/((m_train_end-3+1)*sum_xi2-sum_xi*sum_xi);

			//std::cerr<<"e^"<<a<<"+n*"<<b<<std::endl;

			double alpha = -b;
			final_ab.first = a;
			final_ab.second = a;

			//std::cerr<<nom<<" "<<den<<" "<<alpha<<std::endl;
			return alpha;//-this->_cost_fn(&this->train_set, &angles[0]).mean;
		}, num_params);
		//std::cerr<<"e^"<<final_ab.first<<"+n*"<<final_ab.second<<std::endl;

		/*std::vector<double> initial_params;
		std::mt19937 gen(0); //rd() instead of 0 - seed
		std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
		for(int i = 0; i < num_params/2; ++i){
			double param1 = dis(gen);
			double param2 = dis(gen);
			//std::cerr<<param1<<" "<<param2<<std::endl;
			initial_params.push_back(param1);
			initial_params.push_back(param2);
		}*/
		std::vector<double> initial_params {0.331,0.636, 0.645, 0.534, 0.731, 0.463, 0.837, 0.360, 1.009, 0.259, 1.126, 0.139};

		std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
		std::vector<double> upperBounds(initial_params.size(), 3.141592654);

		logd("QAOA starting optimization", this->loglevel);
		FastVQA::OptResult result = this->qaoaOptions->optimizer->optimize(f, initial_params, this->qaoaOptions->ftol, this->qaoaOptions->max_iters, lowerBounds, upperBounds);
		logd("QAOA finishing optimization", this->loglevel);

		xs.push_back(m_train_end);
		ys_train.push_back(result.first.first);
		//std::cerr<<"alpha: "<< result.first.first <<"\n";
		//std::cerr<<"num_iters: "<<iteration_i<<std::endl;
		if(m_train_end == 9){
			std::cerr<<"Final angles: ";
			for(auto &a: result.first.second){
				std::cerr<<a<<" ";
			}
		}
		//std::cerr<<"\n"<<nlopt_res_to_str2(result.second)<<std::endl;

		//EVALUATE TEST DATASET
		std::vector<double> final_angles(result.first.second);

		double nom=0, den=0;

		/*for(auto &dim: test_dataset){
			double overlap = this->_cost_fn(dim, &final_angles[0]);
			nom += log(overlap) * dim[0].m;
			den += dim[0].m * dim[0].m;
		}*/
		//double  alpha = -nom/den;
		//std::cerr<<"Simple alpha: "<< -nom/den << std::endl;


		double sum_yi=0;
		double sum_xi=0;
		double sum_xi2=0;
		double sum_xi_yi=0;
		for(auto &dim: test_dataset){

			if(dim[0].m > m_train_end)
				continue;

			double overlap = this->_cost_fn(dim, &final_angles[0]);
			//overlaps.push_back(log(overlap));
			sum_yi+=log(overlap);
			sum_xi+=dim[0].m;
			sum_xi2+=dim[0].m * dim[0].m;
			sum_xi_yi+=log(overlap)*dim[0].m;
		}
		double a=0,b=0;
		a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/((m_train_end-3+1)*sum_xi2-sum_xi*sum_xi);
		b=((m_train_end-3+1)*sum_xi_yi-sum_xi*sum_yi)/((m_train_end-3+1)*sum_xi2-sum_xi*sum_xi);
		//std::cerr<<pow(2.71828, a)<<"e^n*"<<b<<std::endl;
		//std::cerr<<"e^"<<a<<"+n*"<<b<<std::endl;
		double alpha = -b;
		ys_test.push_back(alpha);
		//std::cerr<<"Test alpha: "<<alpha<<std::endl;
	}

	std::cerr<<"xs=[";
	int i=0;
	for(auto &x: xs){
		std::cerr<<x;
		if(i++<xs.size())
			std::cerr<<",";
	}
	std::cerr<<"]"<<std::endl;

	std::cerr<<"ys_train=[";
	i=0;
	for(auto &x: ys_train){
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

double G1::_cost_fn(std::vector<AlphaMinimizationExperimentInstance> dataset, const double *angles, bool use_database){
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
