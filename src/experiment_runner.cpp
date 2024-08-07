#include "averaged_svp.h"
#include "experiment_runner.h"
#include "io/logger.h"
#include "hermite_factor.h"
#include "run_paper_experiment.h"
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

	FastVQA::Qaoa qaoa_instance;

	std::vector<std::vector<std::tuple<int /*m*/, double, double, double, double>>> plot(this->p_end - this->p_start + 1);

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
			instance.volume = l.getVolume();

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
				loge("CmQaoaExperiment: Unimplemented, more than 1 solution marked");

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

			for(int p = this->p_start; p <= this->p_end; ++p){
				this->qaoaOptions->p = p;

				logw("Running m="+s(m)+" p="+s(p)+" qs="+s(qs)+" index="+s(i)+"     Total Qubits: "+s(m*qs), this->loglevel);


				if(m != 3 || p!= 3 || qs !=2 || i != 1)
					continue;
				//this->qaoaOptions->p = 30;


				FastVQA::ExperimentBuffer buffer, cm_buffer;
				buffer.storeQuregPtr = true;
				cm_buffer.storeQuregPtr = true;

				double zero_overlap_qaoa = 0, zero_overlap_cm_qaoa = 0;
				double sv_overlap_qaoa = 0, sv_overlap_cm_qaoa = 0;

				qaoa_instance.run_qaoa(&buffer, &instance.h, this->qaoaOptions);

				for(auto &sol: instance.zero_solutions){
					long long int index = sol.index;
					zero_overlap_qaoa    += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
					//zero_overlap_cm_qaoa += cm_buffer.stateVector->stateVec.real[index]*cm_buffer.stateVector->stateVec.real[index]+cm_buffer.stateVector->stateVec.imag[index]*cm_buffer.stateVector->stateVec.imag[index];
				}

				for(auto &sol: instance.sv_solutions){
					long long int index = sol.index;
					sv_overlap_qaoa    += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
					//sv_overlap_cm_qaoa += cm_buffer.stateVector->stateVec.real[index]*cm_buffer.stateVector->stateVec.real[index]+cm_buffer.stateVector->stateVec.imag[index]*cm_buffer.stateVector->stateVec.imag[index];
				}

				calculate_hermite_factor(m, pow(instance.volume, 1./(double)instance.m), buffer.stateVector, &refEnergies);

				qaoa_instance.run_cm_qaoa(&cm_buffer, &instance.h, this->qaoaOptions, zero_index);

				for(auto &sol: instance.zero_solutions){
					long long int index = sol.index;
					//zero_overlap_qaoa    += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
					zero_overlap_cm_qaoa += cm_buffer.stateVector->stateVec.real[index]*cm_buffer.stateVector->stateVec.real[index]+cm_buffer.stateVector->stateVec.imag[index]*cm_buffer.stateVector->stateVec.imag[index];
				}

				for(auto &sol: instance.sv_solutions){
					long long int index = sol.index;
					//sv_overlap_qaoa    += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
					sv_overlap_cm_qaoa += cm_buffer.stateVector->stateVec.real[index]*cm_buffer.stateVector->stateVec.real[index]+cm_buffer.stateVector->stateVec.imag[index]*cm_buffer.stateVector->stateVec.imag[index];
				}

				std::cerr<<"zero_overlap_qaoa: "    << zero_overlap_qaoa<<std::endl;
				std::cerr<<"zero_overlap_cm_qaoa: " << zero_overlap_cm_qaoa<<" / rg="<<1/(pow(2,instance.h.nbQubits))<<std::endl;
				std::cerr<<"sv_overlap_qaoa: "    << sv_overlap_qaoa<<std::endl;
				std::cerr<<"sv_overlap_cm_qaoa: " << sv_overlap_cm_qaoa<<std::endl;

				plot[p - this->p_start].push_back(std::tuple<int, double, double, double, double>(m, zero_overlap_qaoa, zero_overlap_cm_qaoa, sv_overlap_qaoa, sv_overlap_cm_qaoa));

				//PRINT
				/*int max_index=999999999; double max_val=-1;
				std::cerr<<"min: "<<min<<std::endl;
				for(int j = 0; j < buffer.stateVector->numAmpsTotal; ++j){
					double val;
					for(int k = 0; k < refEnergies.size(); ++k){
						if(refEnergies[k].index == j){
							val = refEnergies[k].value;
							break;
						}
					}

					double prob = buffer.stateVector->stateVec.real[j]*buffer.stateVector->stateVec.real[j]+buffer.stateVector->stateVec.imag[j]*buffer.stateVector->stateVec.imag[j];
					if(prob > max_val){
						max_index = j;
						max_val = prob;
					}

					std::cerr<<j<<": "<<prob<<"    "<<val<<"\n";
				}
				std::cerr<<std::endl;
				std::cerr<<max_val<<" "<<max_index<<std::endl;
				//PRINT

				throw;*/


			}

		}

	}

	for(int p = this->p_start; p <= this->p_end; ++p){
		std::cerr<<p<<std::endl;
		for(auto &t: plot[p - this->p_start]){
			std::cerr<<std::get<0>(t)<<" "<<std::get<1>(t)<<" "<<std::get<2>(t)<<" "<<std::get<3>(t)<<" "<<std::get<4>(t)<<" "<<std::endl;
		}
	}

}

CmQaoaExperiment::Cost CmQaoaExperiment::_cost_fn(CmQaoaExperiment::Instance*, bool use_database){
	CmQaoaExperiment::Cost cost;



	return cost;
}


AngleResultsExperiment::AngleResultsExperiment(int loglevel, int m_end, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database, int seed){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->database = database;

	logi("p="+std::to_string(this->qaoaOptions->p), this->loglevel);

	this->logfile.open("log.txt");
	this->angleAnalysisLog.open("angleAnalysis.txt");

	this->m_end = m_end;
	this->seed = seed;

	//if(this->qaoaOptions->p != 2)
	//	throw_runtime_error("Angleres is only for p=2!");

}

std::vector<AngleResultsExperiment::Instance> AngleResultsExperiment::_generate_dataset(int n, int m, bool penalise){

	std::vector<AngleResultsExperiment::Instance> dataset;

	bool new_way=true;
	long long int num_instances;

	/*std::vector<Lattice*> lattices;
	initialize_paper_experiment("qary_4_20", lattices, m, -1);
	int num_lattices = lattices.size();
	logi("Dataset qary_4_20 succesfully loaded", 0);
*/
	std::vector<HamiltonianWrapper> gramian_wrappers;

	if(new_way){
		GeneratorParam param(q, n, m, true, 97, this->max_num_instances); //q, n, m, shuffle, seed, cutoff
		gramian_wrappers = generateFromEvalDecomposition(param);//generateQaryUniform(param);

		/*num_instances = q;
		for(int i = 0; i < (m-n); ++i){ //pow
			num_instances *= q;
			if(max_num_instances < num_instances){
				num_instances = max_num_instances;
				break;
			}
		}

		if(max_num_instances < num_instances)
			num_instances = max_num_instances;*/

		num_instances = gramian_wrappers.size();

	}else{
		std::vector<Lattice*> lattices;
		initialize_paper_experiment("qary_4_20", lattices, m, -1);
		int num_lattices = lattices.size();
		//logi("Dataset qary_4_20 succesfully loaded", 0);
		for(auto &l: lattices){

			fplll::ZZ_mat<mpz_t>* zz_mat = l->get_current_lattice();
			int nrows = l->n_rows;
			int ncols = l->n_cols;

			Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix(nrows, ncols);

			for(int yy = 0; yy < nrows; ++yy)
				for(int xx = 0; xx < ncols; ++xx)
					matrix(yy,xx)=(*zz_mat)(yy,xx).get_d();

			HamiltonianWrapper HW = HamiltonianWrapper(matrix*matrix.transpose(),"");
			gramian_wrappers.push_back(HW);

			//std::cerr<<"size: "<<HW.hamiltonian.rows()<<std::endl;
		}

		num_instances = lattices.size();

	}

	int nbQubits = -1;

	//logw("Saving eigensace which is not needed and very costly");
	//logw("!!!Random guess is now after the penalization!!!");

	for(int i = 0; i < num_instances; ++i){

		AngleResultsExperiment:Instance instance;

		Lattice l(gramian_wrappers[i].hamiltonian, gramian_wrappers[i].name);
		if(penalise)
			mapOptions->penalty = 5*l.getSquaredLengthOfFirstBasisVector(); //penalty set to length of first vector squared

//std::cerr<<gramian_wrappers[i].hamiltonian<<std::endl;
		instance.h = l.getHamiltonian(mapOptions);

		/*if(i < 10)
			std::cerr<<instance.h.getMatrixRepresentation2(true)<<std::endl;
		else
			throw;*/

		if(nbQubits < 0)
			nbQubits = instance.h.nbQubits;
		else if(nbQubits != instance.h.nbQubits)
			loge("Instances with different number of qubits found");

		qaoaOptions->accelerator->initialize(&instance.h);
		qaoaOptions->accelerator->options.createQuregAtEachInilization = false;
		long long int numAmpsTotal = qaoaOptions->accelerator->getQuregPtr()->numAmpsTotal;
		FastVQA::RefEnergies refEnergies = qaoaOptions->accelerator->getEigenspace();//delete
		qreal min = QREAL_MAX;//refEnergies[0].value;

		for(long long int j = 0; j < numAmpsTotal; ++j){

			/*if(i < 10)
						std::cerr<<refEnergies[j].index<<" "<<refEnergies[j].value<<std::endl;
					else
						throw;*/

			if(refEnergies[j].value == min)
				instance.sv_solutions.push_back(FastVQA::RefEnergy(min, refEnergies[j].index, false));
			else if(refEnergies[j].value > 0 && refEnergies[j].value < min){
				instance.sv_solutions.clear();
				min = refEnergies[j].value;
				instance.sv_solutions.push_back(FastVQA::RefEnergy(min, refEnergies[j].index, false));
			}
		}

		instance.zero_solutions = qaoaOptions->accelerator->getSolutions();
		if(instance.zero_solutions.size() > 1)
			loge("CmQaoaExperiment: Unimplemented, more than 1 solution marked");

		for(auto &sol: instance.zero_solutions){
			if(sol.value != 0)
				throw_runtime_error("CmQaoaExperiment: Something else than 0 state marked as a solution");
		}
		//long long int zero_index = instance.zero_solutions[0].index;

		/*qreal min_energy = instance.solutions[0].value;
		int num_sols_with_min_energy = 0;
		for(const auto &sol: instance.solutions){
			if(sol.value < min_energy){
				min_energy = sol.value;
				num_sols_with_min_energy = 1;
			}else if(sol.value == min_energy)
				num_sols_with_min_energy++;
		}
		instance.min_energy = min_energy;*/

		//THIS CHOOSES WHICH RANDOM GUESS IS BEING USED
		instance.random_guess = (qreal)(1./pow(2, nbQubits)) * instance.sv_solutions.size();
		//instance.random_guess = l.get_random_guess_one_vect() * instance.h.custom_solutions.size();

		//std::cerr<<nbQubits<<" "<<instance.solutions.size()<<" "<<instance.random_guess<<std::endl;

		instance.volume	= l.getVolume();
		instance.sv1Squared= l.getSquaredLengthOfFirstBasisVector();
		instance.q = q;
		instance.m = m;
		instance.n = n;

		/*if(i==0 || i == 1){

			std::cerr<<"solutions\n";
			for(const auto &sol: instance.sv_solutions){
				std::cerr<<sol.value<<" "<<sol.index<<std::endl;
			}
			std::cerr<<"zero solutions\n";
			for(const auto &sol: instance.zero_solutions){
				std::cerr<<sol.value<<" "<<sol.index<<std::endl;
			}
			std::cerr<<"SS: "<<instance.sv1Squared<<std::endl;
		}*/

		dataset.push_back(instance);
	}

	//Need to destroy qureg, because next time experiments with different number of qubits will be run
	qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
	qaoaOptions->accelerator->finalize();

	//logi("Experiment dataset generated", this->loglevel);
	return dataset;
}

void AngleResultsExperiment::run_qaoa_with_optimizer(){

	bool penalise;

	std::string meta_data;
	std::string alphas_output="alphas=";
	std::string z_alphas_output="zero_alphas=";


	std::string python_output="ms=[";
	for(int m = this->m_start; m <= this->m_end; ++m){
		if(m>this->m_start)
			python_output+=", ";
		python_output+=std::to_string(m);
	}
	python_output+="]\nys=";

	std::cerr<<"q="<<q<<" p="<<this->qaoaOptions->p<<" #="<<max_num_instances<<std::endl;
	//std::cerr<<"Will take around "<<30*max_num_instances<<" seconds"<<std::endl;

	this->mapOptions->penalty = 0;

	for(int index = 0; index < 2; ++index){

		if(index == 0){
			std::cerr<<std::endl<<"optCM-QAOA"<<std::endl;
			python_output+="{\"CMQAOA\": [";
			alphas_output+="{\"CMQAOA\": ";
			z_alphas_output+="{\"CMQAOA\": ";
			meta_data="optCMQAOA";
			penalise=false;
		}else if(index == 1){
			std::cerr<<std::endl<<"optQAOA penalty=0"<<std::endl;
			python_output+=", \"QAOA non_pen\": [";
			alphas_output+=", \"QAOA non_pen\": ";
			z_alphas_output+=", \"QAOA non_pen\": ";
			meta_data="optQAOAnonpen";
			penalise=false;
		}else{
			std::cerr<<std::endl<<"optQAOA penalised"<<std::endl;
			python_output+=", \"QAOA penalised\": [";
			alphas_output+=", \"QAOA penalised\": ";
			z_alphas_output+=", \"QAOA penalised\": ";
			meta_data="optQAOApen";
			penalise=true;
		}


		std::map<std::pair<int, int>, double> mean_map;
		std::map<std::pair<int, int>, double> stdev_map;
		std::map<std::pair<int, int>, double> num_sols_map;

		const int colWidth = 20;

		std::vector<int> n_list;

		n_list.push_back(3);
		loge("Changed n from 3 to 8, m_start from 4 to 9, m_end from 20 to 15!");
		//for(int i = 1; i < m_end; ++i)
		//	n_list.push_back(i);



		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout << std::setw(colWidth) << std::internal << "alpha";
		std::cout << std::setw(colWidth) << std::internal << "alpha_ext";
		std::cout<<std::endl;

		for(int n : n_list){
			std::cout << std::setw(3) << std::internal << n << "   ";
			double nom=0, den=0, z_nom=0, z_den=0;
			double sum_yi=0;
			double sum_xi=0;
			double sum_xi2=0;
			double sum_xi_yi=0;
			double z_sum_yi=0;
			double z_sum_xi=0;
			double z_sum_xi2=0;
			double z_sum_xi_yi=0;
			double counter=0;
			for(int m = this->m_start; m <= m_end; ++m){

			/*	if(m == 4 || m == 5 || m == 6){
					loge("Skipping m=4 or 5 or 6");
					continue;
				}*/

				logfile << m << std::endl << std::flush;
				//this->qaoaOptions->p = 6/*m*(8./3.)+((-26./3.))*/;///this->angles.size()/2;
				logi("p="+std::to_string(this->qaoaOptions->p), this->loglevel);

				if(n >= m){
					std::cout << std::setw(colWidth) << std::internal << "x";
					continue;
				}

				std::vector<AngleResultsExperiment::Instance> dataset = this->_generate_dataset(n, m, penalise);

				Cost cost;
				if(index == 0)
					cost = this->_cost_fn(&dataset, &this->angles_cmqaoa[0], meta_data, /*true*/false, this->seed);
				else if(index == 1)
					cost = this->_cost_fn(&dataset, &this->angles_optqaoa[0], meta_data, /*true*/false, this->seed);
				else
					throw_runtime_error("Not implemented");

				double mean = cost.mean;
				double stdev = cost.stdev;
				double mean_zero = cost.mean_zero;
				double num_sols = cost.mean_num_of_sols;

				double overlap = mean;
				nom += log2(overlap) * m;
				den += m * m;

				sum_yi+=log2(overlap);
				sum_xi+=m;
				sum_xi2+=m * m;
				sum_xi_yi+=log2(overlap)*m;
				counter++;

				double zero_overlap = mean_zero;
				z_nom += log2(zero_overlap) * m;
				z_den += m * m;

				z_sum_yi+=log2(zero_overlap);
				z_sum_xi+=m;
				z_sum_xi2+=m * m;
				z_sum_xi_yi+=log2(zero_overlap)*m;

				mean_map.emplace(std::pair<int, int>(n,m), mean);
				stdev_map.emplace(std::pair<int, int>(n,m), stdev);
				num_sols_map.emplace(std::pair<int, int>(n,m), num_sols);
				if(m>this->m_start)
					python_output+=", ";
				if(index == 0)
					python_output+="("+std::to_string(mean)+", "+std::to_string(stdev)+")";
				else
					python_output+="("+std::to_string(mean)+", "+std::to_string(mean_zero)+", "+std::to_string(stdev)+")";
				std::cout << std::setw(colWidth) << std::internal << mean<<"/"<<zero_overlap << "/"<<stdev/* << "/" << stdev */<< std::flush;
			}

			//double  alpha = -nom/den;
			std::cout << std::setw(colWidth) << std::internal << -nom/den << "/z="<< -z_nom/z_den << std::flush;

			double a=0,b=0,z_a=0,z_b=0;
			a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(counter*sum_xi2-sum_xi*sum_xi);
			b=(counter*sum_xi_yi-sum_xi*sum_yi)/(counter*sum_xi2-sum_xi*sum_xi);
			z_a=(z_sum_yi*z_sum_xi2-z_sum_xi*z_sum_xi_yi)/(counter*z_sum_xi2-z_sum_xi*z_sum_xi);
			z_b=(counter*z_sum_xi_yi-z_sum_xi*z_sum_yi)/(counter*z_sum_xi2-z_sum_xi*z_sum_xi);
			std::cout << std::setw(colWidth) << std::internal << "2^"<<a<<"+n*"<<b << " z="<< "2^"<<z_a<<"+n*"<<z_b <<std::flush;

			alphas_output += std::to_string(b);
			z_alphas_output += std::to_string(z_b);

			//std::cerr<<"Simple alpha: "<< -nom/den << std::endl;
			python_output+="]";
			std::cout<<std::endl;
		}

		std::cout<<std::endl<<std::endl;
		std::cout<<"   Standard deviations:"<<std::endl;
		std::cout << " n \\ m";
		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout<<std::endl;

		for(int n : n_list){
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
		std::cout<<std::endl<<std::endl;
		std::cout<<"   Num sols"<<std::endl;
		std::cout << " n \\ m";
		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout<<std::endl;

		for(int n : n_list){
			std::cout << std::setw(3) << std::internal << n << "   ";
			for(int m = this->m_start; m <= m_end; ++m){
				if(n >= m){
					std::cout << std::setw(colWidth) << std::internal << "x";
					continue;
				}
				std::cout << std::setw(colWidth) << std::internal << num_sols_map[std::pair<int, int>(n,m)];
			}
			std::cout<<std::endl;
		}
	}

	python_output+="}";
	alphas_output+="}";

	std::cerr<<std::endl<<python_output<<std::endl<<alphas_output<<std::endl;

}

void AngleResultsExperiment::run(){

	bool penalise;

	std::string meta_data;
	std::string alphas_output="alphas=";
	std::string z_alphas_output="zero_alphas=";

	assert(this->angles_optqaoa.size() == this->angles_cmqaoa.size());
	this->qaoaOptions->p = this->angles_optqaoa.size()/2;

	std::string python_output="ms=[";
	for(int m = this->m_start; m <= this->m_end; ++m){
		if(m>this->m_start)
			python_output+=", ";
		python_output+=std::to_string(m);
	}
	python_output+="]\nys=";

	std::cerr<<"q="<<q<<" p="<<this->qaoaOptions->p<<" #="<<max_num_instances<<std::endl;
	//std::cerr<<"Will take around "<<30*max_num_instances<<" seconds"<<std::endl;

	this->mapOptions->penalty = 0;

	for(int index = 0; index < 2; ++index){

		if(index == 0){
			std::cerr<<std::endl<<"CM-QAOA"<<std::endl;
			python_output+="{\"CMQAOA\": [";
			alphas_output+="{\"CMQAOA\": ";
			z_alphas_output+="{\"CMQAOA\": ";
			meta_data="";
			penalise=false;
		}else if(index == 1){
			std::cerr<<std::endl<<"QAOA penalty=0"<<std::endl;
			python_output+=", \"QAOA non_pen\": [";
			alphas_output+=", \"QAOA non_pen\": ";
			z_alphas_output+=", \"QAOA non_pen\": ";
			meta_data="QAOAnonpen";
			penalise=false;
		}else{
			std::cerr<<std::endl<<"QAOA penalised"<<std::endl;
			python_output+=", \"QAOA penalised\": [";
			alphas_output+=", \"QAOA penalised\": ";
			z_alphas_output+=", \"QAOA penalised\": ";
			meta_data="QAOApen";
			penalise=true;
		}


		std::map<std::pair<int, int>, double> mean_map;
		std::map<std::pair<int, int>, double> stdev_map;
		std::map<std::pair<int, int>, double> num_sols_map;

		const int colWidth = 20;

		std::vector<int> n_list;

		n_list.push_back(3);
		loge("Changed n from 3 to 8, m_start from 4 to 9, m_end from 20 to 15!");
		//for(int i = 1; i < m_end; ++i)
		//	n_list.push_back(i);



		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout << std::setw(colWidth) << std::internal << "alpha";
		std::cout << std::setw(colWidth) << std::internal << "alpha_ext";
		std::cout<<std::endl;

		for(int n : n_list){
			std::cout << std::setw(3) << std::internal << n << "   ";
			double nom=0, den=0, z_nom=0, z_den=0;
			double sum_yi=0;
			double sum_xi=0;
			double sum_xi2=0;
			double sum_xi_yi=0;
			double z_sum_yi=0;
			double z_sum_xi=0;
			double z_sum_xi2=0;
			double z_sum_xi_yi=0;
			double counter=0;
			for(int m = this->m_start; m <= m_end; ++m){
				if(n >= m){
					std::cout << std::setw(colWidth) << std::internal << "x";
					continue;
				}

				std::vector<AngleResultsExperiment::Instance> dataset = this->_generate_dataset(n, m, penalise);

				Cost cost;
				if(index == 0)
					cost = this->_cost_fn(&dataset, &this->angles_cmqaoa[0], meta_data, /*true*/false, this->seed);
				else if(index == 1)
					cost = this->_cost_fn(&dataset, &this->angles_optqaoa[0], meta_data, /*true*/false, this->seed);
				else
					throw_runtime_error("Not implemented");

				double mean = cost.mean;
				double stdev = cost.stdev;
				double mean_zero = cost.mean_zero;
				double num_sols = cost.mean_num_of_sols;

				double overlap = mean;
				nom += log2(overlap) * m;
				den += m * m;

				sum_yi+=log2(overlap);
				sum_xi+=m;
				sum_xi2+=m * m;
				sum_xi_yi+=log2(overlap)*m;
				counter++;

				double zero_overlap = mean_zero;
				z_nom += log2(zero_overlap) * m;
				z_den += m * m;

				z_sum_yi+=log2(zero_overlap);
				z_sum_xi+=m;
				z_sum_xi2+=m * m;
				z_sum_xi_yi+=log2(zero_overlap)*m;

				mean_map.emplace(std::pair<int, int>(n,m), mean);
				stdev_map.emplace(std::pair<int, int>(n,m), stdev);
				num_sols_map.emplace(std::pair<int, int>(n,m), num_sols);
				if(m>this->m_start)
					python_output+=", ";
				if(index == 0)
					python_output+="("+std::to_string(mean)+", "+std::to_string(stdev)+")";
				else
					python_output+="("+std::to_string(mean)+", "+std::to_string(mean_zero)+", "+std::to_string(stdev)+")";
				std::cout << std::setw(colWidth) << std::internal << mean/* << "/" << stdev */<< std::flush;
			}

			//double  alpha = -nom/den;
			std::cout << std::setw(colWidth) << std::internal << -nom/den << "/z="<< -z_nom/z_den << std::flush;

			double a=0,b=0,z_a=0,z_b=0;
			a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(counter*sum_xi2-sum_xi*sum_xi);
			b=(counter*sum_xi_yi-sum_xi*sum_yi)/(counter*sum_xi2-sum_xi*sum_xi);
			z_a=(z_sum_yi*z_sum_xi2-z_sum_xi*z_sum_xi_yi)/(counter*z_sum_xi2-z_sum_xi*z_sum_xi);
			z_b=(counter*z_sum_xi_yi-z_sum_xi*z_sum_yi)/(counter*z_sum_xi2-z_sum_xi*z_sum_xi);
			std::cout << std::setw(colWidth) << std::internal << "2^"<<a<<"+n*"<<b << " z="<< "2^"<<z_a<<"+n*"<<z_b <<std::flush;

			alphas_output += std::to_string(b);
			z_alphas_output += std::to_string(z_b);

			//std::cerr<<"Simple alpha: "<< -nom/den << std::endl;
			python_output+="]";
			std::cout<<std::endl;
		}

		std::cout<<std::endl<<std::endl;
		std::cout<<"   Standard deviations:"<<std::endl;
		std::cout << " n \\ m";
		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout<<std::endl;

		for(int n : n_list){
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
		std::cout<<std::endl<<std::endl;
		std::cout<<"   Num sols"<<std::endl;
		std::cout << " n \\ m";
		for(int m = this->m_start; m <= this->m_end; ++m)
			std::cout << std::setw(colWidth) << std::internal << m;
		std::cout<<std::endl;

		for(int n : n_list){
			std::cout << std::setw(3) << std::internal << n << "   ";
			for(int m = this->m_start; m <= m_end; ++m){
				if(n >= m){
					std::cout << std::setw(colWidth) << std::internal << "x";
					continue;
				}
				std::cout << std::setw(colWidth) << std::internal << num_sols_map[std::pair<int, int>(n,m)];
			}
			std::cout<<std::endl;
		}
	}
	python_output+="}";
	alphas_output+="}";

	std::cerr<<std::endl<<python_output<<std::endl<<alphas_output<<std::endl;

}

AngleSearchExperiment::AngleSearchExperiment(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	mapOptions->penalty = 0;

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
		//mapOptions->penalty = l.getSquaredLengthOfFirstBasisVector(); //penalty set to length of first vector squared
		instance.h = l.getHamiltonian(mapOptions);
		if(nbQubits < 0)
			nbQubits = instance.h.nbQubits;
		else if(nbQubits != instance.h.nbQubits)
			loge("Instances with different number of qubits found");
		qaoaOptions->accelerator->initialize(&instance.h);
		qaoaOptions->accelerator->options.createQuregAtEachInilization = false;
		long long int numAmpsTotal = qaoaOptions->accelerator->getQuregPtr()->numAmpsTotal;
		FastVQA::RefEnergies refEnergies = qaoaOptions->accelerator->getEigenspace();//delete
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

		instance.zero_solutions = qaoaOptions->accelerator->getSolutions();
		if(instance.zero_solutions.size() != 1)
			loge("CmQaoaExperiment: Unimplemented, more than 1 solution marked");

		for(auto &sol: instance.zero_solutions){
			if(sol.value != 0)
				throw_runtime_error("CmQaoaExperiment: Something else than 0 state marked as a solution");
		}
		long long int zero_index = instance.zero_solutions[0].index;

		/*qreal min_energy = instance.solutions[0].value;
		int num_sols_with_min_energy = 0;
		for(const auto &sol: instance.solutions){
			if(sol.value < min_energy){
				min_energy = sol.value;
				num_sols_with_min_energy = 1;
			}else if(sol.value == min_energy)
				num_sols_with_min_energy++;
		}
		instance.min_energy = min_energy;*/


		//THIS CHOOSES WHICH RANDOM GUESS IS BEING USED
		instance.random_guess = (qreal)(1./pow(2, nbQubits)) * instance.sv_solutions.size();
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

std::pair<double, double> AngleExperimentBase::try_many_starts(std::string meta_data, Instance* instance, FastVQA::Qaoa* qaoa_instance, int seed){

	const int num_starts = /*5000*/5000;//(instance->h.nbQubits * (163) - (633));/*20*/;


	FastVQA::ExperimentBuffer buffer;
	buffer.storeQuregPtr = true;
	if(!seeded){
		gen19937.seed(seed); //rd() instead of 0 - seed
		logi("Resetting seed to " + std::to_string(seed));
		seeded=true;
	}

	std::vector<double> best_params;

	double max_ground_state_overlap=0, zero_overlap_global=0;

	int i;
	for(i = 0; i < num_starts/*max_ground_state_overlap < 0.16*/; ++i){

		if(i % 10 == 0)
			logfile << "Angle set " << i << " / " << num_starts << std::endl << std::flush;

		this->qaoaOptions->initial_params.clear();
		std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
		for(int i = 0; i < this->qaoaOptions->p; ++i){
			double param1 = dis(gen19937);
			double param2 = dis(gen19937);
			this->qaoaOptions->initial_params.push_back(param1);
			this->qaoaOptions->initial_params.push_back(param2);
		}

		double ground_state_overlap = 0, zero_overlap = 0;
		if(meta_data == "optCMQAOA"){
			qaoa_instance->run_cm_qaoa(&buffer, &instance->h, this->qaoaOptions, instance->zero_solutions[0].index);
		}else if(meta_data == "optQAOAnonpen"){
			qaoa_instance->run_qaoa(&buffer, &instance->h, this->qaoaOptions);
		}

		for(auto &sol: instance->sv_solutions){
			long long int index   = sol.index;
			//std::cerr<<"y"<<sol.index<<" "<<sol.value<<std::endl;
			ground_state_overlap += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
		}

		for(auto &sol: instance->zero_solutions){
			long long int index   = sol.index;
			//std::cerr<<"x"<<sol.index<<" "<<buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index]<<std::endl;
			zero_overlap += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
		}

		//logfile << ground_state_overlap << std::endl << std::flush;


		if(max_ground_state_overlap < ground_state_overlap){
			max_ground_state_overlap = ground_state_overlap;
			zero_overlap_global = zero_overlap;
			best_params = this->qaoaOptions->initial_params;

			//logfile << ground_state_overlap << " after #=" << i << std::endl << std::flush;

		}

	}

	logfile << i << " start points were needed " << std::endl << std::flush;


/*	angleAnalysisLog << "[";
	int j = 0;
	for(auto &param: best_params){
		angleAnalysisLog << param;
		if(j++ < best_params.size()-1)
			angleAnalysisLog << ",";
	}
	angleAnalysisLog << "],"<<std::flush;*/


	return std::pair<double, double>(max_ground_state_overlap, zero_overlap_global);

}

AngleExperimentBase::Cost AngleExperimentBase::_cost_fn(std::vector<Instance>* dataset, const double *angles, std::string meta_data, bool use_database, int seed){

	std::vector<int> num_sols;


	int i = 0;

	double mean;
	double stdev;
	double mean_zero;
	FastVQA::Qaoa qaoa_instance;

	std::vector<double> gs_overlaps, zero_overlaps;

	/*loge("Overriding angles");
	bool print=false;//std::cerr<<angles[0] <<" " <<angles[1]<<"\n";
	if(angles[0] == 0.01 && angles[1]>0.3899 && angles[1]<0.3901){
		angles[0] = 2;angles[1] = 0;
		print=true;
	}*/

	//this->logfile << "Angles 1/"





	for(auto &instance: (*dataset)){

		if((meta_data == "optQAOAnonpen" || meta_data == "optCMQAOA") && i >= 5){
			loge("Breaking after 5 instances");
			break;
		}

		logfile << "i=" << i << std::endl << std::flush;

		if(instance.h.nbQubits != this->qaoaOptions->accelerator->getNumQubitsInQureg()){
			this->qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
			this->qaoaOptions->accelerator->finalize();
		}

		FastVQA::ExperimentBuffer buffer;
		buffer.storeQuregPtr = true;

		double ground_state_overlap = 0, zero_overlap = 0;
		if(use_database){throw;
			Database::DatasetRow output_row;
			this->database->getOrCalculate_qary_with_fixed_angles(&buffer, angles, 6, &instance.h, &output_row, this->qaoaOptions, &qaoa_instance);
			for(auto &sol: instance.sv_solutions){
				long long int index = sol.index;
				ground_state_overlap+=output_row.finalStateVectorMap[index].second;
			}
		}else{
			//qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);

			if(meta_data == "QAOAnonpen"){
				qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
			}else if(meta_data == "QAOApen"){
				qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
			}else if(meta_data == "optCMQAOA"){

				/*for(int angle_index = 0; angle_index < num_starts; ++angle_index){
					//this->qaoaOptions->initial_params.size();
					std::mt19937 gen(0); //rd() instead of 0 - seed
					std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
					for(int i = 0; i < this->qaoaOptions->p; ++i){
						double param1 = dis(gen);
						double param2 = dis(gen);
						this->qaoaOptions->initial_params.push_back(param1/2);
						this->qaoaOptions->initial_params.push_back(param2);
						//buffer->initial_params.push_back(std::pair<std::string, double>("alpha_"+std::to_string(i),param1));
						//buffer->initial_params.push_back(std::pair<std::string, double>("beta_"+std::to_string(i),param2));
					}
				}*/
				std::pair<double, double> res = this->try_many_starts(meta_data, &instance, &qaoa_instance, seed);
				ground_state_overlap = res.first;
				zero_overlap = res.second;

				qaoa_instance.run_cm_qaoa(&buffer, &instance.h, this->qaoaOptions, instance.zero_solutions[0].index);
			}else if(meta_data == "optQAOAnonpen"){
				std::pair<double, double> res = this->try_many_starts(meta_data, &instance, &qaoa_instance, seed);
				ground_state_overlap = res.first;
				zero_overlap = res.second;
				//qaoa_instance.run_qaoa(&buffer, &instance.h, this->qaoaOptions);
			}else if(meta_data == "optQAOApen"){
				throw_runtime_error("unimplemented");
			}
			//else CM_QAOA
			else
				qaoa_instance.run_cm_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles, instance.zero_solutions[0].index);
			//qaoa_instance.run_cm_qaoa(&buffer, &instance.h, this->qaoaOptions, instance.zero_solutions[0].index);

			/*for(auto &f: buffer.finalParams){
				std::cerr<<f<<" ";
			}std::cerr<<std::endl;*/


				/*for(int jj = 0; jj < 16; ++jj){
					long long int index   = jj;
					std::cerr<<jj<<" "<<buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index]<<std::endl;
				}
				std::cerr<<std::endl;*/
			if(meta_data != "optCMQAOA" && meta_data != "optQAOAnonpen"){
				for(auto &sol: instance.sv_solutions){
					long long int index   = sol.index;
					//std::cerr<<"y"<<sol.index<<" "<<sol.value<<std::endl;
					ground_state_overlap += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
				}

				for(auto &sol: instance.zero_solutions){
					long long int index   = sol.index;
					//std::cerr<<"x"<<sol.index<<" "<<buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index]<<std::endl;
					zero_overlap += buffer.stateVector->stateVec.real[index]*buffer.stateVector->stateVec.real[index]+buffer.stateVector->stateVec.imag[index]*buffer.stateVector->stateVec.imag[index];
				}
			}

		}

		/*if(num_sols.size() == 0){
			std::cerr<<"x"<<instance.h.custom_solutions.size()<<" "<<instance.sv_solutions.size()<<std::endl;
		}*/
		num_sols.push_back(instance.h.custom_solutions.size());
		qreal improvement_ratio = ground_state_overlap;//zero_overlap;/* / instance.random_guess*/;
		//std::cerr<<improvement_ratio<<std::endl;
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

		/*if(i == 0){
			std::cerr<<instance.h.getMatrixRepresentation2(true)<<std::endl;
			std::cerr<<"improvement ratio: "<< improvement_ratio<<std::endl<<std::endl;
			for(int j = 0; j < 12; ++j)
				std::cerr<<angles[j]<<" ";
			std::cerr<<std::endl;
		}loge("throw here");throw;*/

		gs_overlaps.push_back(improvement_ratio);
		zero_overlaps.push_back(zero_overlap/* / 	(qreal)(1./pow(2, instance.h.nbQubits)) * instance.zero_solutions.size()*/);

		//std::cerr<<ground_state_overlap<<" "<<instance.random_guess<<"\n";
		i++;

		this->qaoaOptions->accelerator->options.createQuregAtEachInilization = false;

	}
	/*std::cerr<<std::endl;std::cerr<<std::endl;
	for(auto &gg: gs_overlaps){
		std::cerr<<gg<<" ";
	}std::cerr<<std::endl;std::cerr<<std::endl;*/

	double sum = std::accumulate(gs_overlaps.begin(), gs_overlaps.end(), 0.0);
	mean = sum / gs_overlaps.size();

	long double accum = 0.0;
	std::for_each (std::begin(gs_overlaps), std::end(gs_overlaps), [&](const double d) {
	    accum += (d - mean) * (d - mean);
	});
	stdev = sqrt(accum / (gs_overlaps.size()-1));

	//double sq_sum = std::inner_product(gs_overlaps.begin(), gs_overlaps.end(), gs_overlaps.begin(), 0.0);
	//stdev = std::sqrt(sq_sum / gs_overlaps.size() - mean * mean);

	double sum_z = std::accumulate(zero_overlaps.begin(), zero_overlaps.end(), 0.0);
	mean_zero = sum_z / zero_overlaps.size();

	//mean = median(gs_overlaps, gs_overlaps.size());

	this->qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
	this->qaoaOptions->accelerator->finalize();


	double sum_mean_num_of_sols = std::accumulate(num_sols.begin(), num_sols.end(), 0.0);
	double mean_num_of_sols = sum_mean_num_of_sols / num_sols.size();

	//std::cerr<<mean<<"  "<<stdev<<std::endl;


	return Cost(mean, stdev, mean_zero, mean_num_of_sols);

}

void AngleSearchExperiment::run(){

	run_cobyla();
	//run_p6_test();
	//run_p2();
	return;

	if(this->qaoaOptions->p == 1)
		run_p1();
	else if(this->qaoaOptions->p == 2)
		run_cobyla();
		//run_p2_full_bruteforce();
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
							cost = this->_cost_fn(&this->train_set, angles, "");
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

	double mean_threshold = 1.2;
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
					cost = this->_cost_fn(&this->train_set, angles, "", false);
					//final_plot_means[x].push_back(cost.mean);
					//final_plot_stdevs[x].push_back(cost.stdev);

					if(cost.mean < debug_lowest_mean){
						debug_lowest_mean = cost.mean;
						lgamma1=gamma1;
						lbeta1=beta1;
						lgamma2=gamma2;
						lbeta2=beta2;
					}

					if(cost.mean / 1./(pow(2, m)) >= mean_threshold && cost.stdev <= stdev_threshold){
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


void AngleSearchExperiment::run_cobyla(){

	//if(this->num_params != 4)
	//	throw_runtime_error("Unimplemented for other depth than 2");

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
		return -this->_cost_fn(&this->train_set, &angles[0], "").mean;
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

	double mean_threshold = 1.232;//1.2;
	double stdev_threshold = 2000000;//0.019;

	double beta_min = 0;//pi/16;
	double beta_max = pi;//pi;//;pi/8;

	double gamma_min = 0;//0;
	double gamma_max = pi;//2*pi;

	double range_beta = beta_max - beta_min;
	double range_gamma = gamma_max - gamma_min;
	double incr_beta = 0.01;///0.12;
	double incr_gamma = 0.01;//0.04;///0.01;

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

	double debug_lowest_mean=-1;
	double lgamma, lbeta;

	double i = 0;
	AngleSearchExperiment::Cost cost;
	for(double gamma = gamma_min; gamma < gamma_max; gamma+=incr_gamma){
		for(double beta = beta_min; beta < beta_max; beta+=incr_beta){
			bar.tick();
			angles[0] = gamma;
			angles[1] = beta;

			cost = this->_cost_fn(&this->train_set, angles, "", false);
			logd("NOT USING DATABASE", this->loglevel);
			final_plot_means[x].push_back(cost.mean);
			final_plot_stdevs[x].push_back(cost.stdev);

			//loge("ahoj");

			if(cost.mean > debug_lowest_mean){
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

1.24181 2.78 0.79

 *
 */

	std::vector<std::tuple<double, double>> first_round_angles;
	first_round_angles.push_back(std::tuple<double, double>(2.78, 0.79));
	/*first_round_angles.push_back(std::tuple<double, double>(0.4, 0.48));
	first_round_angles.push_back(std::tuple<double, double>(0.4, 0.6));*/


	double *angles = (double*) malloc(4 * sizeof(double));
	for(int j = 0; j < this->num_params; ++j){
		angles[j]=0;//dis(gen));
	}

		double mean_threshold = 1.26;
		double stdev_threshold = 1000;
		double beta_min = 0;
		double beta_max = pi/2;

		double gamma_min = 0;
		double gamma_max = 2*pi;

		double range_beta = beta_max - beta_min;
		double range_gamma = gamma_max - gamma_min;
		double incr_beta = 0.01;//0.12;
		double incr_gamma = 0.01;//0.01;

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
		double debug_lowest_mean = -1, lbeta1, lbeta2, lgamma1, lgamma2;
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
					cost = this->_cost_fn(&this->train_set, angles, "");
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

				if(cost.mean > debug_lowest_mean){
					debug_lowest_mean = cost.mean;
					lgamma1=angles[0];
					lgamma2=angles[2];
					lbeta1=angles[1];
					lbeta2=angles[3];
				}

			}
			x++;
		}

		std::cerr<<debug_lowest_mean <<" "<<lgamma1<<" "<<lbeta1<<" "<<lgamma2<<" "<<lbeta2<<std::endl;


			/*double mx, my;int  mb;
			 * while(true){
				gp.getMouse(mx, my, mb, "");
				printf("You pressed mouse button %d at x=%f y=%f\n", mb, mx, my);
				//printf("%f", final_plot[(int)mx][(int)my]);
			}*/


	free(angles);

}

void AngleSearchExperiment::run_p6_test(){

	struct Angl{
		double g1,b1,g2,b2,g3,b3,g4,b4,g5,b5,g6,b6;
		double train_mean;
		double train_std;
		Angl(double g1, double b1, double g2, double b2, double g3,double b3, double g4, double b4, double g5, double b5, double g6,double b6, double train_mean, double train_std){
			this->g1=g1;this->b1=b1;
			this->g2=g2;this->b2=b2;
			this->g3=g3;this->b3=b3;
			this->g4=g4;this->b4=b4;
			this->g5=g5;this->b5=b5;
			this->g6=g6;this->b6=b6;

			this->train_mean=train_mean;this->train_std=train_std;
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
	second_round_angles.push_back(Angl(2.2287, 2.13607 ,2.23057 ,2.17114, 0.736168, -0.775112, -1.17933, -2.84297, -1.32729 ,-0.37258 ,1.90813 ,-0.247229, 1.16567, 0));

	for(auto &setting : second_round_angles){
		double angles[12];
		angles[0] = setting.g1;
		angles[1] = setting.b1;
		angles[2] = setting.g2;
		angles[3] = setting.b2;
		angles[4] = setting.g3;
		angles[5] = setting.b3;
		angles[6] = setting.g4;
		angles[7] = setting.b4;
		angles[8] = setting.g5;
		angles[9] = setting.b5;
		angles[10] = setting.g6;
		angles[11] = setting.b6;
		AngleSearchExperiment::Cost cost = this->_cost_fn(&this->test_set, angles, "");

		std::cerr<<"Train/Test mean="<<setting.train_mean<<"/"<<cost.mean<<"   std="<<setting.train_std<<"/"<<cost.stdev<<"\n";

		//loge("Break from second_roung_angles");
		//break;

	}


}


inline double AlphaMinimizationExperiment::strategy_alpha_c(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset, std::vector<double> angles, std::string meta_data){

	double nom=0, den=0;

	/*for(auto &dim: train_dataset){
		double overlap = this->_cost_fn(dim, &angles[0]);
		//std::cerr<<overlap<<std::endl;
		nom += log2(overlap) * dim[0].m;
		den += dim[0].m * dim[0].m;
		//overlaps.push_back(this->_cost_fn(dim, &angles[0]));
	}
	double alpha = -nom/den;
	final_ab.first=0;
	final_ab.second=alpha;*/


	//std::vector<double> overlaps;


	double sum_yi=0;
	double sum_xi=0;
	double sum_xi2=0;
	double sum_xi_yi=0;
	double alpha_calc_dataset_size=0;
	for(auto &dim: train_dataset){
		//overlaps.push_back(log2(overlap));
		//if(dim[0].m > 14){
			double overlap = this->_cost_fn(dim, &angles[0], meta_data);
			sum_yi+=log2(overlap);
			sum_xi+=dim[0].m;
			sum_xi2+=dim[0].m * dim[0].m;
			sum_xi_yi+=log2(overlap)*dim[0].m;
			alpha_calc_dataset_size++;
		//}
 	}
	double a=0,b=0;
	a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(alpha_calc_dataset_size*sum_xi2-sum_xi*sum_xi);
	b=(alpha_calc_dataset_size*sum_xi_yi-sum_xi*sum_yi)/(alpha_calc_dataset_size*sum_xi2-sum_xi*sum_xi);
	std::cerr<<"2^"<<a<<"+n*"<<b<<std::endl;

	double alpha = -b;
	//final_ab.first = a;
	//final_ab.second = a;

	return /*alpha*/a;
}

inline double AlphaMinimizationExperiment::strategy_inv_diff(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset, std::vector<double> angles, std::string meta_data){

	//double alpha = strategy_alpha_c(train_dataset, angles, meta_data);

	double den=0;
	for(auto &dim: train_dataset){

		double overlap = this->_cost_fn(dim, &angles[0], meta_data);
		double diff = overlap * pow(2, dim[0].m)-1;//overlap / pow(2, -dim[0].m);
		//if(diff < 0)
		//	diff = 0;
		den += diff * diff;
	}

	double res = 1./den;
	std::cerr << res << std::endl;
	return res;

}


AlphaMinimizationExperiment::AlphaMinimizationExperiment(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){
	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->p = qaoaOptions->p;
	this->mapOptions->penalty = 0;
}


void AlphaMinimizationExperiment::run(){

	//std::string meta_data = "fixedQAOA";
	std::string meta_data = "fixedCMQAOA";

	std::cerr<<meta_data<<std::endl;

	//int p = /*6*/7;

	int q = 97;
	int n = 3;
	int m_start = 4;
	int m_end = 10;

	bool new_way=/*false*/true;
	loge("new_way="+std::to_string(new_way));

	int max_num_instances = 100;//1000;
	double test_ratio = 0;//.2;

	int num_params = this->p*2;

	this->qaoaOptions->p = this->p;
	logi("p="+std::to_string(this->qaoaOptions->p), this->loglevel);

	std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset;
	std::vector<std::vector<AlphaMinimizationExperimentInstance>> test_dataset;

	for(int m = m_start; m <= m_end; ++m){

		std::vector<AlphaMinimizationExperimentInstance> m_train_instances;
		std::vector<AlphaMinimizationExperimentInstance> m_test_instances;

		long long int num_instances;


		std::vector<HamiltonianWrapper> gramian_wrappers;
		if(new_way){/*
			GeneratorParam param(q, n, m, true, 97, num_instances); //q, n, m, shuffle, seed, cutoff
			gramian_wrappers = generateQaryUniform(param);

			num_instances = q;
			for(int i = 0; i < (m-n); ++i){ //pow
				num_instances *= q;
				if(max_num_instances < num_instances){
					num_instances = max_num_instances;
					break;
				}
			}

			if(max_num_instances < num_instances)
				num_instances = max_num_instances;
				*/
			GeneratorParam param(q, n, m, true, 97, 100); //q, n, m, shuffle, seed, cutoff
			gramian_wrappers = generateFromEvalDecomposition(param);//generateQaryUniform(param);

					/*num_instances = q;
					for(int i = 0; i < (m-n); ++i){ //pow
						num_instances *= q;
						if(max_num_instances < num_instances){
							num_instances = max_num_instances;
							break;
						}
					}

					if(max_num_instances < num_instances)
						num_instances = max_num_instances;*/

					num_instances = gramian_wrappers.size();

		}
		else{
				std::vector<Lattice*> lattices;
				initialize_paper_experiment("qary_4_20", lattices, m, -1);
				int num_lattices = lattices.size();
				//logi("Dataset qary_4_20 succesfully loaded", 0);
				for(auto &l: lattices){

					fplll::ZZ_mat<mpz_t>* zz_mat = l->get_current_lattice();
					int nrows = l->n_rows;
					int ncols = l->n_cols;

					Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix(nrows, ncols);

					for(int yy = 0; yy < nrows; ++yy)
						for(int xx = 0; xx < ncols; ++xx)
							matrix(yy,xx)=(*zz_mat)(yy,xx).get_d();

					HamiltonianWrapper HW = HamiltonianWrapper(matrix*matrix.transpose(),"");
					gramian_wrappers.push_back(HW);

					//std::cerr<<"size: "<<HW.hamiltonian.rows()<<std::endl;
				}

				num_instances = lattices.size();

			}
		int num_test_instances = num_instances * test_ratio;
		int num_train_instances = num_instances - num_test_instances;
		logi("m="+std::to_string(m)+" Num_instances="+std::to_string(num_instances)+" Test ratio="+std::to_string(test_ratio)+" "+"train_size="+std::to_string(num_train_instances)+" "+"test_size="+std::to_string(num_test_instances), this->loglevel);

		if(num_test_instances == 0)
			loge("Zero test instances!");


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
				loge("CmQaoaExperiment: Unimplemented, more than 1 solution marked");

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

	this->qaoaOptions->ftol = 1e-10;
	this->qaoaOptions->max_iters = 1000;

	ProgressBar bar{
		option::BarWidth{50},
		option::MaxProgress{this->qaoaOptions->max_iters},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{"Running Angle Search Experiment p="+std::to_string(p)+" with COBYLA"},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::ForegroundColor{Color::yellow},
		option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};

	unsigned int iteration_i = 0;
loge("cost value is sum_yi not alpha!");
	//std::pair<double, double> final_ab;
	FastVQA::OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
		iteration_i++;
		bar.tick();
		std::vector<double> angles(x);
		//std::cerr<<this->_cost_fn(&this->train_set, &angles[0]).mean<<std::endl;







		return strategy_inv_diff(train_dataset, angles, meta_data);
		//return strategy_alpha_c(train_dataset, angles, meta_data);








		//double cost = this->_cost_fn(train_dataset[train_dataset.size()-1], &angles[0], meta_data);
		//std::cerr<<cost<<std::endl;
		//return -cost;


		//std::cerr<<nom<<" "<<den<<" "<<alpha<<std::endl;



		//return alpha;//-sum_yi/alpha_calc_dataset_size;//alpha;//-this->_cost_fn(&this->train_set, &angles[0]).mean;
	}, num_params);

	std::vector<double> initial_params;
	std::mt19937 gen(0); //rd() instead of 0 - seed
	std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
	for(int i = 0; i < num_params/2; ++i){
		double param1 = /*pi/4.;*/dis(gen);
		double param2 = /*pi/8.;*/dis(gen);
		std::cerr<<param1<<" "<<param2<<std::endl;
		initial_params.push_back(param1);
		initial_params.push_back(param2);
	}

	std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
	std::vector<double> upperBounds(initial_params.size(), 3.141592654);

	logd("QAOA starting optimization", this->loglevel);
	FastVQA::OptResult result = this->qaoaOptions->optimizer->optimize(f, initial_params, this->qaoaOptions->ftol, this->qaoaOptions->max_iters, lowerBounds, upperBounds);
	logd("QAOA finishing optimization", this->loglevel);

	//std::cerr<<".   2^"<<result.first.first<<"n+"<<result.second<<std::endl;
	std::cerr<<"alpha: "<< result.first.first <<"\n";
	std::cerr<<"num_iters: "<<iteration_i<<std::endl;
	if(meta_data == "fixedQAOA")
		std::cerr<<"const std::vector<double> angles_optqaoa{";
	else if(meta_data == "fixedCMQAOA")
		std::cerr<<"const std::vector<double> angles_cmqaoa{";
	else
		throw_runtime_error("Not implemented");
	bool start=true;
	for(auto &a: result.first.second){
		if(!start)
			std::cerr<<", ";
		std::cerr<<std::setprecision (15)<<a;
		start=false;
	}std::cerr<<"};";
	std::cerr<<"\n"<<nlopt_res_to_str(result.second)<<std::endl;

	//EVALUATE TRAIN DATASET
	std::vector<double> final_angles(result.first.second);

	double sum_yi=0;
	double sum_xi=0;
	double sum_xi2=0;
	double sum_xi_yi=0;
	for(auto &dim: train_dataset){
		double overlap = this->_cost_fn(dim, &final_angles[0], meta_data);
		//overlaps.push_back(log2(overlap));
		sum_yi+=log2(overlap);
		sum_xi+=dim[0].m;
		sum_xi2+=dim[0].m * dim[0].m;
		sum_xi_yi+=log2(overlap)*dim[0].m;
	}
	double a=0,b=0;

	a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(train_dataset.size()*sum_xi2-sum_xi*sum_xi);
	b=(train_dataset.size()*sum_xi_yi-sum_xi*sum_yi)/(train_dataset.size()*sum_xi2-sum_xi*sum_xi);
	//std::cerr<<pow(2.71828, a)<<"2^n*"<<b<<std::endl;
	std::cerr<<"2^"<<a<<"+n*"<<b<<std::endl;
	double alpha = -b;

	std::cerr<<"Test alpha: "<<alpha<<std::endl;

}
double AlphaMinimizationExperiment::_cost_fn(std::vector<AlphaMinimizationExperimentInstance> dataset, const double *angles, std::string meta_data, bool use_database){
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

			/*if(i >= 5)
				continue;*/

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
				if(meta_data == "fixedQAOA")
					qaoa_instance.run_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles);
				else if(meta_data == "fixedCMQAOA")
					qaoa_instance.run_cm_qaoa_fixed_angles(&buffer, &instance.h, this->qaoaOptions, angles, instance.zero_solutions[0].index);
				else
					throw_runtime_error("Not implemented");
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

			/*if(i == 0){
				std::cerr<<instance.h.getMatrixRepresentation2(true)<<std::endl;
				std::cerr<<"improvement ratio: "<< improvement_ratio<<std::endl<<std::endl;
				for(int j = 0; j < 12; ++j)
					std::cerr<<angles[j]<<" ";
				std::cerr<<std::endl;
			}*/

			gs_overlaps.push_back(improvement_ratio);
			//std::cerr<<ground_state_overlap<<" "<<instance.random_guess<<"\n";
			i++;

			this->qaoaOptions->accelerator->options.createQuregAtEachInilization = false;

		}

		/*std::cerr<<std::endl;std::cerr<<std::endl;
		for(auto &gg: gs_overlaps){
			std::cerr<<gg<<" ";
		}std::cerr<<std::endl;std::cerr<<std::endl;*/

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
