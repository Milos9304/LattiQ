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
#include <sstream>

//delete
//#include <chrono>


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


template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 15)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

/*CmQaoaExperiment::CmQaoaExperiment(FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database, int loglevel){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->database = database;
}*/

/*void CmQaoaExperiment::run(){

	int nbQubits_acc = -1;
	this->mapOptions->penalty = 0;

	FastVQA::Qaoa qaoa_instance;

	std::vector<std::vector<std::tuple<int , double, double, double, double>>> plot(this->p_end - this->p_start + 1);

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
}*/


AngleResultsExperiment::AngleResultsExperiment(int loglevel, int m_start, int m_end, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database, int seed, bool use_database_to_load_dataset){

	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->database = database;

	logi("p="+std::to_string(this->qaoaOptions->p), this->loglevel);

	this->logfile.open("log.txt");
	this->angleAnalysisLog.open("angleAnalysis.txt");

	this->m_start = m_start;
	this->m_end = m_end;
	this->seed = seed;

	this->use_database_to_load_dataset = use_database_to_load_dataset;

	//if(this->qaoaOptions->p != 2)
	//	throw_runtime_error("Angleres is only for p=2!");

}

std::vector<AngleExperimentBase::Instance> AngleExperimentBase::_generate_dataset(int n, int m, bool penalise){

	//logfile << "dataset gen start" <<std::endl<<std::flush;

	std::vector<AngleExperimentBase::Instance> dataset;

	bool new_way=true;
	long long int num_instances;

	/*std::vector<Lattice*> lattices;
	initialize_paper_experiment("qary_4_20", lattices, m, -1);
	int num_lattices = lattices.size();
	logi("Dataset qary_4_20 succesfully loaded", 0);
*/
	std::vector<HamiltonianWrapper> gramian_wrappers;

	if(new_way){
		GeneratorParam param(this->q, n, m, true, 97, this->max_num_instances); //q, n, m, shuffle, seed, cutoff
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

//auto t_start = std::chrono::high_resolution_clock::now();


//std::cerr<<gramian_wrappers[i].hamiltonian<<std::endl;

		l.setSolutionsToZeroVector();
		instance.h = l.getHamiltonian(mapOptions);
//auto t_end = std::chrono::high_resolution_clock::now();
//double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
//logfile << "del include, " << elapsed_time_ms <<std::endl<<std::flush;


		/*if(i < 10)
			std::cerr<<instance.h.getMatrixRepresentation2(true)<<std::endl;
		else
			throw;*/

		if(nbQubits < 0)
			nbQubits = instance.h.nbQubits;
		else if(nbQubits != instance.h.nbQubits)
			loge("Instances with different number of qubits found");

		if(this->use_database_to_load_dataset){

			FastVQA::Accelerator::DiagonalOpDuplicate diagonalOpDuplicate;
			std::vector<long double> real;
			bool found = this->database->getDataset(nbQubits, i, this->mapOptions->num_qbits_per_x, &real);

			if(found){
				diagonalOpDuplicate.numQubits = nbQubits;
				diagonalOpDuplicate.real = real;
				qaoaOptions->accelerator->initialize(&instance.h, true, &diagonalOpDuplicate);
			}else{
				diagonalOpDuplicate.numQubits = -1; //flag that not found
				qaoaOptions->accelerator->initialize(&instance.h, true, &diagonalOpDuplicate);

				this->database->insertDataset(m, i, mapOptions->num_qbits_per_x, diagonalOpDuplicate.real);

			}

			instance.diagOpDuplicate = diagonalOpDuplicate;

		}else{
			qaoaOptions->accelerator->initialize(&instance.h, false, nullptr);
		}

		//t_end = std::chrono::high_resolution_clock::now();
		//elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		//logfile << "del include2, " << elapsed_time_ms <<std::endl<<std::flush;

		qaoaOptions->accelerator->options.createQuregAtEachInilization = false;
		long long int numAmpsTotal = qaoaOptions->accelerator->getQuregPtr()->numAmpsTotal;
		FastVQA::RefEnergies refEnergies = qaoaOptions->accelerator->getEigenspace();//delete
		qreal min = QREAL_MAX;//refEnergies[0].value;

		//t_end = std::chrono::high_resolution_clock::now();
		//elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		//logfile << "del include3, " << elapsed_time_ms <<std::endl<<std::flush;

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

		//t_end = std::chrono::high_resolution_clock::now();
		//elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		//logfile << "del include4, " << elapsed_time_ms <<std::endl<<std::flush;

		instance.zero_solutions = qaoaOptions->accelerator->getSolutions();
		if(instance.zero_solutions.size() > 1){
			loge("CmQaoaExperiment: Unimplemented, more than 1 solution marked");
		}

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
		//t_end = std::chrono::high_resolution_clock::now();
		//elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		//logfile << "del include10, " << elapsed_time_ms <<std::endl<<std::flush;

	}

	//Need to destroy qureg, because next time experiments with different number of qubits will be run
	qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
	qaoaOptions->accelerator->finalize();

	//logi("Experiment dataset generated", this->loglevel);
	//logfile << "dataset gen complete" <<std::endl<<std::flush;
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
					cost = this->_cost_fn(&dataset, {}, meta_data, /*true*/false, this->seed);
				else if(index == 1)
					cost = this->_cost_fn(&dataset, {}, meta_data, /*true*/false, this->seed);
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

	if(optAngles[this->qaoaOptions->p].cm_meta_data == "nan" || optAngles[this->qaoaOptions->p].qaoa_meta_data == "nan")
		throw_runtime_error("Input in opt angles is missing");

	std::vector<double> angles_cmqaoa = optAngles[this->qaoaOptions->p].cm_angles;

	std::vector<double> angles_optqaoa = optAngles[this->qaoaOptions->p].qaoa_angles;

	assert(angles_optqaoa.size() == angles_cmqaoa.size());
	assert(this->qaoaOptions->p == angles_optqaoa.size()/2);

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
			std::cerr<<std::endl<<"  CM-QAOA"<<std::endl<<std::endl;
			python_output+="{\"CMQAOA\": [";
			alphas_output+="{\"CMQAOA\": ";
			z_alphas_output+="{\"CMQAOA\": ";
			meta_data="";
			penalise=false;
		}else if(index == 1){
			std::cerr<<std::endl<<"  QAOA penalty=0"<<std::endl<<std::endl;
			python_output+=", \"QAOA non_pen\": [";
			alphas_output+=", \"QAOA non_pen\": ";
			z_alphas_output+=", \"QAOA non_pen\": ";
			meta_data="QAOAnonpen";
			penalise=false;
		}else{
			std::cerr<<std::endl<<"  QAOA penalised"<<std::endl<<std::endl;
			python_output+=", \"QAOA penalised\": [";
			alphas_output+=", \"QAOA penalised\": ";
			z_alphas_output+=", \"QAOA penalised\": ";
			meta_data="QAOApen";
			penalise=true;
		}


		std::map<int, double> mean_map;
		std::map<int, double> stdev_map;
		std::map<int, double> num_sols_map;

		const int colWidth = 20;

		//std::cout << std::setw(colWidth) << std::internal << "alpha";
		//std::cout << std::setw(colWidth) << std::internal << "alpha_ext";
		//std::cout<<std::endl;

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

			std::vector<AngleResultsExperiment::Instance> dataset = this->_generate_dataset(3, m, penalise); //3 is ignored

			Cost cost;
			if(index == 0)
				cost = this->_cost_fn(&dataset, &angles_cmqaoa[0], meta_data, /*true*/false, this->seed);
			else if(index == 1)
				cost = this->_cost_fn(&dataset, &angles_optqaoa[0], meta_data, /*true*/false, this->seed);
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

			mean_map.emplace(m, mean);
			stdev_map.emplace(m, stdev);
			num_sols_map.emplace(m, num_sols);
			if(m > this->m_start)
				python_output+=", ";
			if(index == 0)
				python_output+="("+to_string_with_precision(mean)+", "+to_string_with_precision(stdev)+")";
			else
				python_output+="("+to_string_with_precision(mean)+", "+to_string_with_precision(mean_zero)+", "+to_string_with_precision(stdev)+")";

			std::stringstream ss;
			ss << m <<":   " << std::fixed << std::setprecision(15) << mean << "/" << stdev << "/" << zero_overlap << std::endl;
			logi(ss.str());

		}

		//double  alpha = -nom/den;
		std::cout << "alpha:" << std::setw(colWidth-4) << std::internal << nom/den << "/alpha_zero="<< z_nom/z_den << std::endl;

		double a=0,b=0,z_a=0,z_b=0;
		a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(counter*sum_xi2-sum_xi*sum_xi);
		b=(counter*sum_xi_yi-sum_xi*sum_yi)/(counter*sum_xi2-sum_xi*sum_xi);
		z_a=(z_sum_yi*z_sum_xi2-z_sum_xi*z_sum_xi_yi)/(counter*z_sum_xi2-z_sum_xi*z_sum_xi);
		z_b=(counter*z_sum_xi_yi-z_sum_xi*z_sum_yi)/(counter*z_sum_xi2-z_sum_xi*z_sum_xi);

		std::cout << "alpha_ext:" << std::setw(colWidth-14) << std::internal << "2^"<<a<<"+n*"<<b<< " alpha_ext_zero="<< "2^"<<z_a<<"+n*"<<z_b <<std::endl;

		alphas_output += std::to_string(b);
		z_alphas_output += std::to_string(z_b);

		//std::cerr<<"Simple alpha: "<< -nom/den << std::endl;
		python_output+="]";
		std::cout<<std::endl;

		std::cout<<std::endl<<std::endl;
		std::cout<<"   Num sols"<<std::endl;

		for(int m = this->m_start; m <= m_end; ++m){
			std::cout << std::setw(colWidth) << std::internal << to_string_with_precision(num_sols_map[m], 2);
		}
		std::cout<<std::endl;

	}
	python_output+="}";
	alphas_output+="}";

	std::cerr<<std::endl<<python_output<<std::endl<<alphas_output<<std::endl;

}

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

	logfile << i << " random start points were used" << std::endl << std::flush;


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

		if(instance.diagOpDuplicate.numQubits != -1)
			this->qaoaOptions->diagOpDuplicatePtr = &instance.diagOpDuplicate;

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
		this->qaoaOptions->diagOpDuplicatePtr = nullptr;


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


inline double AlphaMinimizationExperiment::strategy_alpha_c(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset,
		std::vector<double> angles, std::string meta_data, std::string* optimized_by){


	bool trivial_alpha=true;

	if(trivial_alpha){

		*optimized_by = "strategy_alpha_trivial";


		double nom = 0, den = 0;

		for(auto &dim: train_dataset){
			//overlaps.push_back(log2(overlap));
			//if(dim[0].m > 14){
				double overlap = this->_cost_fn(dim, &angles[0], meta_data);
				nom += log2(overlap) * dim[0].m;
				den += dim[0].m * dim[0].m;
			//}
		}

		double alpha = -nom/den;
		//final_ab.first = a;
		//final_ab.second = a;

		return alpha;



	}else{

		*optimized_by = "strategy_alpha_c";



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
		//std::cerr<<"2^"<<a<<"+n*"<<b<<std::endl;

		double alpha = -b;
		//final_ab.first = a;
		//final_ab.second = a;

		return alpha;
	}
}

inline double AlphaMinimizationExperiment::strategy_random_alpha_c(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset,
		std::vector<double> angles, std::string meta_data, std::string* optimized_by){

	*optimized_by = "strategy_random_alpha";


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


	int probability_dim = 80;
	int probability_instance = 80;

	int i = 0;
	double sum_yi=0;
	double sum_xi=0;
	double sum_xi2=0;
	double sum_xi_yi=0;
	double alpha_calc_dataset_size=0;

	for(auto &dim: train_dataset){
		i++;
		int p_num = rand() % 100 + 1;  //Generate random number 1 to 100
		if (p_num > probability_dim && (alpha_calc_dataset_size >= 2 || i < train_dataset.size()-1))
			continue;

		double overlap = this->_cost_fn(dim, &angles[0], meta_data, probability_instance);
		//std::cerr<<overlap<<std::endl;
		nom += log2(overlap) * dim[0].m;
		den += dim[0].m * dim[0].m;
		//overlaps.push_back(this->_cost_fn(dim, &angles[0]));
	}
		double alpha = -nom/den;
		std::cerr<<"2^n*-"<<alpha<<std::endl;

		//final_ab.first=0;
		//final_ab.second=alpha;

		/*		for(auto &dim: train_dataset){

					i++;

					int p_num = rand() % 100 + 1;  //Generate random number 1 to 100
					if (p_num > probability_dim && (alpha_calc_dataset_size >= 2 || i < train_dataset.size()-1))
						continue;

					double overlap = this->_cost_fn(dim, &angles[0], meta_data, probability_instance);

					//overlaps.push_back(log2(overlap));
					//if(dim[0].m > 14){
						sum_yi+=log2(overlap);
						sum_xi+=dim[0].m;
						sum_xi2+=dim[0].m * dim[0].m;
						sum_xi_yi+=log2(overlap)*dim[0].m;
						alpha_calc_dataset_size++;
					//}
				}
				double a=0,b=0;
				//std::cerr<<alpha_calc_dataset_size<<"   "<<alpha_calc_dataset_size*sum_xi2-sum_xi*sum_xi<<std::endl;
				a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(alpha_calc_dataset_size*sum_xi2-sum_xi*sum_xi);
				b=(alpha_calc_dataset_size*sum_xi_yi-sum_xi*sum_yi)/(alpha_calc_dataset_size*sum_xi2-sum_xi*sum_xi);
				std::cerr<<"2^"<<a<<"+n*"<<b<<std::endl;

				double alpha = -b;
				//final_ab.first = a;
				//final_ab.second = a;

		*/

	return alpha;
}

inline double AlphaMinimizationExperiment::strategy_random_inv_diff(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset,
		std::vector<double> angles, std::string meta_data, std::string* optimized_by){

	*optimized_by = "strategy_random_inv_diff";

	//double alpha = strategy_alpha_c(train_dataset, angles, meta_data);

	//else
	//	tails++; //This is tail


	int probability_dim = 50;

	int probability_instance = 50;

	double res = 0;
	int num_dims = 0;
	double sum_dims = 0;
	int i = 0;
	for(auto &dim: train_dataset){

			i++;

			int p_num = rand() % 100 + 1;  //Generate random number 1 to 100
			if (p_num > probability_dim && (num_dims > 0 || i < train_dataset.size()-1))
				continue;

			double overlap = this->_cost_fn(dim, &angles[0], meta_data, probability_instance);
			double diff = /*pow(2, dim[0].m)*/2-overlap;//overlap / pow(2, -dim[0].m);
			//if(diff < 0)
			//	diff = 0;
			res += diff * diff;
			num_dims++;
			sum_dims += pow(2, dim[0].m);
	}
	res*=0.25/ /*sum_dims*/num_dims;
	std::cerr << res << std::endl;
	return res;

}

inline double AlphaMinimizationExperiment::strategy_inv_diff(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset,
		std::vector<double> angles, std::string meta_data, std::string* optimized_by){

	*optimized_by = "strategy_inv_diff";

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
	//std::cerr << res << std::endl;
	return res;

}


AlphaMinimizationExperiment::AlphaMinimizationExperiment(int loglevel, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database, int seed){
	this->loglevel = loglevel;
	this->qaoaOptions = qaoaOptions;
	this->mapOptions = mapOptions;
	this->p = qaoaOptions->p;
	this->mapOptions->penalty = 0;
	this->database=database;
	this->seed = seed;

	srand(0);

}


void AlphaMinimizationExperiment::run(bool use_database_to_load_dataset){

	bool append_previous_angles = false; //initialize with prev angles padded with 2 zeros

	this->qaoaOptions->ftol = 1e-12;
	this->qaoaOptions->max_iters = 2000; //1000

	std::string meta_data;
	std::stringstream output;

	std::string space = "				";
	output << "		optAngle(\n"<<space<<this->p<<",\n"<<space<<"//";

	int q = 97;
	int m_start = 4;
	int m_end = 10;

	bool new_way=true;

	int max_num_instances = 100;//1000;
	double test_ratio = 0;//.2;

	int num_params = this->p*2;
	/*f(use_previous_angles){
		num_params = true;
		loge("use_previous_angles = true;");
	}*/

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
			GeneratorParam param(q, 3, m, true, 97, 100); //q, n, m, shuffle, seed, cutoff; n is ignored
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
			instance.m = m;
			if(nbQubits_acc < 0)
				nbQubits_acc = instance.h.nbQubits;
			else if(nbQubits_acc != instance.h.nbQubits){
				//Need to destroy qureg, because next time experiments with different number of qubits will be run
				qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
				qaoaOptions->accelerator->finalize();
			}

			if(use_database_to_load_dataset){

				FastVQA::Accelerator::DiagonalOpDuplicate diagonalOpDuplicate;
				std::vector<long double> real;
				bool found = this->database->getDataset(instance.h.nbQubits, counter, this->mapOptions->num_qbits_per_x, &real);

				if(found){
					diagonalOpDuplicate.numQubits = instance.h.nbQubits;
					diagonalOpDuplicate.real = real;
					qaoaOptions->accelerator->initialize(&instance.h, true, &diagonalOpDuplicate);
				}else{
					diagonalOpDuplicate.numQubits = -1; //flag that not found
					qaoaOptions->accelerator->initialize(&instance.h, true, &diagonalOpDuplicate);

					this->database->insertDataset(m, counter, mapOptions->num_qbits_per_x, diagonalOpDuplicate.real);
				}

				instance.diagOpDuplicate = diagonalOpDuplicate;

			}else{
				qaoaOptions->accelerator->initialize(&instance.h, false, nullptr);
			}

			//qaoaOptions->accelerator->initialize(&instance.h);

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

	for(int indexx = 1; indexx < 2; ++indexx){ //<2

		if(indexx == 0){
					meta_data = "fixedCMQAOA";
					output<<"CM\n"<<space<<"{";
				}
				else{
					meta_data = "fixedQAOA";
					output<<space<<"{";
				}


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

				std::string optimized_by;
				//std::pair<double, double> final_ab;
				FastVQA::OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {
					iteration_i++;
					bar.tick();
					std::vector<double> angles(x);
					//std::cerr<<this->_cost_fn(&this->train_set, &angles[0]).mean<<std::endl;






					//return strategy_inv_diff(train_dataset, angles, meta_data);

					//return strategy_random_inv_diff(train_dataset, angles, meta_data);
					if(indexx == 0) //CM-QAOA
						return strategy_alpha_c(train_dataset, angles, meta_data, &optimized_by);//strategy_random_alpha_c(train_dataset, angles, meta_data, &optimized_by);						//return strategy_inv_diff(train_dataset, angles, meta_data);
					else if(indexx == 1) //QAOA
						return strategy_alpha_c(train_dataset, angles, meta_data, &optimized_by);						//strategy_inv_diff(train_dataset, angles, meta_data);
					else{
						throw_runtime_error("Not implemented conditional case");
						return 0.;
					}



					//double cost = this->_cost_fn(train_dataset[train_dataset.size()-1], &angles[0], meta_data);
					//std::cerr<<cost<<std::endl;
					//return -cost;
\

					//std::cerr<<nom<<" "<<den<<" "<<alpha<<std::endl;



					//return alpha;//-sum_yi/alpha_calc_dataset_size;//alpha;//-this->_cost_fn(&this->train_set, &angles[0]).mean;
				}, num_params);

				std::vector<double> initial_params;

				FastVQA::OptResult best_result;

				for(int i_rand_angles = 0; i_rand_angles < 100; i_rand_angles++){

					initial_params.clear();

					logi("Random angles "+std::to_string(i_rand_angles)+"/100");

					if(append_previous_angles && AngleResultsExperiment::optAngles[p-1].initialized == true){

						loge("Append previous angles");

						if(indexx == 0){
							for(int i = 0; i < num_params/2-1; ++i){
								double param1 = AngleResultsExperiment::optAngles[p-1].cm_angles[2*i];
								double param2 = AngleResultsExperiment::optAngles[p-1].cm_angles[2*i+1];
								std::cerr<<param1<<" "<<param2<<std::endl;
								initial_params.push_back(param1);
								initial_params.push_back(param2);
							}
							initial_params.push_back(0);
							initial_params.push_back(0);
							std::cerr<<0<<" "<<0<<std::endl;

						}else if(indexx == 1){
							for(int i = 0; i < num_params/2-1; ++i){
								double param1 = AngleResultsExperiment::optAngles[p-1].qaoa_angles[2*i];
								double param2 = AngleResultsExperiment::optAngles[p-1].qaoa_angles[2*i+1];
								std::cerr<<param1<<" "<<param2<<std::endl;
								initial_params.push_back(param1);
								initial_params.push_back(param2);
							}
							initial_params.push_back(0);
							initial_params.push_back(0);
							std::cerr<<0<<" "<<0<<std::endl;
						}else{
							throw_runtime_error("Unimplemented condition case");
						}


					}else{

						loge("append_previous_angles = false");

						std::mt19937 gen(i_rand_angles/*this->seed*/); //rd() instead of 0 - seed
						std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);
						for(int i = 0; i < num_params/2; ++i){
							double param1 = /*pi/4.;*/dis(gen);
							double param2 = /*pi/8.;*/dis(gen);
							std::cerr<<param1<<" "<<param2<<std::endl;
							initial_params.push_back(param1);
							initial_params.push_back(param2);
						}
					}

					std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
					std::vector<double> upperBounds(initial_params.size(), 3.141592654);

					bar.set_progress(0);
					logd("QAOA starting optimization", this->loglevel);
					FastVQA::OptResult result = this->qaoaOptions->optimizer->optimize(f, initial_params, this->qaoaOptions->ftol, this->qaoaOptions->max_iters, lowerBounds, upperBounds);
					logd("QAOA finishing optimization", this->loglevel);


					std::cerr<<"min: "<<result.first.first<<std::endl;
					if(i_rand_angles == 0){
						best_result = result;
						std::cerr<<"new min: "<<result.first.first<<std::endl;
					}else if(result.first.first < best_result.first.first){
						best_result = result;
						std::cerr<<"new min: "<<result.first.first<<std::endl;
					}

				}FastVQA::OptResult result = best_result;


				//std::cerr<<".   2^"<<result.first.first<<"n+"<<result.second<<std::endl;
				std::cerr<<"cost_f min: "<< result.first.first <<"\n";
				std::cerr<<"num_iters: "<<iteration_i<<std::endl;
				/*if(meta_data == "fixedCMQAOA"){
					//std::cerr<<"const std::vector<double> angles_cmqaoa{";

				}
				else if(meta_data == "fixedQAOA"){
					//std::cerr<<"const std::vector<double> angles_optqaoa{";

				}
				else
					throw_runtime_error("Not implemented");*/
				bool start=true;
				int jj = 0;
				for(auto &a: result.first.second){
					if(!start){
						output<<", ";
						if(((jj++)-1) % 3 == 0)
							output<<"\n"<<space;
					}
					output<<std::setprecision (15)<<a;
					start=false;
				}output<<"},\n";
				//std::cerr<<"\n"<<nlopt_res_to_str(result.second)<<std::endl;

				//EVALUATE TRAIN DATASET
				std::vector<double> final_angles(result.first.second);

				double nom=0, den=0;
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

					//for simple alpha
					nom += log2(overlap) * dim[0].m;
					den += dim[0].m * dim[0].m;

				}
				double a=0,b=0;

				a=(sum_yi*sum_xi2-sum_xi*sum_xi_yi)/(train_dataset.size()*sum_xi2-sum_xi*sum_xi);
				b=(train_dataset.size()*sum_xi_yi-sum_xi*sum_yi)/(train_dataset.size()*sum_xi2-sum_xi*sum_xi);
				//std::cerr<<pow(2.71828, a)<<"2^n*"<<b<<std::endl;
				std::cerr<<"2^"<<a<<"+n*"<<b<<std::endl;
				output<<space<<a<<",\n"<<space<<b<<",\n"<<space<<nom/den<<",\n"<<space;
				if(meta_data == "fixedCMQAOA"){
					output<<"\"CM: optimized by "<<optimized_by<<", "<<nlopt_res_to_str(result.second)<<", num_iters: "<<iteration_i<<"\",\n"<<space<<"//QAOA\n";
				}else if(meta_data == "fixedQAOA"){
					output<<"\"QAOA: optimized by "<<optimized_by<<", "<<nlopt_res_to_str(result.second)<<", num_iters: "<<iteration_i<<"\"\n		)";
				}

				double alpha = -b;

	}

	std::cout<<output.str()<<std::endl;



}
double AlphaMinimizationExperiment::_cost_fn(std::vector<AlphaMinimizationExperimentInstance> dataset, const double *angles, std::string meta_data, int probability100){
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


			if(probability100 < 100){
				int p_num = rand() % 100 + 1;  //Generate random number 1 to 100
				if (p_num > probability100 && (i < dataset.size()-2 || gs_overlaps.size() > 0))
					continue;
			}


			if(instance.h.nbQubits != this->qaoaOptions->accelerator->getNumQubitsInQureg()){
				this->qaoaOptions->accelerator->options.createQuregAtEachInilization = true;
				this->qaoaOptions->accelerator->finalize();
			}

			if(instance.diagOpDuplicate.numQubits != -1)
				this->qaoaOptions->diagOpDuplicatePtr = &instance.diagOpDuplicate;

			FastVQA::ExperimentBuffer buffer;
			buffer.storeQuregPtr = true;

			double ground_state_overlap = 0;
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
			this->qaoaOptions->diagOpDuplicatePtr = nullptr;

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
