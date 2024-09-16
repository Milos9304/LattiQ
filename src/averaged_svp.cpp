#include "averaged_svp.h"
#include "lattice/lattice.h"
#include "experiment_runner.h"
#include "test_runner.h"
#include "io/sql_io.h"

#include "popl.hpp"

using namespace popl;

int main(int ac, char** av){

	int seed = 1997;

	OptionParser op("Allowed options");
	auto help_option     		= op.add<Switch>("h", "help", "produce help message");
	auto n_opt		     		= op.add<Value<int>>("n", "", "", 1);
	auto m_opt		    		= op.add<Value<int>>("m", "", "", 3);
	auto log_level     		  	= op.add<Value<int>>("", "loglevel", "0 - debug, 1 - info, 2 - warning, 3 - error", 1);
	auto print_hml      		= op.add<Switch>("", "print_hml", "print calculated hamiltonian expression");
	auto niters        	 		= op.add<Value<int>>("i", "iters", "max num of iterations", 1000);
	auto qubits_per_x    		= op.add<Value<int>>("q", "", "qubits per x for uniform assignment (-1 means disabled)", -1);
	auto absolute_bound	 		= op.add<Value<int>>("e", "", "exponent bound for each coefficient, i.e. |x_i|<=2^e (-1 means disabled)", -1);
	auto qaoadepth      		= op.add<Value<int>>("p", "depth", "qaoa depth", 1);
	auto penalty         		= op.add<Value<int>>("l", "penalty", "penalty", 100);
	auto save_eigenspace  		= op.add<Switch>("", "espace", "save eigenspace to file");
	auto param_experiment 		= op.add<Value<int>>("", "paramexp", "experiment with n different random initial parameters (0 for disabled)", 0);
	auto angle_results	  		= op.add<Switch>("", "angleres", "results of constant angles experiment");
	auto angle_results_opt 		= op.add<Switch>("", "angleresopt", "results of opt experiment");
	auto test_variable_subst 	= op.add<Switch>("", "testsubst", "test variable substitution");
	auto performance_calc   	= op.add<Switch>("", "performance", "calculate performance");
	//auto cmqaoa					= op.add<Switch>("", "cm", "run cmqaoa experiment");
	auto alphaminim				= op.add<Switch>("", "alpha", "run alpha minimization experiment");
	auto database_info		   	= op.add<Switch>("", "dinfo", "get database info");
	auto g1					   	= op.add<Switch>("", "g1", "generate graph 1");
	auto g2					   	= op.add<Switch>("", "g2", "generate graph 2");
	auto seed_opt	     		= op.add<Value<int>>("s", "seed", "Seed", 0);
	auto m_start	     		= op.add<Value<int>>("", "mstart", "m_start", 4);
	auto m_end		     		= op.add<Value<int>>("", "mend", "m_end", 20);
	auto aqc_pqc	     		= op.add<Switch>("", "aqcpqc", "aqcpqc");


	op.parse(ac, av);
	if (help_option->is_set()){
		std::cout << op << "\n";
		return 0;
	}

	loge("MaxCoeff normalisation in getHamiltonian!");

	/*Database::DatasetRow row;
	row.type = "qary";
	row.q = 7;
	row.n = 1;
	row.m = 4;
	row.p = 1;
	row.penalty = 2;
	//database.write(&row);
	std::cerr<<database.contains_qary(7,1,4,1,0,true);

	return 0;*/

	int loglevel = log_level->value();

	int n = n_opt->value();
	int m = m_opt->value();

	if( (!performance_calc->is_set()) && ((qubits_per_x->value() == -1 && absolute_bound->value()==-1) || (qubits_per_x->value() != -1 && absolute_bound->value() !=-1)) ){
		throw_runtime_error("Exactly one of 'qubits_per_x' or 'absolute_bound' must be set.");
	}

	FastVQA::AcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";
	acceleratorOptions.log_level = /*performance_calc->is_set() ? 3 :*/ log_level->value();

	FastVQA::NLOptimizer optimizer;

	FastVQA::Accelerator accelerator(acceleratorOptions);

	FastVQA::QAOAOptions qaoaOptions;
	qaoaOptions.log_level = log_level->value();
	qaoaOptions.max_iters = niters->value();
	qaoaOptions.optimizer = &optimizer;
	qaoaOptions.accelerator = &accelerator;
	qaoaOptions.nbSamples_calcVarAssignment=1000;
	qaoaOptions.p = qaoadepth->value();
	qaoaOptions.ftol = 10e-12;
	long long int max_iters = 0;
	//DiagonalHamiltonian h;
	//calculateAverage(n, &h);

	MapOptions mapOptions;
	mapOptions.verbose = print_hml->is_set();
	mapOptions.num_qbits_per_x = qubits_per_x->value();
	mapOptions.absolute_bound = absolute_bound->value();
	//mapOptions.pen_mode = MapOptions::penalty_all;
	mapOptions.bin_map = penalty->value() > 0 ? MapOptions::zeta_omega_exact : MapOptions::naive_overapprox;
	mapOptions.penalty = penalty->value();
	mapOptions.loglevel = log_level->value();

	//GeneratorParam param(n);
	//std::vector<DiagonalHamiltonian> gramiams = generateDiagonalExtensive(param);


	//for(auto &hw: gramian_wrappers){
	//	std::cerr<<hw.hamiltonian<<"\n\n";
	//}
	//return 0;

	ExperimentSetup experimentSetup;
	if(performance_calc->is_set()){

		const std::string database_file = "../experiments/database.db";

		if(database_info->is_set()){
			Database::print_sqlite_info(database_file);
			return 0;
		}

		Database database(database_file, Database::DATABASE_QARY_PERFORMANCE);
		test_execution_time(&qaoaOptions, &database);
		return 0;

	}else if(test_variable_subst->is_set()){
		test_variable_substitution(&mapOptions);
		return 0;
	}else if(alphaminim->is_set()){
		logi("Running alpha minimization experiment", loglevel);

		const std::string database_file = "../experiments/database_eigengen.db";
		if(database_info->is_set()){
			Database::print_sqlite_info(database_file);
			return 0;
		}

		Database database(database_file, Database::DATABASE_EIGENGEN_DATASET);
		AlphaMinimizationExperiment alphaMinimExp(loglevel, &qaoaOptions, &mapOptions, &database, seed_opt->value());
		alphaMinimExp.run(true);

		return 0;
	}else if(angle_results->is_set() || angle_results_opt->is_set()){

		const std::string database_file = "../experiments/database_eigengen.db";
		if(database_info->is_set()){
			Database::print_sqlite_info(database_file);
			return 0;
		}

		//qaoaOptions.p = 6;


		Database database(database_file, Database::DATABASE_EIGENGEN_DATASET);
		AngleResultsExperiment angleResultsExp(loglevel, m_start->value(), m_end->value(), &qaoaOptions, &mapOptions, &database, seed_opt->value(), true);

		if(angle_results_opt->is_set())
			angleResultsExp.run_qaoa_with_optimizer();
		else
			angleResultsExp.run();

		return 0;
	}else if(aqc_pqc->is_set()){

		const std::string database_file = "../experiments/database_eigengen_aqcpqc.db";
		if(database_info->is_set()){
			Database::print_sqlite_info(database_file);
			return 0;
		}

		Database database(database_file, Database::DATABASE_EIGENGEN_AQCPQC_DATASET);
		AqcPqcExperiment aqcPqcExperiment(loglevel, m_start->value(), m_end->value(), &qaoaOptions, &mapOptions, &database, seed_opt->value(), true);
		aqcPqcExperiment.run();

		return 0;
	}else if(g1->is_set()){

		//const std::string database_file = "../experiments/database_g1.db";
		//if(database_info->is_set()){
		//	Database::print_sqlite_info(database_file);
		//	return 0;
		//}

		//Database database(database_file, Database::DATABASE_ANGLERES);
		G1 g1(loglevel, &qaoaOptions, &mapOptions/*, &database*/);
		g1.run();

		return 0;
	}else if(g2->is_set()){

		G2 g2(loglevel, &qaoaOptions, &mapOptions/*, &database*/);
		g2.run();

		return 0;
	}/*else if(cmqaoa->is_set()){
		const std::string database_file = "../experiments/database_cmqaoa.db";
		if(database_info->is_set()){
			Database::print_sqlite_info(database_file);
			return 0;
		}
		Database database(database_file, Database::DATABASE_CM_QAOA);
		CmQaoaExperiment cmQaoaExperiment(&qaoaOptions, &mapOptions, &database, loglevel);
		cmQaoaExperiment.run();
		return 0;

	}*/else if(param_experiment->value() > 0){
		logi("Running manyParams experiment", loglevel);
		experimentSetup.experiment_type = "manyParams";
		experimentSetup.num_rand_params=param_experiment->value();
	}else{
		experimentSetup.experiment_type = "singleInstance";
	}

	GeneratorParam param(7,n,m);
	std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);//generateQaryUniformFPLLLWay(param); //generateQaryUniform(param);

	std::vector<double> hit_rates;

	int counter=0;
	for(auto w : gramian_wrappers){

		if(log_level->value() < 2){
			std::cerr<<100*float(counter++)/gramian_wrappers.size()<<"\%\n";
		}

		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> G = w.hamiltonian;
		std::string name = w.name;

		Lattice l(G, name);

		//std::cerr<<"G"<<G<<"\n";
		FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);
		int nbQubits=h.nbQubits;
		//std::cerr<<"qs: "<<nbQubits<<"\n";

		//l.outputGramianToFile(name);
		//h.to_ising_file(name);

		accelerator.initialize(&h);
		FastVQA::RefEnergies solutions = accelerator.getSolutions();

		//std::cerr<<"q="<<w.K<<"\n\n";
		//std::cerr<<"G="<<G<<"\n";
		qreal energy = solutions[0].value;
		for(const auto &sol: solutions){
			if(sol.value < energy)
				energy = sol.value;
		}
		for(auto &sol: solutions){

			//energy index
			if(energy != sol.value)
				logw("An outcome with different energy marked as a solution! TODO: This was expected. Need to change code for RefEnergies.isConsideredSolution");
			long long int index = sol.index;
			//std::cerr<<index<<"    "<<energy<<"\n";
			VectorInt solVectFromAcc = l.quboToXvector(index, nbQubits);
			Eigen::Vector<int, Eigen::Dynamic> solVect(solVectFromAcc.size());
			for(int i = 0; i < solVectFromAcc.size(); ++i){
				solVect[i] = solVectFromAcc[i].get_si();
				//std::cerr<<solVectFromAcc[i]<<" ";
			}

			Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> r = solVect.transpose() * (G * solVect);
			assert(r.rows() == 1 && r.cols() == 1);
			if(r(0,0) != energy){
				std::stringstream ss;
				ss << "Problem with energy evaluation! x="<<solVect.transpose()<<"\nG=\n"<<G<<"\nx^TGx="<<r(0,0)<<" while the expected minimum is "<<energy;
				//throw_runtime_error(ss.str());
				std::cerr<<ss.str()<<"\n";
			}else{
				//std::cerr<<"E: "<<energy<<" SOLUTION FOUND: "<<solVect.transpose()<<"\n";
			}
			//auto p = solVectT*G;
			//std::cerr<<"\n"<<r;//*solVect<<"\n";
			//std::cerr<<"\ni: "<<index<<"\n";
			//std::cerr<<"e: "<<energy<<"\n";

			//break;
		}

		experimentSetup.qaoaOptions = &qaoaOptions;
		experimentSetup.hamiltonian = &h;
		experimentSetup.minimum_energy = energy;
		experiment_runner(&experimentSetup, name, loglevel);

		/*FastVQA::ExperimentBuffer buffer;
		qaoa_instance.run_qaoa(&buffer, &h, &qaoaOptions);
		if(buffer.opt_val != energy){
			logw("Solution returned has incorrect energy");
		}
		for(auto &param: buffer.initial_params){
			std::cerr<<std::get<0>(param)<<" "<<std::get<1>(param)<<"\n";
		}

		hit_rates.push_back(buffer.getTotalHitRate());*/

		if(save_eigenspace->is_set()){
			saveEigenspaceToFile(name+"_espace", accelerator.getEigenspace());
		}
		break;
	}

	double sum=0;
	for(auto &hr : hit_rates){
		sum+=hr;
	}
	std::cerr<<sum/hit_rates.size()<<"\n";

	return 0;
}
