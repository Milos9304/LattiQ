#include "lattiq.h"

#include "run_paper_experiment.h"

#include "popl.hpp"
#include "FastVQA/fastVQA.h"

using namespace popl;

std::string extract_id(std::string const& s){
    std::string::size_type pos = s.find('_');
    if (pos != std::string::npos)
        return s.substr(0, pos);
    else
    	return s;
}

int main(int ac, char** av){

	int seed = 1997;
	int loglevel = 0;

	fastVQA::AcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";

	OptionParser op("Allowed options");
	auto help_option     = op.add<Switch>("h", "help", "produce help message");
	auto log_level       = op.add<Value<int>>("", "loglevel", "0 - debug, 1 - info, 2 - warning, 3 - error", 0);
	auto logExpecStd	 = op.add<Switch>("e", "", "log interdmediate expectation values to standard output");
	auto qaoa 		     = op.add<Switch>("", "qaoa", "run qaoa algorithm");
	auto vqe 		     = op.add<Switch>("", "vqe", "run vqe algorithm");
	auto ansatz_name     = op.add<Value<std::string>>("a", "ansatz","Ry_CNOT_all2all_Ry/qaoa/EfficientSU2", "Ry_CNOT_all2all_Ry");
	auto seed_option 	 = op.add<Value<int>>("", "seed", "seed for the experiments", seed);
	//auto enumeration     = op.add<Switch>("", "enum", "enumerate all qubo configurations");
	//auto config 	     = op.add<Value<std::string>>("", "config", "config file location", "");
	auto lattice_file    = op.add<Value<std::string>>("", "lattice", "lattice file location", "");
	auto niters          = op.add<Value<int>>("i", "iters", "max num of iterations", 1000);
	auto nbSamples 		 = op.add<Value<int>>("n", "nbSamples", "number of samples in var assigmnent", 1024);
	//auto save_hml        = op.add<Value<std::string>>("", "savehml", "save hamiltonian to file", "");
	//auto load_hml        = op.add<Value<std::string>>("", "loadhml", "save hamiltonian to file", "");
	auto qubits_per_x    = op.add<Value<int>>("q", "", "qubits per x for uniform assignment", 1);
	auto overlap_trick   = op.add<Switch>("o", "", "perform overlap trick");
	auto overlap_penalty = op.add<Value<int>>("p", "", "overlap penalty", 0);
	auto lll_preprocess  = op.add<Switch>("", "lll", "perform LLL preprocessing on the lattice");
	auto second_eigval   = op.add<Switch>("", "second", "pick second lowest energy");

	auto initial_alpha   = op.add<Value<double>>("x", "alpha", "initial alpha value", 1);
	auto linear_alpha    = op.add<Switch>("l", "linear", "linear alpha");
	auto final_alpha     = op.add<Value<double>>("f", "final_alpha", "final alpha value", 0.5);
	auto max_alpha_iters = op.add<Value<int>>("m", "max_alpha_iters", "max alpha iters", 1000);

	auto instance_select = op.add<Value<int>>("s", "iselect", "Select i-th instance in the dataset. If unset (=-1), all instances are being run.", -1);
	auto dataset_name    = op.add<Value<std::string>>("d", "datasetName","dataset name in ./experiments folder", "qary_25_50");
	auto rank_select 	 = op.add<Value<int>>("r", "", "rank truncation for paperexp", 0);
	auto circ_dir_prefix = op.add<Value<std::string>>("c", "circ-dir-prefix", "", "../experiment_files");

	//auto save_ansatz	 = op.add<Switch>("", "saveAnsatz", "save ansatz files");
	//auto load_ansatz	 = op.add<Switch>("", "loadAnsatz", "load ansatz files");

	//auto save_interm  = op.add<Value<std::string>>("", "si", "save intermediate results (for specific experiments only)", "");
	//auto load_interm  = op.add<Value<std::string>>("", "li", "load intermediate results (for specific experiments only)", "");

	auto statsfile_prefix = op.add<Value<std::string>>("", "prefix", "statsfile prefix", "rank");

	op.parse(ac, av);

	if (help_option->is_set()){
		std::cout << op << "\n";
		return 0;
	}

	if(log_level->value() < 0 || log_level->value() > 3){
		throw_runtime_error("Invalid loglevel value.");
	}
	loglevel = log_level->value();

	seed = seed_option->value();
	logi("Using seed " + std::to_string(seed), loglevel);

	int num_lattices;
	std::vector<Lattice*> lattices;

	acceleratorOptions.log_level = loglevel;
	acceleratorOptions.samples_cut_ratio = initial_alpha->value();
	acceleratorOptions.final_alpha = final_alpha->value();
	acceleratorOptions.max_alpha_iters = max_alpha_iters->value();
	if(linear_alpha->is_set()){
		acceleratorOptions.alpha_f = "linear";
		logi("Linearly increasing cvar_alpha from " +
				std::to_string(acceleratorOptions.samples_cut_ratio) + " to " +
				std::to_string(acceleratorOptions.final_alpha) + " in " +
				std::to_string(acceleratorOptions.max_alpha_iters) + " iters", loglevel);
	}else{
		acceleratorOptions.alpha_f = "constant";
		if(acceleratorOptions.samples_cut_ratio == 1)
			logi("Non cvar version", loglevel);
		else
			logi("Constant cvar_alpha= " + std::to_string(acceleratorOptions.samples_cut_ratio), loglevel);
	}

	initialize_paper_experiment(dataset_name->value(), lattices, rank_select->value(), instance_select->value());
	num_lattices = lattices.size();
	logi("Dataset " + dataset_name->value() + " succesfully loaded", loglevel);

	fastVQA::Accelerator accelerator(std::make_shared<fastVQA::AcceleratorOptions>(acceleratorOptions));
	accelerator.env = createQuESTEnv();

	fastVQA::NLOptimizer optimizer;

	fastVQA::QAOAOptions *qaoaOptions;
	fastVQA::VQEOptions *vqeOptions;

	fastVQA::VQAOptions *vqaOptions = qaoa->is_set() ? new fastVQA::QAOAOptions() : (true? new fastVQA::VQEOptions() : new fastVQA::VQAOptions());
	vqaOptions->log_level = log_level->value();
	vqaOptions->max_iters = niters->is_set() ? (niters->value() == 0 ? 5000 : niters->value()): 5000;
	vqaOptions->detailed_log_freq = 0;//50;
	vqaOptions->optimizer = &optimizer;
	vqaOptions->accelerator = &accelerator;
	vqaOptions->logEnergies = true;
	vqaOptions->expectationToStandardOutput = logExpecStd->is_set();
	//vqaOptions->calcVarAssignment = true;
	//vqaOptions->provideHamiltonian = true;
	//vqaOptions->saveIntermediate = save_interm->is_set() ? (save_interm->value() == "" ? false : true) : false;
	//vqaOptions->s_intermediateName = vqaOptions->saveIntermediate ? save_interm->value() : "";
	//vqaOptions->loadIntermediate = load_interm->is_set() ? (load_interm->value() == "" ? false : true) : false;
	//vqaOptions->l_intermediateName = vqaOptions->loadIntermediate ? load_interm->value() : "";
	vqaOptions->overlap_trick = overlap_trick->is_set();
	vqaOptions->nbSamples_calcVarAssignment = nbSamples->value();
	//vqaOptions->save_ansatz = save_ansatz->is_set();
	//vqaOptions->load_ansatz = load_ansatz->is_set();
	vqaOptions->ansatz_name = ansatz_name->value();

	logd("VQAOptions set", loglevel);

	if(qaoa->is_set()){
		qaoaOptions = static_cast<fastVQA::QAOAOptions*>(vqaOptions);
		qaoaOptions->simplifiedSimulation = true;
		qaoaOptions->extendedParametrizedMode = false;//true;
	}else if(vqe->is_set()){
		vqeOptions = static_cast<fastVQA::VQEOptions*>(vqaOptions);
	}

	MapOptions* mapOptions = new MapOptions();
	mapOptions->verbose = log_level->value() < 1;
	mapOptions->num_qbits_per_x = qubits_per_x->value();
	mapOptions->penalty=overlap_penalty->value();
	if(ansatz_name->value() == "qaoa")
		mapOptions->pen_mode = MapOptions::penalty_all;
	else
		mapOptions->pen_mode = MapOptions::no_hml_penalization;

	int inst_counter = 0;
	for(auto &lattice : lattices){

		logi("Running " + lattice->name);

		int new_id = std::stoi(extract_id(lattice->name));

		if(qaoa->is_set()){
			throw_runtime_error("TODO: QAOA not implemented yet");
		}else if(vqe->is_set()){
			fastVQA::Vqe vqe_instance;
			fastVQA::ExperimentBuffer buffer;
		}

		fastVQA::Hamiltonian hamiltonian = lattice->getHamiltonian(mapOptions);
		vqeOptions->zero_reference_states = lattice->getZeroReferenceStates();

		fastVQA::Vqe vqe_instance;
		fastVQA::ExperimentBuffer buffer;

		vqe_instance.run_vqe(&buffer, &hamiltonian, vqeOptions);

		inst_counter++;
	}


	return 0;
}
