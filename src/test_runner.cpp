#include "FastVQA/fastVQA.h"

#include "lattice/lattice.h"
#include "test_runner.h"
#include "averaged_svp.h"

void test_variable_substitution(MapOptions* mapOptions){

	Lattice l("variable_test", mapOptions);

}

void test_execution_time(FastVQA::QAOAOptions* qaoaOptions){

	int penalty = 0;

	MapOptions mapOptions;
	mapOptions.verbose = false;//true;
	//mapOptions.num_qbits_per_x = qubits_per_x->value();
	//mapOptions.absolute_bound = absolute_bound->value();
	mapOptions.pen_mode = MapOptions::penalty_all;
	mapOptions.bin_map = penalty > 0 ? MapOptions::zeta_omega_exact : MapOptions::naive_overapprox;
	mapOptions.penalty = penalty;

	qaoaOptions->log_level = 3;
	qaoaOptions->max_iters = 10000;

	//int m_max = 10;
	int num_samples = 2;
	bool shuffle = false;

		GeneratorParam param(7,1,2, shuffle, 0, num_samples);
		std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

		for(int qs=1; qs<20; ++qs){

			mapOptions.num_qbits_per_x = qs;

			for(auto w : gramian_wrappers){

				Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> G = w.hamiltonian;
				std::string name = w.name;

				Lattice l(G, name);

				//std::cerr<<"G"<<G<<"\n";
				FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);
				int nbQubits=h.nbQubits;

				logi("Running experiment with " + std::to_string(nbQubits) + " qubits");

				qaoaOptions->accelerator->initialize(&h);
				FastVQA::Qaoa qaoa_instance;
				FastVQA::ExperimentBuffer buffer;
				buffer.storeQuregPtr = false;
				qaoa_instance.run_qaoa(&buffer, &h, qaoaOptions);

				std::cerr<<buffer.num_iters<<"\n";
				std::cerr<<buffer.opt_message<<"\n";

			}


	}

}
