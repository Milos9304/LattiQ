#include "popl.hpp"
#include "FastVQA/fastVQA.h"
#include "lattiq.h"

#include <mpi.h>

using namespace popl;

int main(int ac, char** av){

	int seed = 1997;

	int rank, numRanks;

	MPI_Init(&ac,&av);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0){

		OptionParser op("Allowed options");
		auto help_option     = op.add<Switch>("h", "help", "produce help message");
		auto qaoa 		     = op.add<Switch>("", "qaoa", "run qaoa algorithm");
		auto vqe 		     = op.add<Switch>("", "vqe", "run vqe algorithm");
		auto ansatz_name     = op.add<Value<std::string>>("a", "ansatz","Ry_CNOT_all2all_Ry/qaoa/EfficientSU2", "Ry_CNOT_all2all_Ry");
		auto seed_option 	 = op.add<Value<int>>("", "seed", "seed for the experiments", seed);
		//auto enumeration     = op.add<Switch>("", "enum", "enumerate all qubo configurations");
		auto config 	     = op.add<Value<std::string>>("", "config", "config file location", "");
		auto lattice_file    = op.add<Value<std::string>>("", "lattice", "lattice file location", "");
		auto niters          = op.add<Value<int>>("i", "iters", "max num of iterations", 1000);
		auto nbSamples 		 = op.add<Value<int>>("n", "nbSamples", "number of samples in var assigmnent", 1024);
		auto save_hml        = op.add<Value<std::string>>("", "savehml", "save hamiltonian to file", "");
		auto load_hml        = op.add<Value<std::string>>("", "loadhml", "save hamiltonian to file", "");
		auto debug           = op.add<Switch>("d", "debug", "print debug messages");
		auto qubits_per_x    = op.add<Value<int>>("q", "", "qubits per x", 1);
		auto overlap_trick   = op.add<Switch>("o", "", "perform overlap trick");
		auto overlap_penalty = op.add<Value<int>>("p", "", "overlap penalty", 0);
		auto lll_preprocess  = op.add<Switch>("", "lll", "perform LLL preprocessing on the lattice");
		auto second_eigval   = op.add<Switch>("", "second", "pick second lowest energy");

		auto initial_alpha   = op.add<Value<double>>("x", "alpha", "initial alpha value", 1);
		auto linear_alpha    = op.add<Switch>("l", "linear", "linear alpha");
		auto final_alpha     = op.add<Value<double>>("f", "final_alpha", "final alpha value", 0.5);
		auto max_alpha_iters = op.add<Value<int>>("m", "max_alpha_iters", "max alpha iters", 1000);

		auto paper_exp		 = op.add<Switch>("e", "paperexp", "perform experiment as in the paper");
		auto rank_reduce 	 = op.add<Value<int>>("r", "", "rank truncation for paperexp", 0);
		auto circ_dir_prefix = op.add<Value<std::string>>("c", "circ-dir-prefix", "", "../experiment_files");

		auto save_ansatz	 = op.add<Switch>("", "saveAnsatz", "save ansatz files");
		auto load_ansatz	 = op.add<Switch>("", "loadAnsatz", "load ansatz files");

		auto save_interm  = op.add<Value<std::string>>("", "si", "save intermediate results (for specific experiments only)", "");
		auto load_interm  = op.add<Value<std::string>>("", "li", "load intermediate results (for specific experiments only)", "");

		auto statsfile_prefix = op.add<Value<std::string>>("", "prefix", "statsfile prefix", "rank");

		op.parse(ac, av);

		seed = seed_option->value();

	}

	return 0;
}
