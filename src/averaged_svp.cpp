#include "averaged_svp.h"
#include "lattice/lattice.h"

#include "FastVQA/fastVQA.h"
#include "popl.hpp"

using namespace popl;

int n = 4; //dim of basis

int main(int ac, char** av){

	int seed = 1997;
	int loglevel = 0;

	OptionParser op("Allowed options");
	auto help_option     = op.add<Switch>("h", "help", "produce help message");
	auto log_level       = op.add<Value<int>>("", "loglevel", "0 - debug, 1 - info, 2 - warning, 3 - error", 1);
	auto print_hml       = op.add<Switch>("", "print_hml", "print calculated hamiltonian expression");
	auto niters          = op.add<Value<int>>("i", "iters", "max num of iterations", 1000);
	auto qubits_per_x    = op.add<Value<int>>("q", "", "qubits per x for uniform assignment", 1);

	op.parse(ac, av);
	if (help_option->is_set()){
		std::cout << op << "\n";
		return 0;
	}

	FastVQA::AcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";

	FastVQA::NLOptimizer optimizer;

	FastVQA::Accelerator accelerator(acceleratorOptions);

	FastVQA::QAOAOptions qaoaOptions;
	qaoaOptions.log_level = log_level->value();
	qaoaOptions.max_iters = niters->value();
	qaoaOptions.optimizer = &optimizer;
	qaoaOptions.accelerator = &accelerator;
	//DiagonalHamiltonian h;
	//calculateAverage(n, &h);

	MapOptions mapOptions;
	mapOptions.verbose = print_hml->is_set();
	mapOptions.num_qbits_per_x = qubits_per_x->value();
	mapOptions.penalty = 100;

	GeneratorParam param(n);
	std::vector<DiagonalHamiltonian> gramiams = generateDiagonalExtensive(param);

	for(auto G : gramiams){

		FastVQA::Qaoa qaoa_instance;
		Lattice l(G, "name", true);
		FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);

	}


	return 0;
}
