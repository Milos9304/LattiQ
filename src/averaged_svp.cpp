#include "averaged_svp.h"
#include "lattice/lattice.h"

#include "popl.hpp"

using namespace popl;

int main(int ac, char** av){

	int seed = 1997;
	int loglevel = 0;

	OptionParser op("Allowed options");
	auto help_option     = op.add<Switch>("h", "help", "produce help message");
	auto n_opt		     = op.add<Value<int>>("n", "", "", 3);
	auto m_opt		     = op.add<Value<int>>("m", "", "", 3);
	auto log_level       = op.add<Value<int>>("", "loglevel", "0 - debug, 1 - info, 2 - warning, 3 - error", 1);
	auto print_hml       = op.add<Switch>("", "print_hml", "print calculated hamiltonian expression");
	auto niters          = op.add<Value<int>>("i", "iters", "max num of iterations", 1000);
	auto qubits_per_x    = op.add<Value<int>>("q", "", "qubits per x for uniform assignment (-1 means disabled)", -1);
	auto absolute_bound  = op.add<Value<int>>("e", "", "exponent bound for each coefficient, i.e. |x_i|<=2^e (-1 means disabled)", -1);
	auto penalty         = op.add<Value<int>>("p", "penalty", "penalty", 100);
	auto save_eigenspace = op.add<Switch>("", "espace", "save eigenspace to file");

	op.parse(ac, av);
	if (help_option->is_set()){
		std::cout << op << "\n";
		return 0;
	}

	int n = n_opt->value();
	int m = m_opt->value();

	if((qubits_per_x->value() == -1 && absolute_bound->value()==-1) || (qubits_per_x->value() != -1 && absolute_bound->value() !=-1)){
		throw_runtime_error("Exactly one of 'qubits_per_x' or 'absolute_bound' must be set.");
	}

	FastVQA::AcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";
	acceleratorOptions.log_level = log_level->value();

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
	mapOptions.absolute_bound = absolute_bound->value();
	mapOptions.pen_mode = MapOptions::penalty_all;
	mapOptions.bin_map = MapOptions::zeta_omega_exact;
	mapOptions.penalty = penalty->value();

	//GeneratorParam param(n);
	//std::vector<DiagonalHamiltonian> gramiams = generateDiagonalExtensive(param);
	n=1;m=3;
	GeneratorParam param(n,m);
	std::vector<Hamiltonian> gramiams = generateQaryUniform(param);

	for(auto G : gramiams){

		FastVQA::Qaoa qaoa_instance;
		Lattice l(G, "name");
		std::cerr<<"G"<<G<<"\n";
		FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);
		std::cerr<<"NEEEW\n";
		accelerator.initialize(&h);

		if(save_eigenspace->is_set()){
			saveEigenspaceToFile("name", accelerator.getEigenspace());

		}

		break;

		/*FastVQA::Qaoa qaoa_instance;
		Lattice l(G, "name");
		FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);
		std::cerr<<"NEEEW\n";
		accelerator.initialize(&h);
		break;*/
	}


	return 0;
}
