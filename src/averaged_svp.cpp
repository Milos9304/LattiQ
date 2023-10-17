#include "averaged_svp.h"

#include "FastVQA/fastVQA.h"
#include "popl.hpp"

using namespace popl;

int n = 3; //dim of basis

int main(int ac, char** av){

	OptionParser op("Allowed options");
	auto help_option     = op.add<Switch>("h", "help", "produce help message");

	op.parse(ac, av);
	if (help_option->is_set()){
		std::cout << op << "\n";
		return 0;
	}

	DiagonalHamiltonian h;
	calculateAverage(n, &h);

	std::cerr<<h;

	return 0;
}
