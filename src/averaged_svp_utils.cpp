#include "averaged_svp.h"

InstanceGenerator generateDiagonalUniform = [](GeneratorParam param){
    std::vector<DiagonalHamiltonian> res;

    for(int i = 0; i < param.num_instances; ++i){
    	DiagonalHamiltonian G;
    	G.setZero(param.n);//resize(param.n);
        G(0) = param.lambda1;
        G(1) = param.lambda2;

    	res.push_back(G);
    }
    return res;
};

void calculateAverage(int n, DiagonalHamiltonian* h){

	GeneratorParam param(n);
	std::vector<DiagonalHamiltonian> res = generateDiagonalUniform(param);
	*h = res[0];

}
