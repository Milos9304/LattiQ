#include "averaged_svp.h"
//#include <iostream>

InstanceGenerator generateDiagonalUniform = [](GeneratorParam param){

	throw_runtime_error("Not implemented");

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

int myPow (int x, int p) {
  int i = 1;
  for (int j = 1; j <= p; j++)  i *= x;
  return i;
}
extern InstanceGenerator generateDiagonalExtensive= [](GeneratorParam param){
    std::vector<DiagonalHamiltonian> res;

    //Calculate A
		int delta = param.lambda_ub-param.lambda3_lb+1;
		int A_size = myPow(delta,(param.n - 2));
		int** A = (int**) malloc(A_size * sizeof(int*));
		int** A_temp = (int**) malloc(A_size * sizeof(int*));
		for(int i = 0; i < A_size; ++i){
			A[i] = (int*) malloc((param.n-2) * sizeof(int));
			A_temp[i] = (int*) malloc((param.n-2) * sizeof(int));
		}
		for(int i = param.lambda3_lb; i <= param.lambda_ub; ++i){
			A[i-param.lambda3_lb][0]=i;
		}
		for(int i = 0; i < param.n-3; ++i){
			for(int j = 0; j < myPow(delta,(i+1)); ++j){
				for(int z = 0; z < delta; ++z){
					for(int x = 0; x < (i+1); ++x)
						A_temp[j*delta+z][x]=A[j][x];
					A_temp[j*delta+z][i+1]=z+param.lambda3_lb;
				}
			}
			for(int j = 0; j < myPow(delta,(i+2)); ++j){
				for(int z = 0; z < (i+2); ++z){
					A[j][z]=A_temp[j][z];
				}
			}
		}

	    for(int i = 0; i < A_size; ++i){
	    	free(A_temp[i]);
	    }free(A_temp);
	//A calculated

	for(int i = 0; i < A_size; ++i){
		DiagonalHamiltonian G(param.n);
		G(0)=param.lambda1;
		G(1)=param.lambda2;
		for(int j = 0; j < param.n-2; ++j){
			//std::cerr<<A[i][j-2]<<" ";
			G(j+2)=A[i][j];
		}
		res.push_back(G);
		//std::cerr<<"\n";
	}

    for(int i = 0; i < A_size; ++i)
    	free(A[i]);
    free(A);

    return res;
};


void calculateAverage(int n, DiagonalHamiltonian* h){

	//GeneratorParam param(n);
	//std::vector<DiagonalHamiltonian> res = generateDiagonalExtensive(param);
	//for(auto r : res)
	//	std::cerr << r.transpose() <<"\n";

	//*h = res[0];

}
