#include "averaged_svp.h"
#include <fstream>

//#include <iostream>

int myPow (int x, int p) {
  int i = 1;
  for (int j = 1; j <= p; j++)  i *= x;
  return i;
}

void saveEigenspaceToFile(std::string filename, FastVQA::RefEnergies eigenspace){
	std::ofstream f (filename);
	if (f.is_open()){
		for(auto &p: eigenspace){
			qreal energy = std::get<0>(p);
			long long int index = std::get<1>(p);
			f << index << " " << energy << "\n";
		}
	}
	else loge("Unable to create file " + filename);
}

Eigen::Vector<int, Eigen::Dynamic> convertFromBaseTo(int length, unsigned long long int convertFrom, int convertTo) {

	Eigen::Vector<int, Eigen::Dynamic> answer(length);

	int i = 0;
    while (i < length){
    	int digit = convertFrom % convertTo;
        answer(length-(++i))=digit;
        convertFrom /= convertTo;
    }

    return answer;
}

Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> randomVectors(int size, int lb, int ub, int cutoff=-1, int seed=0){

	int delta = ub-lb+1;
	int A_size = myPow(delta, size);

	if(cutoff < 0 || A_size <= cutoff){

		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> A(A_size, size);

		//int** A = (int**) malloc(A_size * sizeof(int*));
		int** A_temp = (int**) malloc(A_size * sizeof(int*));
		for(int i = 0; i < A_size; ++i){
			//A[i] = (int*) malloc(size * sizeof(int));
			A_temp[i] = (int*) malloc(size * sizeof(int));
		}

		for(int i = lb; i <= ub; ++i){
			A(i-lb,0)=i;
		}

		for(int i = 0; i < size-1; ++i){
			for(int j = 0; j < myPow(delta,(i+1)); ++j){
				for(int z = 0; z < delta; ++z){
					for(int x = 0; x < (i+1); ++x)
						A_temp[j*delta+z][x]=A(j,x);
					A_temp[j*delta+z][i+1]=z+lb;
				}
			}
			for(int j = 0; j < myPow(delta,(i+2)); ++j){
				for(int z = 0; z < (i+2); ++z){
					A(j,z)=A_temp[j][z];
				}
			}
		}

		for(int i = 0; i < A_size; ++i){
			free(A_temp[i]);
		}free(A_temp);
		return A;
	}else{ //return at most cutoff number of instances

		if(lb != 0)
			throw_runtime_error("LB=0 unimplemented");

		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> A(cutoff, size);

		int lower = 0;
		int upper = A_size-1;

		// Create and seed the random number generator
		auto gen = std::mt19937(seed);
		auto dist = std::uniform_int_distribution<>(lower, upper);

		bool colision;
		for(int i = 0; i < cutoff; ++i){
			colision = false;
			Eigen::Vector<int, Eigen::Dynamic> rand_vect = convertFromBaseTo(size, dist(gen), delta);
			for(int j = 0; j < i; ++j){
				if(A.row(j).isApprox(rand_vect.transpose())){
					i--;colision=true;break;
				}
			}
			if(!colision)
				A.row(i) = rand_vect;
			if(A(i,0)>100){
				std::cerr<<rand_vect<<std::endl;throw;}

		}
		return A;
	}
}

InstanceGenerator generateQaryUniformFPLLLWay = [](GeneratorParam param){
	if(param.__diagonal)
			throw_runtime_error("Invalid Generator param instance");
		std::vector<HamiltonianWrapper> res;

	int n = param.n;
	int m = param.m;

	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> A=randomVectors((m-n)*n, 0, param.q);

	int counter = 0;
	for(auto row: A.rowwise()){
		Hamiltonian GT(m,m);
		int i, j;
		  /*if (c != r || k > r)
		  {
			FPLLL_ABORT("gen_qary called on an ill-formed matrix");
			return;
		  }*/

		  for (i = 0; i < m - n; i++)
			for (j = 0; j < m - n; j++)
			  GT(m-i-1,m-j-1) = 0;

		  for (i = 0; i < m - n; i++)
			  GT(m-i-1,m-i-1) = 1;

		  for (i = 0; i < m - n; i++)
			for (j = m - n; j < m; j++){
				//G(i,j).randm(q);
				//std::cerr<<i*n+j-(m-n)<<" ";
				GT(m-i-1,m-j-1)=row(i*n+j-(m-n));
			}//std::cerr<<"\n";

		  for (i = m - n; i < m; i++)
			for (j = 0; j < m; j++)
				GT(m-i-1,m-j-1) = 0;

		  for (i = m - n; i < m; i++)
			  GT(m-i-1,m-i-1) = param.q;
		//ZZ_mat<mpz_t> m;
		//m.resize(d, d);
		 std::cerr<<"GT:"<<GT<<"\n";
		 std::cerr<<"G:"<<GT.transpose()<<"\n";
		 auto H=GT*GT.transpose();
		 std::cerr<<"H:"<<H<<"\n\n";
		 HamiltonianWrapper HW = HamiltonianWrapper(GT*GT.transpose(),std::to_string(param.q)+"-ary_"+std::to_string(n)+"x"+std::to_string(m)+"_"+std::to_string(counter));
		  HW.K=GT; //DELETE
		 res.push_back(HW);
		 counter++;
	}
	return res;

};


InstanceGenerator generateQaryUniform = [](GeneratorParam param){
	if(param.__diagonal)
		throw_runtime_error("Invalid Generator param instance");
	std::vector<HamiltonianWrapper> res;

	int n = param.n;
	int m = param.m;

	auto qid = param.q*param.q*Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n);
	auto id = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Identity(m-n,m-n);

	Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> A;

	if(param.cutoff < 0){

		A = randomVectors(n*(m-n), 0, param.q-1, param.cutoff, param.seed);

		for(auto row: A.rowwise()){
			std::cerr<<row<<std::endl;
		}throw;

		if(param.shuffle){
			std::mt19937 rand(param.seed);
			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permX(A.rows());
			permX.setIdentity();
			std::shuffle(permX.indices().data(), permX.indices().data()+permX.indices().size(), rand);
			A = permX * A;   //shuffle row wise
		}

	}else{
		A = randomVectors(n*(m-n), 0, param.q-1, param.cutoff, param.seed);
	}

	int i = 0;
	for(auto row: A.rowwise()){
		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> K = row.reshaped(n,(m-n));
		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> KT=K.transpose();
		Hamiltonian G(m,m);
		G.block(0,0,n,n)=qid;
		G.block(0,n,n,m-n)=K;
		G.block(n,0,m-n,n)=KT;
		G.block(n,n,m-n,m-n)=id+KT*K;
		HamiltonianWrapper HW = HamiltonianWrapper(G,std::to_string(param.q)+"-ary_"+std::to_string(n)+"x"+std::to_string(m)+"_"+std::to_string(i));
		HW.K=K;
		res.push_back(HW);
		i++;
	}

	return res;
};

DiagonalInstanceGenerator generateDiagonalUniform = [](GeneratorParam param){

	if(!param.__diagonal)
		throw_runtime_error("Invalid Generator param instance");

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

/*extern */DiagonalInstanceGenerator generateDiagonalExtensive= [](GeneratorParam param){

	if(!param.__diagonal)
		throw_runtime_error("Invalid Generator param instance");

    std::vector<DiagonalHamiltonian> res;

    //Calculate A
    	std::cerr<<"REWRITE THIS AS RANDOMVECTORS FUNC\n";
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
