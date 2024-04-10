#include "hermite_factor.h"
#include "io/logger.h"

double calculate_hermite_factor(int n, double volume_1n, std::shared_ptr<Qureg> stateVect, FastVQA::RefEnergies* refEnergies){

	std::vector<std::pair<int, double>> pis;
	int temp_len=(*refEnergies)[0].value;
	double temp_prob=0;
	double prevVal = (*refEnergies)[0].value;
	for(long long int k = 0; k < refEnergies->size(); ++k){

		if((*refEnergies)[k].value < prevVal)
			throw_runtime_error("RefEnergies were not sorted by value!");

		if((*refEnergies)[k].value == prevVal){
			//temp_len = (*refEnergies)[k].value;

			long long int j = (*refEnergies)[k].index;
			temp_prob += stateVect->stateVec.real[j]*stateVect->stateVec.real[j]+stateVect->stateVec.imag[j]*stateVect->stateVec.imag[j];
		}
		else if((*refEnergies)[k].value > prevVal){
			pis.push_back(std::pair<int, double>(temp_len, temp_prob));
			long long int j = (*refEnergies)[k].index;
			temp_prob = stateVect->stateVec.real[j]*stateVect->stateVec.real[j]+stateVect->stateVec.imag[j]*stateVect->stateVec.imag[j];;
			temp_len=(*refEnergies)[k].value;
			prevVal = (*refEnergies)[k].value;
		}

		std::cerr<<(*refEnergies)[k].index<<": "<< (*refEnergies)[k].value << std::endl;
	}pis.push_back(std::pair<int, double>(temp_len, temp_prob)); //add the last entry

	if(pis[0].first != 0)
		throw_runtime_error("First successive minima should be zero minima!");

	int k = 50;
	double expectation=0;
	double prob_sum=0;
	double p0=pis[0].second;
	bool zero_minima=true;
	double inner_sum=1-p0;
	for(auto &suc_minima: pis){

		int len = suc_minima.first;
		double pi = suc_minima.second;
		prob_sum += pi;

		if(zero_minima){
			zero_minima=false;
			continue;
		}

		inner_sum -= pi;
		if(inner_sum < 0){
			if(inner_sum < -0.00000000001 || len != pis[pis.size()-1].first)
				logw("inner_sum="+std::to_string(inner_sum)+" < 0, setting to 0", 2);
			inner_sum = 0;
		}

		expectation += len*(pow(inner_sum+p0+pi, k)-pow(inner_sum, k)-pow(p0, k));

	}
	if(prob_sum < 0.9999)
		loge("Probabilities not summing to one!");

	std::cerr<<"expectation: "<<expectation<<"   sqrt="<<sqrt(expectation)<<std::endl;
	std::cerr<<"volume_1n: "<<volume_1n<<std::endl;
	std::cerr<<"hermite factor: "<<pow((sqrt(expectation)/volume_1n),1./(double)n)<<"^n"<<std::endl;

	throw;
	/*for(long long int j = 0; j < pis.size; ++j){

		double prob = stateVect->stateVec.real[j]*stateVect->stateVec.real[j]+stateVect->stateVec.imag[j]*stateVect->stateVec.imag[j];

	}*/


	return 0;
}
