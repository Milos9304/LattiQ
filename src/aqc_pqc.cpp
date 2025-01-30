#include "experiment_runner.h"
#include <algorithm>

void AqcPqcExperiment::run(FastVQA::AqcPqcAcceleratorOptions* options, int num_instances){

	/*
	 * DEFAULTS
	 *
	const int loglevel = 1;
	const int round_decimals = 5; //-1 undefined
	const int opt_strategy = 0;	  //0=trivially, 1=rank_reduction
	const int num_steps = 20;
	const int ansatz_depth = 1;
	const double xtol = 10e-5;
	const double catol = 0.0002;
	const bool classical_esolver_compare = false;
	const bool outputLogToFile = true;
	const bool checkHessian = true;
	const bool printGroundStateOverlap = true;
	const bool print_eps = false;
	const int eval_limit_step = 600; //max iterations per step*/




	FastVQA::AqcPqcAcceleratorOptions acceleratorOptions = *options;

	/*
	 * DEFAULTS
	 *
	 * acceleratorOptions.log_level = loglevel;
	acceleratorOptions.logFileName = "aqc_pqc_log.txt";
	acceleratorOptions.roundDecimalPlaces = round_decimals;
	acceleratorOptions.optStrategy = opt_strategy;
	acceleratorOptions.accelerator_type = "quest";
	acceleratorOptions.nbSteps = num_steps;
	acceleratorOptions.ansatz_name = "Ry_Cz_nn_Ry";//"Ry_CNOT_nn_Rz_CNOT_Rz"
	acceleratorOptions.ansatz_depth = ansatz_depth;
	acceleratorOptions.xtol = xtol;
	acceleratorOptions.catol = catol;
	acceleratorOptions.compareWithClassicalEigenSolver = classical_esolver_compare;
	acceleratorOptions.outputLogToFile = outputLogToFile;
	acceleratorOptions.checkHessian = checkHessian;
	acceleratorOptions.printGroundStateOverlap = printGroundStateOverlap;
	acceleratorOptions.printEpsilons = print_eps;
	acceleratorOptions.eval_limit_step = eval_limit_step;
	acceleratorOptions.initialGroundState = FastVQA::InitialGroundState::PlusState;*/



	//this->mapOptions->penalty = 0;

	std::vector<int> num_iters;
	std::vector<double> final_overlaps, first_excited_overlaps;

	this->max_num_instances = num_instances;//80;

	for(int m = m_start; m <= m_end; ++m){

		loge("m="+std::to_string(m));

		//std::vector<Instance> dataset = _generate_dataset(3, m, true); //3 is arbitrary, true is for penalise
		//this->max_num_instances = 1;
		//loge("Max num instances is 1 instead of 100");
		std::vector<Instance> dataset = _generate_dataset(1, m, this->aqcpqc_penalised);

		int i = 0;

		if(this->aqcpqc_penalised){
		   std::ofstream myfile("dim="+std::to_string(m)+"_steps="+std::to_string(options->nbSteps));
		   //for(auto &o : final_overlaps){
		    //	   myfile<<o/*<<","*/<<std::endl;
		   //}
           myfile.close();
		}else{
			std::ofstream myfile("zero_dim="+std::to_string(m)+"_steps="+std::to_string(options->nbSteps));
			//for(auto &o : final_overlaps){
			//	   myfile<<o/*<<","*/<<std::endl;
			//}
			myfile.close();

			std::ofstream myfile2("sv_dim="+std::to_string(m)+"_steps="+std::to_string(options->nbSteps));
			//for(auto &o : first_excited_overlaps){
			//		myfile2<<o/*<<","*/<<std::endl;
			//}
            myfile2.close();
		}

		for(auto &instance: dataset){

			if(i > /*80*/num_instances && num_instances > 0){
				logw("aqcpqc.cpp breaking after "+std::to_string(num_instances)+"instances");
				break;
			}

			if(num_instances == -1)
				logi("Instance " + std::to_string(i+1) + "/"+std::to_string(dataset.size()));
			else
				logi("Instance " + std::to_string(i+1) + "/"+std::to_string(num_instances));


			logi(std::to_string(instance.h.nbQubits) + " qubits");

			std::vector<long long int> solutions;
			for(auto &sol: instance.sv_solutions){
				solutions.push_back(sol.index);
			}
			if(this->aqcpqc_penalised){
				acceleratorOptions.solutions = solutions;
			}
			else{
				std::vector<long long int> zero_solutions;
				for(auto &sol: instance.zero_solutions){
					zero_solutions.push_back(sol.index);
				}

				acceleratorOptions.solutions = zero_solutions;
				acceleratorOptions.first_excited_states = solutions;
			}

			FastVQA::AqcPqcAccelerator accelerator(acceleratorOptions);

			/*logw("Num sols: " + std::to_string(acceleratorOptions.solutions.size()));
			for(const auto &sol: solutions){
				std::cerr<<"solution index="<<sol<<" with value="<<instance.sv1Squared<<std::endl;
			}

			logw("SV1Squared: " + std::to_string(instance.sv1Squared));
			for(const auto &sol: instance.zero_solutions){
				std::cerr<<"zero:"<<sol.value<<" "<<sol.index<<std::endl;
			}*/

			FastVQA::PauliHamiltonian h0(instance.h.nbQubits);


			Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> matrix = instance.h.getMatrixRepresentation2(true);
			std::vector<double> energies;
			/*for(int j = 0; j < matrix.cols(); ++j){
				std::cerr<<j<<":  "<<matrix(j,j)<<std::endl;
				energies.push_back(matrix(j,j));
			}

			std::sort(energies.begin(), energies.end());
			for(auto &a:energies)
				std::cerr<<a<<" ";*/

			h0.initializeSumMinusSigmaXHamiltonian();

			FastVQA::AqcPqcAcceleratorResult result;
			accelerator.initialize(&h0, &instance.h);
			accelerator.run(&result);

			if(this->aqcpqc_penalised){
				std::cerr<<"Overlap = "<<result.final_state_overlap<<std::endl;
				final_overlaps.push_back(result.final_state_overlap);
			}
			else{
				std::cerr<<"Zero overlap = "<<result.final_state_overlap<<std::endl;
				std::cerr<<"SV overlap = "<<result.first_exc_state_overlap<<std::endl;
				final_overlaps.push_back(result.final_state_overlap);
				first_excited_overlaps.push_back(result.first_exc_state_overlap);
			}

			if(this->aqcpqc_penalised){
                        	std::ofstream myfile("dim="+std::to_string(m)+"_steps="+std::to_string(options->nbSteps), std::ios::app);
                        	//for(auto &o : final_overlaps){
                	                myfile<<final_overlaps[final_overlaps.size()-1]/*<<","*/<<std::endl;
        	                //}
                            myfile.close();
	                }else{
                        	std::ofstream myfile("zero_dim="+std::to_string(m)+"_steps="+std::to_string(options->nbSteps), std::ios::app);
                        	//for(auto &o : final_overlaps){
                                	myfile<<final_overlaps[final_overlaps.size()-1]/*<<","*/<<std::endl;
                        	//}
            				myfile.close();

                        	std::ofstream myfile2("sv_dim="+std::to_string(m)+"_steps="+std::to_string(options->nbSteps), std::ios::app);
            //            	for(auto &o : first_excited_overlaps){
                                	myfile2<<final_overlaps[final_overlaps.size()-1]/*<<","*/<<std::endl;
                        	//}
                            myfile2.close();
                	}

			i++;
		
		}


		logw("Breaking after first m");
		break;
	}

}
