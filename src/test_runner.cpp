#include "FastVQA/fastVQA.h"
#include "lattice/lattice.h"
#include "test_runner.h"
#include "averaged_svp.h"

std::string _experiment_output_directory = "../experiments";

void test_variable_substitution(MapOptions* mapOptions){

	Lattice l("variable_test", mapOptions);

}

void test_execution_time(FastVQA::QAOAOptions* qaoaOptions, Database* database){

	int penalty = 0;

	MapOptions mapOptions;
	mapOptions.verbose = false;//true;
	//mapOptions.num_qbits_per_x = qubits_per_x->value();
	//mapOptions.absolute_bound = absolute_bound->value();
	//mapOptions.pen_mode = MapOptions::penalty_all;
	mapOptions.bin_map = penalty > 0 ? MapOptions::zeta_omega_exact : MapOptions::naive_overapprox;
	mapOptions.penalty = penalty;

	qaoaOptions->accelerator->options.exclude_zero_state = true;
	qaoaOptions->log_level = 3;
	qaoaOptions->ftol = 10e-6;
	qaoaOptions->max_iters = 4000;

	//int m_max = 10;
	int num_samples = 10;
	bool shuffle = false;

		GeneratorParam param(7,1,2, shuffle, 0, num_samples);
		std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

		for(int qs=1; qs<=11; ++qs){

			mapOptions.num_qbits_per_x = qs;

			int counter=0;
			for(auto w : gramian_wrappers){
				for(int p = 1; p <= 6; ++p){
					for(int odd=0; odd < 2; ++odd){

						if(qs == 1 && odd)
							continue;

						mapOptions.__minus_one_qubit_firstvar=odd;
						qaoaOptions->p = p;

						//std::string filename = _experiment_output_directory+"/performance_experiment/test_interim_energies_q"+std::to_string(qs*2+(odd?(-1):0))+"_"+std::to_string(counter)+"_p"+std::to_string(p);
						//std::ifstream f(filename);

						/*if(f.good())
							continue;
						*/
						//std::ofstream interimEnergiesFile(filename);

						Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> G = w.hamiltonian;
						std::string name = w.name;

						Lattice l(G, name);

						//std::cerr<<"G"<<G<<"\n";
						FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);
						int nbQubits=h.nbQubits;

						logi("Running experiment " + std::to_string(counter) + " with " + std::to_string(nbQubits) + " qubits and p="+std::to_string(p));

						FastVQA::Qaoa qaoa_instance;
						FastVQA::ExperimentBuffer buffer;
						buffer.storeQuregPtr = false;

						//time_t start_time = time(0);
						//qaoa_instance.run_qaoa(&buffer, &h, qaoaOptions);
						//time_t end_time = time(0);

						//int duration_s = difftime(end_time,start_time);

						Database::DatasetRow row;
						bool found = database->getOrCalculate_qary(param.q, param.n, param.m, qaoaOptions->p,
								counter, nbQubits, penalty > 0 ? true : false,
										&l, &h, &row, qaoaOptions, &mapOptions);

						/*if(found){
							std::cerr<<row.q<<"\n";
							std::cerr<<row.n<<"\n";
							std::cerr<<row.m<<"\n";
							std::cerr<<row.p<<"\n";
							std::cerr<<row.index<<"\n";
							std::cerr<<row.num_qs<<"\n";
							std::cerr<<row.penalty<<"\n";
							std::cerr<<row.volume<<"\n";
							std::cerr<<row.sv1Squared<<"\n";
							std::cerr<<row.degeneracy<<"\n";
							std::cerr<<row.duration_s<<"\n";
							std::cerr<<row.iters<<" iters\n";
							for(auto&a: row.initAngles)
								std::cerr<<a<<" ";
							std::cerr<<"\n\n";
							for(auto&a: row.finalStateVectorMap)
								std::cerr<<a.first<<" "<<a.second<<"\n";
							std::cerr<<"\n\n";
							for(auto&a: row.intermediateAngles){
								std::cerr<<"{";
								for(auto&aa: a)
									std::cerr<<aa<<" ";
								std::cerr<<"}";
							}
							std::cerr<<"\n\n";
							for(auto&a: row.intermediateEnergies)
								std::cerr<<a<<" ";
							std::cerr<<"\n\n";
							for(auto&a: row.finalAngles)
								std::cerr<<a<<" ";
							std::cerr<<"\n\n";
							std::cerr<<row.probSv1<<" probability\n";
							std::cerr<<row.opt_res<<"\n";
							std::cerr<<row.comment<<"\n";
							return;
						}*/

						/*row.type = "qary";
						row.q = param.q;
						row.n = param.n;
						row.m = param.m;
						row.p = qaoaOptions->p;
						row.penalty = mapOptions.penalty;
						database->write(&row);*/




						//std::cerr<<buffer.num_iters<<"\n";
						//std::cerr<<buffer.opt_message<<"\n";
						//interimEnergiesFile << "Duration=" << duration_s <<"\n";
						//for(double e: buffer.intermediateEnergies)
						//	interimEnergiesFile << e << " ";
						//interimEnergiesFile.close();
					}

				}
				counter++;
			}


	}

}
