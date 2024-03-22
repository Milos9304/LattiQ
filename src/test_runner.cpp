#include "FastVQA/fastVQA.h"
#include "lattice/lattice.h"
#include "test_runner.h"
#include "averaged_svp.h"

std::string _experiment_output_directory = "../experiments";

void test_variable_substitution(MapOptions* mapOptions){

	Lattice l("variable_test", mapOptions);

}

struct PPerformanceData{
	int dim_min;
	int dim_max;

	std::map<int, std::vector<double>> data;
	bool hasQ(int q){
		if(data.count(q) == 1)
			return true;
		return false;
	}
	std::pair<double, double> getMeanAndStdev(int dim){
		double sum = std::accumulate(data[dim].begin(), data[dim].end(), 0.0);
		double mean = sum / data[dim].size();
		double sq_sum = std::inner_product(data[dim].begin(), data[dim].end(), data[dim].begin(), 0.0);
		double stdev = std::sqrt(sq_sum / data[dim].size() - mean * mean);

		return std::pair<double, double>(mean, stdev);
	}

};

void test_execution_time(FastVQA::QAOAOptions* qaoaOptions, Database* database){

	int penalty = 0;
	const int p_min = 1;
	const int p_max = 6;

	const int qs_max = 22; //even number!

	MapOptions mapOptions;
	mapOptions.verbose = false;//true;
	//mapOptions.num_qbits_per_x = qubits_per_x->value();
	//mapOptions.absolute_bound = absolute_bound->value();
	//mapOptions.pen_mode = MapOptions::penalty_all;
	mapOptions.bin_map = penalty > 0 ? MapOptions::zeta_omega_exact : MapOptions::naive_overapprox;
	mapOptions.penalty = penalty;
	mapOptions.loglevel = qaoaOptions->log_level;

	qaoaOptions->accelerator->options.exclude_zero_state = true;
	qaoaOptions->log_level = 3;
	qaoaOptions->ftol = 10e-6;
	qaoaOptions->max_iters = 4000;

	//int m_max = 10;
	int num_samples = 10;
	bool shuffle = false;

	GeneratorParam param(7,1,2, shuffle, 0, num_samples);
	std::vector<HamiltonianWrapper> gramian_wrappers = generateQaryUniform(param);

	PPerformanceData pPerformanceData[p_max]; //index by p=qaoaDepth

	ProgressBar bar{
			option::BarWidth{50},
			option::MaxProgress{qs_max/2*gramian_wrappers.size()*(p_max-p_min+1)*2-qs_max/2},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{"Getting Sv1Mean probability"},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::ForegroundColor{Color::yellow},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		};

	for(int qs=1; qs<=qs_max/2; ++qs){

		mapOptions.num_qbits_per_x = qs;

		int counter=0;
		for(auto w : gramian_wrappers){
			for(int p = p_min; p <= p_max; ++p){
				for(int odd=0; odd < 2; ++odd){

					if(qs == 1 && odd)
						continue;
					bar.tick();

					/*mapOptions.__minus_one_qubit_firstvar=odd;
					qaoaOptions->p = p;

					Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> G = w.hamiltonian;
					std::string name = w.name;

					Lattice l(G, name);

					FastVQA::PauliHamiltonian h = l.getHamiltonian(&mapOptions);
					int nbQubits=h.nbQubits;

					logi("Running experiment " + std::to_string(counter) + " with " + std::to_string(nbQubits) + " qubits and p="+std::to_string(p), qaoaOptions->log_level);

					FastVQA::Qaoa qaoa_instance;
					FastVQA::ExperimentBuffer buffer;
					buffer.storeQuregPtr = false;

					Database::DatasetRow row;
					bool found = database->getOrCalculate_qary(param.q, param.n, param.m, qaoaOptions->p,
							counter, nbQubits, penalty > 0 ? true : false,
									&l, &h, &row, qaoaOptions, &mapOptions);

					pPerformanceData[p-1].data[nbQubits].push_back(row.probSv1);*/

					double prob = database->getSv1Probability(param.q, param.n, param.m, p,
							param.m*qs-odd, counter);
					//std::cerr<<prob<<std::endl;

					pPerformanceData[p-1].data[param.m*qs-odd].push_back(prob);

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

	//std::cout<< "ps: " << p_min << " - " << p_max  << std::endl;
	//std::cout<< "qs: " << 1     << " - " << qs_max << std::endl;
	std::cout<<"{";
	for(int p = p_min; p <= p_max; ++p){
		std::cout <<","<<p<<": {";
		for(int q = 1; q <= 30; ++q){
			if(pPerformanceData[p-1].hasQ(q)){
				std::pair<double, double> meanStdev = pPerformanceData[p-1].getMeanAndStdev(q);
				std::cout<<","<<q<<": ["<<meanStdev.first<<","<<meanStdev.second<<"]";
			}
		}std::cerr<<"}";
	}
	std::cout<<"}";

}
