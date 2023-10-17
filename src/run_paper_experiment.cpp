#include "run_paper_experiment.h"

#include <fstream>
#include <utility>
#include "io/logger.h"

std::shared_ptr<SolutionDataset> __read_experiment_file(std::string experiment_name, int instance_select){

	std::ifstream expFile("../experiments/"+experiment_name+".csv");
	std::ifstream matrixFile("../experiments/"+experiment_name+"_matrices.csv");
	if(expFile.fail() || matrixFile.fail()){
		loge("Failed opening experiments file");
	}
	std::shared_ptr<SolutionDataset> solutionDataset;//(num_ranks, rank_min, dim);

	std::string line;
	int inst_counter = 0;

	int header_index = 0;
	int num_ranks, rank_min, num_instances, dim;

	while(std::getline(expFile, line)){

		if(header_index == 0){
			if(line.substr(0,10) != "num_ranks ")
				throw_runtime_error("Invalid file header format");

			std::istringstream s(line);
			std::string temp, int_val;
			s >> temp >> int_val;

			try{num_ranks = std::stoi(int_val);}
			catch(...){throw_runtime_error("Invalid file header format");}

			header_index++;
			continue;
		}else if(header_index == 1){
			if(line.substr(0,9) != "rank_min ")
				throw_runtime_error("Invalid file header format");

			std::istringstream s(line);
			std::string temp, int_val;
			s >> temp >> int_val;

			try{rank_min = std::stoi(int_val);}
			catch(...){throw_runtime_error("Invalid file header format");}

			header_index++;
			continue;
		}else if(header_index == 2){
			if(line.substr(0,14) != "num_instances ")
				throw_runtime_error("Invalid file header format");

			std::istringstream s(line);
			std::string temp, int_val;
			s >> temp >> int_val;

			try{num_instances = std::stoi(int_val);}
			catch(...){throw_runtime_error("Invalid file header format");}

			solutionDataset = std::make_shared<SolutionDataset>(num_ranks, rank_min, dim);
			header_index++;
			continue;
		}else if(header_index == 3){
			if(line.substr(0,4) != "dim ")
				throw_runtime_error("Invalid file header format");

			std::istringstream s(line);
			std::string temp, int_val;
			s >> temp >> int_val;

			try{dim = std::stoi(int_val);}
			catch(...){throw_runtime_error("Invalid file header format");}

			solutionDataset = std::make_shared<SolutionDataset>(num_ranks, rank_min, dim);
			header_index++;
			continue;
		}

		if(instance_select != -1 && instance_select != inst_counter){
			inst_counter++;
			continue;
		}

		std::istringstream s(line);
		std::string field;

		while (getline(s, field,';')){ //iterate over different ranks

			Solution sol;

			//loge(s.str());
			std::string solution = field.substr(1, field.size()-3);

			int i = 0;
			while(solution[i++]!=',');

			double svLength = std::stod(solution.substr(0,i-1));
			std::string rest = solution.substr(i+2, solution.size()-1);

			std::vector<int> vectOfCoeffs;
			std::stringstream ss(rest);
		    for (std::string i; ss >> i;) {
		    	if(i[i.size()-1]==',')
		    		i=i.substr(0,i.size()-1);
		    	vectOfCoeffs.push_back(std::stoi(i));
		        while(ss.peek() == ' ')
		            ss.ignore();
		    }

		    sol.lattice_id = inst_counter;
		    sol.rank = vectOfCoeffs.size();
		    sol.svLength = svLength;
		    sol.coeffs = vectOfCoeffs;
		    solutionDataset->addDataset(sol);

		}
		inst_counter++;
	}

	int matrix_counter = 0;
	while(std::getline(matrixFile, line)){

		if(instance_select != -1 && instance_select != matrix_counter){
			matrix_counter++;
			continue;
		}

		MatrixInt lattice;
		int rows = solutionDataset->rank_min+solutionDataset->num_ranks;
		int cols = solutionDataset->dim;
		lattice.resize(rows, cols);

		std::istringstream s(line);
		std::string field;

		int row_index = 0;
		while (getline(s, field,';')){ //iterate over martix rows

			//loge(s.str());
			//logw(field);

			//std::vector<int> vectOfCoeffs;
			std::stringstream ss(field.substr(1,field.size()-2));
			int col_index = 0;
			for (std::string i; ss >> i;) {
				if(i[i.size()-1]==',')
					i=i.substr(0,i.size()-1);
				//vectOfCoeffs.push_back(std::stoi(i));
				//loge(i);
				lattice(row_index, col_index++) = std::stoi(i);
				while(ss.peek() == ' ')
					ss.ignore();
			}
			row_index++;

			//lattice(matrix_counter, col_index++) = n;
			//return;

		}
		solutionDataset->addLattice(lattice);
		matrix_counter++;
	}

	expFile.close();
	matrixFile.close();

	return solutionDataset;
}

void initialize_paper_experiment(std::string experiment_name, std::vector<Lattice*> &lattices, int rank_select, int instance_select){

	std::shared_ptr<SolutionDataset> solutionDataset = __read_experiment_file(experiment_name, instance_select);
	std::pair<std::vector<MatrixInt>, std::vector<Solution>> dataset = solutionDataset->getMatricexAndDataset();
	std::vector<MatrixInt> matrices = std::get<0>(dataset);
	std::vector<Solution> solutions = std::get<1>(dataset);
	int i = 0;
	i = rank_select - solutionDataset->rank_min;

	if(rank_select == 0){
		logw("-r option missing. Selecting rank " + std::to_string(solutionDataset->rank_min));
		i = 0;
	}

	if(i < 0 || i >= solutionDataset->num_ranks){
		throw_runtime_error("Invalid rank_reduce value.");
	}

	logi("Following instances have been loaded (i_rank) and experiments will run in this order:");
	for(auto &m: matrices){
		Lattice* new_lattice = new Lattice(m, std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank));
		new_lattice->reduce_rank(solutions[i].rank);
		new_lattice->svLenSquared = solutions[i].svLength;
		lattices.push_back(new_lattice);
		std::cout << (std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank)) << " ";

		i+=solutionDataset->num_ranks;

		//logw("Loading only one lattice");
		//break;

	}std::cout << "\n";

}
