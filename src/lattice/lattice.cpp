/*
 * lattice.cpp
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */


#include "lattice.h"

bool pen_initialized = false;

void Lattice::generate_qubo(bool print){

	expression_qubo = new FastVQA::Expression(*expression_penalized);
	expression_qubo->name = "QUBO";

	int qubit = 0;
	for(auto &var : expression_qubo->getVariables()){

		if(var->id < 0) //id
			continue;

		std::pair<int, std::string> z = expression_qubo -> addZ(qubit);
		//qubo_to_bin_map.emplace(z.second, var);

		std::map<int, mpq_class> subs_expr; //id, coeff

		//(1-z)/2
		subs_expr.emplace(-1, 0.5);
		subs_expr.emplace(z.first, -0.5);
		expression_qubo->substitute(var->id, subs_expr);

		qbit_to_varId_map.emplace(qubit, var->id);

		qubit++;

	}

	if(print)
		expression_qubo -> print();

}

void Lattice::penalize_expr(int penalty, MapOptions::penalty_mode mode, bool print){

	expression_penalized = new FastVQA::Expression(*expression_bin);
	expression_penalized->name = "expression_penalized";

	if(mode == MapOptions::penalty_all){

		std::vector<FastVQA::Var*>::iterator xn_m1_it;

		expression_penalized -> addConstant(penalty);
		std::vector<FastVQA::Var*> variables = expression_penalized->getVariables();
		std::vector<FastVQA::Var*> variables_to_be_penalized;

		if(variables[0]->id != -1){
			loge("Error! id not the first val");
			return;
		}

		int counter = 0;
		int num_penalized_vars = 0;
		for(std::vector<FastVQA::Var*>::iterator it = variables.begin() + 1;
				it != variables.end(); ++it){

			if((*it)->extra_information=="P0" || (*it)->extra_information=="P1"){
				variables_to_be_penalized.push_back((*it));
				num_penalized_vars++;
			}
		}

		for(std::vector<FastVQA::Var*>::iterator it = variables_to_be_penalized.begin()/* + 1*/;
				it != variables_to_be_penalized.end(); ++it){

			//if((*it)->extra_information!="P0" || (*it)->extra_information!="P1")
			//	continue;

			int z_id = expression_penalized -> addBinaryVar("z_"+(*it)->name);

			if(counter == num_penalized_vars-1)
				zn_id = z_id;
			else if(counter == num_penalized_vars-2){
				zn_m1_id = z_id;
				xn_m1_it = it;
			}
			if((*it)->extra_information=="P0"){
				expression_penalized -> addNewTerm((*it)->id, z_id, -penalty);
			}else{ //(*it)->extra_information[0]=="P1"
				expression_penalized -> addNewTerm(-1, z_id, -penalty);
				expression_penalized -> addNewTerm((*it)->id, z_id, penalty);
			}

			for(std::vector<FastVQA::Var*>::iterator it2 = it+1; it2 != variables_to_be_penalized.end(); it2++){
				if((*it)->extra_information=="P0"){
					expression_penalized -> addNewTerm((*it2)->id, z_id, penalty);
				}else{ //(*it)->extra_information[0]=="P1"
					expression_penalized -> addTermCoeff(-1, z_id, penalty);
					expression_penalized -> addNewTerm((*it2)->id, z_id, -penalty);
				}
			}

		counter++;
		}

		//add z_n=1, z_n-1=x_n-1
		expression_penalized->substituteVarToDouble(zn_id, 1);
		std::map<int, mpq_class> subs_expr;//id, coeff

		if(num_penalized_vars > 1){
			xn_m1_id = (*xn_m1_it)->id;
			if((*xn_m1_it)->extra_information=="P0")
				subs_expr.emplace(xn_m1_id, 1);
			else{ //(*xn_m1_it)->extra_information=="P1")
				subs_expr.emplace(-1, 1);
				subs_expr.emplace(xn_m1_id, -1);
			}
			expression_penalized->substitute(zn_m1_id, subs_expr);
		}
		//std::cout << "subs " << z1_id << " c" << 1 << "\n";
		//std::cout << "subs " << z2_id << " " << (*x2_it)->id << "\n";
	}else if(mode == MapOptions::no_hml_penalization){
		//do nothing as no penalty qubits needed
	}

	if(print)
		expression_penalized->print();

}


void Lattice::__single_variable_test(MapOptions* options){
	init_x(options->x_mode, options->num_qbits_per_x, options->absolute_bound, true, true);
	init_expr_bin(options->bin_map, true);

}


void Lattice::init_x(MapOptions::x_init_mode mode, int num_qbits_per_x, int absolute_bound, bool print, bool testing_single_var){

	Z_NR<mpz_t> coeff;

	auto addVar = [&](int k){
		if(mode == MapOptions::x_symmetric){

			if(num_qbits_per_x == 1){
					int id = expression_int->addBinaryVar("x"+std::to_string(k));
					x_ids.push_back(id);
			}else if(num_qbits_per_x > 1){

				int lb = -pow(2, num_qbits_per_x)/ 2 + 1;
				int ub = 1-lb;

				int id = expression_int->addIntegerVar("x"+std::to_string(k), lb, ub);
				x_ids.push_back(id);
			}
			else if(absolute_bound != -1){
				int lb = -absolute_bound;
				int ub = absolute_bound;

				int id = expression_int->addIntegerVar("x"+std::to_string(k), lb, ub);
				x_ids.push_back(id);
			}else{
				throw_runtime_error("bounds undefined");
			}
		}
	};

	if(testing_single_var){
		addVar(0);
		expression_int->addNewTerm(-1, expression_int->getId("x"+std::to_string(0)), 1);
		expression_int->print();
		return;
	}


	for(int i = 0; i < n_rows; ++i){

		if(gramian){
			if(gramian_diag){
				mpz_class c(diagonalGramian(i));
				if(c!=0){
					addVar(i);
					expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(i)), c); // G_ii*x_i^2
					//uncomment below only for test purposes
					//std::cerr<<"CHANGE THIS\n";
					//expression_int->addNewTerm(-1, expression_int->getId("x"+std::to_string(i)), c);
				}
			}
			else{
				mpz_class c(nonDiagGramian(i,i));
				if(c!=0){
					addVar(i);
					expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(i)), c);
				}
				for(int j = 0; j < i; ++j){
					mpz_class c(nonDiagGramian(i,j));
					if(c!=0){
						try {expression_int->getId("x"+std::to_string(i));}
						catch (const std::out_of_range& e) {addVar(i);}
						try {expression_int->getId("x"+std::to_string(j));}
						catch (const std::out_of_range& e) {addVar(j);}
						expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(j)), 2*c); //2*G_ij*xi
					}
				}
			}
		}else{
			gso_current->get_int_gram(coeff, i, i);
			mpz_class c(coeff.get_data());
			if(c!=0){
				addVar(i);
				expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(i)), c/*.coeff(i, i)*/); // G_ii*x_i^2
			}
			for(int j = 0; j < i; ++j){
				mpz_class c(gso_current->get_int_gram(coeff, i, j).get_data());
				//std::cout << "i" << i << " j" << j << " c"<<c << "\n";
				if(c!=0){
					try {expression_int->getId("x"+std::to_string(i));}
					catch (const std::out_of_range& e) {addVar(i);}
					try {expression_int->getId("x"+std::to_string(j));}
					catch (const std::out_of_range& e) {addVar(j);}
					expression_int->addNewTerm(expression_int->getId("x"+std::to_string(i)), expression_int->getId("x"+std::to_string(j)), 2*c/*.coeff(i, j)*/); //2*G_ij*xi
				}
			}
		}
	}

	if(print)
		expression_int->print();

}

void Lattice::init_expr_bin(MapOptions::bin_mapping mapping, bool print){

	expression_bin = new FastVQA::Expression(*expression_int);
	expression_bin->name = "expression_bin";

	for(auto &var : expression_bin->getVariables()){

		if(var->id < 0) //identity coeff
			continue;

		int lb = var -> lb;
		int ub = var -> ub;

		std::string name = var -> name;

		std::map<int, mpq_class> subs_expr; //id, coeff

		if(var->isBinary() || mapping == MapOptions::naive_overapprox){

			subs_expr.emplace(-1, lb); //set lb to identity coeff
			for(int i = 0; i < ceil(log2(ub-lb+1)); ++i){

				int id = expression_bin->addBinaryVar(name + "_b"+std::to_string(i), "P0"); //P0 means it is the one to be penalized for being zero
				subs_expr.emplace(id, pow(2, i));

				varId_to_zero_ref_map[id]=0; //all zero as no constant is involved in naive_overapprox
			}

		}else if(mapping == MapOptions::zeta_omega_overapprox || mapping == MapOptions::zeta_omega_exact){
			subs_expr.emplace(-1, lb); //set lb to identity coeff
			int id = expression_bin->addBinaryVar(name + "_zeta", "P1"); //P means it is the one to be penalized for being one
			subs_expr.emplace(id, -lb);
			id = expression_bin->addBinaryVar(name + "_omega");
			subs_expr.emplace(id, -lb+1);
			for(int i = 0; i <= (mapping == MapOptions::zeta_omega_overapprox ? floor(log2((-lb)-1)) : floor(log2((-lb)-1))-1); ++i){
				int id = expression_bin->addBinaryVar(name + "_b"+std::to_string(i));
				subs_expr.emplace(id, pow(2, i));
			}
			if(mapping == MapOptions::zeta_omega_exact && (-lb)>1){
				int id = expression_bin->addBinaryVar(name + "_b"+std::to_string(int(floor(log2((-lb)-1)))));
				subs_expr.emplace(id, -lb-pow(2,floor(log2((-lb)-1))));
			}
		}else{
			loge("NOT IMPLEMENTED");
			throw;
		}

		expression_bin->substitute(var->id, subs_expr);
		int_to_bin_map.emplace(var->id, subs_expr);
	}

	if(print)
		expression_bin->print();
}

void Lattice::calcHamiltonian(MapOptions* options, bool print){

		if(!gso_current_initialized){

			ZZ_mat<mpz_t> blank;

			gso_current = new MatGSO<Z_NR<mpz_t>, FP_NR<double>>(current_lattice, blank, blank, GSO_INT_GRAM);
			gso_current->update_gso();

			gso_current_initialized = true;
		}

		if(!x_initialized){
			init_x(options->x_mode, options->num_qbits_per_x, options->absolute_bound, print);
			x_initialized = true;
		}

		if(!bin_initialized){
			init_expr_bin(options->bin_map, print);
			bin_initialized = true;
		}

		if(!pen_initialized){

			if(gramian){
				logi("Penalty="+std::to_string(options->penalty), print ? 0 : 3);
			}else{
				int first_vect_len = 0;

				for(int i = 0; i < n_cols; ++i){
					first_vect_len+=current_lattice.matrix[0][i].get_si()*current_lattice.matrix[0][i].get_si();
				}
				options->penalty = 5 * first_vect_len;
				logi("Overriding penalty with 5*firstVectLenSq="+std::to_string(options->penalty), print ? 0 : 3);
			}

			penalize_expr(options->penalty, options->pen_mode, print);
			pen_initialized = true;
		}

		if(!qubo_generated){
			generate_qubo(print);
			qubo_generated = true;
		}
}



/*xacc::quantum::PauliOperator Lattice::getHamiltonian(MapOptions* options){

	calcHamiltonian(options, options->verbose);

	std::map<int, std::pair<std::string, std::complex<double>>> operators;

	for(auto &term : expression_qubo->polynomial){

		std::pair<int, int> vars = term.first;
		if(vars.first == -1){

			if(vars.second == -1){
				//operators.emplace(expression_qubo->getQubit(var), std::pair<std::string, std::complex<double>>("Z", std::complex<double>(term.second.get_d(),0)));
			}else{
				operators.emplace(expression_qubo->getQubit(vars.second), std::pair<std::string, std::complex<double>>("Z", std::complex<double>(term.second.get_d(),0)));
			}

		}else{
			//operators.emplace(expression_qubo->getQubit(vars.second), std::pair<std::string, std::complex<double>>("Z", std::complex<double>(term.second.get_d(),0)));
		}


	}

	operators.emplace(0, std::pair<std::string, std::complex<double>>("Z", std::complex<double>(5,0)));
	operators.emplace(1, std::pair<std::string, std::complex<double>>("Z", std::complex<double>(6,0)));

	hamiltonian = xacc::quantum::PauliOperator(operators);

	return hamiltonian;

}*/

/*std::string Lattice::toHamiltonianString(){
	if(!qubo_generated){
		loge("Hamiltonian referenced but not yet generated!");
		return "";
	}
	return expression_qubo->expression_line_print();

}

std::string Lattice::toHamiltonianString(MapOptions* options){

	calcHamiltonian(options, options->verbose);
	return expression_qubo->expression_line_print();

}*/

std::vector<long long unsigned int> Lattice::getZeroReferenceStates(){
	if(!qubo_generated){
		throw_runtime_error("Hamiltonian referenced but not yet generated!");
	}

	logw("TODO: Be careful when returning zero reference states. There can be more than one.");

	long long unsigned int zero_ref_state = 0;

	//q_n-1 ... q_0

	for(int i = 0; i < this->getNumQubits(); ++i){
		int var_id = qbit_to_varId_map[i];
		zero_ref_state += varId_to_zero_ref_map[var_id] == 0 ? 0 : (1<<i);
	}

	return std::vector<long long unsigned int>{zero_ref_state};
}

void Lattice::reduce_rank(int reduced_rank){

	if(reduced_rank < 0 || reduced_rank > this -> current_lattice.get_rows()){
		loge("Invalid reduced rank");
		return;
	}

	this -> n_rows = reduced_rank;//lattice.get_rows();
}

