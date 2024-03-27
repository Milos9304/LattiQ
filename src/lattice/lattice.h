/*
 * lattice.h
 *
 *  Created on: Apr 2, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICE_H_
#define SRC_LATTICE_H_

#include <gmpxx.h>
#include "fplll.h"
#include "../io/logger.h"
#include "FastVQA/symbolic_manipulation.h"
#include "FastVQA/pauliHamiltonian.h"
#include <vector>

typedef fplll::ZZ_mat<mpz_t> MatrixInt;
typedef std::vector<mpz_class> VectorInt;
typedef Eigen::Vector<int, Eigen::Dynamic> DiagonalHamiltonian;

class MapOptions{

	public:

		bool verbose;
		int loglevel;

		/*
		 * x_symmetric: x_i \ in [-2^(B-1) + 1 , 2^(B-1)]
		 *
		 * */
		enum x_init_mode { x_symmetric };
		enum bin_mapping { naive_overapprox, zeta_omega_overapprox, zeta_omega_exact };

		/*
		 * penalty_all: apply penalty
		 * no_hml_penalization
		 *
		 *
		 */
		//enum penalty_mode { penalty_all, no_hml_penalization };

		x_init_mode x_mode;
		bin_mapping bin_map;
		//penalty_mode pen_mode;

		int penalty;
		int num_qbits_per_x;
		int absolute_bound;

		bool __minus_one_qubit_firstvar=false; //only for testing purposes!

		MapOptions(x_init_mode x_mode=x_symmetric,
				bin_mapping bin_map=naive_overapprox,
				//penalty_mode pen_mode=penalty_all,
				int penalty_val=1000,
				int num_qbits_per_x=1,
				int absolute_bound=-1,
				bool verbose=false){

			if(num_qbits_per_x != -1 && bin_map != naive_overapprox){
				throw_runtime_error("Not implemented");
			}
			if(absolute_bound != -1 && bin_map != zeta_omega_overapprox && bin_map != zeta_omega_exact){
				throw_runtime_error("Not implemented");
			}

			this->x_mode = x_mode;
			this->bin_map = bin_map;
			//this->pen_mode = pen_mode;
			this->penalty = penalty_val;
			this->num_qbits_per_x = num_qbits_per_x;
			this->absolute_bound = absolute_bound;
			this->verbose = verbose;
		}
};

class Lattice {

	public:

	std::string name;
	int n_rows, n_cols;

	MatrixInt* lll_transformation;
	bool lll_preprocessed;

	//long long int svLenSquared;

		Lattice(std::string s, MapOptions* mapOptions){
			if(s!="variable_test")
				throw_runtime_error("Invalid lattice initialization");

     		this -> expression_int = new FastVQA::Expression("expression_int");
			this -> __single_variable_test(mapOptions);
		}

		Lattice(MatrixInt lattice, std::string name = "", int reduced_rank=0){ // @suppress("Class members should be properly initialized")

			this -> n_rows = lattice.get_rows();
			this -> n_cols = lattice.get_cols();

			this -> name = name;

			this->gramian = false;
			if(reduced_rank != 0){
				this -> n_rows = reduced_rank;//lattice.get_rows();
				if(reduced_rank < 0 || reduced_rank > lattice.get_rows()){
					loge("Invalid reduced rank");
					return;
				}
			}

			this -> lll_preprocessed = false;
			this -> orig_lattice = lattice;
			this -> orig_lattice_transposed = MatrixInt(lattice);
			this -> orig_lattice_transposed.transpose();
			this -> current_lattice = MatrixInt(lattice);
			this -> orig_gh_sq = calculate_gh_squared(&orig_lattice);

			loge("firstVectorLengthSquared calculation not yet implemented");
			this->firstVectorLengthSquared=0;

			//if(lattice.get_rows()/*.rows()*/ != lattice.get_cols()/*.cols()*/){
			//	loge("Non-square lattice not supported");
			//	return;
			//}

     		this -> expression_int = new FastVQA::Expression("expression_int");
		}

		Lattice(DiagonalHamiltonian gramian, std::string name = ""){

			this -> gramian = true;
			this -> gramian_diag = true;
			this -> diagonalGramian = gramian;
			this -> n_rows = gramian.size();
			this -> n_cols = 1;
			this -> name = name;

			gso_current_initialized = true;
			this->firstVectorLengthSquared=gramian(0,0);

			this -> expression_int = new FastVQA::Expression("expression_int");
		};

		Lattice(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> gramian, std::string name = ""){

			this -> gramian = true;
			this -> gramian_diag = false;
			this -> nonDiagGramian = gramian;
			this -> n_rows = gramian.rows();
			this -> n_cols = gramian.cols();
			this -> name = name;

			gso_current_initialized = true;
			this->firstVectorLengthSquared=gramian(0,0);

			this -> expression_int = new FastVQA::Expression("expression_int");
		};

		void reduce_rank(int reduced_rank);

		//decode qubo optimal config to x config
		VectorInt quboToXvector(long long int measurement, int nbQubits);
		VectorInt quboToXvector(std::string measurement);
		VectorInt quboToXvector(bool* measurement, int n);

		int getNumQubits(){return expression_qubo->getIdMapSize()-1;}
		double getVolume();

		long long int calculateSVLength();

		mpq_class get_orig_gh(){return orig_gh_sq;}

		MatrixInt* get_orig_lattice(){ return &orig_lattice; }
		MatrixInt* get_orig_lattice_transposed(){ return &orig_lattice_transposed; }
		MatrixInt* get_current_lattice(){ return &current_lattice; }

		FastVQA::PauliHamiltonian getHamiltonian(MapOptions* options);
		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> getMatrix(){
			if(gramian && gramian_diag)
				return diagonalGramian;
			if(gramian && !gramian_diag)
				return nonDiagGramian;
			throw_runtime_error("Non implemented");
			//return orig_lattice;
		}

		/*std::string toHamiltonianString();
		std::string toHamiltonianString(MapOptions* options);*/

		//classical state that corresponds to zero eigen-value of Hamiltonian
		std::vector<long long unsigned int> getZeroReferenceStates();

		void outputGramianToFile(std::string filename);

		int getSquaredLengthOfFirstBasisVector(){
			return this->firstVectorLengthSquared;
		}

		long long int firstVectorLengthSquared=0;

		qreal get_random_guess_one_vect(){
			if(bin_initialized)
				return random_guess_one_vect;
			throw_runtime_error("random_guess_one_vect_prior_penalization cannot be returned at this point");return -1;
		}

	private:

		bool gramian;
		bool gramian_diag;
		DiagonalHamiltonian diagonalGramian;
		Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> nonDiagGramian;

		mpq_class orig_gh_sq; //gaussian heuristics
		MatrixInt orig_lattice, orig_lattice_transposed, current_lattice;

		bool gso_current_initialized = false, gso_orig_initialized = false;
		fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>>* gso_current, *gso_orig;

		FastVQA::Expression *expression_int, *expression_bin, *expression_penalized, *expression_qubo;
		std::vector<std::map<FastVQA::Var*, int>> solutions;
		bool solutions_calculated = false;
		mpq_class solutions_length_squared=-1;
		void calculate_solutions(bool allow_for_zero_ground_state, bool print=false);
		void _bruteForceSolutions(int n, std::map<FastVQA::Var*, int> *varBoolMap, int i, bool allow_for_zero_ground_state);

		//std::map<std::string, Var*> qubo_to_bin_map;
		std::map<int, int> qbit_to_varId_map;

		bool x_initialized = false;
		void init_x(MapOptions::x_init_mode mode, int num_qbits_per_x, int exponent_bound, bool print=false, bool testing_single_var=false, bool __minus_one_qubit_firstvar=false);

		std::vector<int> x_ids;
		std::map<int, std::map<int, mpq_class>> int_to_bin_map;

		// x_i->c+sum(c_i*x'_i);
		// specifies values for x'_i s.t. x_i=0
		std::map<int, int> varId_to_zero_ref_map;
		bool bin_initialized = false;
		void init_expr_bin(MapOptions::bin_mapping mapping, bool print=false);

		qreal random_guess_one_vect;

		int zn_id=0, zn_m1_id=0, xn_m1_id;
		bool pen_initialized = false;
		void penalize_expr(int penalty, /*MapOptions::penalty_mode mode, */bool print=false);

		bool qubo_generated = false;
		void generate_qubo(bool print=false);

		void calcHamiltonian(MapOptions* options, bool print);

		mpq_class calculate_gh_squared(MatrixInt* lattice);

		void __single_variable_test(MapOptions*);


};


#endif /* SRC_LATTICE_H_ */
