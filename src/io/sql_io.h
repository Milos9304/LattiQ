/*
 * saveProgress.hpp
 *
 *  Created on: Mar 06, 2024
 *      Author: Milos Prokop
 */

#ifndef SRC_IO_SQL_IO_HPP_
#define SRC_IO_SQL_IO_HPP_

#include "SQLiteCpp/SQLiteCpp.h"
#include "../averaged_svp.h"
#include "../lattice/lattice.h"

#include "logger.h"

typedef std::vector<std::pair<qreal, double>> FinalStateVectorMap;

class Database{
public:

	enum DATABASE_TYPE{ DATABASE_QARY_PERFORMANCE, DATABASE_ANGLERES };

	struct DatasetRow{
		std::string type;
		int q;
		int n;
		int m;
		int p;
		int index;
		int num_qs;
		int penalty;
		double volume;
		int sv1Squared;
		int degeneracy;
		int duration_s;
		int iters;
		std::vector<double> initAngles;
		FinalStateVectorMap finalStateVectorMap;
		std::vector<std::vector<double>> intermediateAngles;
		std::vector<double> intermediateEnergies;
		std::vector<double> finalAngles;
		double probSv1;
		std::string opt_res;
		std::string comment;
	};

	static void print_sqlite_info(std::string filename);

	Database(std::string filename, DATABASE_TYPE, int loglevel=0);
	bool contains_qary(int q, int n, int m, int p, int index, int num_qs, bool penaltyUsed);
	double getSv1Probability(int q, int n, int m, int p, int num_qs, int index);
	bool getOrCalculate_qary(int q, int n, int m, int p, int index, int num_qs, bool penaltyUsed, Lattice *l, FastVQA::PauliHamiltonian* h, DatasetRow* output_row, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions);
	bool getOrCalculate_qary_with_fixed_angles(	FastVQA::ExperimentBuffer* buffer,const double *angles, int angle_size, FastVQA::PauliHamiltonian* h, DatasetRow* row, FastVQA::QAOAOptions* qaoaOptions, FastVQA::Qaoa* qaoa_instance);

private:
	SQLite::Database* db;
	int loglevel = 0;
	void write(DatasetRow* row, DATABASE_TYPE database_type);

};

#endif /* SRC_IO_SQL_IO_HPP_ */
