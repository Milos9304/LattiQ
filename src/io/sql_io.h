/*
 * saveProgress.hpp
 *
 *  Created on: Mar 06, 2024
 *      Author: Milos Prokop
 */

#ifndef SRC_IO_SQL_IO_HPP_
#define SRC_IO_SQL_IO_HPP_

#include "SQLiteCpp/SQLiteCpp.h"
#include "logger.h"

typedef std::vector<std::tuple<long, double>> FinalStateVectorMap;

class Database{
public:

	struct DatasetRow{
		std::string type;
		int q;
		int n;
		int m;
		int p;
		int penalty;
		double volume;
		int sv1Squared;
		int degeneracy;
		int num_qs;
		std::vector<double> initAngles;
		FinalStateVectorMap finalStateVectorMap;
		std::vector<std::vector<double>> intermediateAngles;
		std::vector<double> intermediateEnergies;
		std::vector<double> finalAngles;
		double probSv1;
		std::string comment;
	};

	static void print_sqlite_info(std::string filename);

	Database(std::string filename, int loglevel=0);
	void write(DatasetRow* row);
	bool contains_qary(int q, int n, int m, int p, bool penaltyUsed);

private:
	SQLite::Database* db;
	int loglevel = 0;
};

#endif /* SRC_IO_SQL_IO_HPP_ */
