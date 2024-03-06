#include "sql_io.h"
#include <iostream>

#define s(X) std::to_string(X)

bool Database::contains_qary(int q, int n, int m, int p, bool penaltyUsed){

	bool found = false;

	try{
		//const std::string value = db->execAndGet("SELECT EXISTS ( SELECT q FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND penaltyBool="+s(penaltyUsed)+") LIMIT 1)");
		//const std::string value = db->execAndGet("SELECT rowid, q, n FROM qary WHERE (q=" +s(q)+/*" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND penaltyBool="+s(penaltyUsed)+*/")");
		SQLite::Statement query(*db, "SELECT q FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND " + (!penaltyUsed ? "NOT " : "")+ "penaltyBool)");
		while (query.executeStep()){
			found = true;
			break;
		}

	}catch (std::exception& e){
		throw_runtime_error(e.what());
	}

	return found;
}


void Database::write(DatasetRow* row){

	if(row->type != "qary")
		throw_runtime_error("Not implemented");

	bool penaltyBool = row->penalty == 0 ? false : true;

	try{
		SQLite::Statement query(*db, "INSERT or IGNORE INTO " + row->type + " VALUES ("+s(row->q)+", "+s(row->n)+", "+s(row->m)+", "+s(row->p)+", "+s(penaltyBool)+", "+s(row->penalty)+", ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
		//query.bind(1, blob, size);
		query.exec();
	}catch (std::exception& e){
		throw_runtime_error(e.what());
	}

}

Database::Database(std::string filename, int loglevel){

	this->loglevel = loglevel;

	try{// Open a database file in create/write mode
		db = new SQLite::Database(filename, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
		logd("SQLite database file '" + db->getFilename() + "' opened successfully\n", loglevel);

		db->exec("CREATE TABLE IF NOT EXISTS qary (q INTEGER, n INTEGER, m INTEGER, p INTEGER, penaltyBool BOOL, penalty INTEGER, volume REAL, "
				"sv1Squared INTEGER, degeneracy INTEGER, num_qs INTEGER, initAngles BLOB, finalStateVectorMap BLOB, "
				"intermediateAngles BLOB, intermediateEnergies BLOB, finalAngles BLOB, probSv1 REAL, comment TEXT, "
				"PRIMARY KEY (q, n, m, p, penaltyBool))");

	}catch (std::exception& e){
		throw_runtime_error(e.what());
	}
}

void Database::print_sqlite_info(std::string filename){
	std::cout << "SQlite3 version " << SQLite::VERSION << " (" << SQLite::getLibVersion() << ")" << std::endl;
	std::cout << "SQliteC++ version " << SQLITECPP_VERSION << std::endl;

	try
	    {
	       const SQLite::Header header = SQLite::Database::getHeaderInfo(filename);

	       // Print values for all header fields
	       // Official documentation for fields can be found here: https://www.sqlite.org/fileformat.html#the_database_header
	        std::cout << "Magic header string: " << header.headerStr << std::endl;
	        std::cout << "Page size bytes: " << header.pageSizeBytes << std::endl;
	        std::cout << "File format write version: " << (int)header.fileFormatWriteVersion << std::endl;
	        std::cout << "File format read version: " << (int)header.fileFormatReadVersion << std::endl;
	        std::cout << "Reserved space bytes: " << (int)header.reservedSpaceBytes << std::endl;
	        std::cout << "Max embedded payload fraction " << (int)header.maxEmbeddedPayloadFrac << std::endl;
	        std::cout << "Min embedded payload fraction: " << (int)header.minEmbeddedPayloadFrac << std::endl;
	        std::cout << "Leaf payload fraction: " << (int)header.leafPayloadFrac << std::endl;
	        std::cout << "File change counter: " << header.fileChangeCounter << std::endl;
	        std::cout << "Database size pages: " << header.databaseSizePages << std::endl;
	        std::cout << "First freelist trunk page: " << header.firstFreelistTrunkPage << std::endl;
	        std::cout << "Total freelist trunk pages: " << header.totalFreelistPages << std::endl;
	        std::cout << "Schema cookie: " << header.schemaCookie << std::endl;
	        std::cout << "Schema format number: " << header.schemaFormatNumber << std::endl;
	        std::cout << "Default page cache size bytes: " << header.defaultPageCacheSizeBytes << std::endl;
	        std::cout << "Largest B tree page number: " << header.largestBTreePageNumber << std::endl;
	        std::cout << "Database text encoding: " << header.databaseTextEncoding << std::endl;
	        std::cout << "User version: " << header.userVersion << std::endl;
	        std::cout << "Incremental vaccum mode: " << header.incrementalVaccumMode << std::endl;
	        std::cout << "Application ID: " << header.applicationId << std::endl;
	        std::cout << "Version valid for: " << header.versionValidFor << std::endl;
	        std::cout << "SQLite version: " << header.sqliteVersion << std::endl;
	    }
	    catch (std::exception& e)
	    {
	        std::cout << "SQLite exception: " << e.what() << std::endl;
	        //return EXIT_FAILURE; // unexpected error : exit the example program
	    }
}
