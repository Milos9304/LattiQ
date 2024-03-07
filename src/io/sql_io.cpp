#include "sql_io.h"
#include <iostream>

#define s(X) std::to_string(X)

bool Database::getOrCalculate_qary(int q, int n, int m, int p, int index, int num_qs, bool penaltyUsed, Lattice *l, FastVQA::PauliHamiltonian* h, DatasetRow* output_row, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){

	assert(qaoaOptions->p == p);

	bool found = false;
	try{
		SQLite::Statement query(*db, "SELECT * FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND indexx="+s(index)+" AND num_qs="+s(num_qs)+" AND " + (!penaltyUsed ? "NOT " : "")+ "penaltyBool)");
		std::cerr<<"SELECT * FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND indexx="+s(index)+" AND num_qs="+s(num_qs)+" AND " + (!penaltyUsed ? "NOT " : "")+ "penaltyBool)";
		while (query.executeStep()){
			output_row->q 						= query.getColumn("q").getInt();
			output_row->n 						= query.getColumn("n").getInt();
			output_row->m 						= query.getColumn("m").getInt();
			output_row->p 						= query.getColumn("p").getInt();
			output_row->index 					= query.getColumn("indexx").getInt();
			output_row->num_qs 					= query.getColumn("num_qs").getInt();
			output_row->penalty					= query.getColumn("penalty").getInt();
			output_row->volume					= query.getColumn("volume").getInt();
			output_row->sv1Squared				= query.getColumn("sv1Squared").getInt();
			output_row->degeneracy				= query.getColumn("degeneracy").getInt();
			output_row->duration_s				= query.getColumn("duration_s").getInt();
			output_row->iters					= query.getColumn("iters").getInt();
			/*output_row->initAngles			= query.getColumn("initAngles");
			output_row->finalStateVectorMap		= query.getColumn("finalStateVectorMap");
			output_row->intermediateAngles		= query.getColumn("intermediateAngles");
			output_row->intermediateEnergies	= query.getColumn("intermediateEnergies");
			output_row->finalAngles				= query.getColumn("finalAngles");
			output_row->probSv1					= query.getColumn("probSv1");*/
			output_row->opt_res					= query.getColumn("opt_res").getText();
			/*output_row->comment					= query.getColumn("comment");*/
			found = true;
			break;
		}
	}catch (std::exception& e){
		throw_runtime_error(e.what());
	}

	if(found)
		return true;

	FastVQA::Qaoa qaoa_instance;
	FastVQA::ExperimentBuffer buffer;
	buffer.storeQuregPtr = false;

	time_t start_time = time(0);
	qaoa_instance.run_qaoa(&buffer, h, qaoaOptions);
	time_t end_time = time(0);
	int duration_s = difftime(end_time,start_time);

	output_row->type					= "qary";
	output_row->q 						= q;
	output_row->n 						= n;
	output_row->m 						= m;
	output_row->p 						= p;
	output_row->index 					= index;
	output_row->num_qs 					= num_qs;
	output_row->penalty					= mapOptions -> penalty;
	output_row->volume					= l->getVolume();
	output_row->sv1Squared				= l->getSquaredLengthOfFirstBasisVector();
	output_row->degeneracy				= 0;
	output_row->duration_s				= duration_s;
	output_row->iters					= buffer.intermediateEnergies.size();


	output_row->opt_res					= buffer.opt_message;
	this->write(output_row);

	return false;
}


bool Database::contains_qary(int q, int n, int m, int p, int index, int num_qs, bool penaltyUsed){

	bool found = false;

	try{
		SQLite::Statement query(*db, "SELECT q FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND indexx="+s(index)+" AND num_qs="+s(num_qs)+" AND " + (!penaltyUsed ? "NOT " : "")+ "penaltyBool)");
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
		SQLite::Statement query(*db, "INSERT or IGNORE INTO " + row->type + " VALUES ("+
				s(row->q)+", "+s(row->n)+", "+s(row->m)+", "+s(row->p)+", "+s(row->index)+
				", "+s(row->num_qs)+", "+s(penaltyBool)+", "+s(row->penalty)+", "+s(row->volume)+
				", "+s(row->sv1Squared)+", "+s(row->degeneracy)+", "+s(row->duration_s)+", "+s(row->iters)+", ?, ?, ?, ?, ?, ?, \""+row->opt_res+"\", ?)");
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

		db->exec("CREATE TABLE IF NOT EXISTS qary (q INTEGER, n INTEGER, m INTEGER, p INTEGER, indexx INTEGER, num_qs INTEGER, penaltyBool BOOL, penalty INTEGER, volume REAL, "
				"sv1Squared INTEGER, degeneracy INTEGER, duration_s INTEGER, iters INTEGER, initAngles BLOB, finalStateVectorMap BLOB, "
				"intermediateAngles BLOB, intermediateEnergies BLOB, finalAngles BLOB, probSv1 REAL, opt_res TEXT, comment TEXT, "
				"PRIMARY KEY (q, n, m, p, indexx, num_qs, penaltyBool))");

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
