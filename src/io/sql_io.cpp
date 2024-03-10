#include "sql_io.h"
#include <iostream>
/*#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>*/
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#define s(X) std::to_string(X)

inline unsigned int to_uint(char ch)
{
    // EDIT: multi-cast fix as per David Hammen's comment
    return static_cast<unsigned int>(static_cast<unsigned char>(ch));
}


//https://stackoverflow.com/questions/51230764/serialization-deserialization-of-a-vector-of-integers-in-c
/*template<typename POD>
std::ostream& serialize(std::ostream& os, std::vector<POD> const& v)
{
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only serialize POD types with this function");

    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(POD));
    return os;
}

template<typename POD>
std::istream& deserialize(std::istream& is, std::vector<POD>& v)
{
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
        "Can only deserialize POD types with this function");

    decltype(v.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(POD));
    return is;
}*/


inline std::stringstream query_unbind(std::string col_name, SQLite::Statement* query){

	std::stringstream ss;
	SQLite::Column colInitAngles		= query->getColumn(col_name.c_str());
	const void* const blob_init_angles	= colInitAngles.getBlob();
	const size_t size_init_angles   	= colInitAngles.getBytes();
	std::string initAngles_str((const char* const)blob_init_angles, size_init_angles+1);
	ss<<initAngles_str;
	//std::cerr<<size_init_angles;
	//std::cerr<<ss.str();
	//boost::archive::binary_iarchive ia(ss);
	return ss;
	//ia = new boost::archive::binary_iarchive(ss);
	/*std::cerr<<"h";
	ia >> output_row->initAngles;

	for(auto &a: output_row->initAngles)
		std::cerr<<a<<" ";*/

}

bool Database::getOrCalculate_qary(int q, int n, int m, int p, int index, int num_qs, bool penaltyUsed, Lattice *l, FastVQA::PauliHamiltonian* h, DatasetRow* output_row, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions){

	assert(qaoaOptions->p == p);

	bool found = false;
	try{
		SQLite::Statement query(*db, "SELECT * FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND indexx="+s(index)+" AND num_qs="+s(num_qs)+" AND " + (!penaltyUsed ? "NOT " : "")+ "penaltyBool)");
		//std::cerr<<"SELECT * FROM qary WHERE (q=" +s(q)+" AND n="+s(n)+" AND m="+s(m)+" AND p="+s(p)+" AND indexx="+s(index)+" AND num_qs="+s(num_qs)+" AND " + (!penaltyUsed ? "NOT " : "")+ "penaltyBool)";
		while (query.executeStep()){
			output_row->q 						= query.getColumn("q").getInt();
			output_row->n 						= query.getColumn("n").getInt();
			output_row->m 						= query.getColumn("m").getInt();
			output_row->p 						= query.getColumn("p").getInt();
			output_row->index 					= query.getColumn("indexx").getInt();
			output_row->num_qs 					= query.getColumn("num_qs").getInt();
			output_row->penalty					= query.getColumn("penalty").getInt();
			output_row->volume					= query.getColumn("volume").getDouble();
			output_row->sv1Squared				= query.getColumn("sv1Squared").getInt();
			output_row->degeneracy				= query.getColumn("degeneracy").getInt();
			output_row->duration_s				= query.getColumn("duration_s").getInt();
			output_row->iters					= query.getColumn("iters").getInt();

			//"initAngles"

			/*std::stringstream ss;
			SQLite::Column colInitAngles		= query.getColumn("initAngles");
			const void* const blob_init_angles	= colInitAngles.getBlob();
			const size_t size_init_angles   	= colInitAngles.getBytes();
			std::string initAngles_str((const char* const)blob_init_angles, size_init_angles+1);
			ss<<initAngles_str;
			std::cerr<<size_init_angles;
			//std::cerr<<ss.str();
			boost::archive::binary_iarchive ia(ss);*/
			std::stringstream ss = query_unbind("initAngles", &query);
			boost::archive::binary_iarchive ia(ss);
			ia >> output_row->initAngles;
			//ss.str(std::string());

			ss = query_unbind("finalStateVectorMap", &query);
			boost::archive::binary_iarchive ia2(ss);
			ia2 >> output_row->finalStateVectorMap;

			ss = query_unbind("intermediateAngles", &query);
			boost::archive::binary_iarchive ia3(ss);
			ia3 >> output_row->intermediateAngles;

			ss = query_unbind("intermediateEnergies", &query);
			boost::archive::binary_iarchive ia4(ss);
			ia4 >> output_row->intermediateEnergies;

			ss = query_unbind("finalAngles", &query);
			boost::archive::binary_iarchive ia5(ss);
			ia5 >> output_row->finalAngles;

			output_row->probSv1					= query.getColumn("probSv1").getDouble();
			output_row->opt_res					= query.getColumn("opt_res").getText();
			output_row->comment					= query.getColumn("comment").getText();
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
	buffer.storeQuregPtr = true;

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
	output_row->degeneracy				= buffer.final_solutions.size();
	output_row->duration_s				= duration_s;
	output_row->iters					= buffer.intermediateEnergies.size();

	std::vector<double> init_params;
	for(auto &pair : buffer.initial_params){
		init_params.push_back(std::get<1>(pair));
	}
	output_row->initAngles = init_params;

	std::shared_ptr<Qureg> stateVector = buffer.stateVector;
	FinalStateVectorMap finalStateVectorMap(stateVector->numAmpsTotal);

	//they are sorted in their energy
	FastVQA::RefEnergies eigenSpace = qaoaOptions->accelerator->getEigenspace();
	if(stateVector->numAmpsTotal != pow(2, output_row->num_qs))
		throw_runtime_error("stateVector->numAmpsTotal != pow(2, output_row->num_qs)");

	if(eigenSpace.size() != stateVector->numAmpsTotal)
		throw_runtime_error("eigenSpace.size() != stateVector->numAmpsTotal");

	//here we sort in their indices
	std::sort(eigenSpace.begin(), eigenSpace.end(), [&](FastVQA::RefEnergy i, FastVQA::RefEnergy j){return std::get<1>(i) < std::get<1>(j);});
	//for(auto e: eigenSpace){
	//	std::cerr<<std::get<1>(e) <<" "<<std::get<0>(e)<<"\n";
	//}

	for(long long int i = 0; i < stateVector->numAmpsTotal; ++i){
		finalStateVectorMap[i] = std::pair<qreal, double>(std::get<0>(eigenSpace[i]), stateVector->stateVec.real[i]*stateVector->stateVec.real[i]+stateVector->stateVec.imag[i]*stateVector->stateVec.imag[i]);
	}

	output_row->finalStateVectorMap 	= finalStateVectorMap;
	output_row->intermediateAngles 		= buffer.intermediateAngles;
	output_row->intermediateEnergies 	= buffer.intermediateEnergies;
	output_row->finalAngles				= buffer.finalParams;
	output_row->probSv1					= buffer.getTotalHitRate();
	output_row->opt_res					= buffer.opt_message;
	output_row->comment					= "";
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

template<typename POD>
inline void query_bind(int i, SQLite::Statement* query, std::vector<POD> const& v){

	std::stringstream ss;
	boost::archive::binary_oarchive oa(ss);
	oa << v;

	std::string blob_str = ss.str();
	query->bind(i, (void*)blob_str.c_str(), blob_str.size());
	//query->exec();
}

void Database::write(DatasetRow* row){

	if(row->type != "qary")
		throw_runtime_error("Not implemented");

	bool penaltyBool = row->penalty == 0 ? false : true;

	try{
		SQLite::Statement query(*db, "INSERT or IGNORE INTO " + row->type + " VALUES ("+
				s(row->q)+", "+s(row->n)+", "+s(row->m)+", "+s(row->p)+", "+s(row->index)+
				", "+s(row->num_qs)+", "+s(penaltyBool)+", "+s(row->penalty)+", "+s(row->volume)+
				", "+s(row->sv1Squared)+", "+s(row->degeneracy)+", "+s(row->duration_s)+", "+s(row->iters)+", ?, ?, ?, ?, ?," +s(row->probSv1)+", \""+row->opt_res+"\", \""+row->comment+"\")");

		query_bind(1, &query, row->initAngles);
		query_bind(2, &query, row->finalStateVectorMap);
		query_bind(3, &query, row->intermediateAngles);
		query_bind(4, &query, row->intermediateEnergies);
		query_bind(5, &query, row->finalAngles);
		query.exec();

		/*std::stringstream ss;
		boost::archive::binary_oarchive oa(ss);
		oa << row->initAngles;

		std::string blob_str = ss.str();
		query.bind(1, (void*)blob_str.c_str(), blob_str.size());
		query.exec();

		boost::archive::binary_iarchive ia(ss);
		ia >> row->initAngles;
		for(auto &a:row->initAngles)
			std::cerr<<a<<" ";*/


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
