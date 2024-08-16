/*
 * experiment_runner.h
 *
 *  Created on: Nov 20, 2023
 *      Author: Milos Prokop
 */

#ifndef SRC_EXPERIMENTRUNNER_H_
#define SRC_EXPERIMENTRUNNER_H_

#include "FastVQA/fastVQA.h"
#include "lattice/lattice.h"
#include "io/sql_io.h"
#include <string>
#include "paper_experiments/g1.h"
#include "paper_experiments/g2.h"


class ExperimentSetup{
public:

	std::string experiment_type;
	FastVQA::PauliHamiltonian* hamiltonian;
	FastVQA::QAOAOptions* qaoaOptions;
	qreal minimum_energy;

	//case specific
	int num_rand_params;
};

class AngleExperimentBase{
public:

	FastVQA::QAOAOptions* qaoaOptions;

	int q = 97;

	int m_start; /*9*/;
	int m_end;// = 20;//10;///20; //10;

	int max_num_instances = 100;//3000;

	struct Cost{
			double mean;
			double stdev;
			double mean_zero;
			double mean_num_of_sols;

			Cost(){}

			Cost(double mean, double stdev, double mean_zero, double mean_num_of_sols){
				this->mean = mean;
				this->stdev = stdev;
				this->mean_zero = mean_zero;
				this->mean_num_of_sols=mean_num_of_sols;
			}
		};


	~AngleExperimentBase(){
		this->logfile.close();
		this->angleAnalysisLog.close();
	}
protected:

	Database* database;
	int loglevel = 1;
	MapOptions* mapOptions;

	struct Instance{
			FastVQA::PauliHamiltonian h;
			FastVQA::RefEnergies zero_solutions;
			FastVQA::RefEnergies sv_solutions;
			//FastVQA::RefEnergies eigenspace; //for debug
			qreal min_energy;
			qreal random_guess;

			double volume;
			int sv1Squared;

			int q,m,n;

			FastVQA::Accelerator::DiagonalOpDuplicate diagOpDuplicate;
	};

	std::ofstream logfile, angleAnalysisLog;
	std::mt19937 gen19937;
	bool seeded=false;

	std::vector<Instance> _generate_dataset(int n, int m, bool penalise=false);
	int seed=0;
	bool use_database_to_load_dataset;

	std::pair<double, double> try_many_starts(std::string meta_data, Instance* instance, FastVQA::Qaoa* qaoa_instance, int seed);
	Cost _cost_fn(std::vector<Instance>*, const double *angles, std::string meta_data, bool use_database=false, int seed=0);
};

class AqcPqcExperiment : AngleExperimentBase{

public:
	AqcPqcExperiment(int loglevel, int m_start, int m_end, FastVQA::QAOAOptions* qaoaOptions, MapOptions* mapOptions, Database* database, int seed, bool use_database_to_load_dataset){

		this->loglevel = loglevel;
		this->qaoaOptions = qaoaOptions;
		this->mapOptions = mapOptions;
		this->database = database;
		this->m_start = m_start;
		this->m_end = m_end;
		this->seed = seed;

		this->use_database_to_load_dataset = use_database_to_load_dataset;

	}

	void run();

};

class AngleResultsExperiment : AngleExperimentBase{
private:
	struct optAngle{

		//2**(c+alpha*n)

		int p;
		std::vector<double> cm_angles;
		double cm_c;
		double cm_n;
		std::string cm_meta_data = "nan";

		std::vector<double> qaoa_angles;
		double qaoa_c;
		double qaoa_n;
		std::string qaoa_meta_data = "nan";

		optAngle(int p, const std::vector<double> cm_angles, double cm_c, double cm_n, std::string cm_meta_data, const std::vector<double> qaoa_angles, double qaoa_c, double qaoa_n, std::string qaoa_meta_data){
			this->p = p;
			this->cm_angles = cm_angles;
			this->cm_c=cm_c;
			this->cm_n=cm_n;
			this->cm_meta_data = cm_meta_data;

			this->qaoa_angles = qaoa_angles;
			this->qaoa_c = qaoa_c;
			this->qaoa_n = qaoa_n;
			this->qaoa_meta_data = qaoa_meta_data;
		}

		optAngle(int p){this->p = p;};

	};
public:

	const std::vector<optAngle> optAngles{
		//0 to 4 just to make index of arrays match p
		optAngle(0),
		optAngle(1),
		optAngle(2),
		optAngle(3),
		optAngle(4),
		/*optAngle(
				5,
				//CM
				{0.593466619587856, 2.16406666372329, 2.91680843513582,
				2.90777554501237, 0.776326446551238, -0.709028070636096,
				-1.25692251514367, -2.78536849493982, -1.43170870653588,
				-0.125794565420379},
				-0.353648801738946,
				-0.863370504453672,
				"CM: optimized by diff, MAXEVAL_REACHED, num_iters: 1000",

				//QAOA
				{2.16757160644228, 2.18531279230608, 2.41412345494909,
				2.1843093644902, 0.77104656059519, -0.74637995068122,
				-1.27432665185911, -2.78883772481428, -1.43167967368264,
				-0.169083861050817},
				-0.660683735633075,
				-0.681604944157161,
				"QAOA: optimized by diff, FTOL_REACHED, num_iters: 951"
		),*/
		optAngle(
		                                5,
		                                //CM
		                                {2.16268938746618, 2.90226322313521,
		                                2.23487712130271, 2.19073802435072, 0.860098887796126,
		                                -0.715472691321374, -1.28469680884718, -2.78874181966366,
		                                -1.41465552670908, -0.147919049349518},
		                                0.00959032096969687,
		                                -1.01586546126246,
		                                "CM: optimized by diff, ROUNDOFF_LIMITED, num_iters: 2282",
		                                //QAOA
		                                {0.538465841789183, 2.91914145726405,
		                                2.24183247223159, 3.03404533172038, 0.53282011478364,
		                                -0.809791281103069, -1.55591470883738, -2.84944171187325,
		                                -1.72760678216298, -0.529438854036467},
		                                -0.8696403423757,
		                                -0.744237716560413,
		                                "QAOA: optimized by diff, ROUNDOFF_LIMITED, num_iters: 2180"
		                ),
		optAngle(
				6,
				//CM
				{0.524413709222256, 2.50291231189887,
				2.93878239982546, 2.93574145025635, 0.7554670199668,
				-0.830204156706467, -1.68347094505389, -2.81104037624286,
				-1.55329305113912, -0.244212853734048, 1.95196042386488,
				-0.228645879681772},
				-0.408258795390403,
				-0.86303052301192,
				"CM: optimized by diff, MAXEVAL_REACHED, num_iters: 1000",
				//QAOA
				{2.14754824256986, 2.16836131486609,
				2.12837988684515, 2.18142325302405, 0.785432792698188,
				-0.726508211167727, -1.27266845242185, -2.78530660573082,
				-1.43028367472466, -0.136623665542789, 2.84748587759566,
				-0.126958255879974},
				-0.483670508627367,
				-0.677773130871272,
				"QAOA: optimized by diff, MAXEVAL_REACHED, num_iters: 1000"
		),
		optAngle(7),
		optAngle(8),
		optAngle(9),
		optAngle(
10,
//CM
{2.15415622223479, 2.16308546716135,
2.24903864898094, 2.18184700594222, 0.776373637641152,
0.844345195271818, -1.27212762305763, -2.78525451219963,
-1.42844260038797, -0.140334204675264, 1.96141396258403,
1.44498920786814, -0.673653079195817, 2.11164518069261,
-1.02167004727948, 0.930991353047629, -0.827862873397143,
2.87239055072913, -2.2597427073567, 2.32532677849527},
-0.98589263139679,
-0.860452253113452,
"CM: optimized by diff, MAXEVAL_REACHED, num_iters: 1000",
//QAOA
{0.580242709724772, 2.89704003868955,
2.92041964300038, 2.18269660071074, 2.35155244972844,
-0.725050623624204, -1.2706730138929, -2.51146888788683,
-1.42814647314992, -0.139997673451478, 1.96136734063202,
-0.129099841689273, -0.67369134489236, 2.88866390099968,
0.548914407435724, 0.927984830221724, -0.829945872645291,
3.07282948236154, -1.59870132989752, 2.32523201768276},
-0.543107907573785,
-0.537928642679504,
"QAOA: optimized by diff, MAXEVAL_REACHED, num_iters: 1000"
)

		//optAngle()
	};


	//const std::vector<double> angles{0.4,0.48,5.56,0.28}; //work interesting

	//WORKS GREAT
	//const std::vector<double> angles{2.2287, 2.13607 ,2.23057 ,2.17114, 0.736168, -0.775112, -1.17933, -2.84297, -1.32729 ,-0.37258 ,1.90813 ,-0.247229};
	//const std::vector<double> angles{3.08091, 3.11479, 2.24896, 2.26394, 0.805914, -0.754133, -1.2905, -2.78654, -1.4643, -0.166783, 1.9449, -0.141213};

	//Found by AlphaMinim
	//new_way
	//const std::vector<double> angles{1.84887, 2.15068, 2.24906, 2.18152, 0.776378, -0.71996, -1.26273, -2.78525, -1.42938 ,-0.140317, 1.96114, -0.126679};


	//old way
	//const std::vector<double> angles{2.15445, 2.16252, 2.24904 ,2.18185, 0.776374 ,-0.726451, -1.27213 ,-2.78525 ,-1.42844, -0.140334, 1.96141, -0.125807};

	//Found by AlphaMinim qary_4_20 q=97
	//const std::vector<double> angles{3.09402 ,2.61273, 2.27356, 2.1883 ,0.784678, -0.830595, -1.23473, -2.78089 ,-1.40327, -0.133618, 1.97257 ,-0.103609 };

	//here I fixed p=6
	//const std::vector<double> angles{0.620733, 2.18047, 2.25891 ,2.22203};//, 2.34926, -0.347949, -1.28916 ,-2.78512 ,-1.43107, -0.195526 ,1.96525 ,-0.0549031};
	//const std::vector<double> angles{2.15416 ,2.16309, 2.24904, 2.90166 ,0.776374, -0.726451 ,-1.27213 ,-2.78525, -1.42844 ,-0.140334 ,1.96141 ,-0.125807 };
	//const std::vector<double> angles{2.13774 ,2.17964, 2.23325, 2.15159 ,0.792933, -0.642643 ,-1.2983 ,-2.79203, -1.45799, -0.0792081, 1.96191, -0.086081 };
	//const std::vector<double> angles{0.621138, 3.14106 ,2.47842, 2.17828 ,1.08465 ,0.850709 ,-1.273, -2.78661, 0.145223, -0.140685 ,1.95694 ,-0.127101 };
	//const std::vector<double> angles{0.572054, 2.18904, 2.254 ,2.20524 ,0.788806, 0.0820688 ,0.270073 ,-2.78492, -1.42865, -0.0629617, 1.95564 ,-0.299644};
	//const std::vector<double> angles{0.469226, 2.07037, 2.13885, 2.18051, 0.729963, 1.50839, -1.33214, -2.4887, -1.76529, -1.06572, 1.8024, -1.03263};
	//const std::vector<double> angles{0.469225527665256, 2.0703697424767, 2.13885358583399, 2.1805121131481, 0.729963277183836, 1.50838535672113, -1.33213908612272, -2.48869765417828, -1.76529454649408, -1.0657212157689, 1.80239788470083, -1.03262650444359};
	//const std::vector<double> angles{0.572053699085232, 2.18904039343271, 2.25400248594688, 2.20523665388413, 0.788805944209345, 0.0820687583023222, 0.270073448207576, -2.78491502689493, -1.42864987471102, -0.0629616583429195, 1.95564143867776, -0.299643861461467};
	//const std::vector<double> angles{2.29564492153543, 2.23457510881661, 2.23840220369636, 2.90886358876134, 0.669492961403253, 0.76354956132728, 0.30587570182067, -2.7997137981176, 0.0402611083789431, -0.061731654401902, 2.02276087531338, -0.32073364358472};
	//const std::vector<double> angles{0.586233360039837, 2.09152545293368, 2.24799347012238, 2.90045252659239, 0.775105822307031, -0.726154588790109, -1.27186023022667, -2.51829028197983, -1.43104071313081, -0.141133384711111, 1.96088220572089, -0.103904749531919};






	//const std::vector<double> angles_optqaoa{0.668373713325827, 0.395489639960678, 1.01189488174712, 0.0716902196238647, 2.36312900989193, 0.3910579937851, 2.37481345843973, 0.364219890147565, 0.809596592930284, 0.356048847055572, 0.799885077937716, 0.345458300846393};
	//const std::vector<double> angles_optqaoa{2.15117201255203, 2.16366715121114, 2.24799259022539, 2.18206089715351, 0.776037279297707, -0.723740787022484, 0.299514557193188, -2.78580448703939, -1.43471061039816, -0.14667707028971, 1.95826294997079, -0.127924272009467};

	//p=7
	//const std::vector<double> angles_optqaoa{0.712592176356089, 1.13183206034174, -0.267304592085007, 2.17025670667795, 3.04413825248703, -0.172930821841033, 1.69485996479564, -0.133274769031307, 2.48402745282503, 0.899294329780828, 2.23541397782318, 0.91098670585426, -1.37356405013328, 0.870114051197184};


	//const std::vector<double> angles_cmqaoa{0.668373713325827, 0.395489639960678, 1.01189488174712, 0.0716902196238647, 2.36312900989193, 0.3910579937851, 2.37481345843973, 0.364219890147565, 0.809596592930284, 0.356048847055572, 0.799885077937716, 0.345458300846393};
	//const std::vector<double> angles_cmqaoa{0.0622837466897899, 1.60530356338513, 2.37145183079299, 0.329686755158126, 0.544763061077257, 0.287636536564108, 1.49021949821313, 1.994686824452, 0.685993589190685, 0.295995789511559, 0.716793729178239, 0.278672884908343};


	//const std::vector<double> angles_cmqaoa{2.78475622527346, 2.42485902183718, 2.89407293479389, 2.17258774439195, 1.1832036873172, 0.514674880973178, 0.315284023148172, -2.51582479341127, 0.153056596878831, -0.304275696701542, 1.95172146690934, -0.261807499106509};
	//const std::vector<double> angles_cmqaoa{0.621138037244616, 3.14106105652719, 2.47841824458321, 2.17827658373638, 1.08465084556728, 0.850709306987419, -1.27300211754974, -2.78660735037908, 0.145222608167673, -0.140684892497685, 1.95693664215292, -0.127101424814812};
	//const std::vector<double> angles_cmqaoa{3.03073316671942, 1.70005902905338, 0.984953130900807, 2.96912368985253, 3.04539918603008, 0.807491753884046, 1.45250406962383, 0.149390050574869, 0.588001376393152, 0.932185268727973, 2.26467396419548, 0.968928395222161};


	//p=7
	//const std::vector<double> angles_cmqaoa{1.0398938452859, -0.301279449107924, 0.0132584027202693, 2.14324194192504, 3.11981478542475, 0.248099885485753, 1.62539346673195, 0.267705037628689, 0.811361593108557, 0.977545822278086, 2.91355687724525, 1.0111198658054, -1.36780171714015, 0.977907606632986};

	//test inv diff
	//const std::vector<double> angles_cmqaoa{2.02064224262685, 2.01564522328784, 2.99845395184209, 2.26047337282221, 2.9279783883639, -0.458233018857491, 0.306887385399807, -2.78600959913084, -1.43778920468067, -0.1822121881121, 1.8688398837591, -0.246610441872544, -0.685379479630748, 2.11076624833216};
	//const std::vector<double> angles_cmqaoa{0.60153895632161, 2.29118105158059, 2.92259102919125, 2.89587839310268, 0.723863020973766, -0.707890807050217, -1.25967440672003, -2.7906023907625, -1.41419435751768, -0.187750638872502, 1.95362922597325, -0.0547654081190543};


	//p=5
	/*
	 * num_iters: 1000
const std::vector<double> angles_cmqaoa{0.593466619587856, 2.16406666372329, 2.91680843513582, 2.90777554501237, 0.776326446551238, -0.709028070636096, -1.25692251514367, -2.78536849493982, -1.43170870653588, -0.125794565420379};
MAXEVAL_REACHED
2^-0.353648801738946+n*-0.863370504453672
	 *
	 */

	//p=6
	//best so far, optimized by 1/sum(f-1)
	//const std::vector<double> angles_cmqaoa{0.524413709222256, 2.50291231189887, 2.93878239982546, 2.93574145025635, 0.7554670199668, -0.830204156706467, -1.68347094505389, -2.81104037624286, -1.55329305113912, -0.244212853734048, 1.95196042386488, -0.228645879681772};

	//2nd best, optimized by c in 2^c+xn
	//const std::vector<double> angles_cmqaoa{0.624252283373159, 2.88050061815466, 2.19376394838984, 2.35633892704388, 0.642474779825364, -0.572273387754715, -1.34097902496242, -2.79962116072776, -0.0398835051960656, -0.136707786034164, 2.09058913512657, -0.245809192212849};


	//const std::vector<double> angles_cmqaoa{2.02064224028929, 2.01564522338183, 2.9984539516002, 2.26047337314089, 2.92797838753706, -0.458233019499983, 0.306887386203445, -2.78600959919879, -1.43778919841884, -0.18221218857827, 1.86883987711935, -0.246610442523591, -0.685379478880007, 2.11076624811066};
	//const std::vector<double> angles_cmqaoa{2.66555823675693, -0.344026235752078, -0.384395298158486, 2.47724904479005, 3.06134781841481, 0.80077194180158, 2.79134167551185, 0.161486735808091, 0.73042857913784, 0.563941079933157, 2.16704279224347, 0.85274651397496, -1.36581140094452, 0.94669047269831};


	//optimised by overlap on dim15
	//const std::vector<double> angles{-0.918725, 3.06621, 2.24903, 2.18185, 0.776357, -0.726439, -1.27211 ,-2.78525, -1.42844, -0.140346 ,1.96142, -0.1258 };
	/*
 *  q=97
   Averages:
 n \ m                   4                   5                   6                   7                   8                   9                  10
  1              0.0641671           0.0333538           0.0157872          0.00804493          0.00401294          0.00211669          0.00100611
  2              0.0684507           0.0324452           0.0162824          0.00841373          0.00680977          0.00600239          0.00476806
  3               0.106731           0.0515169           0.0207424           0.0126368            0.005107          0.00519504           0.0044405
  4                      x            0.102938           0.0451122           0.0091706          0.00585934          0.00539514          0.00442476
  5                      x                   x           0.0859569           0.0323786          0.00566087          0.00542296          0.00443021
  6                      x                   x                   x           0.0508255           0.0153818          0.00288475          0.00273762
  7                      x                   x                   x                   x           0.0364387          0.00482137          0.00144237
  8                      x                   x                   x                   x                   x           0.0213782          0.00144802
  9                      x                   x                   x                   x                   x                   x           0.0122817
 *
 */


	//const std::vector<double> angles{0.583336, 2.16313 ,2.24903 ,2.18185};

	AngleResultsExperiment(int loglevel, int m_start, int m_end, FastVQA::QAOAOptions*, MapOptions*, Database*, int seed, bool use_database_to_load_dataset);

	void run_qaoa_with_optimizer();

	void run();

private:
	//std::vector<Instance> _generate_dataset(int n, int m, bool penalise=false);
	int seed=0;
	bool use_database_to_load_dataset;

};

class AlphaMinimizationExperiment{

	public:
		int loglevel = 1;
		int p = 1;

	AlphaMinimizationExperiment(int loglevel, FastVQA::QAOAOptions*, MapOptions*, Database*);

	void run(bool use_database_to_load_dataset);

	private:

	struct AlphaMinimizationExperimentInstance{
		FastVQA::PauliHamiltonian h;
		FastVQA::RefEnergies zero_solutions;
		FastVQA::RefEnergies sv_solutions;

		double volume;

		int num_qubits_per_dim;
		int q,m;

		FastVQA::Accelerator::DiagonalOpDuplicate diagOpDuplicate;
	};


	inline double strategy_alpha_c(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset, std::vector<double> angles, std::string meta_data);
	inline double strategy_inv_diff(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset, std::vector<double> angles, std::string meta_data);

	inline double strategy_random_alpha_c(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset, std::vector<double> angles, std::string meta_data);
	inline double strategy_random_inv_diff(std::vector<std::vector<AlphaMinimizationExperimentInstance>> train_dataset, std::vector<double> angles, std::string meta_data);

	double _cost_fn(std::vector<AlphaMinimizationExperimentInstance>, const double *angles, std::string meta_data, int probability100=100);
	FastVQA::QAOAOptions* qaoaOptions;
	MapOptions* mapOptions;
	Database* database;
};

void experiment_runner(ExperimentSetup*, std::string experiment_name, int loglevel);

#endif /* SRC_EXPERIMENTRUNNER_H_ */
