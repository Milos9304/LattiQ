#include "experiment_runner.h"

const std::vector<AngleResultsExperiment::optAngle> AngleResultsExperiment::optAngles = std::vector<optAngle>{	//0 to 4 just to make index of arrays match p
	optAngle(0),
		optAngle(
			 1,
			 //CM
			 {-1.9818560023213, 2.49080900690849},
			 -1.1606568894798,
			 -0.760314510728059,
			 -0.91360881688577,
			 "CM: optimized by strategy_alpha_c, FTOL_REACHED, num_iters: 13013",
			 //QAOA
			{-2.93822689293179, 0.649554060019249},
                        0.205486995094955,
                        -0.728346781110618,
                        -0.701206989305623,
                        "QAOA: optimized by strategy_alpha_trivial, FTOL_REACHED, num_iters: 17908"
			),
	optAngle(
				2,
				//CM
				{-2.62904184702749, 0.934365389567761, 
				-0.891064087361401, 3.01942666836609},
				-0.566291448854019,
				-0.790916721882354,
				-0.865709932108357,
				"CM: optimized by strategy_alpha_trivial, FTOL_REACHED, num_iters: 75954",
				//QAOA
				/*{3.0790487715258, 2.44047352891003, 
				-2.46598803170775, 0.282795837640577},
				0.425185806166174,
				-0.74683475147189,
				-0.69067813556315,
				"QAOA: optimized by strategy_alpha_trivial, FTOL_REACHED, num_iters: 74343"*/
				{0.422893724621816, 2.89751136368467,
				-0.953980008922711, -2.45539483558622},
				0.173598095358752,
				-0.681497316357708,
				-0.658569266027306,
				"QAOA: optimized by strategy_inv_diff, FTOL_REACHED, num_iters: 81944"
		
		),

		/*this comes from angle append
		optAngle(
		                                2,
		                                //CM
		                                {-1.98185600236407, 2.49080900692214,
		                                -3.67596645898487e-12, -3.67879179490481e-10},
		                                -1.16065595497184,
		                                -0.760314607237011,
		                                -0.913608789969142,
		                                "CM: optimized by strategy_alpha_trivial, FTOL_REACHED, num_iters: 184",
		                                //QAOA
		                                {-2.9382268929318, 0.649554060464992,
		                                1.33536894092762e-11, 1.13900167317228e-09},
		                                0.205487133476458,
		                                -0.728346799063334,
		                                -0.701206988981537,
		                                "QAOA: optimized by strategy_alpha_trivial, FTOL_REACHED, num_iters: 178"
		                )*/
	optAngle(
				3,
				//CM
				/*{1.37405229519096, 2.97982492443362, 
				0.120331074004508, 3.04451969902882, 0.964965837096694, 
				-0.902157465495305},*/ //orig
				
				//{2.81244106643164, -0.851164688746971, -2.73660062175586, -0.814277074521615, 3.04786397656108, -2.54704642885341}, //i=6/100 0.95
				//{0.947062504145166, -2.38894340512297, 1.17210113972164, -0.0329207107072488, 1.22177736147284, -2.48760461968259}, //23/100 0.93
				//{-2.52283141032084, -2.46881701870303, -1.22784214095808, 3.00986295003153, -1.14687745903958, 0.578282979949387}, //24/100 0.905
				//{-2.73056439888557, 0.765343081755736, -0.136395377303467, 0.285642239137269, 1.38937578794491, -3.08357159721367}, //i=31/100 0.92
				//{-2.34395385124363, 2.34219246844936, 2.53796527402623, -2.80946629083266, 3.14159031244396, -2.73528763842277}, //37/100 0.93				
				//{-0.452783874059261, -0.0580396747446087, 3.07514966700645, 0.623947132395434, 2.83254506046648, -0.382216147302585}, //51/100 0.94
				//{1.91502254594642, 0.863659217150407, 3.01057590695428, -0.640787552203429, 2.82087260035087, 2.7719473655417}, //59/100 0.93
				//{2.78431854353105, 0.65208521308778, 2.99201248328994, 0.38886841499274, 3.14139744029428, -0.424551887562068}, //60/100 0.94
				//{0.122554609026554, -2.8571310975287, 1.30884475955703, 0.678293734557697, 2.23375586290952, 0.427026918238707}, //76/100 0.891
				//{-3.02497559112185, 0.935680072552535, -2.9347414147973, 2.96320023709121, 1.69689158441636, 3.09239311983661}, //77/100 0.93
				//{3.13010514150112, -0.600292040677666, 1.87252074479868, 3.03341052432501, -0.275519866443202, 0.225622942121158}, //91/100 0.94

				
				//trying very new
	//{0.122558853838647, -2.85713076166752, 1.30884241780053, 0.678292668594285, 2.23375551404907,	0.42702904320225}, //most recent 0.91
//{1.69335094142114, -2.8571310970736, 2.87964108632242, 0.678293736306142, 2.23375586037836, 0.427026917954336}, //0.897328
	//{0.122554608137382, -2.85713109722187, 1.30884477335296, 0.678293744253285, 2.23375586214312, 0.427026912807376}, 0.917
	// {1.69335092483654, -2.85713111924581, 2.87964110263768, 0.678293760190778, 2.23375585424438, 0.427027007200606}, //0.900
	{0,0,0,0,0,0},
				-0.110263471008016,
				-0.864217587631346,
				-0.878780687575802,
				"CM: optimized by strategy_alpha_trivial, MAXEVAL_REACHED, num_iters: 106788",
				//QAOA
				{2.15716733995899, 2.15884232758313, 
				2.25424333068428, 2.17236853814546, 0.763099959521827, 
				-0.734511748171839},
				-0.0644310433896435,
				-0.667310341316709,
				-0.675820101764399,
				"QAOA: optimized by strategy_inv_diff, FTOL_REACHED, num_iters: 128962"
		),
	optAngle(4),
	optAngle( //newest today 80 80
	                                5,
	                                //CM
	                                {2.12192095086874, 2.16362508301005,
	                                2.92040151806505, 2.90298685961388, 0.77293889858348,
	                                -0.669307433197668, -1.2739339249801, -2.786095148819,
	                                -1.43004934263276, -0.144646908663418},
	                                -0.154599807724666,
	                                -0.918378730087432,
									0,
	                                "CM: optimized by diff, ROUNDOFF_LIMITED, num_iters: 2314",
	                                //QAOA
	                                {2.13567836738189, 2.1596482051349,
	                                2.23987352699944, 2.18642413418565, 1.55343911594894,
	                                -0.751462804672511, -1.27298967244183, -2.76699346253223,
	                                0.142592564243054, -0.144818073999694},
	                                -0.314404252698797,
	                                -0.699699552159603,
									0,
	                                "QAOA: optimized by diff, ROUNDOFF_LIMITED, num_iters: 2329"
	                ),
	/*optAngle(
			5,
			//CM
			{2.15627340276711, 2.16452766247427,
2.91983504050558, 2.90324169432516, 0.781637549571124,
-0.719204499212141, -1.26956910380845, -2.78478709100459,
-1.43139899090179, -0.143913264618946},
-0.0418465244415251,
-0.950104738875151,
"CM: optimized by diff, ROUNDOFF_LIMITED, num_iters: 2221",
			//QAOA
			{2.16757160644228, 2.18531279230608, 2.41412345494909,
			2.1843093644902, 0.77104656059519, -0.74637995068122,
			-1.27432665185911, -2.78883772481428, -1.43167967368264,
			-0.169083861050817},
			-0.660683735633075,
			-0.681604944157161,
			"QAOA: optimized by diff, FTOL_REACHED, num_iters: 951"
	),*/

	/*optAngle(
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
	),*/
					optAngle(
					                                6,
					                                //CM
					                                {2.12192096350434, 2.16362504979583,
					                                2.92040149454857, 2.90298686465215, 0.772938921189123,
					                                -0.669307879245874, -1.27393391511215, -2.78609514932928,
					                                -1.59459524471504, -0.14464694874021, -0.0808453637579762,
					                                9.33605255133508e-09},
					                                -0.231226196869362,
					                                -0.900253110229867,
													0,
					                                "CM: optimized by random alpha, MAXEVAL_REACHED, num_iters: 4000",
					                                //QAOA
													{2.15117201258113, 2.16366715116623,
													2.24799259021588, 2.18206089726526, 0.776037279415087,
													-0.723740787036124, 0.299514557157671, -2.78580448707204,
													-1.43471061122679, -0.146677070147529, 1.9582629509832,
													-0.127924272028654},
													-2.06315843846787,
													-0.535304471708114,
													0,
													"QAOA: optimized by diff, MAXEVAL_REACHED, num_iters: 2000"

					                ),/*
	optAngle(
		6,
		//CM
		{2.12469761917628, 2.89456926114369,
		2.91688607811889, 2.90258593492434, 0.781854702527559,
		-0.578190200500547, -1.28042499897594, -2.78548182008663,
		-1.43255788277925, -0.146776692484744, -0.0020591767175171,
		0.00117419205764549},
		0.0894174017630761,
		-0.920054573815141,
		"CM: optimized by diff, ROUNDOFF_LIMITED, num_iters: 2689",
		 //QAOA
		{2.13567836819851, 2.15964819966764,
		2.2398735449, 2.18642418344552, 1.55343910949933,
		-0.751462688982622, -1.27298966023324, -2.76699345851517,
		1.71338869914112, -0.144818064314646, 1.12036090687682e-08,
		-6.67099428470163e-08},
		-0.666661182708427,
		-0.600106908259479,
		"QAOA: optimized by diff, FTOL_REACHED, num_iters: 1015"
	),*/

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
0,
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
0,
"QAOA: optimized by diff, MAXEVAL_REACHED, num_iters: 1000"
)
};
