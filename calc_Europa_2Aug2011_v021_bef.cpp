// remove int asID
#include <utility>

void calc(g4m::dataStruct &it, g4m::incrementTab *species[8], g4m::ageStruct &cohort, g4m::ageStruct &newCohort,
		  dat &singleCell, griddata2<char> &managedForest, griddata &maiForest, griddata &rotationForest,
		  griddata &rotationForestNew, griddata &thinningForest, griddata &thinningForestNew, 
		  griddata &harvestSWGrid, griddata &harvestRWGrid, int year, double priceC, int asID, griddata2<float> &OforestShGrid)

		  // griddata &harvestGrid
{
	
	float EmissionsCur=0.;
	float EmissionsProductCur= 0.;  
	float EmissionsLitterCur = 0.;  
	float EmissionsSOCCur = 0.;      
	float EmissionsSlashBurnCur = 0.;
	float EmissionsDeadBurnCur = 0.;
	float EmissionsCRootBurnCur = 0.;
	float EmissionsSOCAfforCur = 0.; 
	float EmissionsLitterAfforCur = 0.;
	float EmissionsAfforCur = 0.;
	float PlantPhytHaBmGr = 0.;
	float plantPhytHaBmGr20 = 0.;
	float PlantPhytHaBmGrBef = 0.;
	float plantPhytHaBmGrBef20 = 0.;
	float plantArea20_rel = 0.;
	//  float plantArea20 = 0.;

	//  float PlantPhytHaBmGrNP =0.;
	float PlantPhytHaBmGrG =0.;
	float PlantPhytHaBlGr = 0.;
	//  float PlantPhytHaBlGrNP = 0.;
	float sawLogsO;
	float restLogsO;
	float LogsOLost;
	float LogsO_FC;
	float LogsO_TH;

	float sawLogsN;		
	float restLogsN;	
	float LogsNLost;
	float LogsN_FC;
	float LogsN_TH;

	harvestWood objHarvWood;
	

	pair <double, double> DeforestedWoodM3;			//.first - sawnwood, .second - restwood
	pair <double, double> HarvestedWoodM3;			//.first - sawnwood, .second - restwood
	DeforestedWoodM3.first = 0.;
	DeforestedWoodM3.second = 0.;
	HarvestedWoodM3.first = 0.;
	HarvestedWoodM3.second = 0.;

	float forNPV_CH = 0.;
	float abBiomassO = 0.;
	float harvWood = 0.; // Total current harvested wood in the cell, m3
	float harvWoodLost = 0.; // Total current "lost" wood in the cell, tC (in remote forests)
	float harvWoodNew = 0.; // Total current harvested wood in the cell, m3
	float harvWoodLostNew = 0.; // Total current "lost" wood in the cell, tC (in remote forests)

	float harvWoodPlus = 0.; // harvestable wood (m3/ha) including 50% of residues
	float harvWoodNewPlus = 0.; // harvestable wood (m3/ha) including 50% of residues  
	float defBiomass = 0.;
	float deforSW = 0;
	float deforRW = 0;
	//  decision.setYear(year);
	//MG: setup forest Age
	int Age = (year-byear)/modTimeStep;
	int xi = (it.x);
	int yi = (it.y);
	int Country = (int)it.COUNTRY[2000];
	float X = (it.x)*GridStepLon+GridStepLon/2-180;
	float Y = (it.y)*GridStepLat+GridStepLat/2-90;


	float OforestShare = singleCell.OforestShare;
	float AforestShare = singleCell.AforestShare;
	if (OforestShare < 0.) {OforestShare = 0.;}	       
	if (AforestShare < 0.) {AforestShare = 0.;}	
	//----------Initialise cohorts-----------------

	int rotationTimeCurr = rotationForest.get(xi,yi);
	singleCell.Rotation = rotationTimeCurr;
	if (year == byear){
		cohort.setU(rotationTimeCurr);       
		newCohort.setU(rotationForest.get(xi,yi));       
		rotationForestNew.set(xi,yi,rotationForest.get(xi,yi));
		thinningForestNew.set(xi,yi,thinningForest.get(xi,yi));
	}
	pair<g4m::ageStruct::v, g4m::ageStruct::v> res = cohort.aging(); //MG
	pair<g4m::ageStruct::v, g4m::ageStruct::v> newRes = newCohort.aging();

	float harvAreaO = 0;
	float harvAreaN = 0;
	float realAreaO = 0;
	float realAreaN = 0;

	if (thinningForest.get(xi,yi)<0) {
		if (OforestShare>0){
			harvAreaO = res.second.area;
			realAreaO = cohort.getArea();

			if (realAreaO>0) {
				//sawnWlost = res.second.sw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				//restWlost  = res.second.rw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//sawnThWlost  = res.first.sw/realAreaO/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				//restThWlost = res.first.rw/realAreaO/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
				//bmHlost  = res.second.bm  * harvAreaO/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				//bmThlost  = res.first.bm/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				//harvLlost = (bmHlost +bmThlost  - (sawnWlost  + restWlost  + sawnThWlost  + restThWlost )); // MG: harvest residues for the set (old) forest tC/ha
				objHarvWood.InitData (res, harvAreaO, realAreaO, modTimeStep);
				LogsOLost = objHarvWood.harvL;
				harvWoodLost = objHarvWood.sumSW + objHarvWood.sumRW + LogsOLost; // Total current "lost" wood in the cell, tC (in remote forests)
			}
		}
		if (AforestShare>0){
			harvAreaN = newRes.second.area;
			realAreaN = newCohort.getArea();
			if (realAreaN>0) {
				//sawnWlostNew = newRes.second.sw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				//restWlostNew = newRes.second.rw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//sawnThWlostNew = newRes.first.sw /realAreaN/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				//restThWlostNew = newRes.first.rw /realAreaN/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
				//bmHlostNew = newRes.second.bm * harvAreaN/realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				//bmThlostNew = newRes.first.bm /realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				//harvLlostNew = (bmThlostNew+bmThlostNew - (sawnWlostNew + restWlostNew + sawnThWlostNew + restThWlostNew)); // MG: harvest residues for the planted (new) forest tC/ha				
				objHarvWood.InitData(newRes,harvAreaN, realAreaN, modTimeStep);
				LogsNLost = objHarvWood.harvL;
				harvWoodLostNew = objHarvWood.sumSW + objHarvWood.sumRW + LogsNLost; // Total current "lost" wood in the cell, tC (in remote forests)
			}
		}  
	} else if (thinningForest.get(xi,yi)>0) {
		if (OforestShare>0){
			harvAreaO = res.second.area;
			realAreaO = cohort.getArea();
			if (realAreaO>0) {
				//sawnW = res.second.sw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				//restW = res.second.rw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//sawnThW = res.first.sw/realAreaO/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				//restThW = res.first.rw/realAreaO/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
				//bmH = res.second.bm * harvAreaO/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				//bmTh = res.first.bm/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				//harvL = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) ; // MG: harvest residues for the set (old) forest tC/ha
				
				objHarvWood.InitData(res, harvAreaO, realAreaO, modTimeStep);
				sawLogsO = objHarvWood.sumSW;
				restLogsO = objHarvWood.sumRW;
				LogsO_FC = objHarvWood.fcut;
				LogsO_TH = objHarvWood.thinning;

			}
		}
//		if ((year==2005)&&(Country==224)){cout<<"objHarvO  "<<objHarvWood.sumRW *  OforestShare * singleCell.LandAreaHa <<endl;}
		if (AforestShare>0){
			harvAreaN = newRes.second.area;
			realAreaN = newCohort.getArea();
			if (realAreaN>0) {
				//sawnWnew = newRes.second.sw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				//restWnew = newRes.second.rw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//sawnThWnew = newRes.first.sw /realAreaN/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				//restThWnew = newRes.first.rw /realAreaN/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
				//bmHnew = newRes.second.bm * harvAreaN/realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				//bmThnew = newRes.first.bm /realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				//harvLNew = (bmHnew+bmThnew - (sawnWnew + restWnew + sawnThWnew + restThWnew)); // MG: harvest residues for the planted (new) forest tC/ha
				
				objHarvWood.InitData(newRes,harvAreaN, realAreaN, modTimeStep);
				sawLogsN = objHarvWood.sumSW;
				restLogsN = objHarvWood.sumRW;
				LogsN_FC = objHarvWood.fcut;
				LogsN_TH = objHarvWood.thinning;
			}
		}   
	}

	if (OforestShare>0){
		
		if (realAreaO > 0){ 
			abBiomassO = cohort.getBm();
			harvWood = (sawLogsO + restLogsO) * it.FTIMBER[byear]; // Total current harvested wood in the cell, m3			
			harvWoodPlus = harvWood + resUse * LogsOLost * it.FTIMBER[byear]; // harvestable wood (m3/ha) including 50% of residues
			if (realAreaO < 1.){
				abBiomassO /= realAreaO;
				             
			}
		}
	}
//	cout << "CALC: harvWood	" << harvWood << endl;
	if (AforestShare>0){
		if (realAreaN > 0){
			harvWoodNew = (sawLogsN + restLogsN) * it.FTIMBER[byear]; // Total current harvested wood in the cell, m3			
			harvWoodNewPlus = harvWoodNew + resUse * LogsNLost * it.FTIMBER[byear]; // harvestable wood (m3/ha) including 50% of residues  
		}
	}
	
	
	int rotMAI = 1;
	if (singleCell.forestShare > 0 && it.CABOVEHA[byear] > 0 && maiForest.get(xi,yi)> 0) {
		rotMAI = species[int(it.SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
		if (rotMAI<1) rotMAI=1;
	}	
	if ((singleCell.Rotation <= 1)||(harvWood < 0)) harvWood = 0.;
	// MG: setup price of carbon
	
	if (year > refYear) {
		coeff.PriceC.clear();		
		coeff.PriceC.insert(0, priceC * it.CORRUPTION[byear]);}
	string regprice, regprice0;
	char prstr[5];
	char regstr[3];
	int2str((int)it.POLESREG[byear],regstr);
	{
		regprice = "re"+string(regstr)+"price0";
	}
	regprice0 = "re"+string(regstr)+"price0";
	{
		// ---- Checking NPV of forestry for current harvest------------       
		dima decision(year
			, it.NPP
			, it.SPOPDENS
			, it.SAGRSUIT
			, it.PRICEINDEX
			, coeff.PriceIndexE
			, it.R
			, coeff.PriceC
			, coeff.PlantingCostsR
			, coeff.PriceLandMinR
			, coeff.PriceLandMaxR
			, coeff.MaxRotInter
			, coeff.MinRotInter
			, coeff.decRateL
			, coeff.decRateS
			, it.FRACLONGPROD
			, coeff.baseline
			, it.FTIMBER
			, coeff.PriceTimberMaxR
			, coeff.PriceTimberMinR
			, coeff.FCuptake
			, it.GDP
			, coeff.HarvLoos
			, OforestShare
			, wprice[regprice]
		, wprice[regprice0].g(2000)
			, rotationTimeCurr
			, harvWood + harvWoodLost * it.FTIMBER[byear] );
		
		forNPV_CH = decision.forValExt();
		//    } // MG: use External NPV
		if (forNPV_CH < 0 || forNPV_CH == 0 || forNPV_CH > 0) 
		{forNPV_CH=forNPV_CH;
		}else{forNPV_CH=0;}
		//-------------------------------------------------------------------------------
	}

	float harvMAI = maiForest.get(xi,yi)*it.FTIMBER[byear]*(1-coeff.HarvLoos[byear]); 

	dima decision(year
		, it.NPP
		, it.SPOPDENS
		, it.SAGRSUIT
		, it.PRICEINDEX
		, coeff.PriceIndexE
		, it.R
		, coeff.PriceC
		, coeff.PlantingCostsR
		, coeff.PriceLandMinR
		, coeff.PriceLandMaxR
		, coeff.MaxRotInter
		, coeff.MinRotInter
		, coeff.decRateL
		, coeff.decRateS
		, it.FRACLONGPROD
		, coeff.baseline
		, it.FTIMBER
		, coeff.PriceTimberMaxR
		, coeff.PriceTimberMinR
		, coeff.FCuptake
		, it.GDP
		, coeff.HarvLoos
		, OforestShare
		, wprice[regprice]
	, wprice[regprice0].g(2000)
		, rotMAI
		, harvMAI);

	//Timber price
	float TimberPrice = 0.;  
	TimberPrice = decision.priceTimberComb(); // MG: use Combined G4M + external TimberPrice

	//Forestry Value
	float fval = 0.;  
	fval = decision.forValComb(); //Combined G4M + external

	if (fval < 0 || fval == 0 || fval > 0) 
	{fval=fval;
	}else{fval=0;}

	//Agricultural Value
	float aval = 0.;
	float minFS =1.;
	aval = decision.agrVal() * lprice[regprice].g(year)/lprice[regprice0].g(2000); // MG: use G4M & GLOBIOM prices

	if (cellInteract)
	{
		minFS = OforestShGrid.GetMin(xi,yi); // minimun of forest share in the neighbour cells
		if ((it.POTVEG[2000] > 2) && (singleCell.prevOForShare < minFS))
		{aval = 2 * aval * (1.-0.5*singleCell.prevOForShare);     // higher agriculture value if closer to recently deforested places
		}
		else{
			aval = 2 * aval * (1. - 0.5*minFS); // agriculture is more attractive if neighbor cells are deforested aval = 2 * aval * (1 - minFS);
		}
	}
	bool roadInfo = true;
	 
	if (roadInfo && singleCell.road > 0.0) {aval *= (1.2+singleCell.road*0.0044); }                // presence of road increases the probability of deforestation   

	bool deforInit = true;
	if (deforInit){
		if ((it.POTVEG[2000] ==1 || it.POTVEG[2000] ==2) && (singleCell.deforPrev > 0.00014)) {aval *= 2.;} // higher agriculture value if closer to recently deforested places
	}
	
	//---------------------------------------------------------
	double crop = it.CROP.g(year);
	double builtup = it.BUILTUP.g(year);
	if (year<2010){crop = it.CROP.g(2010);builtup = it.BUILTUP.g(2010);} //CROP.g(2010) == simulation of the dataset as of November 2010 for consistency of the results

	double maxfor = 1. - (builtup + crop);
	if (maxfor < 0.) maxfor = 0;
	float maxaffor1 = 1. - it.GLOBIOM_RESERVED[2000];
	float maxaffor = 1.;
	if (maxfor < maxaffor1) {maxaffor = maxfor;}else{maxaffor=maxaffor1;}
	//  double freeLand = maxfor - OforestShare;
	float spopdens = it.SPOPDENS.g(year);
	//defrorestation speed
	float defShare = 0.;
	//MG: Afforestation speed
	float affShare = 0.;

	float gdp = 1644.; // MG: gdp definition
	if (it.POPDENS.g(year) > 0.) {
		gdp = 200000. * it.GDP.g(year) / it.POPDENS.g(year);
		gdp = gdp * deflator; // Global GDP deflator GDP(1995)/GDP(2000)=8.807/10 (World Bank)
	}

	// Deforestation due to new infrastructure and minimal crop production
	
	if (OforestShare + AforestShare > maxfor){OforestShare = maxfor-AforestShare;}
	if (OforestShare<0) OforestShare=0;
	
	float OforestShareTmp = OforestShare;
	float defIncome = 0.;
	//MG: Deforestation due to agriculture and forestry land competition
	if(it.AGRSUIT.g(year) > 0 && OforestShare > 0) {
		defShare = 0.05/(1. + exp(1.799e+00
			+ 2.200e-01/OforestShare
			+ 1.663e-01/it.AGRSUIT.g(year)
			- 4.029e-02 * it.POPDENS.g(year)
			+ 5.305e-04 * it.POPDENS.g(year) * it.POPDENS.g(year)
			+ 1.282e-04 * gdp ));

		if (deforInit){
			if (it.POTVEG[2000]==1 || it.POTVEG[2000]==2)              
			{
				if (year == byear)       
				{         
					if (defShare > 1.e-17 && singleCell.deforPrev > 1.e-17) {singleCell.deforRateCoeffCell = (singleCell.deforPrev)/defShare; } //To match observed defor rate in the tropics	
					else if (defShare > 1.e-17) {singleCell.deforRateCoeffCell = (0.000139)/defShare;} // less than 0.00014, the threshold of increased agriculture attractiveness
				}else
				{             defShare *= singleCell.deforRateCoeffCell; 
				}
			}
		}    
		defShare *= deforRate_opt[Country-1] ;
		// MG: End of gdp definition
		float defIncome = 0.;
		//Pay if Carbon get's to air (Harvest goes to products)
		
		float pDefIncome = it.CABOVEHA[byear] *           
			(TimberPrice * it.FTIMBER[byear]
		* (1. -coeff.HarvLoos[byear])
			- coeff.PriceC.g(year) * (1 + it.R.g(year))
			* (it.FRACLONGPROD.g(year) * coeff.decRateL.g(year)
			/ (coeff.decRateL.g(year) + it.R.g(year))
			+ it.FRACLONGPROD.g(year) * coeff.decRateS.g(year) //MG: mistake: must be DECRATES but not DECRATEL
			/ (coeff.decRateS.g(year) + it.R.g(year))));
		//*******************************************************************************************************
		
		//*******************************************************************************************************
		float pDefLoss = - coeff.PriceC.g(year) * (1 + it.R.g(year))
			* (it.CLITTERHA[byear]*(0.3*it.DECWOOD.g(year)/(it.DECWOOD.g(year)+it.R.g(year))
			+ 0.7*it.DECHERB.g(year)/(it.DECHERB.g(year)+it.R.g(year))) 
			+ it.CBELOWHA[byear]*0.3*it.DECHERB.g(year)/(it.DECHERB.g(year)+it.R.g(year))
			+ it.SOCHA[byear]*it.DECSOC.g(year)/(it.DECSOC.g(year)+it.R.g(year)));
		//Immediate Pay if deforested (Slash and Burn)
	
		float sDefIncome = it.CABOVEHA[byear] *              
			(TimberPrice * it.FTIMBER[byear]
		* (1. -coeff.HarvLoos[byear])
			- coeff.PriceC.g(year));
		float sDefLoss = - coeff.PriceC.g(year) * (it.CDEADHA[byear]+ it.CBELOWHA.g(year)*0.7);
		defIncome = pDefIncome * (1. - it.SLASHBURN[byear])
			+ sDefIncome * it.SLASHBURN[byear]
		+ pDefLoss + sDefLoss;
		if ((aval + defIncome) > (fval * Hurdle_opt[Country-1])){ 
			OforestShare -= defShare * modTimeStep; // Decrease Forest share

			if (OforestShare > maxaffor) OforestShare = maxaffor;
			if (OforestShare < 0.) OforestShare = 0.;
		}
	} // end for if(it.AGRSUIT.g(year) > 0 && OforestShare > 0)

	//MG: Afforestation
	            
	float afformax = 0.5*it.AFFORMAX[2000];
	if (Country==71 || Country==177) afformax = 1;
	
	if (OforestShareTmp + AforestShare < maxaffor && AforestShare < afformax){ 
		if ((it.POTVEG[2000])<=9 && (it.POTVEG[2000])>0 && (it.NPP[2000])>0)  { // MG: We afforest only places, where potential vegetation is forest and savanna   
			affShare = 0.01/(1+exp(1.+0.1/it.AGRSUIT.g(year)+1/(0.001*gdp)));
			if (affShare < 1./(it.LANDAREA[byear]*1000000.)) affShare = 0.;   // minimun one tree (1m^2)
			affShare *= afforRate_opt[Country-1]; 
			if (aval < (fval + coeff.PriceC[byear]* (((it.SOCHA[byear]*0.4/(8*decision.rotInter()))
				+ 5/(1.053*decision.rotInter()))/it.R[byear])) * Hurdle_opt[Country-1])
			{
				
				if (AforestShare < 0.) {
					AforestShare = 0.;
					affShare = 0.;
				}	    				
				if (OforestShareTmp + AforestShare > maxaffor) AforestShare = maxaffor - OforestShareTmp;
				if (AforestShare < 0.) AforestShare = 0.;				
			}
		}
	}else{ 
		
		OforestShare = maxfor - AforestShare;
		if (OforestShare<0) OforestShare=0;
	}      
	if (OforestShare>singleCell.prevOForShare) OforestShare=singleCell.prevOForShare;    
	if (AforestShare<singleCell.AforestSharePrev) AforestShare=singleCell.AforestSharePrev;     
	  
	singleCell.forestAgeShare[Age] = AforestShare - singleCell.AforestSharePrev;
	if (singleCell.forestAgeShare[Age] > 0.02 * modTimeStep) singleCell.forestAgeShare[Age] = 0.02 * modTimeStep; //Limit affor speed to 2% of cell area (0.5deg) per year Hannes Bottcher personal communication, 2013
	if (singleCell.forestAgeShare[Age] < 0.) {singleCell.forestAgeShare[Age]=0.;}
	AforestShare = singleCell.AforestSharePrev + singleCell.forestAgeShare[Age];
	newCohort.afforest(singleCell.forestAgeShare[Age]); // MG: Afforest


	singleCell.deforestA[Age] = singleCell.prevOForShare - OforestShare;
	if (singleCell.deforestA[Age] > 0.05 * modTimeStep) singleCell.deforestA[Age] = 0.05 * modTimeStep; //Limit defo speed to 5% of cell area (0.5deg) per year (see Kindermann et al. 2007)
	if (singleCell.deforestA[Age] < 0.) {singleCell.deforestA[Age]=0.;}
	if (singleCell.deforestA[Age] > 1.) {singleCell.deforestA[Age]=1.;}
	OforestShare = singleCell.prevOForShare - singleCell.deforestA[Age];

	if (year < 2001) { // MG: We start deforestation after 2000, because we have initial data for 2000 (discussion with Hannes Bottcher 10.09.09)
		OforestShare = singleCell.OforestShare;
		AforestShare = 0.;  
	}

	float defArea = 0;
	if (singleCell.deforestA[Age]>0){
		//MG: Correcting the bug causing relAreaO approaching zero before real deforestation starts
		g4m::ageStruct::v resDefor;
		if (year < 2001) { // MG: We start deforestation after 2000, because we have initial data for 2000 (discussion with Hannes Bottcher 10.09.09)  
			g4m::ageStruct cohortTmp = cohort;
			resDefor = cohortTmp.deforest(singleCell.deforestA[Age],0); //MG: deforest set forest and get deforested biomass. Devide by refForShare to get per ha value
			defArea = resDefor.area;
			if (defArea>0 && realAreaO>0)
			{				
				defBiomass = resDefor.bm / defArea / modTimeStep; // tC/per ha (only stem!)
				deforSW = resDefor.sw / modTimeStep;// / defArea; // tC per cell
				deforRW = resDefor.rw / modTimeStep;/// defArea;  // tC per cell
			}
		}else{   
			resDefor = cohort.deforest(singleCell.deforestA[Age],0); //MG: deforest set forest and get deforested biomass. Devide by refForShare to get per ha value
			defArea = resDefor.area;
			if (defArea>0 && realAreaO>0)
			{				
				defBiomass = resDefor.bm / defArea / modTimeStep; // tC/per ha (only stem!)
				deforSW = resDefor.sw / modTimeStep; // /defArea;// / defArea; // tC per cell
				deforRW = resDefor.rw / modTimeStep; // / defArea;// / defArea;  // tC per cell
			}
		}
	}
	
	singleCell.OforestShare = OforestShare;
	singleCell.AforestShare = AforestShare;  
	singleCell.forestShare = OforestShare + AforestShare;
	double deforestHa = singleCell.deforestA[Age]*singleCell.LandAreaHa; 
	singleCell.deforestHaTot += deforestHa;
	double afforestHa = singleCell.forestAgeShare[Age]*singleCell.LandAreaHa;
	singleCell.afforestHaTot += afforestHa;
	singleCell.ProdLongA[Age] = defBiomass*it.FRACLONGPROD[byear]*(1-coeff.HarvLoos[byear]) * (1. - it.SLASHBURN[byear])*deforestHa;
	singleCell.ProdShortA[Age] = singleCell.ProdLongA[Age];
	singleCell.LitterA[Age] = (it.CLITTERHA[byear] + (BEF(int(it.SPECIESTYPE[byear])-1,defBiomass,it.FTIMBER[2000])-1) * defBiomass) * deforestHa; //MG:BEF: we assume that the non-stem aboveground biomass goes to the littecr pool
	singleCell.FineRootA[Age] = it.CBELOWHA[byear]*0.3*deforestHa;
	singleCell.SOCA[Age] = it.SOCHA[byear]*deforestHa;
	singleCell.SOCaffor[Age] = it.SOCHA[byear] * afforestHa;
	singleCell.LitterAffor[Age] = it.CLITTERHA[byear] * afforestHa;

	//********* Emissions from deforestation *************
	for (int i=0;i<=Age;i++) {
		// Calculate Emissions from deforestation in current cell for current year caused by decomposition
		float decRateL_dec = exp(-1*coeff.decRateL[2000]*modTimeStep*i);
		float decRateS_dec = exp(-1*coeff.decRateS[2000]*modTimeStep*i);
		float decHerb_dec = exp(-1*it.DECHERB[2000]*modTimeStep*i);
		float decWood_dec = exp(-1*it.DECWOOD[2000]*modTimeStep*i);
		float decSOC_dec = exp(-1*it.DECSOC[2000]*modTimeStep*i);
		
		EmissionsProductCur = EmissionsProductCur + singleCell.ProdLongA[i]*(1-decRateL_dec)+singleCell.ProdShortA[i]*(1-decRateS_dec); // tC/modTimeStep

		EmissionsLitterCur = EmissionsLitterCur  
			+ singleCell.LitterA[i]*(0.3*(1-decWood_dec)+0.7*(1-decHerb_dec))
			+ singleCell.FineRootA[i]*(1-decHerb_dec);

		if (singleCell.SOCA[i] > it.SOCHA[2000]*0.6*singleCell.deforestA[i]*singleCell.LandAreaHa) {                                                 
			EmissionsSOCCur = EmissionsSOCCur + singleCell.SOCA[i] * (1-decSOC_dec);
		}  
		// MG: Recalculate carbon pools

		singleCell.ProdLongA[i] = singleCell.ProdLongA[i] * decRateL_dec; 
		singleCell.ProdShortA[i] = singleCell.ProdShortA[i] * decRateS_dec;
		singleCell.LitterA[i] = singleCell.LitterA[i] * 0.3 * decWood_dec
			+ singleCell.LitterA[i] * 0.7 * decHerb_dec;

		singleCell.FineRootA[i] = singleCell.FineRootA[i] * decHerb_dec;

		if (singleCell.SOCA[i]>it.SOCHA[2000]*0.6*singleCell.deforestA[i]*singleCell.LandAreaHa) { 
			singleCell.SOCA[i] = singleCell.SOCA[i] * decSOC_dec;
			
		}
	}
	//Emissions from deforestation in current cell for current year caused by burning  
	EmissionsSlashBurnCur = EmissionsSlashBurnCur + abBiomassO*it.SLASHBURN.g(year)* deforestHa;
	EmissionsDeadBurnCur = EmissionsDeadBurnCur + it.CDEADHA[byear]* deforestHa;
	EmissionsCRootBurnCur = EmissionsCRootBurnCur + it.CBELOWHA.g(year)*0.7 * deforestHa;
	//Emissions for current cell summed over years
	singleCell.EmissionsProduct += EmissionsProductCur;
	singleCell.EmissionsLitter += EmissionsLitterCur;
	singleCell.EmissionsSOC += EmissionsSOCCur;
	singleCell.EmissionsSlashBurn += EmissionsSlashBurnCur;
	singleCell.EmissionsDeadBurn += EmissionsDeadBurnCur;
	singleCell.EmissionsCRootBurn += EmissionsCRootBurnCur;
	//Total emissions in current cell for current year
	EmissionsCur = EmissionsProductCur+EmissionsLitterCur+EmissionsSOCCur+EmissionsSlashBurnCur
		+ EmissionsDeadBurnCur+EmissionsCRootBurnCur;
	//Total emissions in current cell summed over years       
	singleCell.EmissionsTot += EmissionsCur;    
	EmissionsCurCountry[Country] += EmissionsCur;
	//*************** END Emissions from deforestation ****************
	//*************** Afforestation "negative" emissions block ********
	
	for (int ia=0; ia<Age; ia++) {
		if (singleCell.forestAgeShare[ia]>0) {
			float abovePhCur = newCohort.getBm((Age-ia)*modTimeStep); // afforested (stem) biomass of age (Age-ia)*modTimeStep per ha
			float abovePhCurBef = abovePhCur *  // afforested biomass of age (Age-ia)*modTimeStep per ha  * afforShare[ia]
				BEF(int(it.SPECIESTYPE[byear])-1,abovePhCur,it.FTIMBER[2000])*singleCell.forestAgeShare[ia];	  
			abovePhCur *= singleCell.forestAgeShare[ia]; // afforested (stem) biomass per ha * afforShare[ia]
			PlantPhytHaBmGr += abovePhCur;// * modTimeStep; // Stem biomass of planted forest, tC per ha 
			PlantPhytHaBmGrBef += abovePhCurBef; //Aboveground biomass of planted forest, tC per ha 
			if ((Age-ia) * modTimeStep >= 20)  //We track area and phytomass of new forest over 20 y.o.
			{plantPhytHaBmGrBef20 += abovePhCurBef;//afforested biomass per ha * afforShare (forest over 20 y.o.)
			plantArea20_rel += singleCell.forestAgeShare[ia];// relative area of forest over 20 y.o.
			}

			if (singleCell.LitterAffor[ia]<5*singleCell.forestAgeShare[ia]*singleCell.LandAreaHa) {
				float CurEmissionsLitterAfforCur = 0.95 * pow((1.-exp(-0.1*abovePhCurBef/singleCell.forestAgeShare[ia])),3)
					* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep;
				EmissionsLitterAfforCur += CurEmissionsLitterAfforCur;
				singleCell.LitterAffor[ia] += CurEmissionsLitterAfforCur;
			}
			if (singleCell.SOCaffor[ia]<=it.SOCHA[2000]*1.4* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa) {
				if (it.POTVEG[2000]==4||it.POTVEG[2000]==6) {
					EmissionsSOCAfforCur += 0.04 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
					/ (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
						* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep;

					singleCell.SOCaffor[ia] += 0.04 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
					/ (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3)
						* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep;                   
				}
				if (it.POTVEG[2000]==8) {
					EmissionsSOCAfforCur += 0.2 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
					/ (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
						* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep;
					singleCell.SOCaffor[ia] += 0.2 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
					/ (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
						* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep; 
				} else {
					EmissionsSOCAfforCur += 0.35 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
					/ (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
						* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep; 
					singleCell.SOCaffor[ia] += 0.35 * pow((1.-exp(-1.2*singleCell.LitterAffor[ia]
					/ (singleCell.forestAgeShare[ia]*singleCell.LandAreaHa))),3) 
						* singleCell.forestAgeShare[ia]*singleCell.LandAreaHa * modTimeStep;
				}  
			}           // End for if SOCaffor <= ...
			//     }            // End for if relAffor>0
		}             // End for if forestShare > 0
	}               // End for age loop
	
	float coeffBl = 0.22;
	if (it.POTVEG[byear]==1 || it.POTVEG[byear]==2){coeffBl = 0.18;}
	if (it.POTVEG[byear]==6 || it.POTVEG[byear]==7) {coeffBl = 0.25;}
	if (it.POTVEG[byear]==3 || it.POTVEG[byear]==4 || it.POTVEG[byear]==5||(it.POTVEG[byear]==8)) {coeffBl = 0.22;}
	float CurPlantPhytBmGr = 0;
	if (AforestShare > 0)
	{
		PlantPhytHaBlGr = PlantPhytHaBmGrBef * coeffBl; //Belowground phytomass

		if (singleCell.AforestSharePrev > 0)
			CurPlantPhytBmGr=(PlantPhytHaBmGrBef - singleCell.prevPlantPhytHaBmGrBef) * singleCell.LandAreaHa; // Aboveground phytomass increment in the cell, tC/modTimeStep
	}
	float CurPlantPhytBlGr = (PlantPhytHaBlGr-singleCell.prevPlantPhytHaBlGr) * singleCell.LandAreaHa;// Belowground phytomass increment in the cell, tC/modTimeStep 
	// Emissions in current cell summed over years 
	singleCell.EmLitterAffor += EmissionsLitterAfforCur;
	singleCell.EmSOCAffor += EmissionsSOCAfforCur;
	//     Total (negative) emissions from afforestation for current year and current cell
	EmissionsAfforCur = (CurPlantPhytBmGr+CurPlantPhytBlGr) + EmissionsLitterAfforCur+EmissionsSOCAfforCur;
	EmissionsCurAfforCountry[Country] += EmissionsAfforCur;
	//Total emissions in current cell summed over years       
	singleCell.EmissionsAffor += EmissionsAfforCur;
		
	DeforestedWoodM3.first = deforSW * (1.-it.SLASHBURN[byear]) * it.FTIMBER[byear]*singleCell.LandAreaHa;		// deforested sawnwood in the cell, m3
	DeforestedWoodM3.second = deforRW * (1.-it.SLASHBURN[byear]) * it.FTIMBER[byear]*singleCell.LandAreaHa;		// deforested restwood in the cell, m3

	float harvestDfM3Ha = (deforSW+deforRW)*(1.-it.SLASHBURN[byear]) * it.FTIMBER[byear]; // wood obtained at deforestation, m3/ha
	float deforWoodTotM3 = harvestDfM3Ha * singleCell.LandAreaHa; // wood obtained from harvesting (FM) in the cell, m3
	float harvestFmTotM3 = (harvWood*OforestShare+harvWoodNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa;// Harvested new and old wood (from FM) in the cell
	float harvestTotM3 = harvestFmTotM3 + deforWoodTotM3;// total wood (from FM and deforestation) in the cell

	char countrych[4];
	int ii = Country;
	int2str(ii,countrych);
	string countryprice = "re"+string(countrych)+"price0";
	
	// OT:  25.06.15
	HarvestedWoodM3.second = (restLogsO  * OforestShare + restLogsN *  singleCell.AforestSharePrev) * singleCell.LandAreaHa * it.FTIMBER [byear]
							 + DeforestedWoodM3.second;	
	if (CountriesSawnwoodHarvestM3Year.get(Country,year) > wprod_sawnwood[countryprice].g(year)){
	HarvestedWoodM3.second += (sawLogsO *  OforestShare + sawLogsN  *  singleCell.AforestSharePrev)*  singleCell.LandAreaHa * it.FTIMBER [byear] 
						   + DeforestedWoodM3.first;
	}else{
	HarvestedWoodM3.first =   (sawLogsO *  OforestShare + sawLogsN  *  singleCell.AforestSharePrev)*  singleCell.LandAreaHa * it.FTIMBER [byear] 
	  				          +   DeforestedWoodM3.first;
	}


	//HarvestedWoodM3.first =   (sawLogsO *  OforestShare + sawLogsN  *  singleCell.AforestSharePrev)*  singleCell.LandAreaHa * it.FTIMBER [byear] 
	//						 + DeforestedWoodM3.first;   
	//HarvestedWoodM3.second = (restLogsO  * OforestShare + restLogsN *  singleCell.AforestSharePrev) * singleCell.LandAreaHa * it.FTIMBER [byear]
	//						 + DeforestedWoodM3.second;
	 	
	//						 if (year==2005){cout << "HarvestedWoodM3.first" << HarvestedWoodM3.first << endl;}

	harvestSWGrid.set(xi, yi, HarvestedWoodM3.first);
	harvestRWGrid.set(xi, yi, HarvestedWoodM3.second);
	//cout << "harvestSWGrid	" << harvestSWGrid.get << endl;

	float harvestFcM3 = (LogsO_FC + LogsN_FC)* it.FTIMBER[byear]; // current harvested wood from final cut in the cell, m3/ha
	float harvestThM3 =(LogsO_TH + LogsN_TH) * it.FTIMBER[byear]; // current harvested wood from thinning in the cell, m3/ha
	float harvestTotPlusM3 = (harvWoodPlus*OforestShare + harvWoodNewPlus*singleCell.AforestSharePrev+(defBiomass-(1-resUse)*(defBiomass-(deforSW+deforRW))) * (1. - it.SLASHBURN[byear])*singleCell.deforestA[Age]* it.FTIMBER[byear])*singleCell.LandAreaHa; //Total current harvested wood in the cell, m3, including residues (harvest losses)
	
	float CAI_m3ha = ((abBiomassO-singleCell.ObiomassPrev)/modTimeStep + sawLogsO + restLogsO + LogsOLost + harvWoodLost) * it.FTIMBER[byear]; //Current annual increment m3/ha
	CAI_m3ha = CAI_m3ha < 0. ? 0. : CAI_m3ha; // m3/ha/year
	float phytNewIncr = 0.; if (AforestShare > 0) {phytNewIncr = CurPlantPhytBmGr/(AforestShare * singleCell.LandAreaHa);}
	float CAI_new_m3ha = (phytNewIncr + harvWoodLostNew + sawLogsN + restLogsN + LogsNLost) * it.FTIMBER[2000];
	CAI_new_m3ha = CAI_new_m3ha < 0. ? 0. : CAI_new_m3ha;
	// should we account for deforested biomass when estimating CAI???
	// should we account for deforested biomass when estimating harvest?

	
	float FMsink_ab = (abBiomassO * BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[byear]) 
		- singleCell.ObiomassPrev * BEF(int(it.SPECIESTYPE[byear])-1,singleCell.ObiomassPrev,it.FTIMBER[byear]))
		* 3.6666666667*unitConv/modTimeStep; // tCO2/ha/year
	float FMsink_bm = FMsink_ab * (1+coeffBl);  // tCO2/ha/year

	int biomassRot = 1;
	if (thinningForest.get(xi,yi) < 0) {
		biomassRot = species[int(it.SPECIESTYPE[byear])-1]->gU(abBiomassO, maiForest.get(xi,yi))+1;  
	} else {
		biomassRot = species[int(it.SPECIESTYPE[byear])-1]->gUSdTab(abBiomassO, maiForest.get(xi,yi), thinningForest.get(xi,yi))+1;
	}

	if (fff.is_open()) {  
		if (countriesList.find(Country) != countriesList.end()){

			fff<<asID<<"\t"<<year<<"\t"<<Country<<"\t"<<OforestShare<<"\t"<<AforestShare<<"\t"<<it.CABOVEHA[byear]<<"\t"<<abBiomassO<<"\t"<<PlantPhytHaBmGr<<"\t"<<maiForest.get(xi,yi)*it.FTIMBER[byear]<<"\t"<< rotMAI;
			fff<<"\t"<<(LogsO_TH) * it.FTIMBER[byear]<<"\t"<<(LogsO_FC)* it.FTIMBER[byear]<<"\t"<<(LogsN_TH) * it.FTIMBER[byear]<<"\t"<<(LogsN_FC)* it.FTIMBER[byear]<<"\t"<<harvestTotM3<<"\t"<<HarvestedWoodM3.first<<"\t"<<HarvestedWoodM3.second<<"\t"<<harvMAI<<"\t"<<singleCell.rotBiomass;
			fff<<"\t"<< (fval + coeff.PriceC[byear]* (((it.SOCHA[byear]*0.4/(8*decision.rotInter()))+5/(1.053*decision.rotInter()))/it.R[byear])) * Hurdle_opt[Country-1];
			fff<<"\t"<<fval<<"\t"<<aval<<"\t"<<rotationTimeCurr<<"\t"<<rotationForestNew.get(xi,yi)<<"\t"<<TimberPrice;
			fff<<"\t"<<(thinningForest.get(xi,yi));
			fff<<"\t"<<species[int(it.SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), optimMAI);
			fff<<"\t"<<realAreaO;
			fff<<"\t"<<singleCell.deforestA[Age];
			fff<<"\t"<<defArea;
			fff<<"\t"<<BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[byear]);
			
			fff<< endl;
		}
	}
	// Output results for countries

	CountriesNforCover.inc(Country,year,AforestShare * singleCell.LandAreaHa);
	CountriesAfforHaYear.inc(Country,year,afforestHa/modTimeStep); 
	CountriesAfforAccumHa.inc(Country,year,singleCell.afforestHaTot);
	double AforBm = (PlantPhytHaBmGrBef + PlantPhytHaBlGr) * singleCell.LandAreaHa ; // accumulated tot living biomass of planted forest, tC
	CountriesNforTotC.inc(Country,year,AforBm);
	CountriesNfor_stem_C.inc(Country,year,PlantPhytHaBmGr * singleCell.LandAreaHa); //accumulated living stem biomass of planted forest, tC
	CountriesNfor_ab_C.inc(Country,year,PlantPhytHaBmGrBef * singleCell.LandAreaHa);//accumulated living aboveground biomass of planted forest, tC
	CountriesAfforCYear.inc(Country,year,EmissionsAfforCur*3.6666666667*unitConv/modTimeStep);  //????? per year or modTimeStep???
	CountriesAfforCYear_ab.inc(Country,year, CurPlantPhytBmGr * 3.6666666667*unitConv/modTimeStep);   //????? per year 
	CountriesAfforCYear_bl.inc(Country,year, CurPlantPhytBlGr * 3.6666666667*unitConv/modTimeStep);    //????? per year or modTimeStep???
	float EmAffBm = (CurPlantPhytBmGr + CurPlantPhytBlGr) * 3.6666666667*unitConv/modTimeStep; // CO2 sink in forest living biomass of planted forest, mtCO2/year
	CountriesAfforCYear_biom.inc(Country,year, EmAffBm);   //????? per year or modTimeStep???
	CountriesAfforCYear_dom.inc(Country,year, EmissionsLitterAfforCur * 3.6666666667*unitConv/modTimeStep);    
	CountriesAfforCYear_soil.inc(Country,year, EmissionsSOCAfforCur * 3.6666666667*unitConv/modTimeStep);     
	CountriesAfforCover20.inc(Country,year, plantArea20_rel * singleCell.LandAreaHa);   //????? per year or modTimeStep???
	CountriesAfforTotC20.inc(Country,year,plantPhytHaBmGrBef20 * singleCell.LandAreaHa);
	//--------- 
	CountriesOforCover.inc(Country,year,OforestShare * singleCell.LandAreaHa);
	CountriesDeforHaYear.inc(Country,year,deforestHa/modTimeStep);  // ha/year
	CountriesOfor_stem_C.inc(Country,year,abBiomassO*OforestShare*singleCell.LandAreaHa);
	CountriesOfor_ab_C.inc(Country,year,abBiomassO*BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[2000])*OforestShare*singleCell.LandAreaHa);
	double OforBm = (abBiomassO*BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[2000])+it.CBELOWHA[byear]); // tot living biomass of old forest, tC/ha
	CountriesOforC_biom.inc(Country,year, OforBm*OforestShare*singleCell.LandAreaHa);   
	CountriesDeforCYear.inc(Country,year,EmissionsCur*3.6666666667*unitConv/modTimeStep);    // tC/year
	CountriesDeforCYear_ab.inc(Country,year, defBiomass*BEF(int(it.SPECIESTYPE[byear])-1,defBiomass,it.FTIMBER[2000])*deforestHa*3.6666666667*unitConv/modTimeStep);
	CountriesDeforCYear_bl.inc(Country,year,it.CBELOWHA[2000]*deforestHa*3.6666666667*unitConv/modTimeStep);      
	float EmDefBm = (defBiomass*BEF(int(it.SPECIESTYPE[byear])-1,defBiomass,it.FTIMBER[2000])+it.CBELOWHA[2000])*3.6666666667*unitConv; // CO2 lost from living biomass of old forests due to deforestation, mtCO2/hayear
	CountriesDeforCYear_biom.inc(Country,year,EmDefBm*deforestHa);       // tC/ha/year
	CountriesDeforCYear_dom.inc(Country,year, (EmissionsDeadBurnCur+EmissionsLitterCur)*3.6666666667*unitConv/modTimeStep);     
	CountriesDeforCYear_soil.inc(Country,year, EmissionsSOCCur*3.6666666667*unitConv/modTimeStep);     
	//--------- 
	

	CountriesSawnwoodHarvestM3Year.inc(Country, year, HarvestedWoodM3.first);
	CountriesRestwoodHarvestM3Year.inc(Country, year, HarvestedWoodM3.second);	
		
	CountriesWoodHarvestM3Year.inc(Country,year,harvestTotM3);
	CountriesWoodHarvestPlusM3Year.inc(Country,year,harvestTotPlusM3);    
	CountriesWoodHarvestFmM3Year.inc(Country,year,harvestFmTotM3);
	CountriesWoodHarvestDfM3Year.inc(Country,year,deforWoodTotM3);
	CountriesWoodHarvestFc_oldM3Year.inc(Country,year, (LogsO_FC) * it.FTIMBER[byear] *OforestShare*singleCell.LandAreaHa);
	CountriesWoodHarvestTh_oldM3Year.inc(Country,year, (LogsO_TH) * it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa);
	CountriesWoodHarvestFc_newM3Year.inc(Country,year,(LogsN_FC) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa);//corrected 21.08.2013
	CountriesWoodHarvestTh_newM3Year.inc(Country,year,(LogsN_TH) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa);//corrected 21.08.2013     
	CountriesWoodLoosCYear.inc(Country,year,(harvWoodLost*OforestShare+harvWoodLostNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa*3.6666666667*unitConv);
	CountriesHarvLossesYear.inc(Country,year,(LogsOLost*OforestShare+LogsNLost*singleCell.AforestSharePrev)* it.FTIMBER[byear]*singleCell.LandAreaHa);
	//---------
	if (thinningForest.get(xi,yi)>0 & thinningForestNew.get(xi,yi)>0){
		CountriesManagedCount.inc(Country,year,1);
		CountriesManagedForHa.inc(Country,year, (AforestShare+OforestShare) * singleCell.LandAreaHa);
	}else if (thinningForest.get(xi,yi)>0 & thinningForestNew.get(xi,yi)<0){
		CountriesManagedCount.inc(Country,year,1);
		CountriesManagedForHa.inc(Country,year,OforestShare * singleCell.LandAreaHa);
	}
	CountriesMAI.inc(Country,year,maiForest.get(xi,yi)*it.FTIMBER[byear]);
	CountriesCAI.inc(Country,year,CAI_m3ha * OforestShare * singleCell.LandAreaHa);

	CountriesCAI_new.inc(Country,year,CAI_new_m3ha * singleCell.AforestSharePrev*singleCell.LandAreaHa);     
	CountriesFM.inc(Country,year,FMsink_ab*OforestShare*singleCell.LandAreaHa);
	CountriesFMbm.inc(Country,year,FMsink_bm*OforestShare*singleCell.LandAreaHa);

	profit = ((LogsO_FC + LogsO_TH + LogsOLost*resUse)* it.FTIMBER[byear] * TimberPrice - decision.plantingCosts() * cohort.getArea(0)) ; // profit only from harvesting old forest per ha
	CountriesProfit.inc(Country,year,profit);
	if (thinningForest.get (xi,yi)>0) {
		CountryRotation.inc(Country,year,rotationTimeCurr);};

	
	CountryregSawnwoodHarvestM3Year.inc(Country, year, HarvestedWoodM3.first);
	CountryregRestwoodHarvestM3Year.inc(Country, year, HarvestedWoodM3.second);	

	CountryregWoodHarvestM3Year.inc(Country,year,harvestTotM3);
	CountryregWoodHarvestFmM3Year.inc(Country,year,harvestFmTotM3);
	CountryregWoodHarvestDfM3Year.inc(Country,year,deforWoodTotM3);   
//std::cout << "CALC:: CountriesSawnwoodHarvestM3Year " << CountriesSawnwoodHarvestM3Year.get(Country,year) << endl;
 
	
	CountryregWsawnprod.setVal(Country,year,wprod_sawnwood[countryprice].g(year));
	CountryregWrestprod.setVal(Country,year,wprod_restwood[countryprice].g(year));
//std::cout << "CALC::CountryregWsawnprod  "<<CountryregWsawnprod.get(Country,year)<< endl;
	CountryregMaxHarvest.inc(Country,year,deforWoodTotM3 + 
		(maiForest.get(xi,yi)*it.FTIMBER[byear]*(1.-coeff.HarvLoos[byear])*OforestShare+ harvWoodNew * singleCell.AforestSharePrev)*singleCell.LandAreaHa);


	OforestShGrid.set(xi,yi,OforestShare);
	OforestShGrid.update();
	singleCell.prevOForShare = OforestShare;
	singleCell.prevPlantPhytHaBmGr = PlantPhytHaBmGr;//  new forest stem wood, tC/ha
	singleCell.prevPlantPhytHaBmGrBef = PlantPhytHaBmGrBef;//  new forest aboveground biomass, tC/ha
	singleCell.prevPlantPhytHaBmGrBef20 = plantPhytHaBmGrBef20;//  new forest over 20 y.o. aboveground biomass, tC/ha
	singleCell.prevPlantPhytHaBlGr = PlantPhytHaBlGr;// new forest belowground biomass, tC/ha
	singleCell.rotBiomass = biomassRot;
	singleCell.SD = thinningForest.get(xi,yi);
	singleCell.FMsink = FMsink_ab; // FM sink in GgCO2/yr/ha
	singleCell.FMsink_Bm = FMsink_bm;
	singleCell.ObiomassPrev = abBiomassO; // old forest stem wood, tC/ha
	singleCell.DeforWoodM3.first = DeforestedWoodM3.first;
	singleCell.DeforWoodM3.second = DeforestedWoodM3.second;
	singleCell.deforWoodTotM3 = deforWoodTotM3;
	singleCell.CAI = CAI_m3ha;
	if (year > 2000){singleCell.deforestShare = singleCell.deforestA[Age];}else{singleCell.deforestShare = 0.;}
	singleCell.afforestHaYear = afforestHa;
	singleCell.deforestHaYear = deforestHa;
	if (year > 2000){singleCell.deforPrev = singleCell.deforestA[Age];}
	singleCell.harvSawnWood = (sawLogsO + sawLogsN +deforSW * (1.-it.SLASHBURN[byear])) * it.FTIMBER[byear]; 
	singleCell.harvRestWood = (restLogsO + restLogsN +deforRW * (1.-it.SLASHBURN[byear])) * it.FTIMBER[byear];
	singleCell.harvestTot = (harvestFcM3 + harvestThM3) + harvestDfM3Ha; 
	singleCell.harvestFcM3Ha = harvestFcM3;
	singleCell.harvestThM3Ha = harvestThM3;
	singleCell.oforestBm = OforBm;  

	singleCell.emissionsD_Bm = EmDefBm;  

	if (OforestShare>0)
	{
		singleCell.emissionsD_S = EmissionsSOCCur/(OforestShare*singleCell.LandAreaHa)*3.6666666667*unitConv;}else{singleCell.emissionsD_S = 0;
	}
	if (AforestShare * singleCell.LandAreaHa>0)
	{
		singleCell.aforestBm = AforBm/(singleCell.AforestShare * singleCell.LandAreaHa);
		singleCell.emissionsA_Bm = EmAffBm/(singleCell.AforestShare * singleCell.LandAreaHa);                  
		singleCell.emissionsA_S = EmissionsSOCAfforCur/(AforestShare * singleCell.LandAreaHa)*3.6666666667*unitConv;
	}else{
		singleCell.emissionsA_S = 0;
		singleCell.aforestBm = 0;
		singleCell.emissionsA_Bm = 0;      
	}


	singleCell.AforestSharePrev = AforestShare;
	//cout<<"---------------end Calc      ------------"<<endl;
}
