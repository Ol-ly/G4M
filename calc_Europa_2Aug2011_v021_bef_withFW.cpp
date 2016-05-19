// remove int asID

void calc(g4m::dataStruct &it, g4m::incrementTab *species[8], g4m::ageStruct &cohort, g4m::ageStruct &newCohort,
		  dat &singleCell, griddata2<char> &managedForest, griddata &maiForest, griddata &rotationForest,
		  griddata &rotationForestNew, griddata &thinningForest, griddata &thinningForestNew, 
		  griddata &harvestGrid, int year, double priceC, int asID, griddata2<float> &OforestShGrid)

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
	float PlantPhytHaBmGrG =0.;
	float PlantPhytHaBlGr = 0.;
	float sawnW = 0.;
	float restW = 0.;
	float fuelW = 0.;										//1.
	float sawnThW = 0.;
	float restThW = 0.;
	float fuelThW = 0.;									//2.
	float sawnWlost = 0.;
	float restWlost = 0.;
	float fuelWlost = 0.;									//3.
	float sawnThWlost = 0.;
	float restThWlost = 0.;
	float fuelThWlost = 0.;								//4.
	float sawnWnew = 0.;
	float restWnew = 0.;
	float fuelWnew = 0.;									//5.
	float sawnThWnew = 0.;
	float restThWnew = 0.;
	float fuelThWnew = 0.;								//6.
	float sawnWlostNew = 0.;
	float restWlostNew = 0.;
	float fuelWlostNew = 0.;								//7.
	float sawnThWlostNew = 0.;
	float restThWlostNew = 0.;
	float fuelThWlostNew = 0.;							//8.
	float forNPV_CH = 0.;
	float abBiomassO = 0.;
	float harvLlost = 0.;
	float harvLlostNew = 0.;
	float harvL = 0.;
	float harvLNew = 0.;

	float harvWood = 0.; // Total current harvested wood in the cell, m3
	float harvWoodLost = 0.; // Total current "lost" wood in the cell, tC (in remote forests)
	float harvWoodNew = 0.; // Total current harvested wood in the cell, m3
	float harvWoodLostNew = 0.; // Total current "lost" wood in the cell, tC (in remote forests)

	float harvWoodPlus = 0.; // harvestable wood (m3/ha) including 50% of residues
	float harvWoodNewPlus = 0.; // harvestable wood (m3/ha) including 50% of residues  
	float bmHlost = 0.;
	float bmThlost = 0.;
	float bmHlostNew = 0.;
	float bmThlostNew = 0.;
	float bmH = 0.;
	float bmTh = 0.;
	float bmHnew = 0.;
	float bmThnew = 0.;

	float defBiomass = 0.;
	float deforSW = 0;
	float deforRW = 0;
	float deforFW = 0;									//9.
	//  decision.setYear(year);
	//MG: setup forest Age
	int Age = (year-byear)/modTimeStep;
	int xi = (it.x);
	int yi = (it.y);
	int Country = (int)it.COUNTRY[2000];
	float X = (it.x)*GridStepLon+GridStepLon/2-180;		//
	float Y = (it.y)*GridStepLat+GridStepLat/2-90;
	float sawlogs;
	float restlogs;
	//  float forestShare0 = 0;
	//   if (it.FOREST[2000]+(it.CROP[2000])+(it.BUILTUP[2000])>1)
	//      {forestShare0 = (1-(it.CROP[2000]+it.BUILTUP[2000]));}
	//   else {forestShare0 = it.FOREST[2000];}  
	//if (asID==3895) cout<<"OforestShare1\t"<<singleCell.OforestShare<<endl;            
	float OforestShare = singleCell.OforestShare;
	float AforestShare = singleCell.AforestShare;
	if (OforestShare < 0.) {OforestShare = 0.;}	       
	if (AforestShare < 0.) {AforestShare = 0.;}	
	//----------Initialise cohorts-----------------
	int rotationTimeCurr = rotationForest.get(xi,yi);

	//Check condition
	// if ((year==2000)&&(Country==224)&&(thinningForest.get(xi,yi)>0)){
	//  cout << "rotationTimeCurr " << rotationTimeCurr << endl;
	//}

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
		//   if (Country == 33 || Country == 156) {  
		if (OforestShare>0){
			harvAreaO = res.second.area;
			realAreaO = cohort.getArea();
			if (realAreaO>0) {
				sawnWlost = res.second.sw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				restWlost  = res.second.rw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut.			  																												//10.  
				sawnThWlost  = res.first.sw/realAreaO/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				restThWlost = res.first.rw/realAreaO/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
				//11.  
				bmHlost  = res.second.bm  * harvAreaO/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				bmThlost  = res.first.bm/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				harvLlost = (bmHlost +bmThlost  - (sawnWlost  + restWlost  + sawnThWlost  + restThWlost + fuelWlost +fuelThWlost));	//MG: harvest residues for the set (old) forest tC/ha							//12.
			}
		}
		if (AforestShare>0){
			harvAreaN = newRes.second.area;
			realAreaN = newCohort.getArea();
			if (realAreaN>0) {
				sawnWlostNew = newRes.second.sw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				restWlostNew = newRes.second.rw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//13. 
				sawnThWlostNew = newRes.first.sw /realAreaN/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				restThWlostNew = newRes.first.rw /realAreaN/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.         
				//14. 
				bmHlostNew = newRes.second.bm * harvAreaN/realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				bmThlostNew = newRes.first.bm /realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				harvLlostNew = (bmThlostNew+bmThlostNew - (sawnWlostNew + restWlostNew + sawnThWlostNew + restThWlostNew +fuelWlostNew +fuelThWlostNew)); // MG: harvest residues for the planted (new) forest tC/ha			//15.
			}

		}  
	} else if (thinningForest.get(xi,yi)>0) {
		if (OforestShare>0){
			harvAreaO = res.second.area;
			realAreaO = cohort.getArea();
			if (realAreaO>0) {
				sawnW = res.second.sw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				restW = res.second.rw * harvAreaO/realAreaO/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//16. 
				sawnThW = res.first.sw/realAreaO/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				restThW = res.first.rw/realAreaO/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
				//17. 
				bmH = res.second.bm * harvAreaO/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				bmTh = res.first.bm/realAreaO/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				harvL = (bmH+bmTh - (sawnW + restW + fuelW + sawnThW + restThW + fuelThW)) ; // MG: harvest residues for the set (old) forest tC/ha												//18.
			}
		}
		if (AforestShare>0){
			harvAreaN = newRes.second.area;
			realAreaN = newCohort.getArea();
			if (realAreaN>0) {
				sawnWnew = newRes.second.sw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
				restWnew = newRes.second.rw * harvAreaN/realAreaN/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
				//19.
				sawnThWnew = newRes.first.sw /realAreaN/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
				restThWnew = newRes.first.rw /realAreaN/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.            
				//20.
				bmHnew = newRes.second.bm * harvAreaN/realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
				bmThnew = newRes.first.bm /realAreaN/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
				harvLNew = (bmHnew+bmThnew - (sawnWnew + restWnew + fuelWnew + sawnThWnew + restThWnew + fuelThWnew)); // MG: harvest residues for the planted (new) forest tC/ha			//21.
			}
		}   
	}

	if (OforestShare>0){
		if (realAreaO > 0){ 
			abBiomassO = cohort.getBm();
			harvWood = (sawnW + restW + sawnThW + restThW ) * it.FTIMBER[byear]; // Total current harvested wood in the cell, m3												//22.
			harvWoodLost = (sawnWlost + restWlost + fuelWlost + sawnThWlost + restThWlost +fuelThWlost + harvLlost); // Total current "lost" wood in the cell, tC (in remote forests)			//23.
			harvWoodPlus = harvWood + resUse * harvL * it.FTIMBER[byear]; // harvestable wood (m3/ha) including 50% of residues
			if (realAreaO < 1.){
				abBiomassO /= realAreaO;            
			}
		}
	}
	if (AforestShare>0){
		if (realAreaN > 0){
			harvWoodNew = (sawnWnew + restWnew + fuelWnew + sawnThWnew + restThWnew + fuelThWnew) * it.FTIMBER[byear]; // Total current harvested wood in the cell, m3						//24.
			harvWoodLostNew = (sawnWlostNew + restWlostNew + fuelWlostNew + sawnThWlostNew + restThWlostNew + fuelThWlostNew + harvLlostNew); // Total current "lost" wood in the cell, tC (in remote forests)	//25.
			harvWoodNewPlus = harvWoodNew + resUse * harvLNew * it.FTIMBER[byear]; // harvestable wood (m3/ha) including 50% of residues  
		}

	}

	//  if ((year==2000)&&(Country==224)){ cout<<"X "<<X<<"		Y "<<Y<<"		HarvWood "<<harvWood<<endl;};
	// Rotation time fitted to get certain biomass under certain MAI (w/o thinning)
	//  int biomasRot = 0;
	int rotMAI = 1;
	if (singleCell.forestShare > 0 && it.CABOVEHA[byear] > 0 && maiForest.get(xi,yi)> 0) {
		//    biomasRot = species[int(it.SPECIESTYPE[byear])-1]->gU(it.CABOVEHA"][byear], maiForest.get(xi,yi), 1);
		rotMAI = species[int(it.SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
		if (rotMAI<1) rotMAI=1;
	}
	//  if (biomasRot < 0) biomasRot = 0;
	//  if ((biomasRot == 0)||(harvWood < 0)) harvWood = 0.;
	if ((singleCell.Rotation <= 1)||(harvWood < 0)) harvWood = 0.;
	// MG: setup price of carbon
	//  double PriceCi = 0;
	if (year > refYear) {
		coeff.PriceC.clear();
		////  coeff.PriceC.insert(0, PriceCi * LinPrice2020[Age] * it.CORRUPTION"][byear]);
		coeff.PriceC.insert(0, priceC * it.CORRUPTION[byear]);}
	string regprice, regprice0;
	char prstr[5];
	char regstr[3];
	int2str((int)it.POLESREG[byear],regstr);
	//  if (PriceCi!=0) {
	//    int2str(PriceC,prstr);
	//    regprice = "re"+string(regstr)+"price"+ string(prstr);
	//  } else 
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
		//    if (year<=2005) {
		//      forNPV_CH = decision.forVal(); // MG: use internal NPV
		//    } else {
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

	float spopdens = it.SPOPDENS.g(year);
	//deforestation speed
	float defShare = 0.;
	//MG: Afforestation speed
	float affShare = 0.;

	float gdp = 1644.; // MG: gdp definition
	if (it.POPDENS.g(year) > 0.) {
		gdp = 200000. * it.GDP.g(year) / it.POPDENS.g(year);
		gdp = gdp * deflator; // Global GDP deflator GDP(1995)/GDP(2000)=8.807/10 (World Bank)
	}


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
		//  if ((aval + defIncome) > ((fval ) * Hurdle_opt[Country-1] ) && (OforestShare > 0)) {  // MG: adjust the multiplier to account for a forest saving policy in some countries
		if ((aval + defIncome) > (fval * Hurdle_opt[Country-1])){ 
			OforestShare -= defShare * modTimeStep; // Decrease Forest share
			if (OforestShare > maxaffor) OforestShare = maxaffor;
			if (OforestShare < 0.) OforestShare = 0.;
		}
	} // end for if(it.AGRSUIT.g(year) > 0 && OforestShare > 0)
	//cout<<"NPP= "<< it.NPP"][byear]<<"\t aval= "<<aval<<"\t fval= "<<fval<<"\t maxfor= "<<maxfor<<"\t potveg= "<<it.POTVEG"][byear]<<"\t OforestShare= "<<OforestShare <<endl;

	
	
	// OT: added 07.11.14 Exctra+ction of area distribution between age classes within each grid cell
	

	//				***File for Age classes structure
/*	
		if (outfile.is_open()) {
			double currBiom [22][10]={0.};
			double currBiom0 = cohort.Biomass(0)/0.993351;
			for(int i = 1; i<=220; i++){
				for (int m = 0; m<=21; m++){
					for (int n = 0; n<=9; n++){
						
						currBiom[m][n] = cohort.Biomass(i)/0.993351;
						cout <<"currBiom[m][n]  "<< currBiom[m][n]<<endl;
					}
				}
			}
			double biom10Years [22] = {0.};
			for (int m=0; m<=21; m++){
				for (int n=0; n<=9; n++){
					biom10Years [m] +=currBiom[m][n]; 
				}
outfile<<asID<<"\t"<<X<<"\t"<<Y<< "\t"<<m<<"\t"<<biom10Years [m]/10<<endl;
			}
			

*/

	/*	double currBiom [221]={0.};
		double classArea[221]={0.};
		double biomTabVal [221]={0.};
		double biom10Years [22] = {0.};
		for (int n = 21; n<=21; n++){
		for (int i = 0.; i<=220; i++){
			
			currBiom[i] = cohort.Biomass(i)/0.993351;
			double k = i + 0.0;
			classArea[i]= cohort.getArea(k)*singleCell.LandAreaHa*OforestShare/0.993351;
			biomTabVal [i]=cohort.BiomSdTab(i);
			

			outfile<<asID<<"\t"<<X<<"\t"<<Y<< "\t"<<i<<"\t"<<currBiom[i]<<"\t"<<classArea[i]<<"\t"<<biomTabVal[i]<<endl;
			};
		double biom10Years [22] = {0.};
		for (int n = 21; n<=21; n++){
		*/	
	//	};

if ((year==2000)&&(Country==224)) {	
	
	double siteIndex = maiForest.get(xi,yi);
	extractData(cohort, singleCell, siteIndex, X, Y);

		//if((X==25.75)&&(Y==49.75)){
		//double i;
		//double B_0_10 = (cohort.getBm(0)+cohort.getBm(1)+cohort.getBm(2)+cohort.getBm(3)+cohort.getBm(4)+cohort.getBm(5)+cohort.getBm(6)+cohort.getBm(7)+cohort.getBm(8)+cohort.getBm(9)+cohort.getBm(10))/0.993351/10;
		//cout <<"B_0_10   "<<B_0_10<<endl;
		//for (i=0.; i<=220.; i++){
		//	if (outfile3.is_open()) {
		//		outfile3<<asID<<"\t"<<X<<"\t"<<Y<<"\t"<<i<<"\t"<<cohort.getBm(i)/0.993351<<endl; 
		//	};
		//};
	//	};
		
		double CA_0_10 = (cohort.getArea(0)+cohort.getArea(1)+cohort.getArea(2)+cohort.getArea(3)+cohort.getArea(4)+cohort.getArea(5)+cohort.getArea(6)+cohort.getArea(7)+cohort.getArea(8)+cohort.getArea(9)+cohort.getArea(10))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_11_20 = (cohort.getArea(11)+cohort.getArea(12)+cohort.getArea(13)+cohort.getArea(14)+cohort.getArea(15)+cohort.getArea(16)+cohort.getArea(17)+cohort.getArea(18)+cohort.getArea(19)+cohort.getArea(20))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_21_30 = (cohort.getArea(21)+cohort.getArea(22)+cohort.getArea(23)+cohort.getArea(24)+cohort.getArea(25)+cohort.getArea(26)+cohort.getArea(27)+cohort.getArea(28)+cohort.getArea(29)+cohort.getArea(30))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_31_40 = (cohort.getArea(31)+cohort.getArea(32)+cohort.getArea(33)+cohort.getArea(34)+cohort.getArea(35)+cohort.getArea(36)+cohort.getArea(37)+cohort.getArea(38)+cohort.getArea(39)+cohort.getArea(40))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_41_50 = (cohort.getArea(41)+cohort.getArea(42)+cohort.getArea(43)+cohort.getArea(44)+cohort.getArea(45)+cohort.getArea(46)+cohort.getArea(47)+cohort.getArea(48)+cohort.getArea(49)+cohort.getArea(50))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_51_60 = (cohort.getArea(51)+cohort.getArea(52)+cohort.getArea(53)+cohort.getArea(54)+cohort.getArea(55)+cohort.getArea(56)+cohort.getArea(57)+cohort.getArea(58)+cohort.getArea(59)+cohort.getArea(60))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_61_70 = (cohort.getArea(61)+cohort.getArea(62)+cohort.getArea(63)+cohort.getArea(64)+cohort.getArea(65)+cohort.getArea(66)+cohort.getArea(67)+cohort.getArea(68)+cohort.getArea(69)+cohort.getArea(70))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_71_80 = (cohort.getArea(71)+cohort.getArea(72)+cohort.getArea(73)+cohort.getArea(74)+cohort.getArea(75)+cohort.getArea(76)+cohort.getArea(77)+cohort.getArea(78)+cohort.getArea(79)+cohort.getArea(80))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_81_90 = (cohort.getArea(81)+cohort.getArea(82)+cohort.getArea(83)+cohort.getArea(84)+cohort.getArea(85)+cohort.getArea(86)+cohort.getArea(87)+cohort.getArea(88)+cohort.getArea(89)+cohort.getArea(90))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_91_100 = (cohort.getArea(91)+cohort.getArea(92)+cohort.getArea(93)+cohort.getArea(94)+cohort.getArea(95)+cohort.getArea(96)+cohort.getArea(97)+cohort.getArea(98)+cohort.getArea(99)+cohort.getArea(100))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_101_110 = (cohort.getArea(101)+cohort.getArea(102)+cohort.getArea(103)+cohort.getArea(104)+cohort.getArea(105)+cohort.getArea(106)+cohort.getArea(107)+cohort.getArea(108)+cohort.getArea(109)+cohort.getArea(110))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_111_120 = (cohort.getArea(111)+cohort.getArea(112)+cohort.getArea(113)+cohort.getArea(114)+cohort.getArea(115)+cohort.getArea(116)+cohort.getArea(117)+cohort.getArea(118)+cohort.getArea(119)+cohort.getArea(120))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_121_130 = (cohort.getArea(121)+cohort.getArea(122)+cohort.getArea(123)+cohort.getArea(124)+cohort.getArea(125)+cohort.getArea(126)+cohort.getArea(127)+cohort.getArea(128)+cohort.getArea(129)+cohort.getArea(130))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_131_140 = (cohort.getArea(131)+cohort.getArea(132)+cohort.getArea(133)+cohort.getArea(134)+cohort.getArea(135)+cohort.getArea(136)+cohort.getArea(137)+cohort.getArea(138)+cohort.getArea(139)+cohort.getArea(140))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_141_150 = (cohort.getArea(141)+cohort.getArea(142)+cohort.getArea(143)+cohort.getArea(144)+cohort.getArea(145)+cohort.getArea(146)+cohort.getArea(147)+cohort.getArea(148)+cohort.getArea(149)+cohort.getArea(150))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_151_160 = (cohort.getArea(151)+cohort.getArea(152)+cohort.getArea(153)+cohort.getArea(154)+cohort.getArea(155)+cohort.getArea(156)+cohort.getArea(157)+cohort.getArea(158)+cohort.getArea(159)+cohort.getArea(160))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_161_170 = (cohort.getArea(161)+cohort.getArea(162)+cohort.getArea(163)+cohort.getArea(164)+cohort.getArea(165)+cohort.getArea(166)+cohort.getArea(167)+cohort.getArea(168)+cohort.getArea(169)+cohort.getArea(170))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_171_180= (cohort.getArea(171)+cohort.getArea(172)+cohort.getArea(173)+cohort.getArea(174)+cohort.getArea(175)+cohort.getArea(176)+cohort.getArea(177)+cohort.getArea(178)+cohort.getArea(179)+cohort.getArea(180))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_181_190= (cohort.getArea(181)+cohort.getArea(182)+cohort.getArea(183)+cohort.getArea(184)+cohort.getArea(185)+cohort.getArea(186)+cohort.getArea(187)+cohort.getArea(188)+cohort.getArea(189)+cohort.getArea(190))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_191_200= (cohort.getArea(191)+cohort.getArea(192)+cohort.getArea(193)+cohort.getArea(194)+cohort.getArea(195)+cohort.getArea(196)+cohort.getArea(197)+cohort.getArea(198)+cohort.getArea(199)+cohort.getArea(200))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_201_210= (cohort.getArea(201)+cohort.getArea(202)+cohort.getArea(203)+cohort.getArea(204)+cohort.getArea(205)+cohort.getArea(206)+cohort.getArea(207)+cohort.getArea(208)+cohort.getArea(209)+cohort.getArea(210))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		double CA_211_220= (cohort.getArea(211)+cohort.getArea(212)+cohort.getArea(213)+cohort.getArea(214)+cohort.getArea(215)+cohort.getArea(216)+cohort.getArea(217)+cohort.getArea(218)+cohort.getArea(219)+cohort.getArea(220))*(singleCell.LandAreaHa*OforestShare)/0.993351;
		
		

		
		
/*
		double B_0_10 = (cohort.getBm(0)+cohort.getBm(1)+cohort.getBm(2)+cohort.getBm(3)+cohort.getBm(4)+cohort.getBm(5)+cohort.getBm(6)+cohort.getBm(7)+cohort.getBm(8)+cohort.getBm(9)+cohort.getBm(10))/0.993351/10;
		double B_11_20 = (cohort.getBm(11)+cohort.getBm(12)+cohort.getBm(13)+cohort.getBm(14)+cohort.getBm(15)+cohort.getBm(16)+cohort.getBm(17)+cohort.getBm(18)+cohort.getBm(19)+cohort.getBm(20))/0.993351/10;
		double B_21_30 = (cohort.getBm(21)+cohort.getBm(22)+cohort.getBm(23)+cohort.getBm(24)+cohort.getBm(25)+cohort.getBm(26)+cohort.getBm(27)+cohort.getBm(28)+cohort.getBm(29)+cohort.getBm(30))/0.993351/10;
		double B_31_40 = (cohort.getBm(31)+cohort.getBm(32)+cohort.getBm(33)+cohort.getBm(34)+cohort.getBm(35)+cohort.getBm(36)+cohort.getBm(37)+cohort.getBm(38)+cohort.getBm(39)+cohort.getBm(40))/0.993351/10;
		double B_41_50 = (cohort.getBm(41)+cohort.getBm(42)+cohort.getBm(43)+cohort.getBm(44)+cohort.getBm(45)+cohort.getBm(46)+cohort.getBm(47)+cohort.getBm(48)+cohort.getBm(49)+cohort.getBm(50))/0.993351/10;
		double B_51_60 = (cohort.getBm(51)+cohort.getBm(52)+cohort.getBm(53)+cohort.getBm(54)+cohort.getBm(55)+cohort.getBm(56)+cohort.getBm(57)+cohort.getBm(58)+cohort.getBm(59)+cohort.getBm(60))/0.993351/10;
		double B_61_70 = (cohort.getBm(61)+cohort.getBm(62)+cohort.getBm(63)+cohort.getBm(64)+cohort.getBm(65)+cohort.getBm(66)+cohort.getBm(67)+cohort.getBm(68)+cohort.getBm(69)+cohort.getBm(70))/0.993351/10;
		double B_71_80 = (cohort.getBm(71)+cohort.getBm(72)+cohort.getBm(73)+cohort.getBm(74)+cohort.getBm(75)+cohort.getBm(76)+cohort.getBm(77)+cohort.getBm(78)+cohort.getBm(79)+cohort.getBm(80))/0.993351/10;
		double B_81_90 = (cohort.getBm(81)+cohort.getBm(82)+cohort.getBm(83)+cohort.getBm(84)+cohort.getBm(85)+cohort.getBm(86)+cohort.getBm(87)+cohort.getBm(88)+cohort.getBm(89)+cohort.getBm(90))/0.993351/10;
		double B_91_100 = (cohort.getBm(91)+cohort.getBm(92)+cohort.getBm(93)+cohort.getBm(94)+cohort.getBm(95)+cohort.getBm(96)+cohort.getBm(97)+cohort.getBm(98)+cohort.getBm(99)+cohort.getBm(100))/0.993351/10;
		double B_101_110 = (cohort.getBm(101)+cohort.getBm(102)+cohort.getBm(103)+cohort.getBm(104)+cohort.getBm(105)+cohort.getBm(106)+cohort.getBm(107)+cohort.getBm(108)+cohort.getBm(109)+cohort.getBm(110))/0.993351/10;
		double B_111_120 = (cohort.getBm(111)+cohort.getBm(112)+cohort.getBm(113)+cohort.getBm(114)+cohort.getBm(115)+cohort.getBm(116)+cohort.getBm(117)+cohort.getBm(118)+cohort.getBm(119)+cohort.getBm(120))/0.993351/10;
		double B_121_130 = (cohort.getBm(121)+cohort.getBm(122)+cohort.getBm(123)+cohort.getBm(124)+cohort.getBm(125)+cohort.getBm(126)+cohort.getBm(127)+cohort.getBm(128)+cohort.getBm(129)+cohort.getBm(130))/0.993351/10;
		double B_131_140 = (cohort.getBm(131)+cohort.getBm(132)+cohort.getBm(133)+cohort.getBm(134)+cohort.getBm(135)+cohort.getBm(136)+cohort.getBm(137)+cohort.getBm(138)+cohort.getBm(139)+cohort.getBm(140))/0.993351/10;
		double B_141_150 = (cohort.getBm(141)+cohort.getBm(142)+cohort.getBm(143)+cohort.getBm(144)+cohort.getBm(145)+cohort.getBm(146)+cohort.getBm(147)+cohort.getBm(148)+cohort.getBm(149)+cohort.getBm(150))/0.993351/10;
		double B_151_160 = (cohort.getBm(151)+cohort.getBm(152)+cohort.getBm(153)+cohort.getBm(154)+cohort.getBm(155)+cohort.getBm(156)+cohort.getBm(157)+cohort.getBm(158)+cohort.getBm(159)+cohort.getBm(160))/0.993351/10;
		double B_161_170 = (cohort.getBm(161)+cohort.getBm(162)+cohort.getBm(163)+cohort.getBm(164)+cohort.getBm(165)+cohort.getBm(166)+cohort.getBm(167)+cohort.getBm(168)+cohort.getBm(169)+cohort.getBm(170))/0.993351/10;
		double B_171_180= (cohort.getBm(171)+cohort.getBm(172)+cohort.getBm(173)+cohort.getBm(174)+cohort.getBm(175)+cohort.getBm(176)+cohort.getBm(177)+cohort.getBm(178)+cohort.getBm(179)+cohort.getBm(180))/0.993351/10;
		double B_181_190= (cohort.getBm(181)+cohort.getBm(182)+cohort.getBm(183)+cohort.getBm(184)+cohort.getBm(185)+cohort.getBm(186)+cohort.getBm(187)+cohort.getBm(188)+cohort.getBm(189)+cohort.getBm(190))/0.993351/10;
		double B_191_200= (cohort.getBm(191)+cohort.getBm(192)+cohort.getBm(193)+cohort.getBm(194)+cohort.getBm(195)+cohort.getBm(196)+cohort.getBm(197)+cohort.getBm(198)+cohort.getBm(199)+cohort.getBm(200))/0.993351/10;
		double B_201_210= (cohort.getBm(201)+cohort.getBm(202)+cohort.getBm(203)+cohort.getBm(204)+cohort.getBm(205)+cohort.getBm(206)+cohort.getBm(207)+cohort.getBm(208)+cohort.getBm(209)+cohort.getBm(210))/0.993351/10;
		double B_211_220= (cohort.getBm(211)+cohort.getBm(212)+cohort.getBm(213)+cohort.getBm(214)+cohort.getBm(215)+cohort.getBm(216)+cohort.getBm(217)+cohort.getBm(218)+cohort.getBm(219)+cohort.getBm(220))/0.993351/10;
			
*/
		
/*if (outfile.is_open()) {
			
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc01"<<"\t"<<CA_0_10<<endl; //"\t"<<B_0_10<< "\t"<<refStock_0_10<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc02"<<"\t"<<CA_11_20<<endl; //"\t"<<B_11_20<< "\t"<<refStock_11_20<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc03"<<"\t"<<CA_21_30<<endl; //"\t"<<B_21_30<< "\t"<<refStock_21_30<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc04"<<"\t"<<CA_31_40<<endl; //"\t"<<B_31_40<< "\t"<<refStock_31_40<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc05"<<"\t"<<CA_41_50<<endl; //"\t"<<B_41_50<< "\t"<<refStock_41_50<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc06"<<"\t"<<CA_51_60<<endl; //"\t"<<B_51_60<< "\t"<<refStock_51_60<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc07"<<"\t"<<CA_61_70<<endl; //"\t"<<B_61_70<< "\t"<<refStock_61_70<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc08"<<"\t"<<CA_71_80<<endl; //"\t"<<B_71_80<< "\t"<<refStock_71_80<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc09"<<"\t"<<CA_81_90<<endl; //"\t"<<B_81_90<< "\t"<<refStock_81_90<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc10"<<"\t"<<CA_91_100<<endl; //"\t"<<B_91_100<< "\t"<<refStock_91_100<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc11"<<"\t"<<CA_101_110<<endl; //"\t"<<B_101_110<< "\t"<<refStock_101_110<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc12"<<"\t"<<CA_111_120<<endl; //"\t"<<B_111_120<< "\t"<<refStock_111_120<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc13"<<"\t"<<CA_121_130<<endl; //"\t"<<B_121_130<< "\t"<<refStock_121_130<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc14"<<"\t"<<CA_131_140<<endl; //"\t"<<B_131_140<< "\t"<<refStock_131_140<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc15"<<"\t"<<CA_141_150<<endl; //"\t"<<B_141_150<< "\t"<<refStock_141_150<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc16"<<"\t"<<CA_151_160<<endl; //"\t"<<B_151_160<< "\t"<<refStock_151_160<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc17"<<"\t"<<CA_161_170<<endl; //"\t"<<B_161_170<< "\t"<<refStock_161_170<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc18"<<"\t"<<CA_171_180<<endl; //"\t"<<B_171_180<< "\t"<<refStock_171_180<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc19"<<"\t"<<CA_181_190<<endl; //"\t"<<B_181_190<< "\t"<<refStock_181_190<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc20"<<"\t"<<CA_191_200<<endl; //"\t"<<B_191_200<< "\t"<<refStock_191_200<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc21"<<"\t"<<CA_201_210<<endl; //"\t"<<B_201_210<< "\t"<<refStock_201_210<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tc22"<<"\t"<<CA_211_220<<endl; //"\t"<<B_211_220<< "\t"<<refStock_211_220<<endl;

				

		};	   //				***End of file*/
    

	}; //OT: End of extraction 

	
	
	

	//MG: Afforestation
	//   if (OforestShare + AforestShare < maxfor){
	//   if (OforestShareTmp + AforestShare < maxfor){                    
	//   if (OforestShareTmp + AforestShare < maxaffor){                    
	float afformax = 0.5*it.AFFORMAX[2000];
	//   if (Country==71 || Country==177) afformax = it.AFFORMAX[2000];
	if (Country==71 || Country==177) afformax = 1;
	//if (Country==224) cout<<"maxaffor=\t"<<maxaffor<<"\tafformax=\t"<<afformax<<"\tPOTVEG=\t"<<it.POTVEG[2000]<<"\tNPP=\t"<<it.NPP[2000]<<endl;		OT:07.11.14
	//   if (OforestShareTmp + AforestShare < maxaffor && AforestShare < it.AFFORMAX[2000]){ 
	//   if (OforestShareTmp + AforestShare < maxaffor && AforestShare < 0.5*it.AFFORMAX[2000]){ 
	if (OforestShareTmp + AforestShare < maxaffor && AforestShare < afformax){ 
		//      if ((it.POTVEG[2000])<9 && (it.POTVEG[2000])>0 && (it.NPP[2000])>0)  { // MG: We afforest only places, where potential vegetation is forest    
		if ((it.POTVEG[2000])<=9 && (it.POTVEG[2000])>0 && (it.NPP[2000])>0)  { // MG: We afforest only places, where potential vegetation is forest and savanna   
			affShare = 0.01/(1+exp(1.+0.1/it.AGRSUIT.g(year)+1/(0.001*gdp)));
			if (affShare < 1./(it.LANDAREA[byear]*1000000.)) affShare = 0.;   // minimun one tree (1m^2)
			//if (Country==224) cout<<"Country=\t"<<Country<<"\taffshare= "<<affShare<<"\t afforRate_opt= "<<afforRate_opt[Country-1]<<endl;					OT:07.11.14
			affShare *= afforRate_opt[Country-1]; 
			//if (Country==224) cout<<"defshare= "<<defShare<<"\t affShare= "<<affShare<<"\taval=\t"<<aval<<"\tfval*H=\t"<<fval*Hurdle_opt[Country-1]<<endl;		OT:07.11.14
			if (aval < (fval + coeff.PriceC[byear]* (((it.SOCHA[byear]*0.4/(8*decision.rotInter()))
				+ 5/(1.053*decision.rotInter()))/it.R[byear])) * Hurdle_opt[Country-1])
			{
				AforestShare += affShare * modTimeStep;
				//cout<< "AforestShare0= "<<AforestShare<<endl;
				//          if (AforestShare > maxfor) AforestShare = maxfor;
				if (AforestShare < 0.) {
					AforestShare = 0.;
					affShare = 0.;
				}	    
				//              if (OforestShare + AforestShare > maxfor) AforestShare = maxfor - OforestShare;
				//              if (OforestShareTmp + AforestShare > maxfor) AforestShare = maxfor - OforestShareTmp;
				if (OforestShareTmp + AforestShare > maxaffor) AforestShare = maxaffor - OforestShareTmp;
				if (AforestShare < 0.) AforestShare = 0.;
				//if (Country==68){
				//cout<<"asID=\t"<<asID<<"\tAffor OforestShare=\t"<<OforestShare<<"\tprevOForShare=\t"<<singleCell.prevOForShare<<"\tmaxfor=\t"<<maxfor<<"\tAforestShare=\t"<<AforestShare<<"\taffShare=\t"<<affShare<<endl;               
				//}
				//cout<< "AforestShare2= "<<AforestShare<<endl;
			}
		}
	}else{ 
		//if (Country==68){
		//cout<<"asID=\t"<<asID<<"\tElse OforestShare=\t"<<OforestShare<<"\tprevOForShare=\t"<<singleCell.prevOForShare<<"\tmaxfor=\t"<<maxfor<<"\tAforestShare=\t"<<AforestShare<<endl;               
		//}
		OforestShare = maxfor - AforestShare;
		if (OforestShare<0) OforestShare=0;
	}      
	if (OforestShare>singleCell.prevOForShare) OforestShare=singleCell.prevOForShare;    
	if (AforestShare<singleCell.AforestSharePrev) AforestShare=singleCell.AforestSharePrev;         

	singleCell.forestAgeShare[Age] = AforestShare - singleCell.AforestSharePrev;
	//  if (singleCell.forestAgeShare[Age] > 0.05 * modTimeStep) singleCell.forestAgeShare[Age] = 0.05 * modTimeStep; //Limit affor speed to 5% of cell area (0.5deg) per year similar to defor speed (see Kindermann et al. 2007)
	if (singleCell.forestAgeShare[Age] > 0.02 * modTimeStep) singleCell.forestAgeShare[Age] = 0.02 * modTimeStep; //Limit affor speed to 2% of cell area (0.5deg) per year Hannes Bottcher personal communication, 2013
	//  if (singleCell.forestAgeShare[Age] > 0.01 * modTimeStep) singleCell.forestAgeShare[Age] = 0.01 * modTimeStep; //Limit affor speed to 2% of cell area (0.5deg) per year Hannes Bottcher personal communication, 2013
	//  if (singleCell.forestAgeShare[Age] > 0.005 * modTimeStep) singleCell.forestAgeShare[Age] = 0.005 * modTimeStep; //Limit affor speed to 2% of cell area (0.5deg) per year Hannes Bottcher personal communication, 2013
	if (singleCell.forestAgeShare[Age] < 0.) {singleCell.forestAgeShare[Age]=0.;}
	AforestShare = singleCell.AforestSharePrev + singleCell.forestAgeShare[Age];
	//cout<<"singleCell.AforestSharePrev= "<<singleCell.AforestSharePrev<<"\t singleCell.forestAgeShare[Age]= "<<singleCell.forestAgeShare[Age]<<endl;
	newCohort.afforest(singleCell.forestAgeShare[Age]); // MG: Afforest


	singleCell.deforestA[Age] = singleCell.prevOForShare - OforestShare;
	if (singleCell.deforestA[Age] > 0.05 * modTimeStep) singleCell.deforestA[Age] = 0.05 * modTimeStep; //Limit defo speed to 5% of cell area (0.5deg) per year (see Kindermann et al. 2007)
	if (singleCell.deforestA[Age] < 0.) {singleCell.deforestA[Age]=0.;}
	if (singleCell.deforestA[Age] > 1.) {singleCell.deforestA[Age]=1.;}
	OforestShare = singleCell.prevOForShare - singleCell.deforestA[Age];

	if (year < 2001) { // MG: We start deforestation after 2000, because we have initial data for 2000 (discussion with Hannes Bottcher 10.09.09)
		//    if (year < 2000) { // MG: We start deforestation after 2000, because we have initial data for 2000 (discussion with Hannes Bottcher 10.09.09)
		OforestShare = singleCell.OforestShare;
		AforestShare = 0.;  
	}

	float defArea = 0;
	if (singleCell.deforestA[Age]>0){
		//MG: Correcting the bug causing relAreaO approaching zero before real deforestation starts
		g4m::ageStruct::v resDefor;
		if (year < 2001) { // MG: We start deforestation after 2000, because we have initial data for 2000 (discussion with Hannes Bottcher 10.09.09)  
			g4m::ageStruct cohortTmp = cohort;
			resDefor = cohortTmp.deforest(singleCell.deforestA[Age],0); //MG: deforest set forest and get deforested biomass. Devide by refForSha re to get per ha value
			defArea = resDefor.area;
			if (defArea>0 && realAreaO>0)
			{
				//      defBiomass = (resDefor.sw + resDefor.rw + resDefor.hv);///deforestA[Age]; // deforested aboveground biomass per ha // only stemwood !!!

				defBiomass = resDefor.bm / defArea / modTimeStep; // tC/per ha (only stem!)
				deforSW = resDefor.sw / modTimeStep;// / defArea; // tC per cell
				deforRW = resDefor.rw / modTimeStep;/// defArea;  // tC per cell
				//26.
			}
		}else{   
			resDefor = cohort.deforest(singleCell.deforestA[Age],0); //MG: deforest set forest and get deforested biomass. Devide by refForShare to get per ha value
			defArea = resDefor.area;
			if (defArea>0 && realAreaO>0)
			{
				//      defBiomass = (resDefor.enSw + resDefor.enRw + resDefor.hv);///deforestA[Age]; // deforested aboveground biomass per ha // only stemwood !!!
				//      defBiomass = resDefor.bm/modTimeStep; // MG: divide by realAreaO ?????? tC/(ha year)
				//	  defBiomass = resDefor.bm;
				//           defBiomass = resDefor.bm / realAreaO / modTimeStep;
				//	       deforSW = resDefor.sw / realAreaO / modTimeStep;
				//		   deforRW = resDefor.rw / realAreaO / modTimeStep;
				defBiomass = resDefor.bm / defArea / modTimeStep; // tC/per ha (only stem!)
				deforSW = resDefor.sw / modTimeStep; // /defArea;// / defArea; // tC per cell
				deforRW = resDefor.rw / modTimeStep; // / defArea;// / defArea;  // tC per cell
				//27.
			}
		}
	}
	//if (defBiomass>0){
	//cout<<"asID=\t"<<asID<<"\tcountry=\t"<<Country<<"\tyear=\t"<<year<<"\tabBiomassO=\t"<<abBiomassO<<"\tdefBiomass=\t"<<defBiomass<<"\thv=\t"<<resDefor.hv<<endl;
	//}
	singleCell.OforestShare = OforestShare;
	singleCell.AforestShare = AforestShare;  
	singleCell.forestShare = OforestShare + AforestShare;
	double deforestHa = singleCell.deforestA[Age]*singleCell.LandAreaHa; 
	singleCell.deforestHaTot += deforestHa;
	double afforestHa = singleCell.forestAgeShare[Age]*singleCell.LandAreaHa;
	singleCell.afforestHaTot += afforestHa;
	//  singleCell.ProdLongA[Age] = defBiomass*it.FRACLONGPROD[byear]*(1-coeff.HarvLoos[byear])*deforestHa;
	singleCell.ProdLongA[Age] = defBiomass*it.FRACLONGPROD[byear]*(1-coeff.HarvLoos[byear]) * (1. - it.SLASHBURN[byear])*deforestHa;
	singleCell.ProdShortA[Age] = singleCell.ProdLongA[Age];
	//  singleCell.LitterA[Age] = it.CLITTERHA[byear]*deforestHa;
	singleCell.LitterA[Age] = (it.CLITTERHA[byear] + (BEF(int(it.SPECIESTYPE[byear])-1,defBiomass,it.FTIMBER[2000])-1) * defBiomass) * deforestHa; //MG:BEF: we assume that the non-stem aboveground biomass goes to the litter pool
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
			//	    if (singleCell.SOCA[i]<it.SOCHA[2000]*0.6*singleCell.deforestA[i]*singleCell.LandAreaHa) {
			//            singleCell.SOCA[i] = it.SOCHA[2000]*0.6*singleCell.deforestA[i]*singleCell.LandAreaHa; // SOCagr;
			//         }
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
			//cout<<"BEF[sp:"<<int(it.SPECIESTYPE[byear])<<"\tph:"<<abovePhCur*it.FTIMBER[2000]<<"]=\t"<<BEF(int(it.SPECIESTYPE[byear])-1,abovePhCur,it.FTIMBER[2000])<<endl;
			abovePhCur *= singleCell.forestAgeShare[ia]; // afforested (stem) biomass per ha * afforShare[ia]
			PlantPhytHaBmGr += abovePhCur;// * modTimeStep; // Stem biomass of planted forest, tC per ha 
			PlantPhytHaBmGrBef += abovePhCurBef; //Aboveground biomass of planted forest, tC per ha 
			//cout<<"year=\t"<<year<<"\tAge=\t"<<Age<<"\tia=\t"<<ia<<"\tfShare=\t"<<singleCell.forestAgeShare[ia]<<endl;
			if ((Age-ia) * modTimeStep >= 20)  //We track area and phytomass of new forest over 20 y.o.
			{plantPhytHaBmGrBef20 += abovePhCurBef;//afforested biomass per ha * afforShare (forest over 20 y.o.)
			plantArea20_rel += singleCell.forestAgeShare[ia];// relative area of forest over 20 y.o.
			}


			if (singleCell.LitterAffor[ia]<5*singleCell.forestAgeShare[ia]*singleCell.LandAreaHa) {
				//        float CurEmissionsLitterAfforCur = 0.95 * pow((1.-exp(-0.1*abovePhCur/singleCell.forestAgeShare[ia])),3)
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
	//  PlantPhytHaBlGr = PlantPhytHaBmGr * coeffBl;
	float CurPlantPhytBmGr = 0;
	if (AforestShare > 0)
	{
		PlantPhytHaBlGr = PlantPhytHaBmGrBef * coeffBl; //Belowground phytomass

		//  float CurPlantPhytBmGr=(PlantPhytHaBmGr-singleCell.prevPlantPhytHaBmGr) * singleCell.LandAreaHa;

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
	//////////////////// END Afforestation "negative" emissions block                


	float harvestDfM3Ha = (deforSW+deforRW)*(1.-it.SLASHBURN[byear]) * it.FTIMBER[byear]; // wood obtained at deforestation, m3/ha
	float deforWoodTotM3 = harvestDfM3Ha * singleCell.LandAreaHa; // wood obtained from harvesting (FM) in the cell, m3
	float harvestFmTotM3 = (harvWood*OforestShare+harvWoodNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa;// Harvested new and old wood (from FM) in the cell
	float harvestTotM3 = harvestFmTotM3 + deforWoodTotM3;// total wood (from FM and deforestation) in the cell


	harvestGrid.set(xi,yi,harvestTotM3);
	float harvestFcM3 = (sawnW + restW + fuelW + sawnWnew + restWnew)* it.FTIMBER[byear]; // current harvested wood from final cut in the cell, m3/ha																																			//29.
	float harvestThM3 =(sawnThW + restThW + sawnThWnew + restThWnew) * it.FTIMBER[byear]; // current harvested wood from thinning in the cell, m3/ha                                                                                                                                 //30.
	//float harvestTotPlusM3 = (harvWoodPlus*OforestShare + harvWoodNewPlus*singleCell.AforestSharePrev+(defBiomass-(1-resUse)*resDefor.hv) * (1. - it.SLASHBURN[byear])*singleCell.deforestA[Age]* it.FTIMBER[byear])*singleCell.LandAreaHa; //Total current harvested wood in the cell, m3, including residues (harvest losses)
	float harvestTotPlusM3 = (harvWoodPlus*OforestShare + harvWoodNewPlus*singleCell.AforestSharePrev+(defBiomass-(1-resUse)*(defBiomass-(deforSW+deforRW))) * (1. - it.SLASHBURN[byear])*singleCell.deforestA[Age]* it.FTIMBER[byear])*singleCell.LandAreaHa; //Total current harvested wood in the cell, m3, including residues (harvest losses)
	//float CAI_m3ha = (abBiomassO-singleCell.ObiomassPrev + harvWoodLost + sawnW + restW + sawnThW + restThW + harvL) * it.FTIMBER"][byear]; //Current annual increment m3/ha
	//float CAI_m3ha = (abBiomassO-singleCell.ObiomassPrev + sawnW + restW + sawnThW + restThW + harvL) * it.FTIMBER[byear]; //Current annual increment m3/ha
	//float CAI_m3ha = (abBiomassO-singleCell.ObiomassPrev*(1-singleCell.deforestA[Age]) + sawnW + restW + sawnThW + restThW + harvL + harvWoodLost) * it.FTIMBER[byear]; //Current annual increment m3/ha
	float CAI_m3ha = ((abBiomassO-singleCell.ObiomassPrev)/modTimeStep + sawnW + restW + sawnThW + restThW + + harvL + harvWoodLost) * it.FTIMBER[byear]; //Current annual increment m3/ha                                                                                                  //31.
	CAI_m3ha = CAI_m3ha < 0. ? 0. : CAI_m3ha; // m3/ha/year
	float phytNewIncr = 0.; if (AforestShare > 0) {phytNewIncr = CurPlantPhytBmGr/(AforestShare * singleCell.LandAreaHa);}
	float CAI_new_m3ha = (phytNewIncr + harvWoodLostNew + sawnWnew + sawnThWnew + restThWnew + restWnew + fuelWnew + fuelThWnew + harvLNew) * it.FTIMBER[2000];																																			    //32.
	CAI_new_m3ha = CAI_new_m3ha < 0. ? 0. : CAI_new_m3ha;

	float FMsink_ab = (abBiomassO * BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[byear]) 
		- singleCell.ObiomassPrev * BEF(int(it.SPECIESTYPE[byear])-1,singleCell.ObiomassPrev,it.FTIMBER[byear]))
		* 3.6666666667*unitConv/modTimeStep; // tCO2/ha/year
	float FMsink_bm = FMsink_ab * (1+coeffBl);  // tCO2/ha/year

	//if (Country==11 && year>2010){ cout<<"countryCalc = \t"<<Country<<"\t deforWoodTotM3 = \t"<<singleCell.deforWoodTotM3<<endl;       }

	int biomassRot = 1;
	if (thinningForest.get(xi,yi) < 0) {
		biomassRot = species[int(it.SPECIESTYPE[byear])-1]->gU(abBiomassO, maiForest.get(xi,yi))+1;  
	} else {
		biomassRot = species[int(it.SPECIESTYPE[byear])-1]->gUSdTab(abBiomassO, maiForest.get(xi,yi), thinningForest.get(xi,yi))+1;
	}

	if (fff.is_open()) {  
		if (countriesList.find(Country) != countriesList.end()){

			fff<<asID<<"\t"<<X<<"\t"<<Y<<"\t"<<year<<"\t"<<Country<<"\t"<<OforestShare<<"\t"<<singleCell.LandAreaHa<<"\t"<<AforestShare<<"\t"<<it.CABOVEHA[byear]<<"\t"<<abBiomassO<<"\t"<<PlantPhytHaBmGr<<"\t"<<maiForest.get(xi,yi)*it.FTIMBER[byear]<<"\t"<< rotMAI;
			fff<<"\t"<<(sawnThW + restThW) * it.FTIMBER[byear]<<"\t"<<(sawnW + restW)* it.FTIMBER[byear]<<"\t"<<(sawnThWnew + restThWnew) * it.FTIMBER[byear]<<"\t"<<(sawnWnew + restWnew)* it.FTIMBER[byear]<<"\t"<<harvestTotM3<<"\t"<<harvMAI<<"\t"<<singleCell.rotBiomass;			//33.

			fff<<"\t"<< (fval + coeff.PriceC[byear]* (((it.SOCHA[byear]*0.4/(8*decision.rotInter()))+5/(1.053*decision.rotInter()))/it.R[byear])) * Hurdle_opt[Country-1];
			fff<<"\t"<<fval<<"\t"<<aval<<"\t"<<rotationTimeCurr<<"\t"<<rotationForestNew.get(xi,yi)<<"\t"<<TimberPrice;
			fff<<"\t"<<(thinningForest.get(xi,yi));
			fff<<"\t"<<species[int(it.SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), optimMAI);
			fff<<"\t"<<realAreaO;
			fff<<"\t"<<singleCell.deforestA[Age];
			fff<<"\t"<<defArea;
			fff<<"\t"<<BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[byear]);
			fff<<"\t"<<(abBiomassO*BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[2000]));
			fff<< endl;
		}
	}


	//float TH_10;
	

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
	//fff<<"\t"<<               BEF(int(it.SPECIESTYPE[byear])-1,abBiomassO,it.FTIMBER[byear]);
	CountriesOforC_biom.inc(Country,year, OforBm*OforestShare*singleCell.LandAreaHa);   
	CountriesDeforCYear.inc(Country,year,EmissionsCur*3.6666666667*unitConv/modTimeStep);    // tC/year
	CountriesDeforCYear_ab.inc(Country,year, defBiomass*BEF(int(it.SPECIESTYPE[byear])-1,defBiomass,it.FTIMBER[2000])*deforestHa*3.6666666667*unitConv/modTimeStep);
	CountriesDeforCYear_bl.inc(Country,year,it.CBELOWHA[2000]*deforestHa*3.6666666667*unitConv/modTimeStep);      
	float EmDefBm = (defBiomass*BEF(int(it.SPECIESTYPE[byear])-1,defBiomass,it.FTIMBER[2000])+it.CBELOWHA[2000])*3.6666666667*unitConv; // CO2 lost from living biomass of old forests due to deforestation, mtCO2/hayear
	CountriesDeforCYear_biom.inc(Country,year,EmDefBm*deforestHa);       // tC/ha/year
	CountriesDeforCYear_dom.inc(Country,year, (EmissionsDeadBurnCur+EmissionsLitterCur)*3.6666666667*unitConv/modTimeStep);     
	CountriesDeforCYear_soil.inc(Country,year, EmissionsSOCCur*3.6666666667*unitConv/modTimeStep);     
	//--------- 

	CountriesWoodHarvestM3Year.inc(Country,year,harvestTotM3);
	CountriesWoodHarvestPlusM3Year.inc(Country,year,harvestTotPlusM3);    
	CountriesWoodHarvestFmM3Year.inc(Country,year,harvestFmTotM3);
	CountriesWoodHarvestDfM3Year.inc(Country,year,deforWoodTotM3);

	CountriesWoodHarvestFc_oldM3Year.inc(Country,year, (sawnW + restW) * it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa);																			//34.
	CountriesWoodHarvestTh_oldM3Year.inc(Country,year, (sawnThW + restThW) * it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa);																	//35.


	CountriesWoodHarvestFc_newM3Year.inc(Country,year,(sawnWnew + restWnew) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa);						    
	CountriesWoodHarvestTh_newM3Year.inc(Country,year,(sawnThWnew + restThWnew) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa);	


	if((year>1999)&&(Country==224)){
			sawlogs = ((sawnW + sawnThW )* it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa+ (sawnWnew + sawnThWnew)*it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa);
		    restlogs = ((restW + restThW)* it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa+ (restWnew + restThWnew)*it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa);
			
	//float TH_old = (sawnThW + restThW) * it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa;
	//float TH_new = (sawnThWnew + restThWnew) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa;
		//float FinalCutArea = ((sawnW + restW) * it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa+(sawnWnew + restWnew) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa)/harvestFcM3;
		//float ThinningArea = ((sawnThW + restThW) * it.FTIMBER[byear] * OforestShare*singleCell.LandAreaHa +(sawnThWnew + restThWnew) * it.FTIMBER[byear] * singleCell.AforestSharePrev*singleCell.LandAreaHa)/harvestThM3;
	if (outfile1.is_open()) {
			//    if ((year==2000)&&(Country==224)) {
		outfile1<<year<<"\t"<<X<<"\t"<<Y<<"\t"<<sawlogs<<"\t"<<restlogs<< endl;	
		//cout<<"X "<<X<<"	Y "<<Y<<"	year "<<year<<"		thinning "<<TH_old+TH_new<<endl;
	};
					
		//cout<<"X "<<X<<"	Y "<<Y<<"	year "<<year<<"		thinning "<<TH_old+TH_new<<endl; 
	};


	CountriesWoodLoosCYear.inc(Country,year,(harvWoodLost*OforestShare+harvWoodLostNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa*3.6666666667*unitConv);
	//     if (harvestTotM3>0) 
	CountriesHarvLossesYear.inc(Country,year,(harvL*OforestShare+harvLNew*singleCell.AforestSharePrev)* it.FTIMBER[byear]*singleCell.LandAreaHa);
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

	profit = ((sawnW + restW + fuelW + sawnThW + restThW + fuelThW + harvL*resUse)* it.FTIMBER[byear] * TimberPrice - decision.plantingCosts() * cohort.getArea(0)) ; // profit only from harvesting old forest per ha																	//38.
	CountriesProfit.inc(Country,year,profit);

	CountryRotation.inc(Country,year,rotationTimeCurr);

	//OT: Current rotation time extraction

	//if ((year==2000)&&(Country==224)&&(thinningForest.get(xi,yi)>0))
	/*	 if ((year==2000)&&(Country==224))
	{ 
	cout << "		x= "<<X<<"		y= "<<Y<<"OptRotTime= "<< rotMAI << endl;
	CountryRotationMng.inc(Country,year,rotationTimeCurr);
	//cout <<"year="<<year<<"		Country="<<Country<< "		x= "<<X<<"		y= "<<Y<<"		thinningForest= "<< thinningForest.get(xi,yi)<<"		rotationTimeCurr= "<<rotationTimeCurr<<endl;
	}// OT: end of extraction
	*/
	CountryregWoodHarvestM3Year.inc(Country,year,harvestTotM3);
	CountryregWoodHarvestFmM3Year.inc(Country,year,harvestFmTotM3);
	CountryregWoodHarvestDfM3Year.inc(Country,year,deforWoodTotM3);     


	char countrych[4];
	int ii = Country;
	int2str(ii,countrych);
	string countryprice = "re"+string(countrych)+"price0";

	//}
	CountryregWprod.setVal(Country,year,wprod[countryprice].g(year));

	CountryregMaxHarvest.inc(Country,year,deforWoodTotM3 + 
		(maiForest.get(xi,yi)*it.FTIMBER[byear]*(1.-coeff.HarvLoos[byear])*OforestShare+ harvWoodNew * singleCell.AforestSharePrev)*singleCell.LandAreaHa);


	//----------------------------------------------------------------------------
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
	singleCell.deforWoodTotM3 = deforWoodTotM3;
	singleCell.CAI = CAI_m3ha;
	if (year > 2000){singleCell.deforestShare = singleCell.deforestA[Age];}else{singleCell.deforestShare = 0.;}
	singleCell.afforestHaYear = afforestHa;
	singleCell.deforestHaYear = deforestHa;
	if (year > 2000){singleCell.deforPrev = singleCell.deforestA[Age];}
	singleCell.harvestTot = (harvestFcM3 + harvestThM3) + harvestDfM3Ha; 
	singleCell.harvestFcM3Ha = harvestFcM3;
	singleCell.harvestThM3Ha = harvestThM3;
	singleCell.oforestBm = OforBm;  
	singleCell.emissionsD_Bm = EmDefBm;  

	if (OforestShare>0){singleCell.emissionsD_S = EmissionsSOCCur/(OforestShare*singleCell.LandAreaHa)*3.6666666667*unitConv;}else{singleCell.emissionsD_S = 0;}
	if (AforestShare * singleCell.LandAreaHa>0){
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
