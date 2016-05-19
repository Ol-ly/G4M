void initManagedForest(dataDetStruct &data_all, g4m::incrementTab *species[8], 
					   g4m::ffipolm<double> *ffcov, 
					   g4m::ffipolm<double> *ffcoe, 
					   g4m::ffipolm<bool> *ffdov, 
					   g4m::ffipolm<bool> *ffdoe, 
					   datGlobal &dat_all,
					   ageStructVector &cohort_all, ageStructVector &newCohort_all,
					   griddata &maiForest, griddata &thinningForest,
					   griddata2<char> &rotationType, griddata2<char> &managedForest,
					   griddata &rotationForest, griddata &harvestSWGrid,griddata &harvestRWGrid, griddata2<float> &OforestShGrid) 

{
	double sawnWoodHarvest[NumberOfCountries+1];		// array of harvested sawnwood in each country, m3
	double restWoodHarvest[NumberOfCountries+1];		// array of harvested restwood in each country, m3
	double woodLost[NumberOfCountries+1];
	harvestWood objHarvestWood;					// object of harvestWood class

	for (int i=0; i<NumberOfCountries+1; i++){
		sawnWoodHarvest[i]=0.; 
		restWoodHarvest[i]=0.; 
		woodLost[i]=0.;

	} 

	int year =2000;

	//double sawnW = 0.;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
	//double restW = 0.;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
	//double sawnThW = 0.;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
	//double restThW = 0.;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
	//double bmH = 0.;        // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
	//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning
	//double harvRes = 0.;    // MG: usable harvest residues for the set (old) forest tC/ha


	//   Start Zero initialisation
	//------------------------------------------------------------------------------
	dataDetStruct::iterator iter = data_all.begin();
	while (iter != data_all.end()) {
		if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
			if (iter->PROTECT[2000] == 0) {
				int xi = (iter->x);
				int yi = (iter->y);
				int Country0 = (int)iter->COUNTRY[year];
				double LandAreaHa = iter->LANDAREA[year]*100; 
				double forestShare0 = iter->FOREST[byear];
				if (forestShare0<0) forestShare0=0;
				double maxfor = 1. - (iter->BUILTUP[2010] + iter->CROP[2010]);
				if (maxfor < 0.) maxfor = 0;
				double maxaffor1 = 1. - iter->GLOBIOM_RESERVED[2000];
				double maxaffor = 1.;
				if (maxfor < maxaffor1) {maxaffor = maxfor;}else{maxaffor=maxaffor1;}
				if (forestShare0>maxaffor) forestShare0=maxaffor;    

				OforestShGrid.set(xi,yi,forestShare0);
				OforestShGrid.update(); // populate the OforestShGridPrev with forestShare0 data     
				double forestArea0 = LandAreaHa * forestShare0;
				int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
				int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)    
				//no need for that? //       double harvWood=0.; //MG: harvestable wood, m3
				double abBiomassO = 0.;
				double MAI = 0.;
				double defIncome = 0.;

				int rotUnmanaged = 0;
				int rotMAI = 0;
				int rotMaxBm = 0;
				int rotMaxBmTh = 0;
				int rotHarvFin = 0;
				int rotHarvAve = 0;

				if (forestShare0 > 0){
					MAI = iter->MAIE[year]; // Max mean annual increment (tC/ha) of Existing forest (with uniform age structure and managed with rotation length maximazing MAI)
				}else{
					MAI = iter->MAIN[year]; // Max mean annual increment of New forest (with uniform age structure and managed with rotation length maximazing MAI)
				}
				if (MAI < 0.){MAI = 0.;}   
				maiForest.set(xi,yi,MAI);        
				float harvMAI = MAI*iter->FTIMBER[year]*(1-coeff.HarvLoos[year]);
				if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
					if (int(iter->SPECIESTYPE[year])<1) cout<<"speciestypeInit=\t"<<int(iter->SPECIESTYPE[year])<<endl;//<<"\tbiomassRot=\t"<<species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year],MAI)<<"\tbiomassRotTh=\t"<<species[int(iter->SPECIESTYPE[year])-1]->gUt(iter->CABOVEHA[year],MAI)<<endl;
					biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], MAI);        // rotation time to get current biomass (without thinning)          
					biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUt(iter->CABOVEHA[year], MAI);     // rotation time to get current biomass (with thinning)               //          if (iter->FOREST"].g(1990) >0) forFlag = 1.0;
				}
				if (maiForest.get(xi,yi)> 0){
					rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(MAI, optimMAI);
					rotMaxBm = species[int(iter->SPECIESTYPE[year])-1]->gTopt(MAI, optimMaxBm);                        
					rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(MAI, optimMaxBm);
				}

				dima decision(1990
					, iter->NPP
					, iter->SPOPDENS
					, iter->SAGRSUIT
					, iter->PRICEINDEX
					, coeff.PriceIndexE
					, iter->R
					, coeff.PriceC
					, coeff.PlantingCostsR
					, coeff.PriceLandMinR
					, coeff.PriceLandMaxR
					, coeff.MaxRotInter
					, coeff.MinRotInter
					, coeff.decRateL
					, coeff.decRateS
					, iter->FRACLONGPROD
					, coeff.baseline
					, iter->FTIMBER 
					, coeff.PriceTimberMaxR
					, coeff.PriceTimberMinR
					, coeff.FCuptake
					, iter->GDP
					, coeff.HarvLoos
					, forestShare0
					, wprice["regprice"]
				, wprice["regprice0"].g(2000)
					, rotMAI
					, harvMAI); 
				int Rotation = 0;
				double units = iter->FTIMBER[byear];

				if (iter->MANAGEDFLAG.g(1990) >0) {
					thinningForest.set(xi,yi,1.);
					Rotation = biomasRotTh+1;
					if (Rotation < rotMAI) Rotation = rotMAI;
					rotationType.set(xi,yi,11);


					double pDefIncome = iter->CABOVEHA[year] * 
						(decision.priceTimber() * iter->FTIMBER[year] * (1. -coeff.HarvLoos[year]));
					//Immediate Pay if deforested (Slash and Burn)
					double sDefIncome = iter->CABOVEHA[year] *
						(decision.priceTimber() * iter->FTIMBER[year]
					* (1. -coeff.HarvLoos[year]));
					defIncome = pDefIncome * (1. - iter->SLASHBURN[year])
						+ sDefIncome * iter->SLASHBURN[year];


					if (MAI > MAI_CountryUprotect[Country0-1]) {
						managedForest.set(xi,yi,3.);
					} else {
						if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome)) {
							managedForest.set(xi,yi,2.);
						} else {
							managedForest.set(xi,yi,1.);
						}
					}

				} else {
					thinningForest.set(xi,yi,-1.);
					Rotation = biomasRot+1;
					if (Rotation < rotMAI) Rotation = rotMAI;
					rotationType.set(xi,yi,10);

					double pDefIncome = iter->CABOVEHA[year] * 
						(decision.priceTimber() * iter->FTIMBER[year] * (1. -coeff.HarvLoos[year]));
					//Immediate Pay if deforested (Slash and Burn)
					double sDefIncome = iter->CABOVEHA[year] *
						(decision.priceTimber() * iter->FTIMBER[year]
					* (1. -coeff.HarvLoos[year]));
					defIncome = pDefIncome * (1. - iter->SLASHBURN[year])
						+ sDefIncome * iter->SLASHBURN[year];

					if (MAI > MAI_CountryUprotect[Country0-1]) {
						managedForest.set(xi,yi,0.);
					} else {
						if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome)) {
							managedForest.set(xi,yi,-1.);
						} else {
							managedForest.set(xi,yi,-2.);
						}
					}
				}
				rotationForest.set(xi,yi,Rotation);	

			}          // End if (iter->PROTECT[2000] == 0)
		}// Test only some regions
		iter++;
	}            // End while(iter != data_all.end())

	initLoop(0, data_all, species, ffcov,ffcoe, ffdov, ffdoe, cohort_all, newCohort_all, 
		dat_all, maiForest, thinningForest, rotationForest);
	cout << "end of first pass" << endl;
	//
	//
	//
	////******************************************************************************
	////**************************First Pass********************
	////           Init havWood
	////******************************************************************************

	iter = data_all.begin();
	while (iter != data_all.end()) {
		if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions        
			if (iter->PROTECT[2000]==0) {   // We consider only unprotected land

				double abBiomassO = 0.;  

				int asID = iter->asID;
				short int region = iter->COUNTRY[2000];
				char regionch[4];
				int2str(region,regionch);
				string regprice = "re"+string(regionch)+"price0";
				short int Country0 = iter->COUNTRY[2000]; 
				if (Country0==0 || Country0>244) cout<<"Country0=\t"<<Country0<<endl;
				int xi = (iter->x);
				int yi = (iter->y);
				int rotUnmanaged = 0;
				int rotMAI = 0;
				int rotMaxBm = 0;
				int rotMaxBmTh = 0;
				int rotHarvFin = 0;
				int rotHarvAve = 0;
				int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
				int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
				double LandAreaHa = iter->LANDAREA[year]*100;
				float forestShare0 = 0;
				if (iter->FOREST[year]+(iter->CROP[year])+(iter->BUILTUP[year])>1)
				{forestShare0 = (1-(iter->CROP[year]+iter->BUILTUP[year]));}
				else {forestShare0 = iter->FOREST[year];}
				double units = iter->FTIMBER[byear];
				double forestArea0 = LandAreaHa * forestShare0;
				float harvMAI = maiForest.get(xi,yi)*iter->FTIMBER[year]*(1-coeff.HarvLoos[year]);
				double harvSW=0; //harvested sawnwood in 1 Ha site, m3
				double harvRW=0; //harvested restwood in 1 Ha site, m3
				double defIncome = 0.;
				g4m::ageStruct cohortTmp = *cohort_all[asID];
				cohortTmp.setU(rotationForest.get(xi,yi));
				pair<g4m::ageStruct::v, g4m::ageStruct::v> res;
				int Rotation = 0;

				if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
					biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], maiForest.get(xi,yi));        // rotation time to get current biomass (without thinning)
					biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUSdTab(iter->CABOVEHA[year], maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
					if (thinningForest.get(xi,yi) > 0) {
						Rotation = biomasRotTh+1; 
					}  else {
						Rotation = biomasRot+1;
					}
				}
				if (maiForest.get(xi,yi)> 0){
					rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMAI);           
					rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
				}             
				res = cohortTmp.aging();         
				double realAreaO = cohortTmp.getArea();
				objHarvestWood.calcHarvM3(realAreaO, res, modTimeStep, units, biomasRot);
				dima decision(1990
					, iter->NPP
					, iter->SPOPDENS
					, iter->SAGRSUIT
					, iter->PRICEINDEX
					, coeff.PriceIndexE
					, iter->R
					, coeff.PriceC
					, coeff.PlantingCostsR
					, coeff.PriceLandMinR
					, coeff.PriceLandMaxR
					, coeff.MaxRotInter
					, coeff.MinRotInter
					, coeff.decRateL
					, coeff.decRateS
					, iter->FRACLONGPROD
					, coeff.baseline
					, iter->FTIMBER 
					, coeff.PriceTimberMaxR
					, coeff.PriceTimberMinR
					, coeff.FCuptake
					, iter->GDP
					, coeff.HarvLoos
					, forestShare0
					, wprice[regprice]
				, wprice[regprice].g(2000)
					, rotMAI
					, harvMAI); 
				abBiomassO = cohortTmp.getBm();
				double pDefIncome = abBiomassO * 
					(decision.priceTimber() * iter->FTIMBER[year] * (1. -coeff.HarvLoos[year]));
				//Immediate Pay if deforested (Slash and Burn)
				double sDefIncome = abBiomassO *
					(decision.priceTimber() * iter->FTIMBER[year] 
				* (1. -coeff.HarvLoos[year]));
				defIncome = pDefIncome * (1. - iter->SLASHBURN[year])
					+ sDefIncome * iter->SLASHBURN[year];		       
				harvestSWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWS);   
				harvestRWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWR);
				if (managedForest.get(xi,yi) > 0) {

					if (sawnWoodHarvest[Country0-1] < 0.95 * wprod_sawnwood[regprice].g(year)) {

						if (managedForest.get(xi,yi) >= 2) {
							rotationForest.set(xi,yi,rotMAI);		   
							cohortTmp.setU(rotMAI);
							res = cohortTmp.aging();
							double realAreaO = cohortTmp.getArea();
							objHarvestWood.calcHarvM3(realAreaO, res, modTimeStep, units, biomasRot);
							harvestSWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWS);   
							harvestRWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWR);
							sawnWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWS;						
							restWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWR;
						}            
					}
			

					else if (sawnWoodHarvest[Country0-1] > 1.05 * wprod_sawnwood[regprice].g(year)) {

						if (rotationForest.get(xi,yi) < rotMaxBmTh) {
							rotationForest.set(xi,yi, rotMaxBmTh);
							cohortTmp.setU(rotMaxBmTh);
							res = cohortTmp.aging();
							double realAreaO = cohortTmp.getArea();
							objHarvestWood.calcHarvM3(realAreaO, res, modTimeStep, units, biomasRot);
							harvestSWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWS);   
							harvestRWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWR);
							sawnWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWS;
							restWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWR;
						}
					} else { 
						double realAreaO = cohortTmp.getArea();
						objHarvestWood.calcHarvM3(realAreaO, res, modTimeStep, units, biomasRot);
						harvestSWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWS);   
						harvestRWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWR);
						sawnWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWS;
						restWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWR;
					} 

				}else { // unmanaged forests
					res = cohortTmp.aging();
					objHarvestWood.calcHarvM3(realAreaO, res, modTimeStep, units, biomasRot);
					harvestSWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWS);   
					harvestRWGrid.set(xi,yi, forestArea0 * objHarvestWood.harvestWR);
					sawnWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWS;
					restWoodHarvest[Country0-1] += forestArea0 * objHarvestWood.harvestWR;
				} 
			//	cout<<"111111"<<endl;

			} //else {iter++;}  // end for IF unprotected
		}   // Test only some regions
		iter++;
	} //end for WHILE
	//************************End of First Pass************************************


	//******************************************************************************
	//**************************Second Pass********************
	//******************************************************************************



	iter = data_all.begin();
	while (iter != data_all.end()) {
		if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions        
			if (iter->PROTECT[2000]==0) {   // We consider only unprotected land

				double HarvestSWTmp = 0;
				double HarvestRWTmp = 0;
				double newHarvestSWTmp = 0;
				double newHarvestRWTmp = 0;
				double abBiomassO = 0.;  
				int asID = iter->asID;
				short int region = iter->COUNTRY[2000];
				char regionch[4];
				int2str(region,regionch);
				string regprice = "re"+string(regionch)+"price0";
				short int Country0 = iter->COUNTRY[2000]; 
				int xi = (iter->x);
				int yi = (iter->y);

				int rotUnmanaged = 0;
				int rotMAI = 0;
				int rotMaxBm = 0;
				int rotMaxBmTh = 0;
				int rotHarvFin = 0;
				int rotHarvAve = 0;
				int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
				int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)

				double LandAreaHa = iter->LANDAREA[year]*100;
				double forestShare0 = 0;
				if (iter->FOREST[year]+(iter->CROP[year])+(iter->BUILTUP[year])>1)
				{forestShare0 = (1-(iter->CROP[year]+iter->BUILTUP[year]));}
				else {forestShare0 = iter->FOREST[year];}
				double forestArea0 = LandAreaHa * forestShare0;
				float harvMAI = maiForest.get(xi,yi)*iter->FTIMBER[year]*(1-coeff.HarvLoos[year]);
				double harvSW=0; //harvested sawnwood in 1 Ha site, m3
				double harvRW=0; //harvested restwood in 1 Ha site, m3
				double defIncome = 0.;
				double units = iter ->FTIMBER  [byear];

				g4m::ageStruct cohortTmp = *cohort_all[asID];
				cohortTmp.setU(rotationForest.get(xi,yi));

				pair<g4m::ageStruct::v, g4m::ageStruct::v> res;  // MG: results vector for the set (old) forest 
				int Rotation = 0;   
				if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
					biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], maiForest.get(xi,yi));        // rotation time to get current biomass (without thinning)
					biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUSdTab(iter->CABOVEHA[year], maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
					if (thinningForest.get(xi,yi) > 0) {
						Rotation = biomasRotTh+1; 
					}  else {
						Rotation = biomasRot+1;
					}
				}
				if (maiForest.get(xi,yi)> 0){
					rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMAI);                     
					rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
				}         
				if (Rotation < rotMAI) Rotation = rotMAI;           

				res = cohortTmp.aging();
				double realAreaO = cohortTmp.getArea();
				objHarvestWood.calcHarvM3(realAreaO, res, modTimeStep, units, biomasRot);
				dima decision(1990
					, iter->NPP
					, iter->SPOPDENS
					, iter->SAGRSUIT
					, iter->PRICEINDEX
					, coeff.PriceIndexE
					, iter->R
					, coeff.PriceC
					, coeff.PlantingCostsR
					, coeff.PriceLandMinR
					, coeff.PriceLandMaxR
					, coeff.MaxRotInter
					, coeff.MinRotInter
					, coeff.decRateL
					, coeff.decRateS
					, iter->FRACLONGPROD
					, coeff.baseline
					, iter->FTIMBER
					, coeff.PriceTimberMaxR
					, coeff.PriceTimberMinR
					, coeff.FCuptake
					, iter->GDP
					, coeff.HarvLoos
					, forestShare0
					, wprice[regprice]
				, wprice[regprice].g(2000)
					, rotMAI
					, harvMAI); 
				abBiomassO = cohortTmp.getBm();
				double pDefIncome = abBiomassO * 
					(decision.priceTimber() * iter->FTIMBER[year] * (1. -coeff.HarvLoos[year]));
				//Immediate Pay if deforested (Slash and Burn)
				double sDefIncome = abBiomassO *
					(decision.priceTimber() * iter->FTIMBER[year] 
				* (1. -coeff.HarvLoos[year]));
				defIncome = pDefIncome * (1. - iter->SLASHBURN[year])
					+ sDefIncome * iter->SLASHBURN[year];


				if (managedForest.get(xi,yi) > 0) {
					if (sawnWoodHarvest[Country0-1] < 0.95 * wprod_sawnwood[regprice].g(year)) {

						if ((Rotation > rotMAI) && (rotationForest.get(xi,yi) > Rotation)) { 
							HarvestSWTmp = harvestSWGrid.get(xi,yi);
							HarvestRWTmp = harvestRWGrid.get(xi,yi);
							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							res = cohortTmp.aging();
							double realAreaO = cohortTmp.getArea();	
							objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
							newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
							newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
							harvestSWGrid.set(xi,yi,newHarvestSWTmp);
							harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
							restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);

						} else  { 
							HarvestSWTmp = harvestSWGrid.get(xi,yi);
							HarvestRWTmp = harvestRWGrid.get(xi,yi);
							rotationForest.set(xi,yi,rotMAI);		   
							cohortTmp.setU(rotMAI);
							res = cohortTmp.aging();
							double realAreaO = cohortTmp.getArea();
							objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
							newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
							newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
							harvestSWGrid.set(xi,yi,newHarvestSWTmp);
							harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
							restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);
						}   

					} else if (sawnWoodHarvest[Country0-1] > 1.05 * wprod_sawnwood[regprice].g(year)) {

						if (rotationForest.get(xi,yi) < Rotation) { 
							HarvestSWTmp = harvestSWGrid.get(xi,yi);
							HarvestRWTmp = harvestRWGrid.get(xi,yi);
							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							res = cohortTmp.aging();
							double realAreaO = cohortTmp.getArea();
							objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
							newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
							newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
							harvestSWGrid.set(xi,yi,newHarvestSWTmp);
							harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
							restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);

						} else if (rotationForest.get(xi,yi) < rotMaxBmTh) {

							HarvestSWTmp = harvestSWGrid.get(xi,yi);
							HarvestRWTmp = harvestRWGrid.get(xi,yi);
							rotationForest.set(xi,yi, rotMaxBmTh);
							cohortTmp.setU(rotMaxBmTh);
							res = cohortTmp.aging();
							double realAreaO = cohortTmp.getArea();
							objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
							newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
							newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
							harvestSWGrid.set(xi,yi,newHarvestSWTmp);
							harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
							restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);
						}
					}
					if (restWoodHarvest[Country0-1] < 0.95 * wprod_restwood[regprice].g(year)) {
					if (rotationForest.get(xi,yi) != rotMAI)
					{ 
						HarvestSWTmp = harvestSWGrid.get(xi,yi);
						HarvestRWTmp = harvestRWGrid.get(xi,yi);
						rotationForest.set(xi,yi,rotMAI);		   
						cohortTmp.setU(rotMAI);
						res = cohortTmp.aging();          
						double realAreaO = cohortTmp.getArea();
						objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
						newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
						newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
						harvestSWGrid.set(xi,yi,newHarvestSWTmp);
						harvestRWGrid.set(xi,yi,newHarvestRWTmp);
						sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
						restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);
						if (sawnWoodHarvest[Country0-1] > 0.99 * wprod_sawnwood[regprice].g(year)){
							restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp+newHarvestSWTmp-HarvestSWTmp);
							cout<<"1 restWoodHarvest  "<<restWoodHarvest[224]<<endl;
						}
					}
				}
				}
			} //else {iter++;}  // end for IF unprotected
		}   // Test only some regions
		iter++;
	} //end for WHILE
	//************************End of Second Pass************************************


	//******************************************************************************
	//**************************Forth Pass********************
	//******************************************************************************
	//cout << "Start forth pass" << endl;

	iter = data_all.begin();

	while (iter != data_all.end())
	{if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
		if (iter->PROTECT[2000]==0)    // We consider only unprotected land
		{
			double HarvestSWTmp = 0;
			double HarvestRWTmp = 0;
			double newHarvestSWTmp = 0;
			double newHarvestRWTmp = 0;
			double abBiomassO = 0.;  	 
			int asID = iter->asID;
			short int region = iter->COUNTRY[2000];
			char regionch[4];
			int2str(region,regionch);
			string regprice = "re"+string(regionch)+"price0";
			short int Country0 = iter->COUNTRY[2000]; 
			int xi = (iter->x);
			int yi = (iter->y);

			int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
			int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
			double harvSW=0; //harvested sawnwood in 1 Ha site, m3
			double harvRW=0; //harvested restwood in 1 Ha site, m3
			double units =  iter->FTIMBER[byear];

			if (managedForest.get(xi,yi) > 0)
			{ 

				double LandAreaHa = iter->LANDAREA[year]*100;
				double forestShare0 = 0;
				if (iter->FOREST[year]+(iter->CROP[year])+(iter->BUILTUP[year])>1)
				{forestShare0 = (1-(iter->CROP[year]+iter->BUILTUP[year]));}
				else {forestShare0 = iter->FOREST[year];}

				double forestArea0 = LandAreaHa * forestShare0;

				int rotMAI = 0;
				int rotMaxBm = 0;
				int rotMaxBmTh = 0;
				int rotHarvFin = 0;
				int rotHarvAve = 0;

				g4m::ageStruct cohortTmp = *cohort_all[asID];
				cohortTmp.setU(rotationForest.get(xi,yi));
				pair<g4m::ageStruct::v, g4m::ageStruct::v> res; // MG: results vector for the set (old) forest    
				// End of FM initialisation

				int Rotation = 0;

				if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
					biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], maiForest.get(xi,yi));        // rotation time to get current biomass (without thinning)
					biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUSdTab(iter->CABOVEHA[year], maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
					if (thinningForest.get(xi,yi) > 0)
					{Rotation = biomasRotTh+1; 
					}else
					{Rotation = biomasRot+1;}
				}
				if (maiForest.get(xi,yi)> 0){
					rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
					rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
				}          
				if (Rotation < rotMAI) Rotation = rotMAI;   
				res = cohortTmp.aging();
				double realAreaO = cohortTmp.getArea();		  
				objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
				harvSW = objHarvestWood.harvestWS;
				harvRW = objHarvestWood.harvestWR;


				if (sawnWoodHarvest[Country0-1] < 0.95 * wprod_sawnwood[regprice].g(year)) 
				{
					if (rotationForest.get(xi,yi) != rotMAI)
					{ 
						HarvestSWTmp = harvestSWGrid.get(xi,yi);
						HarvestRWTmp = harvestRWGrid.get(xi,yi);
						rotationForest.set(xi,yi,rotMAI);		   
						cohortTmp.setU(rotMAI);
						res = cohortTmp.aging();          
						double realAreaO = cohortTmp.getArea();
						objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
						newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
						newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
						harvestSWGrid.set(xi,yi,newHarvestSWTmp);
						harvestRWGrid.set(xi,yi,newHarvestRWTmp);
						sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
						restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);

					} 		  
					  } 
				else if (sawnWoodHarvest[Country0-1] > 1.05 * wprod_sawnwood[regprice].g(year)) 
				{
					if (rotationForest.get(xi,yi) < rotMaxBmTh) {    
						HarvestSWTmp = harvestSWGrid.get(xi,yi);
						HarvestRWTmp = harvestRWGrid.get(xi,yi);
						rotationForest.set(xi,yi, rotMaxBmTh);
						cohortTmp.setU(rotMaxBmTh);
						res = cohortTmp.aging();
						newHarvestSWTmp = harvSW * forestArea0;
						newHarvestRWTmp = harvRW * forestArea0;
						harvestSWGrid.set(xi,yi,newHarvestSWTmp);
						harvestRWGrid.set(xi,yi,newHarvestRWTmp);
						sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
						restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);
						
					}
					
					}

				if (restWoodHarvest[Country0-1] < 0.95 * wprod_restwood[regprice].g(year)) {
					if (rotationForest.get(xi,yi) != rotMAI)
					{ 
						HarvestSWTmp = harvestSWGrid.get(xi,yi);
						HarvestRWTmp = harvestRWGrid.get(xi,yi);
						rotationForest.set(xi,yi,rotMAI);		   
						cohortTmp.setU(rotMAI);
						res = cohortTmp.aging();          
						double realAreaO = cohortTmp.getArea();
						objHarvestWood.calcHarvM3 (realAreaO, res, modTimeStep, units, biomasRot);
						newHarvestSWTmp = forestArea0 * objHarvestWood.harvestWS;
						newHarvestRWTmp = forestArea0 * objHarvestWood.harvestWR;
						harvestSWGrid.set(xi,yi,newHarvestSWTmp);
						harvestRWGrid.set(xi,yi,newHarvestRWTmp);
						sawnWoodHarvest[Country0-1] += (newHarvestSWTmp-HarvestSWTmp);
						restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp);
						if (sawnWoodHarvest[Country0-1] > 0.99 * wprod_sawnwood[regprice].g(year)){
							restWoodHarvest[Country0-1] += (newHarvestRWTmp-HarvestRWTmp+newHarvestSWTmp-HarvestSWTmp);
							cout<<"2 restWoodHarvest  "<<restWoodHarvest[224]<<endl;
						}
					} 
				}
			}
			cout<<"sawlogs>1.05"<< sawnWoodHarvest[223]<<endl;

		}// else{iter++;}  // end for IF unprotected
	} // Test only some regions   
	iter++;
	} //end for WHILE


	//************************End of Forth Pass************************************
	cout << "End of 4th pass"<<endl;


} //END of void
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void initLoop(int i, dataDetStruct &data_all, g4m::incrementTab *species[8], 
			  g4m::ffipolm<double> *ffcov, 
			  g4m::ffipolm<double> *ffcoe, 
			  g4m::ffipolm<bool> *ffdov, 
			  g4m::ffipolm<bool> *ffdoe, 
			  ageStructVector &cohort_all, 
			  ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
			  griddata &thinningForest, griddata &rotationForest) 

{

	int asID = 0;  // index in the ageStruct vector
	int kk=0;
	numAgeStruct = 0;

	map<int, vector<double> > ageStructData;
	{

		//cout<<"before reading file..."<<endl;
		//    ifstream fin("ageStructData.txt");
		//    string agestructFile = settings.inputPath + "\\" + "ageStructDataJRC_nc.txt";
		///////    string agestructFile = settings.inputPath + "\\" + "ageStructDataJRC_27June2011_nc.txt";    
		//    string agestructFile = settings.inputPath + "\\" + "ageStructDataJRC_30June2011_nc.txt";    
		string agestructFile = settings.inputPath + "\\" + "ageStructDataJRC_IT_26July2011_nc.txt";    //new age structure for Italy (email from Giacomo Grassi from 26 July 2011)
		ifstream fin(agestructFile.c_str());
		if (!fin.is_open()) {
			cout << "Cannot read file!!!! " << agestructFile << endl;
			system("pause");
			exit(0);
		}
		//cout<<"first line reading..."<<endl;    
		{string tmp; getline(fin, tmp);} //Skip first line
		int country;
		vector<double> ageShares;
		//    ageShares.resize(9);
		ageShares.resize(16);
		while(fin.good() && !fin.eof()) {
			fin >> country;
			double sum=0.;
			//      for(int i=0; i<9; ++i) {fin >> ageShares[i]; sum+=ageShares[i];}
			for(int i=0; i<16; ++i) {fin >> ageShares[i]; sum+=ageShares[i];}
			double add=0.;
			//      if(sum <=0.) {sum=0.; add=1./9.;} else {sum = 1./sum;}
			if(sum <=0.) {sum=0.; add=1./16.;} else {sum = 1./sum;}
			//      for(int i=0; i<9; ++i) {ageShares[i] = add + ageShares[i]*sum;}
			for(int i=0; i<16; ++i) {ageShares[i] = add + ageShares[i]*sum;}      
			ageStructData[country] = ageShares;
		}
	}



	//cout<<"start cell loop"<<endl;
	dataDetStruct::iterator iter = data_all.begin();
	while (iter != data_all.end()) {
		if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
			if (countriesList.find(iter->COUNTRY[2000]) != countriesList.end()) { // Test only some countries    
				if (iter->PROTECT[2000]==0) {
					//if (kk>2700) cout<<"we are here"<<endl;
					int xi = (iter->x);
					int yi = (iter->y);
					double X = (iter->x)*GridStepLon+GridStepLon/2-180;
					double Y = (iter->y)*GridStepLat+GridStepLat/2-90;
					int Country = (int)iter->COUNTRY[byear];
					double LandAreaHa = iter->LANDAREA[byear]*100;
					//        int aveRot = 0;                  // Average rotation time from statistics
					int PriceCi = PriceCiS[i];
					int iprice=i+1;
					coeff.PriceC.clear();
					coeff.PriceC.insert(0, PriceCi * iter->CORRUPTION[byear]);
					double harvWood=0.; //MG: harvestable wood, m3
					double harvWoodLost=0.; // wood lost in unmanaged forests, tC/ha
					double harvWoodNew=0.; //MG: harvestable wood, m3 (new forest)
					double harvWoodLostNew=0.; // wood lost in unmanaged forests, tC/ha (new forest)

					int beyear = (int)coeff.bYear;
					double forestShare = iter->FOREST[byear];
					if (forestShare<0) forestShare=0;
					//		double maxfor = 1. - (iter->BUILTUP[byear] + iter->CROP[byear]);
					double maxfor = 1. - (iter->BUILTUP[2010] + iter->CROP[2010]);
					if (maxfor < 0.) maxfor = 0;
					double maxaffor1 = 1. - iter->GLOBIOM_RESERVED[2000];
					double maxaffor = 1.;
					if (maxfor < maxaffor1) {maxaffor = maxfor;}else{maxaffor=maxaffor1;}
					if (forestShare>maxaffor) forestShare=maxaffor;
					double forestArea0 = LandAreaHa * forestShare;
					double OforestShare = forestShare;
					double AforestShare = 0.;              //MG: Actual forest share
					double refForShare = forestShare;      //forest share of ref. year	
					double OfsNoPay = forestShare;         //MG: Forest share (deforested) without payment
					double AfsNoPay = 0.;                  //MG: Forest share (afforested) without payment
					double fsNoPay = OfsNoPay + AfsNoPay;  //MG: Forest share without payment

					hlv.clear();
					//        hlv.insert(7, 0.8*(1-countryLosses[Country-1]));
					hlv.insert(25, 0.8*(1-countryLosses[Country-1]));
					hlv.insert(50, 0.8*(1-0.7*countryLosses[Country-1]));        

					hle.clear();
					//        hle.insert(7, (1-countryLosses[Country-1]));
					hle.insert(25, (1-countryLosses[Country-1]));        
					hle.insert(50, (1-0.7*countryLosses[Country-1]));   

					ffhlv.overwrite(hlv);
					ffhle.overwrite(hle);

					double forFlag = 0.; //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
					double sawnWpot = 0.;
					double restWpot = 0.;
					double sawnThWpot = 0.;
					double restThWpot = 0.;
					double bmHpot = 0.;
					double bmThPot = 0.;  
					double harvResPot = 0.;


					int optimumType = 3;
					int MAIRot = 1;  //MG: optimal rotation time (see above)
					//        int optRotUnmanaged = 0;
					int rotationTimeCurr = 1; 
					int Rotation = 1;
					int rotMaxBm = 1;
					int rotMaxBmTh = 1;
					int biomasRot=1;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
					int biomasRotTh=1;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
					//        if (refForShare >0 && iter->CABOVEHA"][year] > 0 && maiForest.get(xi,yi)> 0) {
					if (iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi)> 0) {
						biomasRot = species[int(iter->SPECIESTYPE[byear])-1]->gU(iter->CABOVEHA[byear], maiForest.get(xi,yi));
						if (biomasRot<1) biomasRot=1;
						//          biomasRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(iter->CABOVEHA"][byear], maiForest.get(xi,yi), 1); // with thinning
						biomasRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(iter->CABOVEHA[byear], maiForest.get(xi,yi), thinningForest.get(xi,yi)); // with thinning
						if(biomasRotTh<1) biomasRotTh=1;
						if (iter->FOREST[byear]>0) forFlag = 1.0;
					}
					if (maiForest.get(xi,yi)> 0){
						MAIRot = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
						rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm);                        
						rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
					}
					rotationTimeCurr = rotationForest.get(xi,yi);
					if (biomasRot < 1) biomasRot = 1;
					if (MAIRot < 1) MAIRot = 1;
					if (thinningForest.get(xi,yi) >0) {
						Rotation = biomasRotTh; 
					} else {
						Rotation = biomasRot;
					}

					g4m::ageStruct *cohort;
					cohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[byear])-1], &ffsws, &ffhlv, &ffhle, ffcov, ffcoe, ffdov, ffdoe, maiForest.get(xi,yi)
						, 0  //What stands the value of u for
						, 1  //Rotation time
						, 0, 0, 0
						//     ,1 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, MAIRot
						,0 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, MAIRot
						, 0  //reference of minrot
						, 0. //Flexibility of stocking degree
						, &ffsdMaxH, &ffsdMinH,1);

					double abBiomass0 = 0.;
					if (iter->POTVEG[2000]<10 && iter->FOREST[2000]>0 && iter->MAIE[2000]>0){
						map<int, vector<double> >::iterator PageStructData = ageStructData.find(Country);
						if((PageStructData != ageStructData.end()) && thinningForest.get(xi,yi)>0 && forFlag > 0) {
							cohort->createNormalForest(321, 0., 1.);

							double ageBreaks[] = {10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,999.};
							double ageSize[] = {11.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
							int ageGroups = 16;


							int oldestAgeGroup = 0;
							int   oag = 0;    
							int ageGroup = 0;

							for(int i=1; i<161; ++i) {
								while(i>ageBreaks[ageGroup] && ageGroup+1 < ageGroups) {++ageGroup;}
								cohort->setArea(i,PageStructData->second[ageGroup]/ageSize[ageGroup]);
								if (PageStructData->second[ageGroup]>0) {oldestAgeGroup = ageGroup;}
							}
							double area = cohort->getArea();
							double biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);

							if(area > 0. && biomass > 0.) {

								//Tune age structure for current cell


								if (biomass < 0.95*iter->CABOVEHA[byear]){                
									int young = 0;
									oag = oldestAgeGroup; // current oldest age group
									while (biomass < 0.95*iter->CABOVEHA[byear] && oag < 30){                          
										if (ageSize[young]>0 && oag>young){   
											for (int i=0;i<10;++i) {
												double areaTmp = cohort->getArea(i+young*10+1);  
												cohort->setArea(i+young*10+1,areaTmp/2);
												cohort->setArea(i+(oag+1)*10+1,areaTmp/2);
												biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);
											}
										}   
										young +=1;
										oag+=1;
									}    

								} else if (biomass > 2 * iter->CABOVEHA[byear]){        //v_24_11                
									int young = 0;
									oag = oldestAgeGroup;          
									while (biomass > 2 * iter->CABOVEHA[byear] && oag >2){              
										if (oag > young){
											if (ageSize[oag]>0 && ageSize[young]>0){      
												for (int i=0;i<ageSize[oag];++i) {
													double areaTmp_oag = cohort->getArea(i+oag*10+1);                    
													double areaTmp_young = cohort->getArea(i+young*10+1);
													cohort->setArea(i+oag*10+1,0.);
													cohort->setArea(i+young*10+1,(areaTmp_oag+areaTmp_young));
													biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);
												}
											}
										} else {
											if (ageSize[oag]>0 && ageSize[oag-1]>0){      
												for (int i=0;i<ageSize[oag];++i) {
													double areaTmp_oag = cohort->getArea(i+oag*10+1);                    
													double areaTmp_young = cohort->getArea(i+(oag-1)*10+1);                   
													cohort->setArea(i+oag*10+1,0.);
													cohort->setArea(i+(oag-1)*10+1,(areaTmp_oag+areaTmp_young));
													biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);
												}
											}       
										}    
										young +=1;
										oag-=1;
									}    
								}  
								double stockingDegree = iter->CABOVEHA[byear]/(cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000])/cohort->getArea());
								if(stockingDegree < 0.) {stockingDegree = 0.;}



								cohort->setStockingdegreeMin(stockingDegree*sdMinCoeff);
								cohort->setStockingdegreeMax(stockingDegree*sdMaxCoeff);
								thinningForest.set(xi,yi, stockingDegree);
								for(int i=0; i<321; ++i) {
									cohort->setBm(i, stockingDegree*cohort->getBm(i));
								}
								cohort->setU(321);
								cohort->aging();
							}


							Rotation = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(cohort->getBm()/cohort->getArea(), maiForest.get(xi,yi), thinningForest.get(xi,yi))+1;
							if (Rotation < MAIRot) Rotation = MAIRot;
							cohort->setU(Rotation);  
						} else if (forFlag > 0) {
							if (Rotation < MAIRot) Rotation = MAIRot;
							cohort->createNormalForest(Rotation, forFlag, thinningForest.get(xi,yi));
							cohort->setStockingdegreeMin(thinningForest.get(xi,yi)*sdMinCoeff);         
							cohort->setStockingdegreeMax(thinningForest.get(xi,yi)*sdMaxCoeff);         
							cohort->setU(Rotation);

						} else {cohort->createNormalForest(1, forFlag, 1);} // MG: create an existing forest with 0 area for consistency of the singleCell structure

						rotationForest.set(xi,yi,Rotation);   
						double abBiomass0 = cohort->getBm(); // modelled biomass at time 0, tC/ha
					} //End for FOREST>0

					g4m::ageStruct *newCohort;
					newCohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[byear])-1], &ffsws, &ffhlv, &ffhle, ffcov, ffcoe, ffdov, ffdoe, maiForest.get(xi,yi)
						, 0  //What stands the value of u for
						, Rotation  //Rotation time
						, 0, 0, 0
						//     ,1 , thinningForest.get(xi,yi)*sdMaxCoef, thinningForest.get(xi,yi)*sdMinCoeff8, 30, MAIRot
						,0 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, MAIRot
						, 0  //reference of minrot
						, 0. //Flexibility of stocking degree
						, &ffsdMaxH, &ffsdMinH,1);

					newCohort->createNormalForest(Rotation,0.,thinningForest.get(xi,yi));

					dat singleCell;
					//************************************
					//**** Element for preparing output for GUI 
					//************************************
					singleCell.simUnit = sMap.getSimu(X,Y);
					singleCell.Rotation = rotationForest.get(xi,yi);
					singleCell.LandAreaHa = LandAreaHa;
					singleCell.potHarvest = 0;
					singleCell.forestShare = forestShare;
					singleCell.OforestShare = OforestShare;
					singleCell.AforestShare = AforestShare;
					singleCell.prevOForShare = refForShare; //MG: Old Forest share in the previous reporting year
					singleCell.prevOForShareRP = OforestShare; //MG: Old Forest share in the previous reporting year
					singleCell.prevAForShareRP = 0.; //MG: New (afforested) Forest share in the previous reporting year
					singleCell.AforestSharePrev = 0.;
					singleCell.savedCarbonPrev = 0.;
					singleCell.gainedCarbonPrev = 0.;
					singleCell.EmissionsTotPrev = 0.;
					singleCell.EmissionsAfforPrev = 0.;
					singleCell.prevPlantPhytHaBmGr = 0.;
					singleCell.prevPlantPhytHaBlGr = 0.;
					singleCell.deforestHaTot=0.;
					singleCell.afforestHaTot=0.;
					singleCell.EmissionsProduct= 0.;  
					singleCell.EmissionsLitter = 0.;  
					singleCell.EmissionsSOC = 0.;      
					singleCell.EmissionsSlashBurn = 0.;
					singleCell.EmissionsDeadBurn = 0.;
					singleCell.EmissionsCRootBurn = 0.;    
					singleCell.EmissionsTot = 0.;     
					singleCell.EmLitterAffor =0.;
					singleCell.EmSOCAffor = 0.; 
					singleCell.EmissionsAffor = 0.;
					singleCell.ObiomassPrev = abBiomass0;
					singleCell.Obiomass0 = abBiomass0;                   // Modelled biomass at time 0, tC/ha
					singleCell.rotBiomass = Rotation;
					singleCell.SD = thinningForest.get(xi,yi);
					for (int i=0; i<110; ++i) {
						singleCell.forestAgeShare[i] = 0.;
						singleCell.BDeadA[i]=0.;
						singleCell.LitterA[i]=0.;
						singleCell.SOCA[i]=0;
						singleCell.ProdLongA[i]=0.;
						singleCell.ProdShortA[i]=0.;
						singleCell.deforestA[i]=0;
						singleCell.FineRootA[i]=0;
						singleCell.LitterAffor[i]=0.;
						singleCell.SOCaffor[i]=0.;
					}
					singleCell.deforWoodTotM3 = 0.;// stem wood obtained from deforestation averaged for last 5 years  
					singleCell.DeforWoodM3.first = 0.;
					singleCell.DeforWoodM3.second = 0.;
					singleCell.deforestShare = 0.;//                 
					singleCell.afforestShare = 0.;//                         
					for (int i=0;i<5;i++) {singleCell.deforWoodArray[i] = 0.;
					singleCell.deforAreaArray[i] = 0.;}
					singleCell.prevReportYear = coeff.bYear;
					singleCell.ireportYear=0;

					singleCell.deforPrev = iter->FORLOSS[2000];
					singleCell.road = iter->ROAD[2000];
					singleCell.deforRateCoeffCell = 1.;
					// saving results to initial vectors
					cohort_all.push_back(cohort);
					newCohort_all.push_back(newCohort);
					dat_all.push_back(singleCell);
					iter->asID = asID;
					asID++;
					numAgeStruct++;
				}     // End for IF unprotected 
			}      // countries
		}      // Test only some regions 
		iter++;
		kk++;
	}       // End of WHILE cells loop
}
