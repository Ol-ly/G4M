#ifndef fm_cpol_cpp
#define fm_cpol_cpp

void fm_cpol(dataDetStruct &data_all, g4m::incrementTab *species[8], ageStructVector &cohort_all, 
			 ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
			 griddata &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
			 griddata &rotationForestNew, griddata &thinningForestNew, griddata2<char> &manageChForest,
			 griddata2<char> &rotationType, griddata &harvestSWGrid, griddata &harvestRWGrid, int year, griddata2<char> &unmanaged, 
			 double fm_hurdle, double &maxDiff, double priceC,
			 set<int> &useChange, griddata &thinningForestO,griddata2<int> &rotationForestO,griddata2<int> &rotationForestNewO) 


{
	double sawnWoodHarvest[NumberOfCountries+1];		// array of harvested sawnwood in each country, m3
	double restWoodHarvest[NumberOfCountries+1];		// array of harvested restwood in each country, m3
	sawnWoodHarvest[0]=0.;
	restWoodHarvest[0]=0.;
	for (int i=1; i<=NumberOfCountries; i++){
		sawnWoodHarvest[i]=0.;
		restWoodHarvest[i]=0.;
		FMs[i]=0.;
	}
	harvestWood objHarvestWood;
	dataDetStruct::iterator iter = data_all.begin();

	while (iter != data_all.end()) {
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions         
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) { // Test only some countries         
				int asID = iter->asID;
				if (iter->PROTECT[2000]==0) {
					int xi = (iter->x);
					int yi = (iter->y);
					float sawnW = 0.;
					float restW = 0.;
					float sawnThW = 0.;
					float restThW = 0.;
					float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
					float sawnWnew = 0.;
					float restWnew = 0.;
					float sawnThWnew = 0.;
					float restThWnew = 0.;
					float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

					float harvestSWTmp = 0;
					float harvestRWTmp = 0;
					float newHarvestSWTmp = 0;
					float newHarvestRWTmp = 0;

					if (thinningForest.get(xi,yi) >0)  
					{
						g4m::ageStruct cohortTmp = *cohort_all[asID];
						g4m::ageStruct cohortTmpNew = *newCohort_all[asID];
						pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

						pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

						float realAreaO = cohortTmp.getArea();
						float realAreaN = cohortTmpNew.getArea();
						objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
						float harvestSO = objHarvestWood.harvestWS;
						float harvestRO = objHarvestWood.harvestWR;
						objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
						float harvestSN = objHarvestWood.harvestWS;
						float harvestRN = objHarvestWood.harvestWR;
						newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
							harvestSN * dat_all[asID].AforestShare) * 
							dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
						newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
							harvestRN * dat_all[asID].AforestShare) * 
							dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;

						sawnWoodHarvest[country] += (newHarvestSWTmp);
						restWoodHarvest[country] += (newHarvestRWTmp);
						harvestSWGrid.set(xi,yi,newHarvestSWTmp);   
						harvestRWGrid.set(xi,yi,newHarvestRWTmp); 
					}else{
						sawnWoodHarvest[country] += dat_all[asID].DeforWoodM3.first;
						restWoodHarvest[country] += dat_all[asID].DeforWoodM3.second;
						harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first); 
						harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
					} 
				}
			}
		}
		iter++;
	}

	if ((fm_hurdle == 1)&&(toAdjust.size()>1)){

		iter = data_all.begin();
		while (iter != data_all.end())
		{
			short int region = iter->POLESREG[2000];
			if (regions.find(region) != regions.end()) { // Test only some regions         
				short int country = iter->COUNTRY[2000];
				if (countriesNoFmcpol.find(country)==countriesNoFmcpol.end()){   // Don't do initial disturbance for countries which cannot produce demanded amount of wood due to lack of forest resourses or very high deforestation
					if (toAdjust.find(country) != toAdjust.end()) {         
						if (iter->PROTECT[2000] == 0)
						{
							int asID = iter->asID;
							char regionch[3];
							int2str(region,regionch);
							string regprice = "re"+string(regionch)+"price0";
							string regprice0 = "re"+string(regionch)+"price0";
							int xi = (iter->x);
							int yi = (iter->y);                    

							if (managedForest.get(xi,yi)>0 && maiForest.get(xi,yi)>0 && dat_all[asID].OforestShare > 0 )
							{

								float rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
								float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];
								float harvestSWTmp = harvestSWGrid.get(xi,yi);  
								float harvestRWTmp = harvestRWGrid.get(xi,yi);
								float forestShare = dat_all[asID].OforestShare;          
								float NPVcw = 0.;   	


								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[byear]);

								{
									g4m::ageStruct cohortTmp = *cohort_all[asID];
									cohortTmp.setU(rotMaxBmTh);
									NPVcw = npv_calc(*iter, cohortTmp, maiV,year,rotMaxBmTh,regprice,regprice0,forestShare,1);
								}                         

								if (NPVcw > NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID])                             
								{
									if (rotationForest.get(xi,yi) < rotMaxBmTh && NPVcw >= 0)
									{
										float rotationTmpD = 1;
										int rotationTmp = 1;

										if (NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID] > 0.){

											rotationTmpD = rotationForest.get(xi,yi) + 10 * (1+((NPVcw-NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID])/(NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID])));


										}else if (NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID] < 0.){
											rotationTmpD = rotationForest.get(xi,yi) + 10 * (1-((NPVcw-NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID])/(NPVbau[(year-refYear-modTimeStep)/modTimeStep][asID])));

										}else{rotationTmpD = rotationForest.get(xi,yi) + 15;}

										if (rotationTmpD>rotMaxBmTh) {rotationTmp = rotMaxBmTh;}else{rotationTmp=int(rotationTmpD);}
										rotationForest.set(xi,yi,rotationTmp);
										cohort_all[asID]->setU(rotationTmp);
										g4m::ageStruct cohortTmp = *cohort_all[asID];
										cohortTmp.setU(rotationTmp);

										g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                 
										int newForAge = newCohort_all[asID]->getActiveAge();
										int rotationNewTmp = rotationForestNew.get(xi,yi);
										if (newForAge > rotationNewTmp){
											newCohort_all[asID]->setU(rotationTmp);
											cohortTmpNew.setU(rotationTmp);
											rotationForestNew.set(xi,yi,rotationTmp);
										}

										pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

										pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

										float realAreaO = cohortTmp.getArea();
										float realAreaN = cohortTmpNew.getArea();

										objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
										float harvestSO = objHarvestWood.harvestWS;
										float harvestRO = objHarvestWood.harvestWR;
										objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
										float harvestSN = objHarvestWood.harvestWS;
										float harvestRN = objHarvestWood.harvestWR;

										float newHarvestSWTmp_cw = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
											harvestSN * dat_all[asID].AforestShare) * 
											dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
										float newHarvestRWTmp_cw = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
											harvestRN * dat_all[asID].AforestShare) * 
											dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;

										sawnWoodHarvest[country] += (newHarvestSWTmp_cw - harvestSWTmp);
										restWoodHarvest[country] += (newHarvestRWTmp_cw - harvestRWTmp);

										harvestSWGrid.set(xi,yi,newHarvestSWTmp_cw);  
										harvestRWGrid.set(xi,yi,newHarvestRWTmp_cw);

										NPVcGrid.set(xi,yi,NPVcw);
									} else 
									{     
										NPVcGrid.set(xi,yi,NPVcw);     
									}

								}else if (rotationForest.get(xi,yi) < rotMaxBmTh)
								{    
									int rotationTmp = rotationForest.get(xi,yi);
									if (coeff.PriceC[2000]<=50){
										rotationTmp += 0.1*coeff.PriceC[2000];
									}else { rotationTmp += 5;  }

									if (rotationTmp>rotMaxBmTh) {rotationTmp = rotMaxBmTh;}

									{   
										rotationForest.set(xi,yi,rotationTmp);
										cohort_all[asID]->setU(rotationTmp);
										g4m::ageStruct cohortTmp = *cohort_all[asID];
										cohortTmp.setU(rotationTmp);

										g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                 
										int newForAge = newCohort_all[asID]->getActiveAge();
										int rotationNewTmp = rotationForestNew.get(xi,yi);
										if (newForAge > rotationNewTmp){
											if (coeff.PriceC[2000]<=50){
												rotationNewTmp += 0.1*coeff.PriceC[2000];
											}else { rotationNewTmp += 5;  }                     
											if (rotationNewTmp>rotMaxBmTh) {rotationNewTmp = rotMaxBmTh;}
											newCohort_all[asID]->setU(rotationNewTmp);
											cohortTmpNew.setU(rotationNewTmp);
											rotationForestNew.set(xi,yi,rotationNewTmp);
										}

										pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

										pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

										float realAreaO = cohortTmp.getArea();
										float realAreaN = cohortTmpNew.getArea();

										objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
										float harvestSO = objHarvestWood.harvestWS;
										float harvestRO = objHarvestWood.harvestWR;
										objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
										float harvestSN = objHarvestWood.harvestWS;
										float harvestRN = objHarvestWood.harvestWR;

										float newHarvestSWTmp_cw = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
											harvestSN * dat_all[asID].AforestShare) * 
											dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
										float newHarvestRWTmp_cw = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
											harvestRN * dat_all[asID].AforestShare) * 
											dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;

										sawnWoodHarvest[country] += (newHarvestSWTmp_cw - harvestSWTmp);
										restWoodHarvest[country] += (newHarvestRWTmp_cw - harvestRWTmp);

										harvestSWGrid.set(xi,yi,newHarvestSWTmp_cw);  
										harvestRWGrid.set(xi,yi,newHarvestRWTmp_cw);                                                                                                                                            

										float NPV50 = npv_calc(*iter, cohortTmp, maiV,year,rotationTmp,regprice,regprice0,forestShare,0);
										NPVcGrid.set(xi,yi,NPV50);
									}                                                                                                                                                                                                                                                                       
								} 
								else 
								{

									NPVcGrid.set(xi,yi,NPVcw);     
								}
							} // if managed
							else if (managedForest.get(xi,yi)<=0 && maiForest.get(xi,yi)>0 && dat_all[asID].OforestShare > 0 )
							{
								float rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm);
								float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

								float forestShare = dat_all[asID].OforestShare;          
								float NPVc50 = 0.;   	

								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[byear]);

								short int used = 0;                                    
								g4m::ageStruct cohortTmp = *cohort_all[asID];

								NPVc50 = npv_calc(*iter, cohortTmp, maiV,year,rotMaxBm,regprice,regprice0,forestShare,0);
								NPVcGrid.set(xi,yi,NPVc50);

							}// else managedForest.get(xi,yi)<=0


						}   // if protected
					}    // if NoFMCpol
				}    // if toAdjust
			}     // if regions 



			iter++;
		}       // while
	}        // if refYear +1 


	//------------------------------------------------------------------------------
	// 
	// -------Zero Adjust thinning if population density changed --------------------
	//
	//------------------------------------------------------------------------------


	if ((year > 2000) && ((year) % 10 == 0) && (fm_hurdle == 1))
	{
		//cout<<"Start Zero Adjust thinning if population density changed"<<endl;
		iter = data_all.begin();

		while (iter != data_all.end())
		{
			short int region = iter->POLESREG[2000];
			if (regions.find(region) != regions.end()) { // Test only some regions         
				short int country = iter->COUNTRY[2000];
				if (toAdjust.find(country) != toAdjust.end()) { // do only for some countries         
					if (iter->PROTECT[2000] == 0)
					{
						int asID = iter->asID;
						char regionch[3];
						int2str(region,regionch);
						string regprice = "re"+string(regionch)+"price0";
						string regprice0 = "re"+string(regionch)+"price0";
						char countrych[4];
						int2str(country,countrych);
						string countryprice = "re"+string(countrych)+"price0";
						country = iter->COUNTRY[2000]; 
						int xi = (iter->x);
						int yi = (iter->y);


						if (sawnWoodHarvest[country] < 0.88 * wprod_sawnwood[countryprice].g(year))   
						{
							if (managedForest.get(xi,yi)<=0)
							{
								if ((iter->POPDENS.g(year) >0) && (iter->GDP.g(year) > 0))
								{
									//float sawnW = 0.;
									//float restW = 0.;
									//float sawnThW = 0.;
									//float restThW = 0.;
									//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
									//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
									//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
									//float sawnWnew = 0.;
									//float restWnew = 0.;
									//float sawnThWnew = 0.;
									//float restThWnew = 0.;
									//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
									//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
									//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha



									//      float harvestTmp = 0.;
									float newHarvestSWTmp = 0.;
									float newHarvestRWTmp = 0.;

									int biomassRot = 0;
									int rotMAI = 0;
									int rotMaxBmTh = 0;
									int Rotation = 0;
									float SD = 1.5;


									float countryHarvestSWTmp = sawnWoodHarvest[country];
									float countryHarvestRWTmp = restWoodHarvest[country];
									float NPVtmp = 0.;
									float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

									g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
									g4m::ageStruct cohortTmp = *cohort_all[asID];
									float harvestSWTmp = harvestSWGrid.get(xi,yi); 
									float harvestRWTmp = harvestRWGrid.get(xi,yi);

									if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
									{
										biomassRot = species[int(iter->SPECIESTYPE[byear])-1]->gU(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
										rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
										rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
									}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
									{
										rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
										rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
										biomassRot = rotMAI;
									}     
									int rotationForestTmp = rotationForest.get(xi,yi);
									int managedForestTmp = managedForest.get(xi,yi);
									float thinningForestTmp = thinningForest.get(xi,yi);
									int rotationForestNewTmp = rotationForestNew.get(xi,yi);
									float thinningForestNewTmp = thinningForestNew.get(xi,yi);

									//---------------------------------------------------
									int newForAge = newCohort_all[asID]->getActiveAge();
									if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
									{

										if (useChange.find(asID)!=useChange.end()){
											managedForest.set(xi,yi,managedForestTmp+3);
											Rotation = rotationForestO.get(xi,yi);
											SD = thinningForestO.get(xi,yi);

										}else if (managedForest.get(xi,yi) == 0)
										{
											managedForest.set(xi,yi,3);
											Rotation = rotMAI;
											rotationType.set(xi,yi,1);
										}else 
										{          managedForest.set(xi,yi,2);
										Rotation = rotMaxBmTh;
										rotationType.set(xi,yi,3);
										}

										rotationForest.set(xi,yi,Rotation);	
										thinningForest.set(xi,yi,SD);
										cohortTmp.setU(Rotation);
										cohortTmp.setStockingdegreeMin(SD*sdMinCoeff);  cohortTmp.setStockingdegreeMax(SD*sdMaxCoeff);
										cohort_all[asID]->setU(Rotation);
										cohort_all[asID]->setStockingdegreeMin(SD*sdMinCoeff);  cohort_all[asID]->setStockingdegreeMax(SD*sdMaxCoeff);


										rotationForestNew.set(xi,yi,Rotation);	
										thinningForestNew.set(xi,yi,SD);	
										cohortTmpNew.setU(Rotation);
										cohortTmpNew.setStockingdegreeMin(SD*sdMinCoeff);  cohortTmpNew.setStockingdegreeMax(SD*sdMaxCoeff);        
										newCohort_all[asID]->setU(Rotation);
										newCohort_all[asID]->setStockingdegreeMin(SD*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(SD*sdMaxCoeff);  

									}else // New forest age <= rotation -> don't change FM for the new forest
									{

										if (useChange.find(asID)!=useChange.end()){
											managedForest.set(xi,yi,managedForestTmp+3);
											Rotation = rotationForestO.get(xi,yi);
											SD = thinningForestO.get(xi,yi);
											rotationForestNew.set(xi,yi,rotationForestNewO.get(xi,yi));	
											thinningForestNew.set(xi,yi,SD);	
											newCohort_all[asID]->setU(rotationForestNewO.get(xi,yi));
											newCohort_all[asID]->setStockingdegreeMin(SD*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(SD*sdMaxCoeff);;  

										}else if (managedForest.get(xi,yi) == 0)
										{
											managedForest.set(xi,yi,3);
											Rotation = rotMAI;
											rotationType.set(xi,yi,1);
										}else 
										{     managedForest.set(xi,yi,2);
										Rotation = rotMaxBmTh;
										rotationType.set(xi,yi,3);
										}

										rotationForest.set(xi,yi,Rotation);	
										thinningForest.set(xi,yi,SD);	
										cohortTmp.setU(Rotation);
										cohortTmp.setStockingdegreeMin(SD*sdMinCoeff);  cohortTmp.setStockingdegreeMax(SD*sdMaxCoeff);
										cohort_all[asID]->setU(Rotation);
										cohort_all[asID]->setStockingdegreeMin(SD*sdMinCoeff);  cohort_all[asID]->setStockingdegreeMax(SD*sdMaxCoeff);
									} // End    else // New forest age < rotation -> don't change FM for the new forest
									pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

									pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

									float realAreaO = cohortTmp.getArea();
									float realAreaN = cohortTmpNew.getArea();
									objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
									float harvestSO = objHarvestWood.harvestWS;
									float harvestRO = objHarvestWood.harvestWR;
									objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
									float harvestSN = objHarvestWood.harvestWS;
									float harvestRN = objHarvestWood.harvestWR;
									newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestSN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
									newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestRN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
									countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
									countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);


									if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {  
										coeff.PriceC.clear();
										coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);
										short int used = 0;
										if (thinningForest.get(xi,yi)>0){used = 1;}
										NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
									}else{NPVtmp = 1e50;};

									if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 

										&& NPVtmp >= NPVcGrid.get(xi,yi)
										){
											sawnWoodHarvest[country] = countryHarvestSWTmp;   
											restWoodHarvest[country] = countryHarvestRWTmp;
											harvestSWGrid.set(xi,yi,newHarvestSWTmp);
											harvestRWGrid.set(xi,yi,newHarvestSWTmp);
											if (useChange.find(asID) == useChange.end()){manageChForest.set(xi,yi,1);}
									}else{ // return old values
										rotationForest.set(xi,yi,rotationForestTmp);	
										thinningForest.set(xi,yi,thinningForestTmp);
										managedForest.set(xi,yi,managedForestTmp);	                   
										cohort_all[asID]->setU(rotationForestTmp);
										cohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   cohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);
										rotationForestNew.set(xi,yi,rotationForestNewTmp);	
										thinningForestNew.set(xi,yi,thinningForestNewTmp);
										//                managedForestNew.set(xi,yi,managedForestNewTmp)	                   
										newCohort_all[asID]->setU(rotationForestNewTmp);
										newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);                 
									}

								}  //End  if ((data["POPDENS"].g(1990) >0) && (data["GDP"].g(1990) > 0))

							} // End  if (managedForest(xi,yi)<=0)

						}else if (sawnWoodHarvest[country] > 1.12 * wprod_sawnwood[countryprice].g(year))  
						{
							if (managedForest.get(xi,yi)>0) 
							{
								if ((iter->POPDENS.g(year) == 0) && (iter->GDP.g(year) == 0))
								{

									/*float sawnW = 0.;
									float restW = 0.;
									float sawnThW = 0.;
									float restThW = 0.;
									float sawnWnew = 0.;
									float restWnew = 0.;
									float sawnThWnew = 0.;
									float restThWnew = 0.;*/
									//      float harvWoodTmp = 0.;
									//      float harvestTmp = 0.;
									float newHarvestSWTmp = 0.;
									float newHarvestRWTmp = 0.;
									int biomassRotTh = 0;
									int rotMAI = 0;
									int rotMaxBm = 0;
									int Rotation = 0;

									//      int managedForestTmp = 0;
									//      double thinningForestTmp = 0;
									//      int rotationForestTmp = 0;
									//      double thinningForestNewTmp = 0;
									//      int managedForestNewTmp = 0;
									//      int rotationForestNewTmp = 0;

									float countryHarvestSWTmp = sawnWoodHarvest[country];
									float countryHarvestRWTmp = restWoodHarvest[country];

									float NPVtmp = 0.;
									float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

									g4m::ageStruct cohortTmp = *cohort_all[asID];
									float harvestSWTmp = harvestSWGrid.get(xi,yi); 
									float harvestRWTmp = harvestRWGrid.get(xi,yi);

									if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
									{
										biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (without thinning)  
										rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
										rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm);  
									}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
									{
										rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
										rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm); 
										biomassRotTh = rotMAI; 
									}     

									int rotationForestTmp = rotationForest.get(xi,yi);
									int managedForestTmp = managedForest.get(xi,yi);
									double thinningForestTmp = thinningForest.get(xi,yi);
									int rotationForestNewTmp = rotationForestNew.get(xi,yi);
									int thinningForestNewTmp = thinningForestNew.get(xi,yi);
									//---------------------------------------------------

									if (managedForest.get(xi,yi) == 2)
									{
										managedForest.set(xi,yi,-1);
										Rotation = rotMaxBm;
										rotationType.set(xi,yi,1);
									}else 
									{     managedForest.set(xi,yi,-2);
									Rotation = rotMaxBm;
									rotationType.set(xi,yi,3);
									}

									rotationForest.set(xi,yi,Rotation);	
									thinningForest.set(xi,yi,-1.);
									//                manageChForest.set(xi,yi,1);	
									cohortTmp.setU(Rotation);
									cohortTmp.setStockingdegreeMin(-1.);  cohortTmp.setStockingdegreeMax(-1.);
									cohort_all[asID]->setU(Rotation);
									cohort_all[asID]->setStockingdegreeMin(-1.);  cohort_all[asID]->setStockingdegreeMax(-1.);

									rotationForestNew.set(xi,yi,Rotation);	
									thinningForestNew.set(xi,yi,-1.);	
									newCohort_all[asID]->setU(Rotation);
									newCohort_all[asID]->setStockingdegreeMin(-1.);  newCohort_all[asID]->setStockingdegreeMax(-1.);  

									countryHarvestSWTmp = countryHarvestSWTmp - harvestSWTmp+dat_all[asID].DeforWoodM3.first;     
									countryHarvestRWTmp = countryHarvestRWTmp - harvestRWTmp+dat_all[asID].DeforWoodM3.second;
									if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {                
										coeff.PriceC.clear();
										coeff.PriceC.insert(0, priceC * iter->CORRUPTION[byear]);
										short int used = 0;
										if (thinningForest.get(xi,yi)>0){used = 1;}
										NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
									}else{NPVtmp = 1e50;};            


									if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 

										&& NPVtmp >= NPVcGrid.get(xi,yi)
										){
											sawnWoodHarvest[country] = countryHarvestSWTmp; 
											restWoodHarvest[country] = countryHarvestRWTmp; 
											harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first);
											harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
											manageChForest.set(xi,yi,1);
									}else{ // return old values
										rotationForest.set(xi,yi,rotationForestTmp);	
										thinningForest.set(xi,yi,thinningForestTmp);
										managedForest.set(xi,yi,managedForestTmp);	                   
										cohort_all[asID]->setU(rotationForestTmp);
										cohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   cohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);
										rotationForestNew.set(xi,yi,rotationForestNewTmp);	
										thinningForestNew.set(xi,yi,thinningForestNewTmp);
										newCohort_all[asID]->setU(rotationForestNewTmp);
										newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);                 
									}

								} // End  if ((data["POPDENS"].g(1990) == 0) && (data["GDP"].g(1990) == 0))

							} // End  if ((managedForest(xi,yi)>0))
						} // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].g(year))   
					} //End for PROTECT == 0
				} // End for some countries 
			} // Test only some regions   
			iter++;
		} // End for WHILE (cell loop) 
		//cout << "Zero pass is finished"<< endl; 
		//if (iter->COUNTRY[2000] == 197){cout<<"woodHarvest_f0= "<<woodHarvest[region-1]<<"\t woodHarvestStat= "<<wprod[countryprice].g(year) <<endl;}
	} // End if year




	//cout<<"Start First pass = adjust rotation time "<<endl;  
	//----First pass = adjust rotation time -------  
	//dataDetStruct::iterator iter = data_all.begin();
	iter = data_all.begin();
	while (iter != data_all.end()) {
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions         
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) {  
				int asID = iter->asID;
				if (iter->PROTECT[2000]==0) {
					char regionch[3];
					int2str(region,regionch);
					string regprice = "re"+string(regionch)+"price0";
					string regprice0 = "re"+string(regionch)+"price0";
					char countrych[4];
					int2str(country,countrych);
					string countryprice = "re"+string(countrych)+"price0";
					country = iter->COUNTRY[2000];     
					int xi = (iter->x);
					int yi = (iter->y);
					short int country = iter->COUNTRY[2000];  
					//float sawnW = 0.;
					//float restW = 0.;
					//float sawnThW = 0.;
					//float restThW = 0.;
					//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
					//float sawnWnew = 0.;
					//float restWnew = 0.;
					//float sawnThWnew = 0.;
					//float restThWnew = 0.;
					//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

					//  double harvWoodTmp = 0.;
					int biomassRotTh = 0;
					int rotMAI = 0;
					int rotMaxBmTh = 0;
					int Rotation = 0;
					//  double harvestTmp = 0;
					float newHarvestSWTmp = 0;
					float newHarvestRWTmp = 0;

					//  int rotationForestTmp = 0;
					//  int rotationForestNewTmp = 0;
					float countryHarvestSWTmp = sawnWoodHarvest[country]; 
					float countryHarvestRWTmp = restWoodHarvest[country]; 

					float NPVtmp = 0.;
					float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

					if (sawnWoodHarvest[country] < 0.99 * wprod_sawnwood[countryprice].g(year))  
					{
						if (managedForest.get(xi,yi)>=2)
						{
							if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
							{
								//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
							{
								//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].prevPlantPhytHaBmGr, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
							}                       


							int rotationForestTmp = rotationForest.get(xi,yi);
							int rotationForestNewTmp = rotationForestNew.get(xi,yi);
							g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                                     
							g4m::ageStruct cohortTmp = *cohort_all[asID];
							int newForAge = newCohort_all[asID]->getActiveAge();         
							float harvestSWTmp = harvestSWGrid.get(xi,yi);
							float harvestRWTmp = harvestRWGrid.get(xi,yi);
							//{cout<<"asID=\t"<<asID<<"\HarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t1<0.99 before"<<endl;}           	                   	                                            
							if (rotMAI < rotationForest.get(xi,yi)) // TEST for AT
							{
								//          Rotation = rotationForest.get(xi,yi) - 10;
								Rotation = rotationForest.get(xi,yi) - 5;          
								if (Rotation < rotMAI) {Rotation = rotMAI;}

								rotationForest.set(xi,yi,Rotation);		   
								cohortTmp.setU(Rotation);
								cohort_all[asID]->setU(Rotation);          
							}else  if (rotMAI > rotationForest.get(xi,yi))
								//       {Rotation = rotationForest.get(xi,yi) + 10;
							{Rotation = rotationForest.get(xi,yi) + 5;
							if (Rotation > rotMAI) {Rotation = rotMAI;}

							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							cohort_all[asID]->setU(Rotation);
							}
							if ((newForAge > biomassRotTh) && (rotMAI < biomassRotTh)) {
								//             Rotation = rotationForestNew.get(xi,yi) - 10;
								Rotation = rotationForestNew.get(xi,yi) - 5;
								if (Rotation < rotMAI) {Rotation = rotMAI;}

								cohortTmpNew.setU(Rotation);
								newCohort_all[asID]->setU(Rotation);    
								rotationForestNew.set(xi,yi,Rotation);
							}
							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

							float realAreaO = cohortTmp.getArea();
							float realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							float harvestSO = objHarvestWood.harvestWS;
							float harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							float harvestSN = objHarvestWood.harvestWS;
							float harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

							if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);       
								short int used = 0;
								if (thinningForest.get(xi,yi)>0){used = 1;}

								//                    NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used)*dat_all[asID].LandAreaHa * dat_all[asID].OforestShare;
								NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
							}else{NPVtmp=1e20;}
							//                if (NPVtmp == 0){NPVtmp +=0.01;}
							//cout<<"country(<0.97)= "<<country<<"\t fm_hurdle= "<<fm_hurdle<<"\t NPVtmp= \t"<<NPVtmp<<"\t NPVbau*fm_hurdle= \t"<<fm_hurdle*NPVbau[year-refYear-1][asID]<<endl;                
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp; 
									restWoodHarvest[country] = countryHarvestRWTmp; 
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								cohort_all[asID]->setU(rotationForestTmp);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								newCohort_all[asID]->setU(rotationForestNewTmp);
							}  

						}  // end for if (managedForest(xi,yi)>=2)          

					} //end for if (woodHarvest[region] <= 0.9 * wprod[region].g(year))
					else if (sawnWoodHarvest[country] > 1.01 * wprod_sawnwood[countryprice].g(year))     
					{
						if ((managedForest.get(xi,yi)>0) && (managedForest.get(xi,yi)<=2))
						{
							if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
							{
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
							{
								//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].prevPlantPhytHaBmGr, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
							}                       
							int rotationForestTmp = rotationForest.get(xi,yi);
							int rotationForestNewTmp = rotationForestNew.get(xi,yi);

							if (rotMaxBmTh > rotationForest.get(xi,yi)) //TEST for AT
							{
								g4m::ageStruct cohortTmp = *cohort_all[asID];
								g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
								int newForAge = newCohort_all[asID]->getActiveAge();  

								float harvestSWTmp = harvestSWGrid.get(xi,yi);
								float harvestRWTmp = harvestRWGrid.get(xi,yi);

								Rotation = rotationForest.get(xi,yi) + 5;
								if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
								rotationForest.set(xi,yi,Rotation);		   
								cohortTmp.setU(Rotation);
								cohort_all[asID]->setU(Rotation);       
								if ((newForAge > biomassRotTh) && (rotMaxBmTh > biomassRotTh)) {
									Rotation = rotationForestNew.get(xi,yi) + 5;
									if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}                          
									cohortTmpNew.setU(Rotation);
									newCohort_all[asID]->setU(Rotation);    
									rotationForestNew.set(xi,yi,Rotation);
								}

								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

								float realAreaO = cohortTmp.getArea();
								float realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
								float harvestSO = objHarvestWood.harvestWS;
								float harvestRO = objHarvestWood.harvestWR;
								objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
								float harvestSN = objHarvestWood.harvestWS;
								float harvestRN = objHarvestWood.harvestWR;
								newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestSN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
								newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestRN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
								countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
								countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

								if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
									coeff.PriceC.clear();
									coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);         
									short int used = 0;
									if (thinningForest.get(xi,yi)>0){used = 1;}

									NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
								}else{NPVtmp=1e20;}
								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
									&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
									){
										sawnWoodHarvest[country] = countryHarvestSWTmp;  
										restWoodHarvest[country] = countryHarvestRWTmp; 
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									cohort_all[asID]->setU(rotationForestTmp);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									newCohort_all[asID]->setU(rotationForestNewTmp);
								}  
								//{cout<<"asID=\t"<<asID<<"\HarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t1>1.01 after"<<endl;}           	                   	                                            
							}


						}  // end for if (managedForest(xi,yi)<=2)       
					}   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].g(year))
				} //End Protect
			}  // End some countries
		}  // Test only some regions 
		iter++;
	} // End While
	//if (iter->COUNTRY[2000] == 197){cout<<"woodHarvest_f1= "<<woodHarvest[region-1]<<"\t woodHarvestStat= "<<wprod[countryprice].g(year) <<endl;}   
	//cout<<" Finished First pass"<<endl;
	// ----- End of First pass


	//----Second pass = adjust rotation time -------  
	//dataDetStruct::iterator iter = data_all.begin();
	iter = data_all.begin();
	while (iter != data_all.end()) {
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions         
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) {  
				if (iter->PROTECT[2000]==0) {
					int asID = iter->asID;
					char regionch[3];
					int2str(region,regionch);
					string regprice = "re"+string(regionch)+"price0";
					string regprice0 = "re"+string(regionch)+"price0";
					char countrych[4];
					int2str(country,countrych);
					string countryprice = "re"+string(countrych)+"price0";
					country = iter->COUNTRY[2000];     
					int xi = (iter->x);
					int yi = (iter->y);

					//float sawnW = 0.;
					//float restW = 0.;
					//float sawnThW = 0.;
					//float restThW = 0.;
					//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
					//float sawnWnew = 0.;
					//float restWnew = 0.;
					//float sawnThWnew = 0.;
					//float restThWnew = 0.;
					//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
					int biomassRotTh = 0;
					int rotMAI = 0;
					int rotMaxBmTh = 0;
					int Rotation = 0;
					//  double harvestTmp = 0;
					float newHarvestSWTmp = 0;
					float newHarvestRWTmp = 0;

					//  int rotationForestTmp = 0;
					//  int rotationForestNewTmp = 0;
					float countryHarvestSWTmp = sawnWoodHarvest[country];  
					float countryHarvestRWTmp = restWoodHarvest[country]; 

					float NPVtmp = 0.;
					float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

					if (sawnWoodHarvest[country] < 0.98 * wprod_sawnwood[countryprice].g(year))  
					{
						if (managedForest.get(xi,yi)>0)
						{
							if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
							{
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
							{
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
							}                       

							int rotationForestTmp = rotationForest.get(xi,yi);
							int rotationForestNewTmp = rotationForestNew.get(xi,yi);
							g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                                    
							g4m::ageStruct cohortTmp = *cohort_all[asID];
							int newForAge = newCohort_all[asID]->getActiveAge();         
							float harvestSWTmp = harvestSWGrid.get(xi,yi);
							float harvestRWTmp = harvestRWGrid.get(xi,yi);
							if (rotMAI < rotationForest.get(xi,yi)) // TEST for AT
							{
								Rotation = rotationForest.get(xi,yi) - 5;          
								if (Rotation < rotMAI) {Rotation = rotMAI;}

								rotationForest.set(xi,yi,Rotation);		   
								cohortTmp.setU(Rotation);
								cohort_all[asID]->setU(Rotation);
							}else if (rotMAI > rotationForest.get(xi,yi))
								//       {Rotation = rotationForest.get(xi,yi) + 10;
							{Rotation = rotationForest.get(xi,yi) + 5;
							if (Rotation > rotMAI) {Rotation = rotMAI;}

							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							cohort_all[asID]->setU(Rotation);        
							}

							if ((newForAge > biomassRotTh) && (rotMAI < biomassRotTh)) {
								//             Rotation = rotationForestNew.get(xi,yi) - 10;
								Rotation = rotationForestNew.get(xi,yi) - 5;
								if (Rotation < rotMAI) {Rotation = rotMAI;}
								cohortTmpNew.setU(Rotation);
								newCohort_all[asID]->setU(Rotation);    
								rotationForestNew.set(xi,yi,Rotation);
							}
							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

							float realAreaO = cohortTmp.getArea();
							float realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							float harvestSO = objHarvestWood.harvestWS;
							float harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							float harvestSN = objHarvestWood.harvestWS;
							float harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
							if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);       
								short int used = 0;
								if (thinningForest.get(xi,yi)>0){used = 1;}

								NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
							}else{NPVtmp=1e20;}
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp; 
									restWoodHarvest[country] = countryHarvestRWTmp;
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								cohort_all[asID]->setU(rotationForestTmp);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								newCohort_all[asID]->setU(rotationForestNewTmp);
							}  
							//        }

							//{cout<<"asID=\t"<<asID<<"\HarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t2<0.98 after"<<endl;}           	                   	                                            
						}  // end for if (managedForest(xi,yi)>=2)          

					} //end for if (woodHarvest[region] <= 0.9 * wprod[region].g(year))

					else if (sawnWoodHarvest[country] > 1.02 * wprod_sawnwood[countryprice].g(year)) 
					{

						if (managedForest.get(xi,yi)>0)
						{

							if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
							{
								//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
							{
								//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].prevPlantPhytHaBmGr, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
								rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
								rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
							}                       

							int rotationForestTmp = rotationForest.get(xi,yi);
							int rotationForestNewTmp = rotationForestNew.get(xi,yi);

							g4m::ageStruct cohortTmp = *cohort_all[asID];
							g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
							int newForAge = newCohort_all[asID]->getActiveAge();  

							float harvestSWTmp = harvestSWGrid.get(xi,yi);
							float harvestRWTmp = harvestRWGrid.get(xi,yi);

							if (rotMaxBmTh > rotationForest.get(xi,yi)) //TEST for AT
							{
								Rotation = rotationForest.get(xi,yi) + 5;
								if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
								rotationForest.set(xi,yi,Rotation);		   
								cohortTmp.setU(Rotation);
								cohort_all[asID]->setU(Rotation);       

								if ((newForAge > biomassRotTh) && (rotMaxBmTh > biomassRotTh)) {
									Rotation = rotationForestNew.get(xi,yi) + 5;
									if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}                          
									cohortTmpNew.setU(Rotation);
									newCohort_all[asID]->setU(Rotation);    
									rotationForestNew.set(xi,yi,Rotation);
								}
								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

								float realAreaO = cohortTmp.getArea();
								float realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
								float harvestSO = objHarvestWood.harvestWS;
								float harvestRO = objHarvestWood.harvestWR;
								objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
								float harvestSN = objHarvestWood.harvestWS;
								float harvestRN = objHarvestWood.harvestWR;
								newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestSN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
								newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestRN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
								countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
								countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

								if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
									coeff.PriceC.clear();
									coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);      
									short int used = 0;
									if (thinningForest.get(xi,yi)>0){used = 1;}

									//                NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used)*dat_all[asID].LandAreaHa * dat_all[asID].OforestShare;
									NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
								}else{NPVtmp = 1e50;};
								//                if (NPVtmp == 0){NPVtmp +=0.01;}                
								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
									&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
									){
										sawnWoodHarvest[country] = countryHarvestSWTmp;   
										restWoodHarvest[country] = countryHarvestRWTmp;   
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
										//{cout<<"asID=\t"<<asID<<"\tnewHarvestTmp=\t"<<newHarvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t2>1.02"<<endl;}           	                   	        
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									cohort_all[asID]->setU(rotationForestTmp);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									newCohort_all[asID]->setU(rotationForestNewTmp);
								}  
							}

							//{cout<<"asID=\t"<<asID<<"\HarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t2>1.02 after"<<endl;}           	                   	                                            

						}  // end for if (managedForest(xi,yi)<=2)       
					}   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].g(year))
				} //End Protect
			}  // End some countries
		}  // Test only some regions 
		iter++;
	} // End While
	//if (iter->COUNTRY[2000] == 197){cout<<"woodHarvest_f1= "<<woodHarvest[region-1]<<"\t woodHarvestStat= "<<wprod[countryprice].g(year) <<endl;}   
	//cout<<" Finished First pass"<<endl;
	// ----- End of Second pass






	////------------------------------------------------------------------------------
	//// 
	//// -------Third pass: Adjust thinning -----------------------------------------------
	////
	////------------------------------------------------------------------------------
	//cout<<"Start Second pass = Adjust thinning"<< endl;
	iter = data_all.begin();
	while (iter != data_all.end())
	{
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions          
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) {  
				if (iter->PROTECT[2000] == 0)
				{
					int xi = (iter->x);
					int yi = (iter->y);
					if (maiForest.get(xi,yi)>0)
					{ 
						int asID = iter->asID;
						char regionch[3];
						int2str(region,regionch);
						string regprice = "re"+string(regionch)+"price0";
						string regprice0 = "re"+string(regionch)+"price0";
						char countrych[4];
						int2str(country,countrych);
						string countryprice = "re"+string(countrych)+"price0";
						country = iter->COUNTRY[2000];     
						//  if (floor(woodHarvest[country]) < 0.95 * wprod[countryprice].g(year))
						//  if (floor(woodHarvest[country]) < 0.94 * wprod[countryprice].g(year))  
						if (floor(sawnWoodHarvest[country]) < 0.97 * wprod_sawnwood[countryprice].g(year))  
						{
							if ((managedForest.get(xi,yi)<=0))
							{
								//float sawnW = 0.;
								//float restW = 0.;
								//float sawnThW = 0.;
								//float restThW = 0.;
								//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
								//float sawnWnew = 0.;
								//float restWnew = 0.;
								//float sawnThWnew = 0.;
								//float restThWnew = 0.;
								//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

								//      float harvWoodTmp = 0.;
								//      float harvestTmp = 0.;
								float newHarvestSWTmp = 0.;
								float newHarvestRWTmp = 0.;

								int biomassRot = 0;
								int biomassRotTh2 = 0;
								int rotMAI = 0;
								int rotMaxBmTh = 0;
								int Rotation = 0;

								//      int managedForestTmp = 0;
								//      double thinningForestTmp = 0;
								//      int rotationForestTmp = 0;
								//      double thinningForestNewTmp = 0;
								//      int managedForestNewTmp = 0;
								//      int rotationForestNewTmp = 0;

								float countryHarvestSWTmp = sawnWoodHarvest[country]; 
								float countryHarvestRWTmp = restWoodHarvest[country]; 

								float NPVtmp = 0.;
								float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];
								float stockingDegree = 1.5;     
								g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
								g4m::ageStruct cohortTmp = *cohort_all[asID];
								float harvestSWTmp = harvestSWGrid.get(xi,yi);
								float harvestRWTmp = harvestRWGrid.get(xi,yi);

								if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
								{
									biomassRotTh2 = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(iter->CABOVEHA[byear], maiForest.get(xi,yi), stockingDegree);     // rotation time to get current biomass (with thinning)            
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
								}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
								{

									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
									biomassRotTh2 = rotMaxBmTh;
								}     

								int rotationForestTmp = rotationForest.get(xi,yi);
								int managedForestTmp = managedForest.get(xi,yi);
								float thinningForestTmp = thinningForest.get(xi,yi);
								int rotationForestNewTmp = rotationForestNew.get(xi,yi);
								float thinningForestNewTmp = thinningForestNew.get(xi,yi);

								//---------------------------------------------------
								int newForAge = newCohort_all[asID]->getActiveAge();  
								if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
								{
									if (useChange.find(asID)!=useChange.end()){
										managedForest.set(xi,yi,managedForestTmp+3);
										Rotation = rotationForestO.get(xi,yi);
										stockingDegree = thinningForestO.get(xi,yi);
									}else if (managedForest.get(xi,yi) == 0)
									{
										managedForest.set(xi,yi,3);
										Rotation = biomassRotTh2+1;
										if (Rotation < rotMAI) {Rotation = rotMAI;}
										rotationType.set(xi,yi,1);
									}else 
									{           managedForest.set(xi,yi,2);
									Rotation = biomassRotTh2+1;
									if (Rotation < rotMAI) {Rotation = rotMAI;}
									//                         if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
									rotationType.set(xi,yi,3);
									}

									rotationForest.set(xi,yi,Rotation);	
									thinningForest.set(xi,yi,stockingDegree);	
									cohortTmp.setU(Rotation);
									cohortTmp.setStockingdegreeMin(stockingDegree*sdMinCoeff);  cohortTmp.setStockingdegreeMax(stockingDegree*sdMaxCoeff);
									cohort_all[asID]->setU(Rotation);
									cohort_all[asID]->setStockingdegreeMin(stockingDegree*sdMinCoeff);  cohort_all[asID]->setStockingdegreeMax(stockingDegree*sdMaxCoeff);

									rotationForestNew.set(xi,yi,Rotation);	
									thinningForestNew.set(xi,yi,stockingDegree);	
									cohortTmpNew.setU(Rotation);
									cohortTmpNew.setStockingdegreeMin(stockingDegree*sdMinCoeff);  cohortTmpNew.setStockingdegreeMax(stockingDegree*sdMaxCoeff);        
									newCohort_all[asID]->setU(Rotation);
									newCohort_all[asID]->setStockingdegreeMin(stockingDegree*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(stockingDegree*sdMaxCoeff);  

								}else{
									if (useChange.find(asID)!=useChange.end()){
										managedForest.set(xi,yi,managedForestTmp+3);
										Rotation = rotationForestO.get(xi,yi);
										stockingDegree = thinningForestO.get(xi,yi);
										rotationForestNew.set(xi,yi,rotationForestNewO.get(xi,yi));	
										thinningForestNew.set(xi,yi,stockingDegree);	
										newCohort_all[asID]->setU(rotationForestNewO.get(xi,yi));
										newCohort_all[asID]->setStockingdegreeMin(stockingDegree*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(stockingDegree*sdMaxCoeff);  
									}else if (managedForest.get(xi,yi) == 0)
									{
										managedForest.set(xi,yi,3);
										Rotation = biomassRotTh2+1;
										if (Rotation < rotMAI) {Rotation = rotMAI;}
										rotationType.set(xi,yi,1);
									}else 
									{     managedForest.set(xi,yi,2);
									Rotation = biomassRotTh2+1;
									if (Rotation < rotMAI) {Rotation = rotMAI;}                         
									//                         if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
									rotationType.set(xi,yi,3);
									}

									rotationForest.set(xi,yi,Rotation);	
									thinningForest.set(xi,yi,stockingDegree);	
									cohortTmp.setU(Rotation);
									cohortTmp.setStockingdegreeMin(stockingDegree*sdMinCoeff);  cohortTmp.setStockingdegreeMax(stockingDegree*sdMaxCoeff);
									cohort_all[asID]->setU(Rotation);
									cohort_all[asID]->setStockingdegreeMin(stockingDegree*sdMinCoeff);  cohort_all[asID]->setStockingdegreeMax(stockingDegree*sdMaxCoeff);
								} // End    else // New forest age < rotation -> don't change FM for the new forest
								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

								float realAreaO = cohortTmp.getArea();
								float realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
								float harvestSO = objHarvestWood.harvestWS;
								float harvestRO = objHarvestWood.harvestWR;
								objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
								float harvestSN = objHarvestWood.harvestWS;
								float harvestRN = objHarvestWood.harvestWR;
								newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestSN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
								newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestRN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
								countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
								countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

								if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) { 
									coeff.PriceC.clear();
									coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);                
									short int used = 0;
									if (thinningForest.get(xi,yi)>0){used = 1;}

									NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
								}else{NPVtmp = 1e50;};
								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 

									&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
									){
										sawnWoodHarvest[country] = countryHarvestSWTmp; 
										restWoodHarvest[country] = countryHarvestRWTmp;
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
										if (useChange.find(asID)== useChange.end()){manageChForest.set(xi,yi,1);}
										//{cout<<"asID=\t"<<asID<<"\tnewHarvestTmp=\t"<<newHarvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t3<0.97"<<endl;}           	                   	        
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									thinningForest.set(xi,yi,thinningForestTmp);
									managedForest.set(xi,yi,managedForestTmp);	                   
									cohort_all[asID]->setU(rotationForestTmp);
									cohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   cohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									thinningForestNew.set(xi,yi,thinningForestNewTmp);
									//                managedForestNew.set(xi,yi,managedForestNewTmp)	                   
									newCohort_all[asID]->setU(rotationForestNewTmp);
									newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);                 
								}

							} // End  if ((managedForest(xi,yi)<=0) & (managedForest(xi,yi)>-2)

						}else if (floor(sawnWoodHarvest[country]) > 1.03 * wprod_sawnwood[countryprice].g(year))   
						{
							if ((managedForest.get(xi,yi)>0) && (managedForest.get(xi,yi)<3))
							{
								//float sawnW = 0.;
								//float restW = 0.;
								//float sawnThW = 0.;
								//float restThW = 0.;
								//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
								//float sawnWnew = 0.;
								//float restWnew = 0.;
								//float sawnThWnew = 0.;
								//float restThWnew = 0.;
								//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha


								//      float harvestTmp = 0.;      
								float newHarvestSWTmp = 0;
								float newHarvestRWTmp = 0;
								int biomassRot = 0;
								int rotMAI = 0;
								int rotMaxBm = 0;
								int Rotation = 0;

								//      int managedForestTmp = 0;
								//      double thinningForestTmp = 0;
								//      int rotationForestTmp = 0;
								//      double thinningForestNewTmp = 0;
								//      int managedForestNewTmp = 0;
								//      int rotationForestNewTmp = 0;

								float countryHarvestSWTmp = sawnWoodHarvest[country];
								float countryHarvestRWTmp = restWoodHarvest[country];

								float NPVtmp = 0.;
								float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];
								g4m::ageStruct cohortTmp = *cohort_all[asID];
								float harvestSWTmp = harvestSWGrid.get(xi,yi); 
								float harvestRWTmp = harvestRWGrid.get(xi,yi);

								//{cout<<"asID=\t"<<asID<<"\HarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t3>1.03 before"<<endl;}           	                   	                                             
								if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
								{
									//          biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (without thinning)  
									//          biomassRot = species[int(iter->SPECIESTYPE[byear])-1]->gU(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), 1);     // rotation time to get current biomass (without thinning)  
									rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
								}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
								{
									//          biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUt(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (without thinning)  
									//          biomassRot = species[int(iter->SPECIESTYPE[byear])-1]->gU(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), 1);     // rotation time to get current biomass (without thinning)  
									rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
								}     

								int rotationForestTmp = rotationForest.get(xi,yi);
								int managedForestTmp = managedForest.get(xi,yi);
								float thinningForestTmp = thinningForest.get(xi,yi);
								int rotationForestNewTmp = rotationForestNew.get(xi,yi);
								//      managedForestNewTmp = managedForestNew.get(xi,yi);
								float thinningForestNewTmp = thinningForestNew.get(xi,yi);

								//---------------------------------------------------
								//           int newForAge = newCohort_all[asID]->getActiveAge();  
								//cout<<"newForAge= "<<newForAge<<endl;
								//           if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
								//            {
								//cout<<"biomassRot = "<< biomassRot<<endl;
								//              g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
								if (managedForest.get(xi,yi) == 2)
								{
									// cout<<"if (managedForest.get(xi,yi) == 2)"<<endl;
									managedForest.set(xi,yi,-1);
									//                         Rotation = biomassRot+1;
									if (rotationForestTmp < rotMaxBm) {
										Rotation = rotMaxBm;
									}else{     
										Rotation = rotationForestTmp;  }
									//                         if (Rotation > rotMaxBm) {Rotation = rotMaxBm;}
									rotationType.set(xi,yi,1);
								}else 
								{     managedForest.set(xi,yi,-2);
								//                         Rotation = biomassRot+1;
								if (rotationForestTmp < rotMaxBm) {
									Rotation = rotMaxBm;
								}else{     
									Rotation = rotationForestTmp;  }
								//                         if (Rotation > rotMaxBm) {Rotation = rotMaxBm;}
								rotationType.set(xi,yi,3);
								//cout<<"if (managedForest.get(xi,yi) == 2) ELSE"<<endl;
								}

								rotationForest.set(xi,yi,Rotation);	
								thinningForest.set(xi,yi,-1.);
								cohortTmp.setU(Rotation);
								cohortTmp.setStockingdegreeMin(-1.);  cohortTmp.setStockingdegreeMax(-1.);
								cohort_all[asID]->setU(Rotation);
								cohort_all[asID]->setStockingdegreeMin(-1.);  cohort_all[asID]->setStockingdegreeMax(-1.);

								rotationForestNew.set(xi,yi,Rotation);	
								thinningForestNew.set(xi,yi,-1.);	

								newCohort_all[asID]->setU(Rotation);
								newCohort_all[asID]->setStockingdegreeMin(-1.);  newCohort_all[asID]->setStockingdegreeMax(-1.);  

								countryHarvestSWTmp = countryHarvestSWTmp-harvestSWTmp+dat_all[asID].DeforWoodM3.first; 
								countryHarvestRWTmp = countryHarvestRWTmp-harvestRWTmp+dat_all[asID].DeforWoodM3.second;

								if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
									coeff.PriceC.clear();
									coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);     
									short int used = 0;
									if (thinningForest.get(xi,yi)>0){used = 1;}

									NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
								}else{NPVtmp = 1e50;};
								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 

									&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
									){
										sawnWoodHarvest[country] = countryHarvestSWTmp; 
										restWoodHarvest[country] = countryHarvestRWTmp;
										harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first);
										harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
										manageChForest.set(xi,yi,1);
										//{cout<<"asID=\t"<<asID<<"\tnewHarvestTmp=\t"<<dat_all[asID].deforWoodTotM3<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t3>1.03"<<endl;}           	                   	        
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									thinningForest.set(xi,yi,thinningForestTmp);
									managedForest.set(xi,yi,managedForestTmp);	                   
									cohort_all[asID]->setU(rotationForestTmp);
									cohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   cohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									thinningForestNew.set(xi,yi,thinningForestNewTmp);
									//                managedForestNew.set(xi,yi,managedForestNewTmp)	                   
									newCohort_all[asID]->setU(rotationForestNewTmp);
									newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);                 
								}
								//         } // End if ((cohort_all[asID]->getArea(0) == 0.)||(cohort_all[asID]->getArea(0) >= 1/400.)){

								//{cout<<"asID=\t"<<asID<<"\HarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t3>1.03 after"<<endl;}           	                   	                                            
							} // End  if ((managedForest(xi,yi)>0) & (managedForest(xi,yi)<3)
						} // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].g(year))   

					} //End for mai>0
				} //End for PROTECT == 0
			}  // End some countries
		} // Test only some regions   
		iter++;
	} // End for WHILE (cell loop) 

	//***********************************************************************************************
	iter = data_all.begin();
	while (iter != data_all.end())
	{
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions          
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) {  
				if (iter->PROTECT[2000] == 0)
				{
					int asID = iter->asID;
					char regionch[3];
					int2str(region,regionch);
					string regprice = "re"+string(regionch)+"price0";
					string regprice0 = "re"+string(regionch)+"price0";
					char countrych[4];
					int2str(country,countrych);
					string countryprice = "re"+string(countrych)+"price0";
					country = iter->COUNTRY[2000];     
					int xi = (iter->x);
					int yi = (iter->y);

					if (floor(sawnWoodHarvest[country]) > 1.05 * wprod_sawnwood[countryprice].g(year))  
					{
						if ((managedForest.get(xi,yi)>0))
						{
							//float sawnW = 0.;
							//float restW = 0.;
							//float sawnThW = 0.;
							//float restThW = 0.;
							//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
							//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
							//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
							//float sawnWnew = 0.;
							//float restWnew = 0.;
							//float sawnThWnew = 0.;
							//float restThWnew = 0.;
							//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
							//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
							//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha


							float newHarvestSWTmp = 0;
							float newHarvestRWTmp = 0;

							int biomassRot = 0;

							int rotMAI = 0;
							int rotMaxBm = 0;
							int Rotation = 0;

							//      int managedForestTmp = 0;
							//      double thinningForestTmp = 0;
							//      int rotationForestTmp = 0;
							//      double thinningForestNewTmp = 0;
							//      int managedForestNewTmp = 0;
							//      int rotationForestNewTmp = 0;

							float countryHarvestSWTmp = sawnWoodHarvest[country];
							float countryHarvestRWTmp = restWoodHarvest[country];

							float NPVtmp = 0.;
							float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];
							g4m::ageStruct cohortTmp = *cohort_all[asID];
							float harvestSWTmp = harvestSWGrid.get(xi,yi); 
							float harvestRWTmp = harvestRWGrid.get(xi,yi); 

							if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
							{
								rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
							}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
							{
								rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
							}     
							int rotationForestTmp = rotationForest.get(xi,yi);
							int managedForestTmp = managedForest.get(xi,yi);
							float thinningForestTmp = thinningForest.get(xi,yi);
							int rotationForestNewTmp = rotationForestNew.get(xi,yi);
							float thinningForestNewTmp = thinningForestNew.get(xi,yi);

							//---------------------------------------------------
							{
								if (managedForest.get(xi,yi) == 3)
								{
									managedForest.set(xi,yi,0);
									if (rotationForestTmp < rotMaxBm) {
										Rotation = rotMaxBm;
									}else{     
										Rotation = rotationForestTmp;  }
									rotationType.set(xi,yi,1);
								}else if (managedForest.get(xi,yi) == 2)
								{
									managedForest.set(xi,yi,-1);
									if (rotationForestTmp < rotMaxBm) {
										Rotation = rotMaxBm;
									}else{     
										Rotation = rotationForestTmp;  }
									rotationType.set(xi,yi,1);
								}else 
								{     managedForest.set(xi,yi,-2);
								if (rotationForestTmp < rotMaxBm) {
									Rotation = rotMaxBm;
								}else{     
									Rotation = rotationForestTmp;  }
								rotationType.set(xi,yi,3);
								}

								rotationForest.set(xi,yi,Rotation);	
								thinningForest.set(xi,yi,-1.);
								cohortTmp.setU(Rotation);
								cohortTmp.setStockingdegreeMin(-1.);  cohortTmp.setStockingdegreeMax(-1.);                
								cohort_all[asID]->setU(Rotation);
								cohort_all[asID]->setStockingdegreeMin(-1.);  cohort_all[asID]->setStockingdegreeMax(-1.);

								rotationForestNew.set(xi,yi,Rotation);	
								thinningForestNew.set(xi,yi,-1.);	

								newCohort_all[asID]->setU(Rotation);
								newCohort_all[asID]->setStockingdegreeMin(-1.);  newCohort_all[asID]->setStockingdegreeMax(-1.);  

							} // End    else // New forest age < rotation -> don't change FM for the new forest

							countryHarvestSWTmp = countryHarvestSWTmp-harvestSWTmp+dat_all[asID].DeforWoodM3.first; 
							countryHarvestRWTmp = countryHarvestRWTmp-harvestRWTmp+dat_all[asID].DeforWoodM3.second;

							if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) { 
								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);   
								short int used = 0;
								if (thinningForest.get(xi,yi)>0){used = 1;}

								NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
							}else{NPVtmp = 1e50;};
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp;  
									sawnWoodHarvest[country] = countryHarvestSWTmp;                                      

									harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first);
									harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
									manageChForest.set(xi,yi,1);
							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								thinningForest.set(xi,yi,thinningForestTmp);
								managedForest.set(xi,yi,managedForestTmp);	                   
								cohort_all[asID]->setU(rotationForestTmp);
								cohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   cohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								thinningForestNew.set(xi,yi,thinningForestNewTmp);
								newCohort_all[asID]->setU(rotationForestNewTmp);
								newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);                 
							}
						} // End  if ((managedForest(xi,yi)>0) & (managedForest(xi,yi)<3)
					} // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].g(year)) 

				} //End for PROTECT == 0
			} //End some countries
		} // Test only some regions   
		iter++;
	} // End for WHILE (cell loop) 


	//*****************************************************************************


	//******************************************************************************
	//**************************Forth Pass********************
	//******************************************************************************
	//cout << "Start forth pass" << endl;
	iter = data_all.begin();
	while (iter != data_all.end()) {
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) {  
				if (iter->PROTECT[2000]==0) {
					int asID = iter->asID;
					char regionch[3];
					int2str(region,regionch);
					string regprice = "re"+string(regionch)+"price0";
					string regprice0 = "re"+string(regionch)+"price0";
					char countrych[4];
					int2str(country,countrych);
					string countryprice = "re"+string(countrych)+"price0";
					country = iter->COUNTRY[2000];     
					int xi = (iter->x);
					int yi = (iter->y);


					//float sawnW = 0.;
					//float restW = 0.;
					//float sawnThW = 0.;
					//float restThW = 0.;
					//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
					//float sawnWnew = 0.;
					//float restWnew = 0.;
					//float sawnThWnew = 0.;
					//float restThWnew = 0.;
					//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha


					int biomassRotTh = 0;
					int biomassRot = 0;  
					int rotMAI = 0;
					int rotMaxBmTh = 0;
					int rotMaxBm = 0;  
					int Rotation = 0;
					//      float harvestTmp = 0;
					float newHarvestSWTmp = 0;
					float newHarvestRWTmp = 0;

					//      int rotationForestTmp = 0;
					//      int rotationForestNewTmp = 0;
					float countryHarvestSWTmp =sawnWoodHarvest[country];
					float countryHarvestRWTmp =restWoodHarvest[country];

					float NPVtmp = 0.;
					float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

					//  if (floor(woodHarvest[country]) < 0.97 * wprod[countryprice].g(year))
					if (floor(sawnWoodHarvest[country]) < 0.99 * wprod_sawnwood[countryprice].g(year))
					{if (managedForest.get(xi,yi)>1)
					{
						if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
						{
							biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
						}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
						{
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							biomassRotTh =  rotMAI;
						}            
						int rotationForestTmp = rotationForest.get(xi,yi);
						int rotationForestNewTmp = rotationForestNew.get(xi,yi);
						g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                                     
						g4m::ageStruct cohortTmp = *cohort_all[asID];
						int newForAge = newCohort_all[asID]->getActiveAge();          
						float harvestSWTmp = harvestSWGrid.get(xi,yi);
						float harvestRWTmp = harvestRWGrid.get(xi,yi);

						if (rotMAI < rotationForest.get(xi,yi)) // TEST for AT
						{
							//          Rotation = rotationForest.get(xi,yi) - 10;
							Rotation = rotationForest.get(xi,yi) - 5;          
							if (Rotation < rotMAI) {Rotation = rotMAI;}

							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							cohort_all[asID]->setU(Rotation);
						}else  if (rotMAI > rotationForest.get(xi,yi))

						{Rotation = rotationForest.get(xi,yi) + 5;
						if (Rotation > rotMAI) {Rotation = rotMAI;}

						rotationForest.set(xi,yi,Rotation);		   
						cohortTmp.setU(Rotation);
						cohort_all[asID]->setU(Rotation);        
						}
						if ((newForAge > biomassRotTh) && (rotMAI < biomassRotTh)) {
							//             Rotation = rotationForestNew.get(xi,yi) - 10;
							Rotation = rotationForestNew.get(xi,yi) - 5;             
							if (Rotation < rotMAI) {Rotation = rotMAI;}
							cohortTmpNew.setU(Rotation);
							newCohort_all[asID]->setU(Rotation);    
							rotationForestNew.set(xi,yi,Rotation);
						}
						pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

						pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

						float realAreaO = cohortTmp.getArea();
						float realAreaN = cohortTmpNew.getArea();
						objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
						float harvestSO = objHarvestWood.harvestWS;
						float harvestRO = objHarvestWood.harvestWR;
						objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
						float harvestSN = objHarvestWood.harvestWS;
						float harvestRN = objHarvestWood.harvestWR;
						newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
							harvestSN * dat_all[asID].AforestShare) * 
							dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
						newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
							harvestRN * dat_all[asID].AforestShare) * 
							dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
						countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
						countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

						if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {  
							coeff.PriceC.clear();
							coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);                
							short int used = 0;
							if (thinningForest.get(xi,yi)>0){used = 1;}

							NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
						}else{NPVtmp = 1e50;};
						//                if (NPVtmp == 0){NPVtmp +=0.01;}                
						if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
							&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
							){
								sawnWoodHarvest[country] = countryHarvestSWTmp;
								restWoodHarvest[country] = countryHarvestRWTmp;                                      

								harvestSWGrid.set(xi,yi,newHarvestSWTmp);
								harvestRWGrid.set(xi,yi,newHarvestRWTmp);
								//{cout<<"asID=\t"<<asID<<"\tnewHarvestTmp=\t"<<newHarvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t4<0.99"<<endl;}           	                   	        
						}else{ // return old values
							rotationForest.set(xi,yi,rotationForestTmp);	
							cohort_all[asID]->setU(rotationForestTmp);
							rotationForestNew.set(xi,yi,rotationForestNewTmp);	
							newCohort_all[asID]->setU(rotationForestNewTmp);
						}        

					}  // end for if (managedForest(xi,yi)>=2)          

					} //end for if (woodHarvest[region] <= 0.9 * wprod[region].g(year))
					else if (floor(sawnWoodHarvest[country]) > 1.01 * wprod_sawnwood[countryprice].g(year))     
					{if ((managedForest.get(xi,yi)>0) && (managedForest.get(xi,yi)<3))
					{
						if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
						{
							biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
						}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
						{
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							biomassRotTh = rotMAI;
						}     

						int rotationForestTmp = rotationForest.get(xi,yi);
						int rotationForestNewTmp = rotationForestNew.get(xi,yi);

						g4m::ageStruct cohortTmp = *cohort_all[asID];
						g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
						int newForAge = newCohort_all[asID]->getActiveAge();  
						float harvestSWTmp = harvestSWGrid.get(xi,yi);
						float harvestRWTmp = harvestRWGrid.get(xi,yi);

						if (rotMaxBmTh > rotationForest.get(xi,yi)) //TEST AT
						{
							//        Rotation = rotationForest.get(xi,yi) + 10;
							Rotation = rotationForest.get(xi,yi) + 5;
							if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}


							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							cohort_all[asID]->setU(Rotation);       

							if ((newForAge > biomassRotTh) && (rotMaxBmTh > biomassRotTh)) {
								Rotation = rotationForestNew.get(xi,yi) + 5;             
								if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
								cohortTmpNew.setU(Rotation);
								newCohort_all[asID]->setU(Rotation);    
								rotationForestNew.set(xi,yi,Rotation);
							}
							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

							float realAreaO = cohortTmp.getArea();
							float realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							float harvestSO = objHarvestWood.harvestWS;
							float harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							float harvestSN = objHarvestWood.harvestWS;
							float harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

							if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {  
								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);      
								short int used = 0;
								if (thinningForest.get(xi,yi)>0){used = 1;}

								NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
							}else{NPVtmp = 1e50;};
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp;    
									restWoodHarvest[country] = countryHarvestRWTmp;                                      

									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);

							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								cohort_all[asID]->setU(rotationForestTmp);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								newCohort_all[asID]->setU(rotationForestNewTmp);
							}        
						}
					}  // end for if (managedForest(xi,yi)<=2)  
					}   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].g(year))
				} // End protect
			} // End some countries
		} // Test only some regions
		iter++;
	} // End While

	//************************End of Forth Pass************************************
	//cout << "End of Forth pass"<<endl;


	//******************************************************************************
	//**************************Fifth Pass********************
	//******************************************************************************
	//cout << "Start fifth pass" << endl;
	iter = data_all.begin();
	while (iter != data_all.end()) {
		short int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions
			short int country = iter->COUNTRY[2000];
			if (toAdjust.find(country) != toAdjust.end()) {  
				if (iter->PROTECT[2000]==0) {
					int asID = iter->asID;
					char regionch[3];
					int2str(region,regionch);
					string regprice = "re"+string(regionch)+"price0";
					string regprice0 = "re"+string(regionch)+"price0";
					char countrych[4];
					int2str(country,countrych);
					string countryprice = "re"+string(countrych)+"price0";
					country = iter->COUNTRY[2000]; 

					int xi = (iter->x);
					int yi = (iter->y);


					//float sawnW = 0.;
					//float restW = 0.;
					//float sawnThW = 0.;
					//float restThW = 0.;
					//float bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
					//float sawnWnew = 0.;
					//float restWnew = 0.;
					//float sawnThWnew = 0.;
					//float restThWnew = 0.;
					//float bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
					//float bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
					//float harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha


					int biomassRotTh = 0;
					int biomassRot = 0;  
					int rotMAI = 0;
					int rotMaxBmTh = 0;
					int rotMaxBm = 0;  
					int Rotation = 0;
					float newHarvestSWTmp = 0;
					float newHarvestRWTmp = 0;
					float countryHarvestSWTmp = sawnWoodHarvest[country];
					float countryHarvestRWTmp = restWoodHarvest[country];


					float NPVtmp = 0.;
					float maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];

					if (floor(sawnWoodHarvest[country]) < 0.98 * wprod_sawnwood[countryprice].g(year))
					{if (managedForest.get(xi,yi)>0)
					{
						if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
						{
							biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         

						}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
						{
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							biomassRotTh = rotMAI;
						}     

						int rotationForestTmp = rotationForest.get(xi,yi);
						int rotationForestNewTmp = rotationForestNew.get(xi,yi);
						g4m::ageStruct cohortTmp = *cohort_all[asID];
						g4m::ageStruct cohortTmpNew = *newCohort_all[asID];            
						int newForAge = newCohort_all[asID]->getActiveAge();          
						float harvestSWTmp = harvestSWGrid.get(xi,yi);
						float harvestRWTmp = harvestRWGrid.get(xi,yi);


						if (rotMAI < rotationForest.get(xi,yi)) // TEST for AT
						{
							//          Rotation = rotationForest.get(xi,yi) - 10;
							Rotation = rotationForest.get(xi,yi) - 5;          
							if (Rotation < rotMAI) {Rotation = rotMAI;}

							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							cohort_all[asID]->setU(Rotation);
						}else  if (rotMAI > rotationForest.get(xi,yi))
							//       {Rotation = rotationForest.get(xi,yi) + 10;
						{Rotation = rotationForest.get(xi,yi) + 5;       
						if (Rotation > rotMAI) {Rotation = rotMAI;}
						else {Rotation = rotMAI;}

						rotationForest.set(xi,yi,Rotation);		   
						cohortTmp.setU(Rotation);
						cohort_all[asID]->setU(Rotation);
						}


						if ((newForAge > biomassRotTh) && (rotMAI < biomassRotTh)) {
							//             Rotation = rotationForestNew.get(xi,yi) - 10;
							Rotation = rotationForestNew.get(xi,yi) - 5;
							if (Rotation < rotMAI) {Rotation = rotMAI;}
							cohortTmpNew.setU(Rotation);
							newCohort_all[asID]->setU(Rotation);    
							rotationForestNew.set(xi,yi,Rotation);
						}
						pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();


						pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

						float realAreaO = cohortTmp.getArea();
						float realAreaN = cohortTmpNew.getArea();
						objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
						float harvestSO = objHarvestWood.harvestWS;
						float harvestRO = objHarvestWood.harvestWR;
						objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
						float harvestSN = objHarvestWood.harvestWS;
						float harvestRN = objHarvestWood.harvestWR;
						newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
							harvestSN * dat_all[asID].AforestShare) * 
							dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
						newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
							harvestRN * dat_all[asID].AforestShare) * 
							dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
						countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
						countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

						if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
							coeff.PriceC.clear();
							coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);                
							short int used = 0;
							if (thinningForest.get(xi,yi)>0){used = 1;}

							NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
						}else{NPVtmp = 1e50;};

						if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
							&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
							){
								sawnWoodHarvest[country] = countryHarvestSWTmp;  
								restWoodHarvest[country] = countryHarvestRWTmp;
								harvestSWGrid.set(xi,yi,newHarvestSWTmp);
								harvestRWGrid.set(xi,yi,newHarvestRWTmp);
						}else{ // return old values
							rotationForest.set(xi,yi,rotationForestTmp);	
							cohort_all[asID]->setU(rotationForestTmp);
							rotationForestNew.set(xi,yi,rotationForestNewTmp);	
							newCohort_all[asID]->setU(rotationForestNewTmp);
						}
					}
					}

					else if (floor(sawnWoodHarvest[country]) > 1.02 * wprod_sawnwood[countryprice].g(year))     
					{if (managedForest.get(xi,yi)>0) 
					{
						if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
						{
							biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
						}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
						{
							rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
							rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
							biomassRotTh = rotMAI;
						}     

						int rotationForestTmp = rotationForest.get(xi,yi);
						int rotationForestNewTmp = rotationForestNew.get(xi,yi);

						g4m::ageStruct cohortTmp = *cohort_all[asID];
						g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
						int newForAge = newCohort_all[asID]->getActiveAge();  
						float harvestSWTmp = harvestSWGrid.get(xi,yi);
						float harvestRWTmp = harvestRWGrid.get(xi,yi);

						if (rotMaxBmTh > rotationForest.get(xi,yi)) //TEST AT
						{
							//        Rotation = rotationForest.get(xi,yi) + 10;
							Rotation = rotationForest.get(xi,yi) + 5;        
							if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
							rotationForest.set(xi,yi,Rotation);		   
							cohortTmp.setU(Rotation);
							cohort_all[asID]->setU(Rotation);       

							if ((newForAge > biomassRotTh) && (rotMaxBmTh > biomassRotTh)) {
								//             Rotation = rotationForestNew.get(xi,yi) + 10;
								Rotation = rotationForestNew.get(xi,yi) + 5;             
								if (Rotation > rotMaxBmTh) {Rotation = rotMaxBmTh;}
								cohortTmpNew.setU(Rotation);
								newCohort_all[asID]->setU(Rotation);    
								rotationForestNew.set(xi,yi,Rotation);
							}
							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();

							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

							float realAreaO = cohortTmp.getArea();
							float realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							float harvestSO = objHarvestWood.harvestWS;
							float harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							float harvestSN = objHarvestWood.harvestWS;
							float harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first;
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second;
							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);

							if (countriesFmcpol.find(country) != countriesFmcpol.end() && (countriesNoFmcpol.find(country)==countriesNoFmcpol.end())) {
								coeff.PriceC.clear();
								coeff.PriceC.insert(0, priceC * iter->CORRUPTION[2000]);                
								short int used = 0;
								if (thinningForest.get(xi,yi)>0){used = 1;}

								NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,Rotation,regprice,regprice0,dat_all[asID].OforestShare,used);
							}else{NPVtmp = 1e50;};
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								&& fm_hurdle*NPVtmp >= NPVcGrid.get(xi,yi)
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp;
									restWoodHarvest[country] = countryHarvestRWTmp;
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);

							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								cohort_all[asID]->setU(rotationForestTmp);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								newCohort_all[asID]->setU(rotationForestNewTmp);
							}        
						}
						//{cout<<"asID=\t"<<asID<<"\tHarvestTmp=\t"<<harvestTmp<<"\tthinning=\t"<<thinningForest.get(xi,yi)<<"\tSD=\t"<<cohort_all[asID]->getStockingdegree()<<"\tRotation=\t"<<rotationForest.get(xi,yi)<<"\tRotPeriod=\t"<<cohort_all[asID]->getRotPeriod()<<"\t5>1.02 before"<<endl;}           	                   	                                            
					}  // end for if (managedForest(xi,yi)<=2)  
					}   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].g(year))
					//if (country==47){cout<<"thinningForest(FM100E)\t"<<thinningForest.get(xi,yi)<<"\tmanagedForest\t"<<int(managedForest.get(xi,yi))<<endl;}
				} // End protect
			} // End some countries
		} // Test only some regions
		iter++;
	} // End While
	// 
	//************************End of Fifth Pass************************************
	//cout << "End of Fifth pass"<<endl;



	if (toAdjust.size()>1){
		for (int i=1; i<=NumberOfCountries; i++){
			char countrych[4];
			int ii = i;
			int2str(ii,countrych);
			string countryprice = "re"+string(countrych)+"price0";
			if (wprod_sawnwood[countryprice].g(year) > 0){

				if ((regions.find(countryRegion[i]) != regions.end()) 
					&& (toAdjust.find(i) != toAdjust.end())  // Test only some countries         
					//                         && i!=114  // Luxenbourg is too small (1 cell) to be considered correctly
					)
				{                                  
					harvDiff[i] = abs(sawnWoodHarvest[i] - wprod_sawnwood[countryprice].g(year))/wprod_sawnwood[countryprice].g(year);                 
					//cout<<"harvDiff["<<i<<"]"<<"\t"<<harvDiff[i]<<endl;
					//                 if (countriesNoFmcpol.find(i)!=countriesNoFmcpol.end()){harvDiff[i] = 0;}  
					cout<<"harvDiff["<<i<<"]="<<"\t"<<harvDiff[i]<<"\t sawnWoodHarvest["<<i<<"]="<<"\t"<<sawnWoodHarvest[i]<<"\t wprod_sawnwood["<<countryprice<<"]"<<year<<"=\t"<<wprod_sawnwood[countryprice].g(year)<<endl;
					if (countriesNoFmcpol.find(i)!=countriesNoFmcpol.end()){harvDiff[i] = 0; cout<<"No FM_Cpol for country =\t"<<i<<endl;}  
				}else{
					harvDiff[i] = 0;}
			}else{
				harvDiff[i] = 0;
			}
			//cout<<"harvDiff["<<i<<"]"<<"\t"<<harvDiff[i]<<endl;
		}
		maxDiff = *max_element(harvDiff+1,harvDiff+NumberOfCountries+1);
		for (int j=0;j<=NumberOfCountries;j++){if (harvDiff[j] == maxDiff){maxDiffCountry = j; 
		break;}}
	}else if (toAdjust.size()==1) {
		set<int>::iterator itt = toAdjust.begin();
		int ii = *itt;
		cout<<"ii=\t"<<ii<<"\t sawnWoodHarvest["<<ii<<"]="<<"\t"<<sawnWoodHarvest[ii]<<endl;          
		char countrych[4];
		int2str(ii,countrych);
		string countryprice = "re"+string(countrych)+"price0";          
		harvDiff[ii] = abs(sawnWoodHarvest[ii] - wprod_sawnwood[countryprice].g(year))/wprod_sawnwood[countryprice].g(year); 
		maxDiff = harvDiff[ii];
	}


}

#endif
