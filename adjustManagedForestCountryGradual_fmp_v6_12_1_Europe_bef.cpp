void adjustManagedForest(dataDetStruct &data_all, g4m::incrementTab *species[8], ageStructVector &cohort_all, 
						 ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
						 griddata &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
						 griddata &rotationForestNew, griddata &thinningForestNew, griddata2<char> &manageChForest,
						 griddata2<char> &rotationType, griddata &harvestSWGrid,griddata &harvestRWGrid, int year, griddata2<char> &unmanaged, double priceC) 


{
	bool harvControl = true; // Additional information to control output of the fm_cpol module
	bool NPV_postcontrol_0 = false; // Control of old forest management NPV at 0 C price: Use only for testing!!!!
	bool NPV_postcontrol = false; // Control of old forest management NPV at non-zero C price
	toAdjust.clear();
	toAdjust = countriesList;     // List of considered countries

	double sawnWoodHarvest[NumberOfCountries+1];		// array of harvested sawnwood in each country, m3
	double restWoodHarvest[NumberOfCountries+1];		// array of harvested restwood in each country, m3
	sawnWoodHarvest[0]=0.;
	restWoodHarvest[0]=0.;
	for (int i=1; i<=NumberOfCountries; i++){
		sawnWoodHarvest[i]=0.;
		restWoodHarvest[i]=0.;
		FMs[i]=0.;

	}
	harvestWood objHarvestWood;							// object of harvestWood class
	if (year == byear+1) {
		dataDetStruct::iterator iter = data_all.begin();
		iter = data_all.begin();
		while (iter != data_all.end()) {
			if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
				short int country = iter->COUNTRY[2000];
				if (countriesList.find(country) != countriesList.end()) { // Test only some countries                  
					if (iter->PROTECT[2000]==0) {
						int xi = (iter->x);
						int yi = (iter->y);
						if (thinningForest.get(xi,yi)<0) {
							unmanaged.set(xi,yi,1);
						}
					}
				}
			}
			iter++;
		}
	}
	//------------------------------------------------------------------------------
	//----------------------- Gradualy adjust thinning -----------------------------
	//------------------------------------------------------------------------------ 
	dataDetStruct::iterator iter = data_all.begin(); 
	while (iter != data_all.end())
	{
		int region = iter->POLESREG[2000];
		if (regions.find(region) != regions.end()) { // Test only some regions   
			short int country = iter->COUNTRY[2000];
			if (countriesList.find(country) != countriesList.end()) { // Test only some countries                  

				int asID = iter->asID;
				int xi = (iter->x);
				int yi = (iter->y);
				double thinningTmp = thinningForest.get(xi,yi);
				if (thinningTmp>1 && manageChForest.get(xi,yi)>0) {thinningTmp -= 0.025; 
				if (thinningTmp<1) thinningTmp = 1;
				thinningForest.set(xi,yi,thinningTmp); cohort_all[asID]->setStockingdegreeMin(thinningTmp*sdMinCoeff);  cohort_all[asID]->setStockingdegreeMax(thinningTmp*sdMaxCoeff);}
				double thinningTmpNew = thinningForestNew.get(xi,yi);
				if (thinningTmpNew>1) {thinningTmpNew -= 0.025; 
				if (thinningTmpNew<1) thinningTmpNew = 1;
				thinningForestNew.set(xi,yi,thinningTmpNew); newCohort_all[asID]->setStockingdegreeMin(thinningTmpNew*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(thinningTmpNew*sdMaxCoeff);} 
			} 
		}
		iter++; 
	}
	//------------------------------------------------------------------------------
	if (fmpol && priceC>0 && year>refYear)
	{
		double maxDiff=0;
		double fm_hurdle=1.;
		double diffMaxDiff = 1.;		
		double maxDiff0[NumberOfCountries+1];		
		countriesNoFmcpol.clear();
		if(year>byear){
			for (int i=1;i<NumberOfCountries+1;i++){
				if (CountryregMaxHarvest.get(i,year-1)<0.9*(CountryregWsawnprod.get(i,year-1)+CountryregWrestprod.get(i,year-1))
					|| CountryregWoodHarvestDfM3Year.get(i,year-1)>(1.05*CountryregWsawnprod.get(i,year-1)+CountryregWrestprod.get(i,year-1)))
				{countriesNoFmcpol.insert(i);}

			}
			set<int>::iterator it; for (it=countriesNoFmcpol.begin();it!=countriesNoFmcpol.end();it++){cout<<"countriesNoFmcpol=\t"<<*it<<endl;}
		}

		set<int> useChange;
		griddata thinningForestO = griddata(ResLongitude,ResLatitude,0.);
		griddata2<int> rotationForestO = griddata2<int>(ResLongitude,ResLatitude,0);
		griddata2<int> rotationForestNewO = griddata2<int>(ResLongitude,ResLatitude,0);

		fm_cpol(data_all, species, cohort_all,
			newCohort_all, dat_all, maiForest, 
			thinningForest, rotationForest, managedForest,
			rotationForestNew, thinningForestNew, manageChForest,
			rotationType, harvestSWGrid, harvestRWGrid, year, unmanaged, fm_hurdle,maxDiff, priceC,
			useChange, thinningForestO,rotationForestO,rotationForestNewO);
//		cout<<"maxDiff = "<<"\t"<<maxDiff<<"\t"<<"fm_hurdle = "<<"\t"<<fm_hurdle<<"\t"<<"diffMaxDiff = "<<"\t"<<diffMaxDiff<<endl;        

		for (int i=0;i<NumberOfCountries;i++){maxDiff0[i] = harvDiff[i];}
		doneList.clear();

		while (doneList.size()< countriesList.size())	{

			while (doneList.find(maxDiffCountry)!=doneList.end()){
				harvDiff[maxDiffCountry]=0;
				maxDiff = *max_element(harvDiff+1,harvDiff+NumberOfCountries+1);
				for (int j=0;j<=NumberOfCountries;j++){if (harvDiff[j] == maxDiff){maxDiffCountry = j; 
				break;}}
			}
			//    if (maxDiff<=0.05) break;
			if (maxDiff<=0.03) break;    

			toAdjust.clear();
			toAdjust.insert(maxDiffCountry);       
			int countryDone = maxDiffCountry;
			//    maxDiffi[0] = maxDiff0[maxDiffCountry];
			int i = 1;
			double maxDiffPrev = maxDiff0[maxDiffCountry];
			fm_hurdle=1;
			//        while (maxDiff > 0.05
			while (maxDiff > 0.03        
				//                               && fm_hurdle > -1 
					&& fm_hurdle > -2      // Try for Ireland 
					//                                          && diffMaxDiff >= 0
					&& (i <= 25 && fm_hurdle >= -1) // Try for Ireland 
					&& (i <= 50) // Try for Ireland 

					) { 

						if (diffMaxDiff <= 0.02){                                                              
							if (fm_hurdle<2.5 && fm_hurdle>0) {fm_hurdle+=0.3;}
							else if (fm_hurdle < 100000){fm_hurdle *= 100000.;}
							else if (fm_hurdle > 100000 && fm_hurdle < 300000){fm_hurdle *= 4.;}      // Try for Ireland
							else if (fm_hurdle == 1000000){fm_hurdle = -1.;}    // Try for Ireland        
							else if (fm_hurdle == -1) {fm_hurdle=-2.;}                      // Try for Ireland        
							else  {fm_hurdle=-2.;} 
						}     

						fm_cpol(data_all, species, cohort_all,
							newCohort_all, dat_all, maiForest, 
							thinningForest, rotationForest, managedForest,
							rotationForestNew, thinningForestNew, manageChForest,
							rotationType, harvestSWGrid, harvestRWGrid, year, unmanaged, fm_hurdle,maxDiff, priceC,
							useChange, thinningForestO,rotationForestO,rotationForestNewO);

						diffMaxDiff = maxDiffPrev - maxDiff;
						maxDiffPrev = maxDiff;

		//				cout<<"maxDiff = "<<"\t"<<maxDiff<<"\t"<<"fm_hurdle_loop = "<<"\t"<<fm_hurdle<<"\t"<<"diffMaxDiff = "<<"\t"<<diffMaxDiff<<endl;        
						i++;
			}

			doneList.insert(countryDone);
		} 
		cout<<"----------------------"<<endl;
		if (NPV_postcontrol){
			iter = data_all.begin();
			while (iter != data_all.end()) {
				short int region = iter->POLESREG[2000];
				if (regions.find(region) != regions.end()) { // Test only some regions         
					short int country = iter->COUNTRY[2000];
					if (countriesList.find(country) != countriesList.end()) {  
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

							coeff.PriceC.clear();
							coeff.PriceC.insert(0, priceC * iter->CORRUPTION[byear]);   
							double NPVtmp = 0.;
							double maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];
							g4m::ageStruct cohortTmp = *cohort_all[asID];
							short int used = 0;
							if (thinningForest.get(xi,yi)>0){used = 1;}
							NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,rotationForest.get(xi,yi),
								regprice,regprice0,dat_all[asID].OforestShare,used);    
	//						cout<<"asID=\t"<<asID<<"\tcountry=\t"<<country<<"\tNPVbau=\t"<<NPVbau[year-refYear-1][asID]<<"\tNPVcur=\t"<<NPVtmp<<endl;
						}
					}
				}
				iter++;
			}
		}// end for NPV_postcontrol
	}
	else
	{
		iter = data_all.begin();
		while (iter != data_all.end()) {
			short int region = iter->POLESREG[2000];
			if (regions.find(region) != regions.end()) { // Test only some regions         
				short int country = iter->COUNTRY[2000];
				if (countriesList.find(country) != countriesList.end()) { // Test only some countries       
					int asID = iter->asID;
					if (iter->PROTECT[2000]==0) {
						int xi = (iter->x);
						int yi = (iter->y);

						double harvestSWTmp = 0;
						double harvestRWTmp = 0;
						double newHarvestSWTmp = 0;
						double newHarvestRWTmp = 0; 

						if (thinningForest.get(xi,yi) >0)  
						{
							g4m::ageStruct cohortTmp = *cohort_all[asID];
							g4m::ageStruct cohortTmpNew = *newCohort_all[asID];
							harvestSWTmp = harvestSWGrid.get(xi,yi);
							harvestRWTmp = harvestRWGrid.get(xi,yi);
							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();               
							pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();                
							double realAreaO = cohortTmp.getArea();
							double realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							double harvestSO = objHarvestWood.harvestWS;
							double harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							double harvestSN = objHarvestWood.harvestWS;
							double harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
							sawnWoodHarvest[country] += (newHarvestSWTmp);
							restWoodHarvest[country] += (newHarvestRWTmp);

							harvestSWGrid.set(xi,yi,newHarvestSWTmp);
							harvestRWGrid.set(xi,yi,newHarvestRWTmp);
						//	if (year==2005)	{cout<<"1 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

						}else{
							sawnWoodHarvest[country] += dat_all[asID].DeforWoodM3.first;
							restWoodHarvest[country] += dat_all[asID].DeforWoodM3.second;
							harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first); 
							harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
							//if (year==2005)	{cout<<"2 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
						}                   
					}
				}
			}
			iter++;
		}

		//------------------------------------------------------------------------------
		// 
		// -------Zero Adjust thinning if population density changed --------------------
		//
		//------------------------------------------------------------------------------
		if ((year > 2000) && ((year) % 10 == 0))
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

							//  if (woodHarvest[country] < 0.9 * wprod[countryprice].g(year))
							//  if (woodHarvest[country] < 0.85 * wprod[countryprice].g(year)) 
//std:: cout << "adjust: wprod=\t"<< wprod_sawnwood[countryprice].g(year)<<endl;
							if (sawnWoodHarvest[country] < 0.88 * wprod_sawnwood[countryprice].g(year))   
							{
								//if (iter->COUNTRY[2000] == 197){cout<<"woodHarvest_0= "<<woodHarvest[region-1]<<"\t 0.9*woodHarvestStat= "<<0.9*wprod[countryprice].g(year) <<endl;}
								if (managedForest.get(xi,yi)<=0)
								{
									if ((iter->POPDENS.g(year) >0) && (iter->GDP.g(year) > 0))
									{
										double newHarvestSWTmp = 0.;
										double newHarvestRWTmp = 0.;
										int biomassRot = 0;
										int rotMAI = 0;
										int rotMaxBmTh = 0;
										int Rotation = 0;
										double SD = 1.5;

										double countryHarvestSWTmp = sawnWoodHarvest[country];
										double countryHarvestRWTmp = restWoodHarvest[country];

										g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
										g4m::ageStruct cohortTmp = *cohort_all[asID];
										double harvestSWTmp = harvestSWGrid.get(xi,yi);
										double harvestRWTmp = harvestRWGrid.get(xi,yi);              


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
										double thinningForestTmp = thinningForest.get(xi,yi);
										int rotationForestNewTmp = rotationForestNew.get(xi,yi);
										double thinningForestNewTmp = thinningForestNew.get(xi,yi);

										//---------------------------------------------------
										int newForAge = newCohort_all[asID]->getActiveAge();
										if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
										{
											if (managedForest.get(xi,yi) == 0)
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
											manageChForest.set(xi,yi,1);	
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
											if (managedForest.get(xi,yi) == 0)
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

										double realAreaO = cohortTmp.getArea();
										double realAreaN = cohortTmpNew.getArea();
										objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
										double harvestSO = objHarvestWood.harvestWS;
										double harvestRO = objHarvestWood.harvestWR;
										objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
										double harvestSN = objHarvestWood.harvestWS;
										double harvestRN = objHarvestWood.harvestWR;
										newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
											harvestSN * dat_all[asID].AforestShare) * 
											dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
										newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
											harvestRN * dat_all[asID].AforestShare) * 
											dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
										countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
										countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
									//	if (year==2005)	{cout<<"3 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

										if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
											){
												sawnWoodHarvest[country] = countryHarvestSWTmp;
												restWoodHarvest[country] = countryHarvestRWTmp;
												harvestSWGrid.set(xi,yi,newHarvestSWTmp);
												harvestRWGrid.set(xi,yi,newHarvestRWTmp);
												manageChForest.set(xi,yi,1);
								//				if (year==2005)	{cout<<"4 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
										}else{ // return old values
											rotationForest.set(xi,yi,rotationForestTmp);	
											thinningForest.set(xi,yi,thinningForestTmp);
											managedForest.set(xi,yi,managedForestTmp);	                   
											cohort_all[asID]->setU(rotationForestTmp);
											cohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   cohort_all[asID]->setStockingdegreeMax(thinningForestTmp*sdMaxCoeff);
											rotationForestNew.set(xi,yi,rotationForestNewTmp);	
											thinningForestNew.set(xi,yi,thinningForestNewTmp);
											newCohort_all[asID]->setU(rotationForestNewTmp);
											newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestNewTmp*sdMaxCoeff);                 
										}
									}  //End  if ((data["POPDENS"].g(1990) >0) && (data["GDP"].g(1990) > 0))
								} // End  if (managedForest(xi,yi)<=0)


							}else if (sawnWoodHarvest[country] > 1.12 * wprod_sawnwood[countryprice].g(year))  
							{
								if (managedForest.get(xi,yi)>0) 
								{
									if ((iter->POPDENS.g(year) == 0) && (iter->GDP.g(year) == 0))
									{
										//double sawnW = 0.;
										//double restW = 0.;
										//double sawnThW = 0.;
										//double restThW = 0.;
										//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
										//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
										//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

										//double sawnWnew = 0.;
										//double restWnew = 0.;
										//double sawnThWnew = 0.;
										//double restThWnew = 0.;
										//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
										//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
										//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

										double newHarvestSWTmp = 0.;
										double newHarvestRWTmp = 0.;
										int biomassRotTh = 0;
										int rotMAI = 0;
										int rotMaxBm = 0;
										int Rotation = 0;
										double countryHarvestSWTmp = sawnWoodHarvest[country];
										double countryHarvestRWTmp = restWoodHarvest[country];

										double harvestSWTmp = harvestSWGrid.get(xi,yi);
										double harvestRWTmp = harvestRWGrid.get(xi,yi);


										if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
										{
											biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (without thinning)  
											rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), 1);
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

										cohort_all[asID]->setU(Rotation);
										cohort_all[asID]->setStockingdegreeMin(-1);  cohort_all[asID]->setStockingdegreeMax(-1);
										rotationForestNew.set(xi,yi,Rotation);	
										thinningForestNew.set(xi,yi,-1.);	

										newCohort_all[asID]->setU(Rotation);
										newCohort_all[asID]->setStockingdegreeMin(-1);  newCohort_all[asID]->setStockingdegreeMax(-1);

										countryHarvestSWTmp = countryHarvestSWTmp - harvestSWTmp+dat_all[asID].DeforWoodM3.first;             
										countryHarvestRWTmp = countryHarvestRWTmp - harvestRWTmp+dat_all[asID].DeforWoodM3.second;
			//							if (year==2005)	{cout<<"5 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
										if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country]
										- wprod_sawnwood[countryprice].g(year))) 

										{
											sawnWoodHarvest[country] = countryHarvestSWTmp;  
											restWoodHarvest[country] = countryHarvestRWTmp;
											harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first);
											harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
				//							if (year==2005)	{cout<<"6 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
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
											newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestNewTmp*sdMaxCoeff);                 
										}


									} // End  if ((data["POPDENS"].g(1990) == 0) && (data["GDP"].g(1990) == 0))
								} // End  if ((managedForest(xi,yi)>0))
							} // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].g(year))   

						} //End for PROTECT == 0
					} // End for some countries 
				} // Test only some regions   
				iter++;
			} // End for WHILE (cell loop) 
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

						//double sawnW = 0.;
						//double restW = 0.;
						//double sawnThW = 0.;
						//double restThW = 0.;
						//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

						//double sawnWnew = 0.;
						//double restWnew = 0.;
						//double sawnThWnew = 0.;
						//double restThWnew = 0.;
						//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha


						int biomassRotTh = 0.;
						int rotMAI = 0.;
						int rotMaxBmTh = 0.;
						int Rotation = 0.;

						double newHarvestSWTmp = 0.;
						double newHarvestRWTmp = 0.;

						double countryHarvestSWTmp = sawnWoodHarvest[country];  
						double countryHarvestRWTmp = restWoodHarvest[country];


						if (sawnWoodHarvest[country] < 0.99 * wprod_sawnwood[countryprice].g(year))  
						{
							if (managedForest.get(xi,yi)>=2)
							{
								if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
								{
									//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
								}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
								{
									//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].prevPlantPhytHaBmGr, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
								}                       


								int rotationForestTmp = rotationForest.get(xi,yi);
								int rotationForestNewTmp = rotationForestNew.get(xi,yi);
								g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                                     
								g4m::ageStruct cohortTmp = *cohort_all[asID];
								int newForAge = newCohort_all[asID]->getActiveAge();         
								double harvestSWTmp = harvestSWGrid.get(xi,yi);
								double harvestRWTmp = harvestRWGrid.get(xi,yi);
								if (rotMAI < rotationForest.get(xi,yi)) // TEST for AT
								{									
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
									Rotation = rotationForestNew.get(xi,yi) - 5;
									if (Rotation < rotMAI) {Rotation = rotMAI;}
									cohortTmpNew.setU(Rotation);
									newCohort_all[asID]->setU(Rotation);    
									rotationForestNew.set(xi,yi,Rotation);
								}
								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmp = cohortTmp.aging();								
								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();								
								double realAreaO = cohortTmp.getArea();
								double realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
								double harvestSO = objHarvestWood.harvestWS;
								double harvestRO = objHarvestWood.harvestWR;
								objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
								double harvestSN = objHarvestWood.harvestWS;
								double harvestRN = objHarvestWood.harvestWR;
								newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestSN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
								newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestRN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
								countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
								countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
	//							if (year==2005)	{cout<<"7 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country]
								- wprod_sawnwood[countryprice].g(year))) 
								{
									sawnWoodHarvest[country] = countryHarvestSWTmp;
									restWoodHarvest[country] = countryHarvestRWTmp;
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
			//						if (year==2005)	{cout<<"8 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									cohort_all[asID]->setU(rotationForestTmp);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									newCohort_all[asID]->setU(rotationForestNewTmp);
								}  
								//        }

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

									double harvestSWTmp = harvestSWGrid.get(xi,yi);
									double harvestRWTmp = harvestRWGrid.get(xi,yi);

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

									double realAreaO = cohortTmp.getArea();
									double realAreaN = cohortTmpNew.getArea();
									objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
									double harvestSO = objHarvestWood.harvestWS;
									double harvestRO = objHarvestWood.harvestWR;
									objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
									double harvestSN = objHarvestWood.harvestWS;
									double harvestRN = objHarvestWood.harvestWR;
									newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestSN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
									newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestRN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
									countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
									countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
		//							if (year==2005)	{cout<<"9 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

									if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country]
									- wprod_sawnwood[countryprice].g(year)))
									{
										sawnWoodHarvest[country] = countryHarvestSWTmp;
										restWoodHarvest[country] = countryHarvestRWTmp;
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
	//									if (year==2005)	{cout<<"10 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
									}else{ // return old values
										rotationForest.set(xi,yi,rotationForestTmp);	
										cohort_all[asID]->setU(rotationForestTmp);
										rotationForestNew.set(xi,yi,rotationForestNewTmp);	
										newCohort_all[asID]->setU(rotationForestNewTmp);
									}  
								}

							}  // end for if (managedForest(xi,yi)<=2)       
						}   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].g(year))
					} //End Protect
				}  // End some countries
			}  // Test only some regions 
			iter++;
		} // End While

		// ----- End of First pass


		//----Second pass = adjust rotation time -------  

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

						//double sawnW = 0.;
						//double restW = 0.;
						//double sawnThW = 0.;
						//double restThW = 0.;
						//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
						//double sawnWnew = 0.;
						//double restWnew = 0.;
						//double sawnThWnew = 0.;
						//double restThWnew = 0.;
						//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

						int biomassRotTh = 0;
						int rotMAI = 0;
						int rotMaxBmTh = 0;
						int Rotation = 0;
						//  double harvestTmp = 0;
						double newHarvestSWTmp = 0;
						double newHarvestRWTmp = 0;
						//  int rotationForestTmp = 0;
						//  int rotationForestNewTmp = 0;
						double countryHarvestSWTmp = sawnWoodHarvest[country];  
						double countryHarvestRWTmp = restWoodHarvest[country];
						//  if (woodHarvest[country] < 0.9 * wprod[countryprice].g(year))
						//  if (woodHarvest[country] < 0.95 * wprod[countryprice].g(year))
						if (sawnWoodHarvest[country] < 0.98 * wprod_sawnwood[countryprice].g(year))  
						{
							if (managedForest.get(xi,yi)>0)
							{
								if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
								{
									//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
								}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
								{
									//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].prevPlantPhytHaBmGr, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
								}                       

								int rotationForestTmp = rotationForest.get(xi,yi);
								int rotationForestNewTmp = rotationForestNew.get(xi,yi);
								g4m::ageStruct cohortTmpNew = *newCohort_all[asID];                                    
								g4m::ageStruct cohortTmp = *cohort_all[asID];
								int newForAge = newCohort_all[asID]->getActiveAge();         
								double harvestSWTmp = harvestSWGrid.get(xi,yi);
								double harvestRWTmp = harvestRWGrid.get(xi,yi);
	//							if (year==2005)	{cout<<"11 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
								if (rotMAI < rotationForest.get(xi,yi)) // TEST for AT
								{
									//          Rotation = rotationForest.get(xi,yi) - 10;
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
								//sawnW = resTmp.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
								//restW = resTmp.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
								//sawnThW = resTmp.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
								//restThW = resTmp.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
								//bmH = resTmp.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//bmTh = resTmp.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
								pair<g4m::ageStruct::v, g4m::ageStruct::v>  resTmpNew = cohortTmpNew.aging();

								double realAreaO = cohortTmp.getArea();
								double realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
								double harvestSO = objHarvestWood.harvestWS;
								double harvestRO = objHarvestWood.harvestWR;
								objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
								double harvestSN = objHarvestWood.harvestWS;
								double harvestRN = objHarvestWood.harvestWR;
								newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestSN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
								newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestRN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
								countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
								countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
	//							if (year==2005)	{cout<<"12 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country]
								- wprod_sawnwood[countryprice].g(year)))
								{
									sawnWoodHarvest[country] = countryHarvestSWTmp;
									restWoodHarvest[country] = countryHarvestRWTmp;
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
			//						if (year==2005)	{cout<<"13 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									cohort_all[asID]->setU(rotationForestTmp);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									newCohort_all[asID]->setU(rotationForestNewTmp);
								}  
								//        }

							}  // end for if (managedForest(xi,yi)>=2)          

						} //end for if (woodHarvest[region] <= 0.9 * wprod[region].g(year))

						else if (sawnWoodHarvest[country] > 1.02 * wprod_sawnwood[countryprice].g(year)) 
						{
							//     if (managedForest.get(xi,yi)<=2)
							if (managedForest.get(xi,yi)>0)
							{
								if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
								{
									//                  biomassRotTh = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(dat_all[asID].ObiomassPrev, maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)  
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);         
								}else if (dat_all[asID].prevPlantPhytHaBmGr > 0)
								{
									rotMAI = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
									rotMaxBmTh = species[int(iter->SPECIESTYPE[byear])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);   
								}                       

								int rotationForestTmp = rotationForest.get(xi,yi);
								int rotationForestNewTmp = rotationForestNew.get(xi,yi);

								g4m::ageStruct cohortTmp = *cohort_all[asID];
								g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
								int newForAge = newCohort_all[asID]->getActiveAge();  

								double harvestSWTmp = harvestSWGrid.get(xi,yi);
								double harvestRWTmp = harvestRWGrid.get(xi,yi);
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

									double realAreaO = cohortTmp.getArea();
									double realAreaN = cohortTmpNew.getArea();
									objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
									double harvestSO = objHarvestWood.harvestWS;
									double harvestRO = objHarvestWood.harvestWR;
									objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
									double harvestSN = objHarvestWood.harvestWS;
									double harvestRN = objHarvestWood.harvestWR;
									newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestSN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
									newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestRN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
									countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
									countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
	//								if (year==2005)	{cout<<"14 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

									if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country]
									- wprod_sawnwood[countryprice].g(year)))
									{
										sawnWoodHarvest[country] = countryHarvestSWTmp;
										restWoodHarvest[country] = countryHarvestRWTmp;
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
	//									if (year==2005)	{cout<<"15 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
									}else{ // return old values
										rotationForest.set(xi,yi,rotationForestTmp);	
										cohort_all[asID]->setU(rotationForestTmp);
										rotationForestNew.set(xi,yi,rotationForestNewTmp);	
										newCohort_all[asID]->setU(rotationForestNewTmp);
									}  
								}

							}  // end for if (managedForest(xi,yi)<=2)       
						}   //end for esle if (woodHarvest[region] >= 1.1 * wprod[region].g(year))
					} //End Protect
				}  // End some countries
			}  // Test only some regions 
			iter++;
		} // End While
		// ----- End of Second pass


		////------------------------------------------------------------------------------
		//// 
		//// -------Third pass: Adjust thinning -----------------------------------------------
		////
		////------------------------------------------------------------------------------
		//cout<<"Start Third pass = Adjust thinning"<< endl;
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

							if (floor(sawnWoodHarvest[country]) < 0.97 * wprod_sawnwood[countryprice].g(year))  
							{
								if ((managedForest.get(xi,yi)<=0))
								{
									//double sawnW = 0.;
									//double restW = 0.;
									//double sawnThW = 0.;
									//double restThW = 0.;
									//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
									//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
									//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
									//double sawnWnew = 0.;
									//double restWnew = 0.;
									//double sawnThWnew = 0.;
									//double restThWnew = 0.;
									//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
									//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
									//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
									double newHarvestSWTmp = 0.;
									double newHarvestRWTmp = 0.;
									int biomassRot = 0;
									int biomassRotTh2 = 0;
									int rotMAI = 0;
									int rotMaxBmTh = 0;
									int Rotation = 0;

									double countryHarvestSWTmp = sawnWoodHarvest[country]; 
									double countryHarvestRWTmp = restWoodHarvest[country];
									double stockingDegree = 1.5;     
									g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
									g4m::ageStruct cohortTmp = *cohort_all[asID];
									double harvestSWTmp = harvestSWGrid.get(xi,yi);              
									double harvestRWTmp = harvestRWGrid.get(xi,yi);
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
									double thinningForestTmp = thinningForest.get(xi,yi);
									int rotationForestNewTmp = rotationForestNew.get(xi,yi);
									double thinningForestNewTmp = thinningForestNew.get(xi,yi);

									//---------------------------------------------------
									int newForAge = newCohort_all[asID]->getActiveAge();  
									if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
									{
										if (managedForest.get(xi,yi) == 0)
										{
											managedForest.set(xi,yi,3);
											Rotation = biomassRotTh2+1;
											if (Rotation < rotMAI) {Rotation = rotMAI;}
											rotationType.set(xi,yi,1);
										}else 
										{           managedForest.set(xi,yi,2);
										Rotation = biomassRotTh2+1;
										if (Rotation < rotMAI) {Rotation = rotMAI;}										
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
										cohortTmpNew.setStockingdegreeMin(stockingDegree*sdMinCoeff); cohortTmpNew.setStockingdegreeMax(stockingDegree*sdMaxCoeff);       
										newCohort_all[asID]->setU(Rotation);
										newCohort_all[asID]->setStockingdegreeMin(stockingDegree*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(stockingDegree*sdMaxCoeff);  

									}else{
										if (managedForest.get(xi,yi) == 0)
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

									double realAreaO = cohortTmp.getArea();
									double realAreaN = cohortTmpNew.getArea();
									objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
									double harvestSO = objHarvestWood.harvestWS;
									double harvestRO = objHarvestWood.harvestWR;
									objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
									double harvestSN = objHarvestWood.harvestWS;
									double harvestRN = objHarvestWood.harvestWR;
									newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestSN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
									newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestRN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
									countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
									countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
		//							if (year==2005)	{cout<<"16 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

									if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country]
									- wprod_sawnwood[countryprice].g(year)))
									{
										sawnWoodHarvest[country] = countryHarvestSWTmp;
										restWoodHarvest[country] = countryHarvestRWTmp;
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
		//								if (year==2005)	{cout<<"17 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
										manageChForest.set(xi,yi,1);
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
										newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestNewTmp*sdMaxCoeff);                 
									}
								} // End  if ((managedForest(xi,yi)<=0) & (managedForest(xi,yi)>-2)


							}else if (floor(sawnWoodHarvest[country]) > 1.03 * wprod_sawnwood[countryprice].g(year))   
							{
								if ((managedForest.get(xi,yi)>0) && (managedForest.get(xi,yi)<3))
								{
									//double sawnW = 0.;
									//double restW = 0.;
									//double sawnThW = 0.;
									//double restThW = 0.;
									//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
									//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
									//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
									//double sawnWnew = 0.;
									//double restWnew = 0.;
									//double sawnThWnew = 0.;
									//double restThWnew = 0.;
									//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
									//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
									//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha																    
									double newHarvestSWTmp = 0;
									double newHarvestRWTmp = 0;
									int biomassRot = 0;
									int rotMAI = 0;
									int rotMaxBm = 0;
									int Rotation = 0;
									double countryHarvestSWTmp = sawnWoodHarvest[country];		
									double countryHarvestRWTmp = restWoodHarvest[country];
									double harvestSWTmp = harvestSWGrid.get(xi,yi);              
									double harvestRWTmp = harvestRWGrid.get(xi,yi); 
									if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
									{

										rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
									}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
									{

										rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
									}     

									int rotationForestTmp = rotationForest.get(xi,yi);
									int managedForestTmp = managedForest.get(xi,yi);
									double thinningForestTmp = thinningForest.get(xi,yi);
									int rotationForestNewTmp = rotationForestNew.get(xi,yi);
									double thinningForestNewTmp = thinningForestNew.get(xi,yi);

									//---------------------------------------------------

									if (managedForest.get(xi,yi) == 2)
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
									cohort_all[asID]->setU(Rotation);
									cohort_all[asID]->setStockingdegreeMin(-1);  cohort_all[asID]->setStockingdegreeMax(-1);

									rotationForestNew.set(xi,yi,Rotation);	
									thinningForestNew.set(xi,yi,-1.);	

									newCohort_all[asID]->setU(Rotation);
									newCohort_all[asID]->setStockingdegreeMin(-1.);  newCohort_all[asID]->setStockingdegreeMax(-1.); 

									countryHarvestSWTmp = countryHarvestSWTmp-harvestSWTmp+dat_all[asID].DeforWoodM3.first;                
									countryHarvestRWTmp = countryHarvestRWTmp-harvestRWTmp+dat_all[asID].DeforWoodM3.second;
	//								if (year==2005)	{cout<<"18 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

									if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
										){
											sawnWoodHarvest[country] = countryHarvestSWTmp;     
											restWoodHarvest[country] = countryHarvestRWTmp;  
											harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first);
											harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
	//										if (year==2005)	{cout<<"19 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
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
										newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestNewTmp*sdMaxCoeff);                 
									}
									//         } // End if ((cohort_all[asID]->getArea(0) == 0.)||(cohort_all[asID]->getArea(0) >= 1/400.)){
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
								//double sawnW = 0.;
								//double restW = 0.;
								//double sawnThW = 0.;
								//double restThW = 0.;
								//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
								//double sawnWnew = 0.;
								//double restWnew = 0.;
								//double sawnThWnew = 0.;
								//double restThWnew = 0.;
								//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
								//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
								//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

								double newHarvestSWTmp = 0;
								double newHarvestRWTmp = 0;
								int biomassRot = 0;

								int rotMAI = 0;
								int rotMaxBm = 0;
								int Rotation = 0;

								double countryHarvestSWTmp = sawnWoodHarvest[country];	
								double countryHarvestRWTmp = restWoodHarvest[country];	
								double harvestSWTmp = harvestSWGrid.get(xi,yi);              
								double harvestRWTmp = harvestRWGrid.get(xi,yi);
								if (dat_all[asID].ObiomassPrev >0 && iter->CABOVEHA[byear] > 0 && maiForest.get(xi,yi) > 0)  
								{

									rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
								}else if (dat_all[asID].prevPlantPhytHaBmGr >0)  
								{

									rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);  
								}     
								int rotationForestTmp = rotationForest.get(xi,yi);
								int managedForestTmp = managedForest.get(xi,yi);
								double thinningForestTmp = thinningForest.get(xi,yi);
								int rotationForestNewTmp = rotationForestNew.get(xi,yi);
								double thinningForestNewTmp = thinningForestNew.get(xi,yi);

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
									cohort_all[asID]->setU(Rotation);
									cohort_all[asID]->setStockingdegreeMin(-1);  cohort_all[asID]->setStockingdegreeMax(-1);

									rotationForestNew.set(xi,yi,Rotation);	
									thinningForestNew.set(xi,yi,-1.);	

									newCohort_all[asID]->setU(Rotation);
									newCohort_all[asID]->setStockingdegreeMin(-1.);  newCohort_all[asID]->setStockingdegreeMax(-1.);  

								} // End    else // New forest age < rotation -> don't change FM for the new forest

								//                countryHarvestTmp += (newHarvestTmp-harvestTmp);
								countryHarvestSWTmp = countryHarvestSWTmp-harvestSWTmp+dat_all[asID].DeforWoodM3.first;                
								countryHarvestRWTmp = countryHarvestRWTmp-harvestRWTmp+dat_all[asID].DeforWoodM3.second;
	//							if (year==2005)	{cout<<"20 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}

								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
									){
										sawnWoodHarvest[country] = countryHarvestSWTmp;     
										restWoodHarvest[country] = countryHarvestRWTmp;  
										harvestSWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.first);
										harvestRWGrid.set(xi,yi,dat_all[asID].DeforWoodM3.second);
									manageChForest.set(xi,yi,1);
	//									if (year==2005)	{cout<<"21 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
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
									newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestNewTmp*sdMaxCoeff);                 
								}
								//         } // End if ((cohort_all[asID]->getArea(0) == 0.)||(cohort_all[asID]->getArea(0) >= 1/400.)){

							} // End  if ((managedForest(xi,yi)>0) & (managedForest(xi,yi)<3)
						} // End   else if (woodHarvest[region-1] >= (1.1 * wprod[region].g(year))   

					} //End for PROTECT == 0
				} //End some countries
			} // Test only some regions   
			iter++;
		} // End for WHILE (cell loop) 

		//cout << "Third pass is finished"<< endl;

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


						//double sawnW = 0.;
						//double restW = 0.;
						//double sawnThW = 0.;
						//double restThW = 0.;
						//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha
						//double sawnWnew = 0.;
						//double restWnew = 0.;
						//double sawnThWnew = 0.;
						//double restThWnew = 0.;
						//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

						double harvSWoodTmp = 0.;
						double harvRWoodTmp = 0.;
						int biomassRotTh = 0;
						int biomassRot = 0;  
						int rotMAI = 0;
						int rotMaxBmTh = 0;
						int rotMaxBm = 0;  
						int Rotation = 0;
						double newHarvestSWTmp = 0;	
						double newHarvestRWTmp = 0;	
						double countryHarvestSWTmp = sawnWoodHarvest[country];
						double countryHarvestRWTmp = restWoodHarvest[country];

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
							double harvestSWTmp = harvestSWGrid.get(xi,yi);
							double harvestRWTmp = harvestRWGrid.get(xi,yi);

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

							double realAreaO = cohortTmp.getArea();
							double realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							double harvestSO = objHarvestWood.harvestWS;
							double harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							double harvestSN = objHarvestWood.harvestWS;
							double harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3

							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
	//						if (year==2005)	{cout<<"22 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp;  
									restWoodHarvest[country] = countryHarvestRWTmp; 
			//						if (year==2005)	{cout<<"23 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								cohort_all[asID]->setU(rotationForestTmp);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								newCohort_all[asID]->setU(rotationForestNewTmp);
							}        
							//        }
							//        } //End if ((cohort_all[asID]->getArea(0) == 0.)||(cohort_all[asID]->getArea(0) >= 1/400.)){
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
							double harvestSWTmp = harvestSWGrid.get(xi,yi);
							double harvestRWTmp = harvestRWGrid.get(xi,yi);

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

								double realAreaO = cohortTmp.getArea();
								double realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
								double harvestSO = objHarvestWood.harvestWS;
								double harvestRO = objHarvestWood.harvestWR;
								objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
								double harvestSN = objHarvestWood.harvestWS;
								double harvestRN = objHarvestWood.harvestWR;
								newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestSN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
								newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
									harvestRN * dat_all[asID].AforestShare) * 
									dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3

								countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
								countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
//								if (year==2005)	{cout<<"24 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
								if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
									){
										sawnWoodHarvest[country] = countryHarvestSWTmp;  
										restWoodHarvest[country] = countryHarvestRWTmp; 
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);
	//									if (year==2005)	{cout<<"25 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
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
		// 
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

						//double sawnW = 0.;
						//double restW = 0.;
						//double sawnThW = 0.;
						//double restThW = 0.;
						//double bmH = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvRes = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

						//double sawnWnew = 0.;
						//double restWnew = 0.;
						//double sawnThWnew = 0.;
						//double restThWnew = 0.;
						//double bmHnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
						//double bmThnew = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
						//double harvResnew = 0.; // MG: usable harvest residues for the set (old) forest tC/ha

						double harvWoodTmp = 0.;
						int biomassRotTh = 0;
						int biomassRot = 0;  
						int rotMAI = 0;
						int rotMaxBmTh = 0;
						int rotMaxBm = 0;  
						int Rotation = 0;
						//  double harvestTmp = 0;
						double newHarvestSWTmp = 0;
						double newHarvestRWTmp = 0;
						//  int rotationForestTmp = 0;
						//  int rotationForestNewTmp = 0;
						double countryHarvestSWTmp = sawnWoodHarvest[country];
						double countryHarvestRWTmp = restWoodHarvest[country];
						//  if (floor(woodHarvest[country]) < 0.9 * wprod[countryprice].g(year))
						//  if (floor(woodHarvest[country]) < 0.95 * wprod[countryprice].g(year))
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
							double harvestSWTmp = harvestSWGrid.get(xi,yi);
							double harvestRWTmp = harvestRWGrid.get(xi,yi);
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

							double realAreaO = cohortTmp.getArea();
							double realAreaN = cohortTmpNew.getArea();
							objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							double harvestSO = objHarvestWood.harvestWS;
							double harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							double harvestSN = objHarvestWood.harvestWS;
							double harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3

							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
	//						if (year==2005)	{cout<<"26 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp;  
									restWoodHarvest[country] = countryHarvestRWTmp; 
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
	//								if (year==2005)	{cout<<"27 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
							}else{ // return old values
								rotationForest.set(xi,yi,rotationForestTmp);	
								cohort_all[asID]->setU(rotationForestTmp);
								rotationForestNew.set(xi,yi,rotationForestNewTmp);	
								newCohort_all[asID]->setU(rotationForestNewTmp);
							}        
							//        }
						}  // end for if (managedForest(xi,yi)>=2)          

						} //end for if (woodHarvest[region] <= 0.9 * wprod[region].g(year))

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
							double harvestSWTmp = harvestSWGrid.get(xi,yi);
							double harvestRWTmp = harvestRWGrid.get(xi,yi);
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

								double realAreaO = cohortTmp.getArea();
								double realAreaN = cohortTmpNew.getArea();
								objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
							double harvestSO = objHarvestWood.harvestWS;
							double harvestRO = objHarvestWood.harvestWR;
							objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
							double harvestSN = objHarvestWood.harvestWS;
							double harvestRN = objHarvestWood.harvestWR;
							newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestSN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
							newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
								harvestRN * dat_all[asID].AforestShare) * 
								dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3

							countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
							countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);
						//	if (year==2005)	{cout<<"28 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
							if (abs(countryHarvestSWTmp - wprod_sawnwood[countryprice].g(year))<(1+tolerance)*abs(sawnWoodHarvest[country] - wprod_sawnwood[countryprice].g(year)) 
								){
									sawnWoodHarvest[country] = countryHarvestSWTmp;  
									restWoodHarvest[country] = countryHarvestRWTmp; 
									harvestSWGrid.set(xi,yi,newHarvestSWTmp);
									harvestRWGrid.set(xi,yi,newHarvestRWTmp);
							//		if (year==2005)	{cout<<"29 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
								}else{ // return old values
									rotationForest.set(xi,yi,rotationForestTmp);	
									cohort_all[asID]->setU(rotationForestTmp);
									rotationForestNew.set(xi,yi,rotationForestNewTmp);	
									newCohort_all[asID]->setU(rotationForestNewTmp);
								}        
							}
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

							if (floor(restWoodHarvest[country]) < 0.99 * wprod_restwood[countryprice].g(year))  
							{
								if ((managedForest.get(xi,yi)<=0))
								{									
									double newHarvestSWTmp = 0.;
									double newHarvestRWTmp = 0.;
									double addcountryHarvestRWTmp=0.;
									int biomassRot = 0;
									int biomassRotTh2 = 0;
									int rotMAI = 0;
									int rotMaxBmTh = 0;
									int Rotation = 0;

									double countryHarvestSWTmp = sawnWoodHarvest[country]; 
									double countryHarvestRWTmp = restWoodHarvest[country];
									double stockingDegree = 1.5;     
									g4m::ageStruct cohortTmpNew = *newCohort_all[asID]; 
									g4m::ageStruct cohortTmp = *cohort_all[asID];
									double harvestSWTmp = harvestSWGrid.get(xi,yi);              
									double harvestRWTmp = harvestRWGrid.get(xi,yi);
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
									double thinningForestTmp = thinningForest.get(xi,yi);
									int rotationForestNewTmp = rotationForestNew.get(xi,yi);
									double thinningForestNewTmp = thinningForestNew.get(xi,yi);

									//---------------------------------------------------
									int newForAge = newCohort_all[asID]->getActiveAge();  
									if (newForAge > biomassRot)  // New forest age > rotation -> change FM for the new forest
									{
										if (managedForest.get(xi,yi) == 0)
										{
											managedForest.set(xi,yi,3);
											Rotation = biomassRotTh2+1;
											if (Rotation < rotMAI) {Rotation = rotMAI;}
											rotationType.set(xi,yi,1);
										}else 
										{           managedForest.set(xi,yi,2);
										Rotation = biomassRotTh2+1;
										if (Rotation < rotMAI) {Rotation = rotMAI;}										
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
										cohortTmpNew.setStockingdegreeMin(stockingDegree*sdMinCoeff); cohortTmpNew.setStockingdegreeMax(stockingDegree*sdMaxCoeff);       
										newCohort_all[asID]->setU(Rotation);
										newCohort_all[asID]->setStockingdegreeMin(stockingDegree*sdMinCoeff);  newCohort_all[asID]->setStockingdegreeMax(stockingDegree*sdMaxCoeff);  

									}else{
										if (managedForest.get(xi,yi) == 0)
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

									double realAreaO = cohortTmp.getArea();
									double realAreaN = cohortTmpNew.getArea();
									objHarvestWood.calculateHarvest(realAreaO, resTmp, modTimeStep);
									double harvestSO = objHarvestWood.harvestWS;
									double harvestRO = objHarvestWood.harvestWR;
									objHarvestWood.calculateHarvest (realAreaN, resTmpNew, modTimeStep);
									double harvestSN = objHarvestWood.harvestWS;
									double harvestRN = objHarvestWood.harvestWR;
									newHarvestSWTmp = (harvestSO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestSN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.first; // Total current harvested wood in the cell, m3
									newHarvestRWTmp = (harvestRO * (dat_all[asID].OforestShare-dat_all[asID].deforestShare) +
										harvestRN * dat_all[asID].AforestShare) * 
										dat_all[asID].LandAreaHa * iter->FTIMBER[byear] + dat_all[asID].DeforWoodM3.second; // Total current harvested wood in the cell, m3
									

									countryHarvestRWTmp += (newHarvestRWTmp-harvestRWTmp);									
									if (countryHarvestSWTmp > wprod_sawnwood[countryprice].g(year)){
									addcountryHarvestRWTmp = newHarvestSWTmp-harvestSWTmp;
									countryHarvestRWTmp += (addcountryHarvestRWTmp);									
									}else
									{
									countryHarvestSWTmp += (newHarvestSWTmp-harvestSWTmp);
									}
								//	if (year==2010)	{cout<<"30 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
									//if (year==2010)	{cout<<"31 restWoodHarvest  "<<restWoodHarvest[224]<<endl;}
									if (abs(countryHarvestRWTmp - wprod_restwood[countryprice].g(year))<(1+tolerance)*abs(restWoodHarvest[country]
									- wprod_restwood[countryprice].g(year)))
									{
										sawnWoodHarvest[country] = countryHarvestSWTmp;
										restWoodHarvest[country] = countryHarvestRWTmp;
										harvestSWGrid.set(xi,yi,newHarvestSWTmp);
										harvestRWGrid.set(xi,yi,newHarvestRWTmp);										
										manageChForest.set(xi,yi,1);
									if (year==2010)	{cout<<"30 sawnWoodHarvest  "<<sawnWoodHarvest[224]<<endl;}
									if (year==2010)	{cout<<"31 restWoodHarvest  "<<restWoodHarvest[224]<<endl;}
										
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
										newCohort_all[asID]->setStockingdegreeMin(thinningForestTmp*sdMinCoeff);   newCohort_all[asID]->setStockingdegreeMax(thinningForestNewTmp*sdMaxCoeff);                 
									}
								} // End  if ((managedForest(xi,yi)<=0) & (managedForest(xi,yi)>-2)
							}

						} //End for mai>0
					} //End for PROTECT == 0
				}  // End some countries
			} // Test only some regions   
			iter++;
		} // End for WHILE (cell loop) 


		if (harvControl){
			for (int i=1; i<=NumberOfCountries; i++){
				char countrych[4];
				int ii = i;
				int2str(ii,countrych);
				string countryprice = "re"+string(countrych)+"price0";
				if (wprod_sawnwood[countryprice].g(year) > 0){

					if ((regions.find(countryRegion[i]) != regions.end()) 
						&& (countriesList.find(i) != countriesList.end())  // Test only some countries         
						//                         && i!=114  // Luxenbourg is too small (1 cell) to be considered correctly
						)
					{                                  
						harvDiff[i] = abs(sawnWoodHarvest[i] - wprod_sawnwood[countryprice].g(year))/wprod_sawnwood[countryprice].g(year);                 
						//cout<<"harvDiff["<<i<<"]"<<"\t"<<harvDiff[i]<<endl;
				//		cout<<"harvDiff["<<i<<"]="<<"\t"<<harvDiff[i]<<"\t sawnWoodHarvest["<<i<<"]="<<"\t"<<sawnWoodHarvest[i]<<"\t wprod["<<countryprice<<"]"<<year<<"=\t"<<wprod_sawnwood[countryprice].g(year)<<endl;
					}else{
						harvDiff[i] = 0;}
				}else{
					harvDiff[i] = 0;
				}
				//cout<<"harvDiff["<<i<<"]"<<"\t"<<harvDiff[i]<<endl;
			}
		}

		if (NPV_postcontrol_0 && fmpol && year>refYear){
			iter = data_all.begin();
			while (iter != data_all.end()) {
				short int region = iter->POLESREG[2000];
				if (regions.find(region) != regions.end()) { // Test only some regions         
					short int country = iter->COUNTRY[2000];
					if (countriesList.find(country) != countriesList.end()) {  
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

							coeff.PriceC.clear();
							coeff.PriceC.insert(0, priceC * iter->CORRUPTION[byear]);   
							double NPVtmp = 0.;
							double maiV = maiForest.get(xi,yi)*iter->FTIMBER[byear];
							g4m::ageStruct cohortTmp = *cohort_all[asID];
							short int used = 0;
							if (thinningForest.get(xi,yi)>0){used = 1;}
							//cout<<"before NPV calc"<<endl;
							NPVtmp = npv_calc(*iter, cohortTmp, maiV,year,rotationForest.get(xi,yi),
								regprice,regprice0,dat_all[asID].OforestShare,used)
								//                                  *dat_all[asID].LandAreaHa * dat_all[asID].OforestShare
								;    

							//profit = ((sawnW + restW)* it.FTIMBER[byear] * TimberPrice - decision.plantingCosts() * cohort.getArea(0)) * OforestShare * singleCell.LandAreaHa; // profit only from harvesting old forest

				//			cout<<"asID=\t"<<asID<<"\tcountry=\t"<<country<<"\tNPVbau=\t"<<NPVbau[year-refYear-1][asID]<<"\tNPVcur=\t"<<NPVtmp<<endl;
						}
					}
				}
				iter++;
			}
		}// end for NPV_postcontrol_0

		//*****************************************************************************
		// Adjusting unmanaged forest rotation length to match reported FM sink
		//*****************************************************************************

		if (adjustFMsink){
			if (year == 2000){
				iter = data_all.begin();
				while (iter != data_all.end()) {
					if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
						short int country = iter->COUNTRY[2000];
						if (countriesList.find(country) != countriesList.end()) { // Test only some countries      
							if (iter->PROTECT[2000]==0) {
								int asID = iter->asID;
								FMs[country-1]+=(cohort_all[asID]->getBm() - dat_all[asID].Obiomass0) *
									(iter->FOREST[2000])*(dat_all[asID].LandAreaHa)*3.6666666667/10000.;  // FM sink averaged over 5 years, GgCO2/yr
							}
						} 
					}
					iter++;
				}

				iter = data_all.begin();
				while (iter != data_all.end()) {
					if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
						short int country = iter->COUNTRY[2000];
						if (countriesList.find(country) != countriesList.end()) { // Test only some countries      
							if (iter->PROTECT[2000]==0) {
								int asID = iter->asID;
								int xi = (iter->x);
								int yi = (iter->y);
								//cout<<"asID= \t"<<asID<<"\t country= \t"<<country<<"\t FMs= \t"<<FMs[country-1]<<"\t Fmstat= \t"<<FM_sink_stat[country-1]<<"\t unmanaged= \t"<<int(unmanaged.get(xi,yi))<<endl;
								if (unmanaged.get(xi,yi)>0){
									int biomassRot = 0;  
									int rotMaxBm = 0;  
									int Rotation = 0;
									if (FMs[country-1] < 0.85*FM_sink_stat[country-1]){
										if (iter->FOREST[1990] >0 && iter->CABOVEHA[1990] > 0 && maiForest.get(xi,yi) > 0)  
										{  
											//                     biomassRot = species[int(iter->SPECIESTYPE[byear])-1]->gU(iter->CABOVEHA[1990], maiForest.get(xi,yi), 1);     // rotation time to get current biomass (without thinning)  
											//                     biomassRot = species[int(iter->SPECIESTYPE[byear])-1]->gU(dat_all[asID].Obiomass0, maiForest.get(xi,yi), 1);     // rotation time to get current biomass (without thinning)  
											rotMaxBm = species[int(iter->SPECIESTYPE[byear])-1]->gTopt(maiForest.get(xi,yi), 2);         

										}     
										//cout<<"rotationForest= \t"<<rotationForest.get(xi,yi)<<"\t rotMaxBm= \t"<< rotMaxBm<<endl;
										if (rotationForest.get(xi,yi) < rotMaxBm)
										{
											float BmTmp = 0.;
											float BmTmpBef = 0.;

											g4m::ageStruct cohortTmp = *cohort_all[asID];
											BmTmp = cohortTmp.getBm();
											BmTmpBef = BmTmp * BEF(int(iter->SPECIESTYPE[byear])-1,BmTmp,iter->FTIMBER[2000]);
											cohort_all[asID]->setU(rotMaxBm);
											rotationForest.set(xi,yi,rotMaxBm);

											cohortTmp.setU(rotMaxBm);
											cohortTmp.aging();
											cohortTmp.aging();
											cohortTmp.aging();
											cohortTmp.aging();                                       
											cohortTmp.aging();   

											cohortTmp.aging();
											cohortTmp.aging();
											cohortTmp.aging();
											cohortTmp.aging();                                       
											cohortTmp.aging();                                 
											//                     FMs[country-1]+=((cohortTmp.getBm()-BmTmp) - (BmTmp-dat_all[asID].Obiomass0)) *
											//                          (iter->FOREST[2000]) * (dat_all[asID].LandAreaHa) * 3.6666666667/10000.;
											FMs[country-1]+=((cohortTmp.getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohortTmp.getBm(),iter->FTIMBER[2000])-BmTmpBef) 
												- (BmTmpBef-dat_all[asID].Obiomass0 * BEF(int(iter->SPECIESTYPE[byear])-1,dat_all[asID].Obiomass0,iter->FTIMBER[2000]))) 
												* (iter->FOREST[2000]) * (dat_all[asID].LandAreaHa) * 3.6666666667/10000.;
										}
									}
								}
							} // Protect
						} // Test some countries
					} // Test some regions
					iter++;
				}// While
			} // if year ==1995
		} //if (adjustFMsink)         
	}
}
