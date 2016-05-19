#include <iostream>
void initManagedForest(dataDetStruct &data_all, g4m::incrementTab *species[8], 
	                   g4m::ffipolm<double> *ffcov, 
                       g4m::ffipolm<double> *ffcoe, 
                       g4m::ffipolm<bool> *ffdov, 
                       g4m::ffipolm<bool> *ffdoe, 
					   datGlobal &dat_all,
                       ageStructVector &cohort_all, ageStructVector &newCohort_all,
                       griddata &maiForest, griddata &thinningForest,
                       griddata2<char> &rotationType, griddata2<char> &managedForest,
                       griddata &rotationForest, griddata &harvestGrid, griddata2<float> &OforestShGrid) 

 {
  double woodHarvest[NumberOfCountries+1];
  double woodLost[NumberOfCountries+1];
//  double woodHarvestStat[NumberOfCountries+1];
//  int managedCount[NumberOfCountries+1];

  for (int i=0; i<NumberOfCountries+1; i++){
	  woodHarvest[i]=0.; 
//cout<<"i=\t"<<i<<"\twoodHarvest=\t"<<woodHarvest[i]<<"\tsize=\t"<<endl;
    woodLost[i]=0.;
//    woodHarvestStat[i]=0.;
//    managedCount[i]=0;
  }

//  int year =1990;
  int year =2000;

  double sawnW = 0.;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
  double restW = 0.;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
  double sawnThW = 0.;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
  double restThW = 0.;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
  double bmH = 0.;        // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
  double bmTh = 0.;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning
  double harvRes = 0.;    // MG: usable harvest residues for the set (old) forest tC/ha
  
////  double defIncome = 0.;
//    //Get optimal rotation time
//    //0 .. Highest average increment
//    //1 . .Maximum average Biomass
//    //2 .. Highest possible age
//    //3 .. Maximum harvest at final cut
//    //4 .. Average Maximum harvest at final cut
//          int optimUnmanaged = 0;
//          int optimMAI = 0;
//          int optimMaxBm = 1;
////          int optimMaxBmTh = 3;
//          int optimHarvFin = 3;
//          int optimHarvAve = 4;
//          int rotUnmanaged = 0;
//          int rotMAI = 0;
//          int rotMaxBm = 0;
//          int rotMaxBmTh = 0;
//          int rotHarvFin = 0;
//          int rotHarvAve = 0;

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
        //double forestShare0 = 0;
        //if (iter->FOREST[year]+(iter->CROP[year])+(iter->BUILTUP[year])>1)
        //   {forestShare0 = (1-(iter->CROP[year]+iter->BUILTUP[year]));}
        //   else {forestShare0 = iter->FOREST[year];}

		double forestShare0 = iter->FOREST[byear];
		if (forestShare0<0) forestShare0=0;
//		double maxfor = 1. - (iter->BUILTUP[year] + iter->CROP[year]);
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
        double harvWood=0.; //MG: harvestable wood, m3
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

//		if (Country0==129) MAI *= 0.9168; // adjustment of Luxembourg MAI, 2 May 2013

        maiForest.set(xi,yi,MAI);
        
        float harvMAI = MAI*iter->FTIMBER[year]*(1-coeff.HarvLoos[year]);


 //       double forFlag = 0.;          //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest

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
//               , Rotation
//               , harvWood );


    int Rotation = 0;

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


//
//    initLoop(0, data_all, species, cohort_all, newCohort_all, dat_all, maiForest, thinningForest, rotationForest);
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
//if (asID>100) cout<<"asID1=\t"<<asID<<endl;
			short int region = iter->COUNTRY[2000];
//cout<<"Country0=\t"<<region<<endl;//"\tmanagedCount=\t"<<managedCount[Country0-1]<<endl;
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
     
           double forestArea0 = LandAreaHa * forestShare0;
          float harvMAI = maiForest.get(xi,yi)*iter->FTIMBER[year]*(1-coeff.HarvLoos[year]);
          double harvWood=0; //MG: harvestable wood, m3
          double defIncome = 0.;

          g4m::ageStruct cohortTmp = *cohort_all[asID];
//cout<<"asID=\t"<<asID<<endl;
          cohortTmp.setU(rotationForest.get(xi,yi));
//          g4m::ageStruct::v res;  // MG: results vector for the set (old) forest 
          pair<g4m::ageStruct::v, g4m::ageStruct::v> res;
          int Rotation = 0;
//          double forFlag = 0.;    //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
          
        if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
          biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], maiForest.get(xi,yi));        // rotation time to get current biomass (without thinning)
          biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUSdTab(iter->CABOVEHA[year], maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
            if (thinningForest.get(xi,yi) > 0) {
              Rotation = biomasRotTh+1; 
            }  else {
              Rotation = biomasRot+1;
            }
//          if (iter->FOREST"].g(1990) >0) forFlag = 1.0;
        }
        if (maiForest.get(xi,yi)> 0){
          rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
//          rotMaxBm = species[int(iter->SPECIESTYPE[year])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm);                        
          rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
        }             
//cout<<"asID2=\t"<<asID<<endl;
          res = cohortTmp.aging();
          //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
          //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
          //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
          //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;

  		  double realAreaO = cohortTmp.getArea();
		  double harvestO = cohortRes(realAreaO,res);
          harvWood = harvestO * iter->FTIMBER[byear];


          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
//cout<<"asID3=\t"<<asID<<endl;
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
//               , Rotation
//               , harvWood );
//cout<<"asID4=\t"<<asID<<endl;
			   abBiomassO = cohortTmp.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * iter->FTIMBER[year] * (1. -coeff.HarvLoos[year]));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * iter->FTIMBER[year] 
		         * (1. -coeff.HarvLoos[year]));
	      defIncome = pDefIncome * (1. - iter->SLASHBURN[year])
		            + sDefIncome * iter->SLASHBURN[year];
		            
          harvestGrid.set(xi,yi,harvWood * forestArea0);   
          if (managedForest.get(xi,yi) > 0) {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}
            if (woodHarvest[Country0-1] < 0.95 * wprod[regprice].g(year)) {
              if (managedForest.get(xi,yi) >= 2) {
//	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,rotMAI);		   
                cohortTmp.setU(rotMAI);
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;

  		         double realAreaO = cohortTmp.getArea();
				 double harvestO = cohortRes(realAreaO,res);
				 harvWood = harvestO * iter->FTIMBER[byear];

                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
			    harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += harvWood * forestArea0;
              }            
            } else if (woodHarvest[Country0-1] > 1.05 * wprod[regprice].g(year)) {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}                   

                if (rotationForest.get(xi,yi) < rotMaxBmTh) {
//                HarvestTmp = harvestGrid.get(xi,yi);
			    rotationForest.set(xi,yi, rotMaxBmTh);
                cohortTmp.setU(rotMaxBmTh);
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
		 	    double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];

                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
			    harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += harvWood * forestArea0;
              }
            } else { // keep biomass rotation
                //res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  				double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];


                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
			    harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += harvWood * forestArea0;
            } 

//            managedCount[Country0-1] +=1;                  
          }else { // unmanaged forests
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  				double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];


                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;

                woodLost[Country0-1] += harvWood * forestArea0;
          }          

//if (iter->asID>100) cout<<"asID2=\t"<<iter->asID<<endl;
    } //else {iter++;}  // end for IF unprotected
//if (iter->asID>100) cout<<"asID3=\t"<<iter->asID<<endl;
   }   // Test only some regions
//if (iter->asID>100) 
//	cout<<"asID4=\t"<<iter->asID<<endl;
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

          double HarvestTmp = 0;
          double newHarvestTmp = 0;
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
          double harvWood=0; //MG: harvestable wood, m3
          double defIncome = 0.;


          g4m::ageStruct cohortTmp = *cohort_all[asID];
          cohortTmp.setU(rotationForest.get(xi,yi));

          pair<g4m::ageStruct::v, g4m::ageStruct::v> res;  // MG: results vector for the set (old) forest 
          int Rotation = 0;
//          double forFlag = 0.;    //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
          
        if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
          biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], maiForest.get(xi,yi));        // rotation time to get current biomass (without thinning)
          biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUSdTab(iter->CABOVEHA[year], maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
            if (thinningForest.get(xi,yi) > 0) {
              Rotation = biomasRotTh+1; 
            }  else {
              Rotation = biomasRot+1;
            }
//          if (iter->FOREST"].g(1990) >0) forFlag = 1.0;
        }
        if (maiForest.get(xi,yi)> 0){
          rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
//          rotMaxBm = species[int(iter->SPECIESTYPE[year])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm);                        
          rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
        }  
        
        if (Rotation < rotMAI) Rotation = rotMAI;           
          
          res = cohortTmp.aging();
          //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
          //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
          //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
          //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  		  double realAreaO = cohortTmp.getArea();
		  double harvestO = cohortRes(realAreaO,res);
          harvWood = harvestO * iter->FTIMBER[byear];
    
          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//          coeff.PriceC.insert(0,0.);
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
//               , Rotation
//               , harvWood );
          abBiomassO = cohortTmp.getBm();
	      double pDefIncome = abBiomassO * 
                   (decision.priceTimber() * iter->FTIMBER[year] * (1. -coeff.HarvLoos[year]));
	      //Immediate Pay if deforested (Slash and Burn)
	      double sDefIncome = abBiomassO *
		           (decision.priceTimber() * iter->FTIMBER[year] 
		         * (1. -coeff.HarvLoos[year]));
	      defIncome = pDefIncome * (1. - iter->SLASHBURN[year])
		            + sDefIncome * iter->SLASHBURN[year];
		            
//          harvestGrid.set(xi,yi,harvWood * forestArea0);   
          if (managedForest.get(xi,yi) > 0) {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}
            if (woodHarvest[Country0-1] < 0.95 * wprod[regprice].g(year)) {
              if ((Rotation > rotMAI) && (rotationForest.get(xi,yi) > Rotation)) { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,Rotation);		   
                cohortTmp.setU(Rotation);
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  				double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];

                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              } else {//if (rotationForest.get(xi,yi) > rotMAI) { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,rotMAI);		   
                cohortTmp.setU(rotMAI);
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  				double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];

                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              }            
            } else if (woodHarvest[Country0-1] > 1.05 * wprod[regprice].g(year)) {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}                   
              if (rotationForest.get(xi,yi) < Rotation) { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,Rotation);		   
                cohortTmp.setU(Rotation);
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  				double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];

                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;

         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              } else if (rotationForest.get(xi,yi) < rotMaxBmTh) {
                HarvestTmp = harvestGrid.get(xi,yi);
			    rotationForest.set(xi,yi, rotMaxBmTh);
                cohortTmp.setU(rotMaxBmTh);
                res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  				double realAreaO = cohortTmp.getArea();
				double harvestO = cohortRes(realAreaO,res);
				harvWood = harvestO * iter->FTIMBER[byear];

                if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;

  			    newHarvestTmp = harvWood * forestArea0;
			    harvestGrid.set(xi,yi,newHarvestTmp);
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              }
            }
          }
//          iter++;
//        } while (((iter-1)->COUNTRY[2000]) == ((iter)->COUNTRY[2000]));   // Check are we in the same country  // end for Within current country
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\twoodHarvestStat= "<<woodHarvestStat[Country0-1] <<endl;}
    } //else {iter++;}  // end for IF unprotected
   }   // Test only some regions
    iter++;
  } //end for WHILE
//************************End of Second Pass************************************

//cout << "end of second pass" << endl;

////
/////////////////////////////////////////////////////
//////
//////                   Third Pass
//////
/////////////////////////////////////////////////////
////
//  iter = data_all.begin();
//
////cout << "Putting data for current cell into conteiner... "<< endl;
//   while (iter != data_all.end())
//   {if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
//	if (iter->PROTECT[2000] == 0)
//	 {
//        double HarvestTmp = 0;
//        double newHarvestTmp = 0;
//        double abBiomassO = 0.;  
//
//          int xi = (iter->x);
//          int yi = (iter->y);
//          int asID = iter->asID;
//          int Country0 = (int)iter->COUNTRY"][year];
//          int rotUnmanaged = 0;
//          int rotMAI = 0;
//          int rotMaxBm = 0;
//          int rotMaxBmTh = 0;
//          int rotHarvFin = 0;
//          int rotHarvAve = 0;
//          int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
//          int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
//          int biomasRotTh2=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning=2)
//
//          double LandAreaHa = iter->LANDAREA"][year]*100;
//          double forestArea0 = LandAreaHa * iter->FOREST"][year];
//          double harvMAI = maiForest.get(xi,yi)*iter->FTIMBER"][year]*(1-coeff.HarvLoos[year]);
//          double harvWood=0; //MG: harvestable wood, m3
//          double defIncome = 0.;
//
//          
//          
//          
//          g4m::ageStruct cohortTmp = *cohort_all[asID];
//          cohortTmp.setU(rotationForest.get(xi,yi));
//
// 
//        g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
//// End of FM initialisation
//
//        double MAI = maiForest.get(xi,yi); //MG: mean annual increment in tC/ha/year
//
//        if (iter->CABOVEHA"][year] > 0 && maiForest.get(xi,yi)> 0) {
//          biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA"][year], MAI, 1);        // rotation time to get current biomass (without thinning)
//          biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUt(iter->CABOVEHA"][year], MAI, thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
//          biomasRotTh2 = species[int(iter->SPECIESTYPE[year])-1]->gUt(iter->CABOVEHA"][year], MAI, 2);     // rotation time to get current biomass (with thinning)     
////          if (iter->FOREST"].g(1990) >0) forFlag = 1.0;
//        }
//        if (maiForest.get(xi,yi)> 0){
//          rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gTopt(MAI, optimMAI);
//          rotMaxBm = species[int(iter->SPECIESTYPE[year])-1]->gTopt(MAI, optimMaxBm);                        
//          rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(MAI, optimMaxBm);
//        }             
//
//        iter->PRICEC"].insert(0,0.);
//        dima decision((int)iter->BYEAR"][year]
//		       , iter->NPP"]
//		       , iter->SPOPDENS"]
//		       , iter->SAGRSUIT"]
//		       , iter->PRICEINDEX"]
//		       , coeff.PriceIndexE
//		       , iter->R"]
//		       , coeff.PriceC
//		       , coeff.PlantingCostsR
//		       , coeff.PriceLandMinR
//		       , coeff.PriceLandMaxR
//		       , coeff.MaxRotInter
//		       , coeff.MinRotInter
//		       , coeff.decRateL
//		       , coeff.decRateS
//		       , iter->FRACLONGPROD"]
//		       , coeff.baseline
//		       , iter->FTIMBER"] 
//		       , coeff.PriceTimberMaxR
//		       , coeff.PriceTimberMinR
//		       , coeff.FCuptake
//		       , iter->GDP"]
//		       , coeff.HarvLoos
//		       , iter->FOREST].g((int)iter->BYEAR][year])
//		       , wprice["regprice"]
//		       , wprice["regprice0"][year]
//               , rotMAI
//               , harvMAI);
//
////  g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
// double Thinning = -1.;
// double ThinningInit = -1.;
// int Rotation = 0;
// int RotationInit = 0;
//
////if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}
//
//if (woodHarvest[Country0-1] < (0.95 * woodHarvestStat[Country0-1])) 
// {
////if (Country0 == 61){cout<<"Inside IF!!???<"<<endl;};
////  if ((managedForest.get(xi,yi)) <= 0 && (managedForest.get(xi,yi) > -2))
//  if ((managedForest.get(xi,yi)) <= 0)
//   {
//                HarvestTmp = harvestGrid.get(xi,yi);
//                ThinningInit = 1.;
//                RotationInit = biomasRotTh2+1;
//                rotationType.set(xi,yi,11);
//                
//                cohortTmp.setStockingdegree(2);
//                thinningForest.set(xi,yi,1.);
//                  
//     if (MAI > MAI_CountryUprotect[Country0-1])
//     {   
//         if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome))
//                  {
//                          managedForest.set(xi,yi,3.);
////                          Rotation = rotMAI;
//                          Rotation = RotationInit;                          
//                          rotationType.set(xi,yi,1);
//                    }else 
//                    {     managedForest.set(xi,yi,2.);
////                          Rotation = rotMaxBmTh;
//                          Rotation = RotationInit;
//                          rotationType.set(xi,yi,2);
//                     }
//       }else
//       {                  managedForest.set(xi,yi,2.);
////                          Rotation = rotMaxBmTh;
//                          Rotation = RotationInit;                          
//                          rotationType.set(xi,yi,3);
//        }
//       rotationForest.set(xi,yi,Rotation);	
//       cohortTmp.setU(Rotation);
//          res = cohortTmp.aging();
//          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
//          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
//          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
//          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
//          harvWood = (sawnW + restW + sawnThW + restThW + resUse*res.hv) * iter->FTIMBER"][year];
//          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//
//   		newHarvestTmp = harvWood * forestArea0;        
//
//        harvestGrid.set(xi,yi,newHarvestTmp);   
//        woodHarvest[Country0-1] += newHarvestTmp;
//        woodLost[Country0-1] -= HarvestTmp; 
//        }           
// }else if (woodHarvest[Country0-1] > (1.05 * woodHarvestStat[Country0-1])) 
// {
////if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}
////if (Country0 == 61){cout<<"Inside IF!!??? >"<<endl;};
//  if (thinningForest.get(xi,yi) > 0)
//   {
//	    		HarvestTmp = harvestGrid.get(xi,yi);
//	    		
//                ThinningInit = -1.;
//                thinningForest.set(xi,yi,-1.);
//                RotationInit = biomasRot+1;
//                rotationType.set(xi,yi,10);
//
//                cohortTmp.setStockingdegree(-1.);                 
//
//      if (MAI > MAI_CountryUprotect[Country0-1])
//       {   if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome))
//           {   
//                managedForest.set(xi,yi,0.);
////                 thinningForest.set(xi,yi,-1.);
//                 Rotation = biomasRot+1;
//                 rotationType.set(xi,yi,1);
//             
//            }else
//            {     
//                 managedForest.set(xi,yi,-1.);
//                  Rotation = biomasRot+1;
//                  rotationType.set(xi,yi,10);
//              
//            }
//       }else 
//       {  if ((decision.forValNC() * Hurdle_opt[Country0-1]) > (decision.agrVal() + defIncome))
//            {          managedForest.set(xi,yi,-1.);
////                   thinningForest.set(xi,yi,-1.);
//                   Rotation = biomasRot+1;
//                   rotationType.set(xi,yi,10);
//               
//             }else
//             {      
//                    managedForest.set(xi,yi,-2.);
////                    Rotation = biomasRot+1;
//                    Rotation = rotMaxBm;
//                    rotationType.set(xi,yi,10);
//               
//              }
//        }
//       rotationForest.set(xi,yi,Rotation);	
//       cohortTmp.setU(Rotation);
////       cohort.setStockingdegree(thinningForest.get(xi,yi));
//          res = cohortTmp.aging();
//          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
//          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
//          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
//          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
//          harvWood = (sawnW + restW + sawnThW + restThW + resUse*res.hv) * iter->FTIMBER"][year];
//          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//   		newHarvestTmp = harvWood * forestArea0;        
//       
//        woodHarvest[Country0-1] -= (HarvestTmp);
//        harvestGrid.set(xi,yi,0);  
//        woodLost[Country0-1] += newHarvestTmp; 
//      
// }
//}
////cout<< managedForest.get(xi,yi)<<endl;
//
//
///*  
//cout <<"country=\t"<<Country0<<"\t woodHarvest=\t"<< woodHarvest[Country0-1];
//cout <<"\t rotation=\t"<< Rotation<<"\t thinning= \t"<<thinningForest.get(xi,yi);
//cout <<"\t sawnThW=\t"<<sawnThW*forestArea0*4<<"\t restThW=\t"<<restThW*forestArea0*4<<"\t";
//cout <<"\t sawnW= \t"<<sawnW*forestArea0*4<<"\t restW= \t"<<restW*forestArea0*4;
//cout << "\t sawnThWha=\t" << sawnThW << "\t restThWha=\t"<<restThW;
//cout <<"\t sawnWha= \t"<<sawnW<<"\t restWha= \t"<<restW<<"\t abbiomassO= \t"<< abBiomassO<<"\t abbiomass= \t"<<iter->BIOMASS"][year]<<"\t";
//cout <<endl;
//*/
////cout<< "managedForestSeting....\n";
//  
//   } //End for PROTECT == 0
//  } // Test only some regions 
//iter++;
// } // End for WHILE (cell loop) 
//
////    initLoop(0, data_all, fi, cohort_all, newCohort_all, dat_all, maiForest, thinningForest, rotationForest);
//    initCohorts(data_all, fi, cohort_all,maiForest, thinningForest, rotationForest); 
//cout << "Third pass is finished"<< endl;




//******************************************************************************
//**************************Forth Pass********************
//******************************************************************************
//cout << "Start forth pass" << endl;

   iter = data_all.begin();

 while (iter != data_all.end())
 {if (regions.find(iter->POLESREG[2000]) != regions.end()) { // Test only some regions
	   if (iter->PROTECT[2000]==0)    // We consider only unprotected land
	   {
      double HarvestTmp = 0;
      double newHarvestTmp = 0;
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
       double harvWood=0; //MG: harvestable wood, m3
  
 
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
//          cohortTmp.setStockingdegree(thinningForest.get(xi,yi));
 
        pair<g4m::ageStruct::v, g4m::ageStruct::v> res; // MG: results vector for the set (old) forest    
// End of FM initialisation



// double Thinning = 1.;
 int Rotation = 0;
 
        if (iter->CABOVEHA[year] > 0 && maiForest.get(xi,yi)> 0) {
          biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA[year], maiForest.get(xi,yi));        // rotation time to get current biomass (without thinning)
          biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUSdTab(iter->CABOVEHA[year], maiForest.get(xi,yi), thinningForest.get(xi,yi));     // rotation time to get current biomass (with thinning)     
          if (thinningForest.get(xi,yi) > 0)
          {Rotation = biomasRotTh+1; 
          }else
          {Rotation = biomasRot+1;}
//          if (iter->FOREST"].g(1990) >0) forFlag = 1.0;
        }
        if (maiForest.get(xi,yi)> 0){
          rotMAI = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMAI);
//          rotMaxBm = species[int(iter->SPECIESTYPE[year])-1]->gTopt(maiForest.get(xi,yi), optimMaxBm);                        
          rotMaxBmTh = species[int(iter->SPECIESTYPE[year])-1]->gToptt(maiForest.get(xi,yi), optimMaxBm);
        }          
        if (Rotation < rotMAI) Rotation = rotMAI;   

          res = cohortTmp.aging();
          //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
          //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
          //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
          //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  		  double realAreaO = cohortTmp.getArea();
		  double harvestO = cohortRes(realAreaO,res);
          harvWood = harvestO * iter->FTIMBER[byear];

          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;


//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"countryB=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1];
//cout<<"\t 0.8harvestStat\t"<< 0.8 * woodHarvestStat[Country0-1]<<"\t 1.2harvestStat\t"<< 1.2 * woodHarvestStat[Country0-1];
//cout<<"\t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotationInit \t"<<Rotation<<endl;      

//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 0.9 * woodHarvestStat= "<<0.9 * woodHarvestStat[Country0-1] <<endl;}
 if (woodHarvest[Country0-1] < 0.95 * wprod[regprice].g(year)) 
{
///if (Country0 == 61){cout<<"Inside IF (4th path) < !!???"<<endl;};
//
//		   if ((Rotation > rotMAI) && (rotationForest.get(xi,yi) > Rotation))
//		    { 
//	    		HarvestTmp = harvestGrid.get(xi,yi);
//                rotationForest.set(xi,yi,Rotation);		   
//                cohort.setU(Rotation);
//#include "growForest.cpp"
//         		newHarvestTmp = harvWood * forestArea0;
//		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
////cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
////cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 0.8harvestStat\t"<< 0.8 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<Rotation<<endl;               
//                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
//
//            } else if (rotationForest.get(xi,yi) > rotMAI)
		   if (rotationForest.get(xi,yi) != rotMAI)
            { 
	    		HarvestTmp = harvestGrid.get(xi,yi);
                rotationForest.set(xi,yi,rotMAI);		   
                cohortTmp.setU(rotMAI);
//---grow forest for one year
//#include "growForest.cpp"
          res = cohortTmp.aging();
          //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
          //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
          //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
          //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
          //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
          //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
          //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
          //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;
  		  double realAreaO = cohortTmp.getArea();
		  double harvestO = cohortRes(realAreaO,res);
          harvWood = harvestO * iter->FTIMBER[byear];

          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;

         		newHarvestTmp = harvWood * forestArea0;
		     	harvestGrid.set(xi,yi,newHarvestTmp);
//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 0.8harvestStat\t"<< 0.8 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<Rotation<<endl;               
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);

              }            
  

		  
 } else if (woodHarvest[Country0-1] > 1.05 * wprod[regprice].g(year)) 
 {
//if (Country0 == 61){cout<<"woodHarvest= "<<woodHarvest[Country0-1]<<"\t 1.1 * woodHarvestStat= "<<1.1 * woodHarvestStat[Country0-1] <<endl;}
//if (Country0 == 61){cout<<"Inside IF (4th path) >!!???"<<endl;};        
//		   if (rotationForest.get(xi,yi) < Rotation)
//		    { 
//	    		HarvestTmp = harvestGrid.get(xi,yi);
//                rotationForest.set(xi,yi,Rotation);		   
//                cohort.setU(Rotation);
//#include "growForest.cpp"
//         		newHarvestTmp = harvWood * forestArea0;
//		     	harvestGrid.set(xi,yi,harvWood * forestArea0);
////cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";		     	
////cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 1.1harvestStat=\t"<< 1.1 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<Rotation<<endl;               
//                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
//             } else if (rotationForest.get(xi,yi) < rotMaxBmTh)
           if (rotationForest.get(xi,yi) < rotMaxBmTh)
			 {    HarvestTmp = harvestGrid.get(xi,yi);
			      rotationForest.set(xi,yi, rotMaxBmTh);
                  cohortTmp.setU(rotMaxBmTh);
//---grow forest for one year
//#include "growForest.cpp"
          res = cohortTmp.aging();
                //sawnW = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                //restW = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                //sawnThW = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                //restThW = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                //bmH = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                //bmTh = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                //harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
                //harvWood = (sawnW + restW + sawnThW + restThW + harvRes) * iter->FTIMBER[year] ;

          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;

  			     newHarvestTmp = harvWood * forestArea0;
			     harvestGrid.set(xi,yi,newHarvestTmp);
//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";			     
//cout<<"country=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t 1.1harvestStat=\t"<< 1.1 * woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotation \t"<<rotationForest.get(xi,yi)<<endl;               
                woodHarvest[Country0-1] += (newHarvestTmp-HarvestTmp);
              }
 }
}

//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"countryE=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t harvestStat=\t"<< woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotationE \t"<<rotationForest.get(xi,yi);               
//cout <<"\t"<<sawnThW*forestArea0*4<<"\t"<<restThW*forestArea0*4<<"\t";
//cout <<sawnW*forestArea0*4<<"\t"<<restW*forestArea0*4;
//cout<<endl;

//cout <<rotationForest.get(xi,yi)<<"\t";
//cout <<sawnThW*forestArea0*4<<"\t"<<restThW*forestArea0*4<<"\t";
//cout <<sawnW*forestArea0*4<<"\t"<<restW*forestArea0*4;// <<"\t"<<abBiomassO*forestArea0<<"\t";
//cout << endl;

//iter++;
//cout<<"iterC="<<((iter-1)->COUNTRY[2000])<<"\t iterC+1="<<((iter)->COUNTRY[2000]) <<endl;
//} while (((iter-1)->COUNTRY[2000]) == ((iter)->COUNTRY[2000]));   // Check are we in the same country  // end for Within current country
//cout<<"iterS="<<((iter-1)->COUNTRY[2000])<<"\t iterCS+1="<<((iter)->COUNTRY[2000]) <<endl;

//cout<<"x \t"<<xi<<"\t y \t"<<yi<<"\t thinning\t"<<thinningForest.get(xi,yi)<<"\t";
//cout<<"countryE=\t"<<Country0<<"\t harvest=\t"<<woodHarvest[Country0-1]<<"\t harvestStat=\t"<< woodHarvestStat[Country0-1]<<"\t \t \t biomassInit\t"<< iter->BIOMASS"][year]<<"\t biomassO\t"<< abBiomassO<<"\t rotationE \t"<<rotationForest.get(xi,yi);               
//cout<<endl;

//cout <<X<<"\t"<<Y<<"\t"<<Country0<<"\t"<<year<<"\t"
//if (Country0 == 61) 
//{cout<<"Country= "<<Country0<<"\t woodHarvestFinal= "<<woodHarvest[Country0-1]<<"\t woodHarvestStat= "<<woodHarvestStat[Country0-1] <<endl;}

}// else{iter++;}  // end for IF unprotected
} // Test only some regions   
iter++;
} //end for WHILE

	
//************************End of Forth Pass************************************
cout << "End of 4th pass"<<endl;



////************* TEST of INIT FM
//for (int i=1;i<=209;i++)
//{
//cout<<"Country= "<<i<<"\t managedCountInitFM = "<< managedCount[i-1]<<"\t woodHarvest = "<< woodHarvest[i-1]<<endl;
//}
//   iter = data_all.begin();
//
//   while (iter != data_all.end())
//   {
//
//  int Country0 = 1;
//  double HarvestTmp = 0;
//  double newHarvestTmp = 0;
////  double forestArea0 = 0.;
//  double abBiomassO = 0.;  	 
//  int xi = 0;
//  int yi = 0;
//  double X = 0.;
//  double Y = 0.;
//
//	   if (iter->PROTECT[2000]==0)    // We consider only unprotected land
//	   {
//
////cout<<Country0<<endl;       
////#include "countrySelect.cpp"
//      if (   (iter->COUNTRY[2000] == 11)
//          || (iter->COUNTRY[2000] == 25)
//          || (iter->COUNTRY[2000] == 33)        
//          || (iter->COUNTRY[2000] == 33)  
//          || (iter->COUNTRY[2000] == 61)
//          || (iter->COUNTRY[2000] == 62)
//          || (iter->COUNTRY[2000] == 69)             
//          || (iter->COUNTRY[2000] == 156)             
//          || (iter->COUNTRY[2000] == 165)             
//          || (iter->COUNTRY[2000] == 179)                     
//          || (iter->COUNTRY[2000] == 197) ) 
//        { // Test only some countries
//
//       int biomasRot=0;  // MG: rotation time fitted to get certain biomass under certain MAI (w/o thinning)
//       int biomasRotTh=0;  // MG: rotation time fitted to get certain biomass under certain MAI (with thinning)
//       double harvWood=0; //MG: harvestable wood, m3
//
//
//    	            xi = (iter->x);
//    	            yi = (iter->y);
//    	            X = (iter->x)*0.5+0.25-180;
//    	            Y = (iter->y)*0.5+0.25-90;
//    	          
//    	   //cout << "Xi = "<< xi <<"   Yi = "<< yi << endl;       
//if (managedForest.get(xi,yi) > 0)
//{ 
//                       
//       Country0 = (int)iter->COUNTRY"][year];
//       double LandAreaHa = iter->LANDAREA"][year]*100;
//       double forestArea0 = LandAreaHa * iter->FOREST"][year];
//
//
//        g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
//        g4m::ipol<double,double> hlv;  //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
//        g4m::ipol<double,double> hle;  //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
//        g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
//        g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)
//
//        sws.insert(10, .0);
//        sws.insert(30, .6);
//        hlv.insert(0, .0);
//        hlv.insert(25, .7);
//        hle.insert(0, .0);
//        hle.insert(25, .7);
//        dbv.insert(0, 2);
//        dbe.insert(0, 3);
// 
//        g4m::ageStruct::v res; // MG: results vector for the set (old) forest    
//// End of FM initialisation
//
//
//
//// double Thinning = 1.;
// int Rotation = 0;
// 
//  double forFlag = 0.; //MG: a forest area for fitting existing forest in the grid: 0-no forest; 1 - 1 ha of normal forest
////  if (iter->FOREST"].g(1990) >0 && iter->BIOMASS"][year] > 0 && maiForest.get(xi,yi) > 0)
//  if (iter->FOREST"].g(1990) >0 && iter->CABOVEHA"][year] > 0 && maiForest.get(xi,yi) > 0)  
//  
//  {
//          biomasRot = species[int(iter->SPECIESTYPE[year])-1]->gU(iter->CABOVEHA"][year], maiForest.get(xi,yi), 1);       // rotation time to get current biomass (without thinning)
//          biomasRotTh = species[int(iter->SPECIESTYPE[year])-1]->gUt(iter->CABOVEHA"][year], maiForest.get(xi,yi), 1);     // rotation time to get current biomass (with thinning)  
//
//          if (thinningForest.get(xi,yi) == 1)
//          {Rotation = biomasRotTh+1; 
//          }else
//          {Rotation = biomasRot+1;
//          }
//
//          forFlag = 1.0;
//  }     
//
////cout << "\t Rotation1 = \t"<<Rotation<<endl;
//
//			g4m::ageStruct cohort(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), 
//			                            Rotation, thinningForest.get(xi,yi),forFlag, 0.75); 
//
//cohort.aging();
//            cohort.setU(rotationForest.get(xi,yi));
//
////---grow forest for one year
////#include "growForest.cpp"
//          res = cohort.aging();
//          sawnW = res.enSw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
//          restW = res.enRw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
//          sawnThW = res.vnSw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
//          restThW = res.vnRw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
//          harvWood = (sawnW + restW + sawnThW + restThW) * iter->FTIMBER"][year];
//          if ((biomasRot ==0)||(harvWood < 0)) harvWood = 0.;
//     CountriesWoodHarvestM3Year.inc(Country0,iter->BYEAR"][year],harvWood*forestArea0);
////     CountriesWoodHarvestM3Year.inc(Country0,iter->BYEAR"][year],harvWood); 
//}// end for if managed
//iter++;
//} else{iter++;}  // end for Country select
//
//} else{iter++;}  // end for IF unprotected
//	   
////int Nc=0;
//
//
//} //end for WHILE
///////////////////******   END of Test of init FM ************************

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

//void initLoop(int i, dataDetStruct &data_all, g4m::incrementTab &fi,  
//              datGlobal &dat_all, griddata &maiForest, 
//              griddata &thinningForest, griddata &rotationForest) 



 {
//cout << "Begin initLoop"<<endl;

//// Deleating cohorts
//    for (int k = 0; k < cohort_all.size();k++) {
//      delete cohort_all[k];
//    }
//    cohort_all.clear();
//
//    for (int k = 0; k < newCohort_all.size();k++) {
//      delete newCohort_all[k];
//    }
//    newCohort_all.clear();
    

//  double woodHarvest[209];
//  double woodLost[209];
//  double woodPot[209];
//
//  for (int i=0; i<=208; i++){
//    woodHarvest[i]=0.; 
//    woodLost[i]=0.;
//    woodPot[i]=0.;
//  }
// 
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
  //  string agestructFile = settings.inputPath + "\\" + "ageStructDataSEF2011_280316ua.txt";    //new age structure for Italy (email from Giacomo Grassi from 26 July 2011)
 string agestructFile = "ageStructDataSEF2011_280316ua.txt";
//	 string agestructFile = settings.inputPath + "\\" +  "ageStructDataJRC_IT_26July2011_nc.txt";
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
//        double maxForInit = 1-iter->BUILTUP.g(beyear)-iter->CROP.g(beyear);
//        if (maxForInit < 0) maxForInit = 0;
//        double forestShare = 0;        //Actual forest share
//          if (iter->FOREST[byear]+(iter->CROP[byear])+(iter->BUILTUP[byear])>1)
//           {forestShare = (1-(iter->CROP[byear]+iter->BUILTUP[byear]));}        
//          if (iter->FOREST[year]+(iter->CROP.g(2010))+(iter->BUILTUP.g(2010))>1) //CROP.g(2010) == simulation of the dataset as of November 2010 
//           {forestShare = (1-(iter->CROP.g(2010)+iter->BUILTUP.g(2010)));}     // for consistency of the results
//		  else {forestShare = iter->FOREST[byear];}

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
//        if (optRotUnmanaged < 0) optRotUnmanaged = 0;
 //       if (aveRot < 0) aveRot = 0;
        if (thinningForest.get(xi,yi) >0) {
          Rotation = biomasRotTh; 
        } else {
          Rotation = biomasRot;
        }

// Existing (old forest)
//        g4m::ageStruct *cohort = new g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),forFlag, 0.75);

     // Initialize age structure
//    g4m::ageStruct *cohort = new g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, 0, 1, thinningForest.get(xi,yi),0, 0.75);
//    g4m::ageStruct *cohort = new g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, 0, 1, thinningForest.get(xi,yi),0, 0.5);    
//    g4m::ageStruct *cohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[year])-1], sws, hlv, hle, dbv, dbe, 0, 0, 0, 1, thinningForest.get(xi,yi),0, 0.75);  
//    ageStruct
//      (incrementTab *it    //Increment table which will be used, the time step width (simulation period legth) of *it will also be used in ageStruct
//       , ffipol<double> *sws //Sawnwood share of harvested wood depending on dbh
//       , ffipol<double> *hlv //1-harvesting losses thinning (Vornutzung) (depending on d) in relation to standing timber (Vorratsfestmeter)
//       , ffipol<double> *hle //1-harvesting losses final felling (Endnutzung) (depending on d) in relation to standing timber (Vorratsfestmeter)
//       , ffipolm<double> *cov //Thinning costs depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
//       , ffipolm<double> *coe //Harvesting costs depending on d and vol
//       , ffipolm<bool> *dov //Do thinning (depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
//       , ffipolm<bool> *doe //Do final felling (depending on d and stocking volume per hectare)
//       , double mai  //mean annual increment in tC stemmwood per hectar and year at increment optimal rotation time
//       , int objOfProd=3 //objective of production: 0..Rotation time in years, 1..ammount of wood which need to beharvested every year, 2..like 1 but the ammount will not be fulfilled if rotation time will be shorter than for Highest average increment; >2 ignore Value of u and claculate it instead  3 .. Highest average increment, 4 . .Maximum average Biomass, 5 .. Highest possible age, 6 .. Maximum harvest at final cut, 7 .. Average Maximum harvest at final cut
//       , double u=0.  //Rotation time if objOfProd 0
//       , double minSw=0. //if objOfProd 1,2 ammount of sawnwood to harvest
//       , double minRw=0. //if objOfProd 1,2 ammount of restwood to harvest
//       , double minHarv=0. //if objOfProd 1,2 ammount of total harvest
//       //Usage of stocking degree: 0..Keep all the time sdMax, 1..alternate between sdMin and sdMax, 3..alternate between sdMinH and sdMaxH
//       , int sdDef=0
//       //Stocking degree: if sd exceeds sdMax do thinning until sdMin. sd > 0 stockingDegree yield table, sd -1 to 0 natural stocking degree
//       //Maybe sdMin and sdMax can be made dependent on h/hmax and MAI
//       , double sdMax=1.
//       , double sdMin=1.
//       , unsigned int maiYears=30  //Years to calculate average mai
//       //Minimal rotation time in years or as share given in minRotRef which needs to be exceedes until final harvests are done
//       , double minRotVal=0.75
//       //meaning of minRotVal value
//       //0..use it as years, 1..minRotVal*u (u>0), 2..*uMaxAvgIncrement, 3..*uMaxAvgBiomass, 4..*uMaxAge, 5..*uMaxHarvest, 6..*uAvgMaxHarvest
//       , int minRotRef=2
//       //how fast should the stoking degree target be reached
//   //0..do only remove caused by stand density  to  1..do only typical removes
//       , double flexSd=0.
//       , ffipol<double> *sdMaxH=NULL//Stockindegree depending on max tree height
//       , ffipol<double> *sdMinH=NULL//Stockindegree depending on max tree height
//       );
                                                                                  
//		if (Country==224){ cout<<  MAIRot<<endl;} 
   g4m::ageStruct *cohort;
//   cohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[byear])-1], &ffsws, &ffhlv, &ffhle, &ffcov, &ffcoe, &ffdov, &ffdoe, maiForest.get(xi,yi)
   cohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[byear])-1], &ffsws, &ffhlv, &ffhle, ffcov, ffcoe, ffdov, ffdoe, maiForest.get(xi,yi)
	   , 0  //What stands the value of u for
     , 1  //Rotation time
     , 0, 0, 0
//     ,1 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, MAIRot
     ,0 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, 1
     , 1  //reference of minrot
     , 0. //Flexibility of stocking degree
     , &ffsdMaxH, &ffsdMinH,1);
   

   	//cohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[byear])-1], &ffsws, &ffhlv, &ffhle, &ffcov, &ffcoe, &ffdov, &ffdoe, 1.5
    //   , 1  //What stands the value of u for
    //   , 100  //Rotation time
    //   , 0, 0, 0
    //   ,0 , 1, 1, 30, MAIRot
    //   , 0  //reference of minrot
    //   , 0. //Flexibility of stocking degree
    //   , &ffsdMaxH, &ffsdMinH
    //   );


 //   cohort = new g4m::ageStruct 
	//(species[0], &ffsws, &ffhlv, &ffhle, &ffcov, &ffcoe, &ffdov, &ffdoe, 1.5
 //      , 0  //What stands the value of u for
 //      , 100  //Rotation time
 //      , 0, 0, 0
 //      ,0 , 1, 1, 30, MAIRot
 //      , 0  //reference of minrot
 //      , 0. //Flexibility of stocking degree
 //      , &ffsdMaxH, &ffsdMinH
 //      );

    double abBiomass0 = 0.;
//    cohort->setMai(maiForest.get(xi,yi));
 if (iter->POTVEG[2000]<10 && iter->FOREST[2000]>0 && iter->MAIE[2000]>0){
  map<int, vector<double> >::iterator PageStructData = ageStructData.find(Country);
  if((PageStructData != ageStructData.end()) && thinningForest.get(xi,yi)>0 && forFlag > 0) {
//  if((PageStructData != ageStructData.end()) && thinningForest.get(xi,yi)>0 && forFlag > 1000) { // Don't use Age Struct Data!!!
//    cohort->createNormalForest(161, 0., 1.);
    cohort->createNormalForest(321, 0., 1.);
    
//    double ageBreaks[] = {10.,20.,40.,60.,80.,100.,120.,140.,999.};
//    double ageSize[] = {11.,10.,20.,20.,20.,20.,20.,20.,20.};
//    int ageGroups = 9;
    double ageBreaks[] = {10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,999.};
    double ageSize[] = {11.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10.};
    int ageGroups = 16;

//    int oldest = 0;
    int oldestAgeGroup = 0;
    int   oag = 0;    
    int ageGroup = 0;
//    for(int i=0; i<161; ++i) {
    for(int i=1; i<161; ++i) {
      while(i>ageBreaks[ageGroup] && ageGroup+1 < ageGroups) {++ageGroup;}
      cohort->setArea(i,PageStructData->second[ageGroup]/ageSize[ageGroup]);
      if (PageStructData->second[ageGroup]>0) {oldestAgeGroup = ageGroup;}
    }
    double area = cohort->getArea();
    double biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);
//cout<<"biomass1=\t"<<biomass<<"\tarea1=\t"<<area<<endl;
    if(area > 0. && biomass > 0.) {

//Tune age structure for current cell
    
//    if (biomass < iter->CABOVEHA"][byear]){
    if (biomass < 0.95*iter->CABOVEHA[byear]){                
       int young = 0;
       oag = oldestAgeGroup; // current oldest age group
//       while (biomass < iter->CABOVEHA"][byear] && oag < 30){
       while (biomass < 0.95*iter->CABOVEHA[byear] && oag < 30){                          
          if (ageSize[young]>0 && oag>young){   
             for (int i=0;i<10;++i) {
//               cohort->setArea(i+young*10,(PageStructData->second[young])/(ageSize[young]*2));
               double areaTmp = cohort->getArea(i+young*10+1);  
               cohort->setArea(i+young*10+1,areaTmp/2);
//               cohort->setArea(i+oag*10+1,(PageStructData->second[young])/(ageSize[young]*2));
               cohort->setArea(i+(oag+1)*10+1,areaTmp/2);
//               cohort->setBm(i+oldestAgeGroup+1, stockingDegree*cohort->getBm(i+oldestAgeGroup+1));
               biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);
               }
           }   
         young +=1;
         oag+=1;
       }    

//     } else if (biomass > 3 * iter->CABOVEHA"][byear]){/
//     } else if (biomass > 2.5 * iter->CABOVEHA"][byear]){            
     } else if (biomass > 2 * iter->CABOVEHA[byear]){        //v_24_11                
        int young = 0;
        oag = oldestAgeGroup;          
//        while (biomass > 3 * iter->CABOVEHA"][byear] && oag >2){
        while (biomass > 2 * iter->CABOVEHA[byear] && oag >2){              
           if (oag > young){
             if (ageSize[oag]>0 && ageSize[young]>0){      
              for (int i=0;i<ageSize[oag];++i) {
               double areaTmp_oag = cohort->getArea(i+oag*10+1);                    
               double areaTmp_young = cohort->getArea(i+young*10+1);
               cohort->setArea(i+oag*10+1,0.);
               cohort->setArea(i+young*10+1,(areaTmp_oag+areaTmp_young));
//               cohort->setBm(i+oldestAgeGroup+1, stockingDegree*cohort->getBm(i+oldestAgeGroup+1));
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
//               cohort->setBm(i+oldestAgeGroup+1, stockingDegree*cohort->getBm(i+oldestAgeGroup+1));
               biomass = cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000]);
               }
              }       
            }    
          young +=1;
          oag-=1;
         }    
//      double stockingDegree = iter->CABOVEHA"][year]/(cohort->getBm()/cohort->getArea());
//cout<< "countryIn = "<<Country<<"\tarea = "<<cohort->getArea()<<"\tbm = "<<cohort->getBm()<<"\tabC= "<<iter->CABOVEHA"][year]<<"\t stockingDegree = "<<stockingDegree<<"\t Th= "<<thinningForest.get(xi,yi)<<"\t oldestAgeGroup= "<<oldestAgeGroup<<"\t oag= "<<oag<<endl;
      }  

      double stockingDegree = iter->CABOVEHA[byear]/(cohort->getBm() * BEF(int(iter->SPECIESTYPE[byear])-1,cohort->getBm(),iter->FTIMBER[2000])/cohort->getArea());
      if(stockingDegree < 0.) {stockingDegree = 0.;}
//	  if (Country==224){
//		  cout<<X<<"\t"<<Y<<"\t"<<stockingDegree<<endl;
//	  }
                 
//         cohort->setArea(i,PageStructData->second[ageGroup]/ageSize[ageGroup]);
//cout<< "countryIn = "<<Country<<"\tarea = "<<cohort->getArea()<<"\tbm = "<<cohort->getBm()<<"\tabC= "<<iter->CABOVEHA"][year]<<"\t Th= "<<thinningForest.get(xi,yi)<<endl;
//cout<< "countryIn = "<<Country<<"\tarea = "<<cohort->getArea()<<"\tbm = "<<cohort->getBm()<<"\tabC= "<<iter->CABOVEHA"][year]<<"\t stockingDegree = "<<stockingDegree<<"\t Th= "<<thinningForest.get(xi,yi)<<"\t oldestAgeGroup= "<<oldestAgeGroup<<"\t oag= "<<oag<<endl;
//      cohort->setStockingdegree(thinningForest.get(xi,yi));
      cohort->setStockingdegreeMin(stockingDegree*sdMinCoeff);
      cohort->setStockingdegreeMax(stockingDegree*sdMaxCoeff);
      thinningForest.set(xi,yi, stockingDegree);
      for(int i=0; i<321; ++i) {
      cohort->setBm(i, stockingDegree*cohort->getBm(i));
      }
	  // Olya T.
	  
	  if (Country==224){
    /*    double CA_0_10 = (cohort->getArea(0)+cohort->getArea(1)+cohort->getArea(2)+cohort->getArea(3)+cohort->getArea(4)+cohort->getArea(5)+cohort->getArea(6)+cohort->getArea(7)+cohort->getArea(8)+cohort->getArea(9)+cohort->getArea(10))*(LandAreaHa*OforestShare)/0.993351;
		double CA_11_20 = (cohort->getArea(11)+cohort->getArea(12)+cohort->getArea(13)+cohort->getArea(14)+cohort->getArea(15)+cohort->getArea(16)+cohort->getArea(17)+cohort->getArea(18)+cohort->getArea(19)+cohort->getArea(20))*( LandAreaHa*OforestShare)/0.993351;
		double CA_21_30 = (cohort->getArea(21)+cohort->getArea(22)+cohort->getArea(23)+cohort->getArea(24)+cohort->getArea(25)+cohort->getArea(26)+cohort->getArea(27)+cohort->getArea(28)+cohort->getArea(29)+cohort->getArea(30))*( LandAreaHa*OforestShare)/0.993351;
		double CA_31_40 = (cohort->getArea(31)+cohort->getArea(32)+cohort->getArea(33)+cohort->getArea(34)+cohort->getArea(35)+cohort->getArea(36)+cohort->getArea(37)+cohort->getArea(38)+cohort->getArea(39)+cohort->getArea(40))*( LandAreaHa*OforestShare)/0.993351;
		double CA_41_50 = (cohort->getArea(41)+cohort->getArea(42)+cohort->getArea(43)+cohort->getArea(44)+cohort->getArea(45)+cohort->getArea(46)+cohort->getArea(47)+cohort->getArea(48)+cohort->getArea(49)+cohort->getArea(50))*( LandAreaHa*OforestShare)/0.993351;
		double CA_51_60 = (cohort->getArea(51)+cohort->getArea(52)+cohort->getArea(53)+cohort->getArea(54)+cohort->getArea(55)+cohort->getArea(56)+cohort->getArea(57)+cohort->getArea(58)+cohort->getArea(59)+cohort->getArea(60))*( LandAreaHa*OforestShare)/0.993351;
		double CA_61_70 = (cohort->getArea(61)+cohort->getArea(62)+cohort->getArea(63)+cohort->getArea(64)+cohort->getArea(65)+cohort->getArea(66)+cohort->getArea(67)+cohort->getArea(68)+cohort->getArea(69)+cohort->getArea(70))*( LandAreaHa*OforestShare)/0.993351;
		double CA_71_80 = (cohort->getArea(71)+cohort->getArea(72)+cohort->getArea(73)+cohort->getArea(74)+cohort->getArea(75)+cohort->getArea(76)+cohort->getArea(77)+cohort->getArea(78)+cohort->getArea(79)+cohort->getArea(80))*( LandAreaHa*OforestShare)/0.993351;
		double CA_81_90 = (cohort->getArea(81)+cohort->getArea(82)+cohort->getArea(83)+cohort->getArea(84)+cohort->getArea(85)+cohort->getArea(86)+cohort->getArea(87)+cohort->getArea(88)+cohort->getArea(89)+cohort->getArea(90))*( LandAreaHa*OforestShare)/0.993351;
		double CA_91_100 = (cohort->getArea(91)+cohort->getArea(92)+cohort->getArea(93)+cohort->getArea(94)+cohort->getArea(95)+cohort->getArea(96)+cohort->getArea(97)+cohort->getArea(98)+cohort->getArea(99)+cohort->getArea(100))*( LandAreaHa*OforestShare)/0.993351;
		double CA_101_110 = (cohort->getArea(101)+cohort->getArea(102)+cohort->getArea(103)+cohort->getArea(104)+cohort->getArea(105)+cohort->getArea(106)+cohort->getArea(107)+cohort->getArea(108)+cohort->getArea(109)+cohort->getArea(110))*( LandAreaHa*OforestShare)/0.993351;
		double CA_111_120 = (cohort->getArea(111)+cohort->getArea(112)+cohort->getArea(113)+cohort->getArea(114)+cohort->getArea(115)+cohort->getArea(116)+cohort->getArea(117)+cohort->getArea(118)+cohort->getArea(119)+cohort->getArea(120))*( LandAreaHa*OforestShare)/0.993351;
		double CA_121_130 = (cohort->getArea(121)+cohort->getArea(122)+cohort->getArea(123)+cohort->getArea(124)+cohort->getArea(125)+cohort->getArea(126)+cohort->getArea(127)+cohort->getArea(128)+cohort->getArea(129)+cohort->getArea(130))*( LandAreaHa*OforestShare)/0.993351;
		double CA_131_140 = (cohort->getArea(131)+cohort->getArea(132)+cohort->getArea(133)+cohort->getArea(134)+cohort->getArea(135)+cohort->getArea(136)+cohort->getArea(137)+cohort->getArea(138)+cohort->getArea(139)+cohort->getArea(140))*( LandAreaHa*OforestShare)/0.993351;
		double CA_141_150 = (cohort->getArea(141)+cohort->getArea(142)+cohort->getArea(143)+cohort->getArea(144)+cohort->getArea(145)+cohort->getArea(146)+cohort->getArea(147)+cohort->getArea(148)+cohort->getArea(149)+cohort->getArea(150))*( LandAreaHa*OforestShare)/0.993351;
		double CA_151_160 = (cohort->getArea(151)+cohort->getArea(152)+cohort->getArea(153)+cohort->getArea(154)+cohort->getArea(155)+cohort->getArea(156)+cohort->getArea(157)+cohort->getArea(158)+cohort->getArea(159)+cohort->getArea(160))*( LandAreaHa*OforestShare)/0.993351;
		double CA_161_170 = (cohort->getArea(161)+cohort->getArea(162)+cohort->getArea(163)+cohort->getArea(164)+cohort->getArea(165)+cohort->getArea(166)+cohort->getArea(167)+cohort->getArea(168)+cohort->getArea(169)+cohort->getArea(170))*( LandAreaHa*OforestShare)/0.993351;
		double CA_171_180= (cohort->getArea(171)+cohort->getArea(172)+cohort->getArea(173)+cohort->getArea(174)+cohort->getArea(175)+cohort->getArea(176)+cohort->getArea(177)+cohort->getArea(178)+cohort->getArea(179)+cohort->getArea(180))*( LandAreaHa*OforestShare)/0.993351;
		double CA_181_190= (cohort->getArea(181)+cohort->getArea(182)+cohort->getArea(183)+cohort->getArea(184)+cohort->getArea(185)+cohort->getArea(186)+cohort->getArea(187)+cohort->getArea(188)+cohort->getArea(189)+cohort->getArea(190))*( LandAreaHa*OforestShare)/0.993351;
		double CA_191_200= (cohort->getArea(191)+cohort->getArea(192)+cohort->getArea(193)+cohort->getArea(194)+cohort->getArea(195)+cohort->getArea(196)+cohort->getArea(197)+cohort->getArea(198)+cohort->getArea(199)+cohort->getArea(200))*( LandAreaHa*OforestShare)/0.993351;
		double CA_201_210= (cohort->getArea(201)+cohort->getArea(202)+cohort->getArea(203)+cohort->getArea(204)+cohort->getArea(205)+cohort->getArea(206)+cohort->getArea(207)+cohort->getArea(208)+cohort->getArea(209)+cohort->getArea(210))*( LandAreaHa*OforestShare)/0.993351;
		double CA_211_220= (cohort->getArea(211)+cohort->getArea(212)+cohort->getArea(213)+cohort->getArea(214)+cohort->getArea(215)+cohort->getArea(216)+cohort->getArea(217)+cohort->getArea(218)+cohort->getArea(219)+cohort->getArea(220))*( LandAreaHa*OforestShare)/0.993351;
		*/

	/*	double B_0_10 = (cohort->getBm(0)+cohort->getBm(1)+cohort->getBm(2)+cohort->getBm(3)+cohort->getBm(4)+cohort->getBm(5)+cohort->getBm(6)+cohort->getBm(7)+cohort->getBm(8)+cohort->getBm(9)+cohort->getBm(10))/0.993351/10;
		double B_11_20 = (cohort->getBm(11)+cohort->getBm(12)+cohort->getBm(13)+cohort->getBm(14)+cohort->getBm(15)+cohort->getBm(16)+cohort->getBm(17)+cohort->getBm(18)+cohort->getBm(19)+cohort->getBm(20))/0.993351/10;
		double B_21_30 = (cohort->getBm(21)+cohort->getBm(22)+cohort->getBm(23)+cohort->getBm(24)+cohort->getBm(25)+cohort->getBm(26)+cohort->getBm(27)+cohort->getBm(28)+cohort->getBm(29)+cohort->getBm(30))/0.993351/10;
		double B_31_40 = (cohort->getBm(31)+cohort->getBm(32)+cohort->getBm(33)+cohort->getBm(34)+cohort->getBm(35)+cohort->getBm(36)+cohort->getBm(37)+cohort->getBm(38)+cohort->getBm(39)+cohort->getBm(40))/0.993351/10;
		double B_41_50 = (cohort->getBm(41)+cohort->getBm(42)+cohort->getBm(43)+cohort->getBm(44)+cohort->getBm(45)+cohort->getBm(46)+cohort->getBm(47)+cohort->getBm(48)+cohort->getBm(49)+cohort->getBm(50))/0.993351/10;
		double B_51_60 = (cohort->getBm(51)+cohort->getBm(52)+cohort->getBm(53)+cohort->getBm(54)+cohort->getBm(55)+cohort->getBm(56)+cohort->getBm(57)+cohort->getBm(58)+cohort->getBm(59)+cohort->getBm(60))/0.993351/10;
		double B_61_70 = (cohort->getBm(61)+cohort->getBm(62)+cohort->getBm(63)+cohort->getBm(64)+cohort->getBm(65)+cohort->getBm(66)+cohort->getBm(67)+cohort->getBm(68)+cohort->getBm(69)+cohort->getBm(70))/0.993351/10;
		double B_71_80 = (cohort->getBm(71)+cohort->getBm(72)+cohort->getBm(73)+cohort->getBm(74)+cohort->getBm(75)+cohort->getBm(76)+cohort->getBm(77)+cohort->getBm(78)+cohort->getBm(79)+cohort->getBm(80))/0.993351/10;
		double B_81_90 = (cohort->getBm(81)+cohort->getBm(82)+cohort->getBm(83)+cohort->getBm(84)+cohort->getBm(85)+cohort->getBm(86)+cohort->getBm(87)+cohort->getBm(88)+cohort->getBm(89)+cohort->getBm(90))/0.993351/10;
		double B_91_100 = (cohort->getBm(91)+cohort->getBm(92)+cohort->getBm(93)+cohort->getBm(94)+cohort->getBm(95)+cohort->getBm(96)+cohort->getBm(97)+cohort->getBm(98)+cohort->getBm(99)+cohort->getBm(100))/0.993351/10;
		double B_101_110 = (cohort->getBm(101)+cohort->getBm(102)+cohort->getBm(103)+cohort->getBm(104)+cohort->getBm(105)+cohort->getBm(106)+cohort->getBm(107)+cohort->getBm(108)+cohort->getBm(109)+cohort->getBm(110))/0.993351/10;
		double B_111_120 = (cohort->getBm(111)+cohort->getBm(112)+cohort->getBm(113)+cohort->getBm(114)+cohort->getBm(115)+cohort->getBm(116)+cohort->getBm(117)+cohort->getBm(118)+cohort->getBm(119)+cohort->getBm(120))/0.993351/10;
		double B_121_130 = (cohort->getBm(121)+cohort->getBm(122)+cohort->getBm(123)+cohort->getBm(124)+cohort->getBm(125)+cohort->getBm(126)+cohort->getBm(127)+cohort->getBm(128)+cohort->getBm(129)+cohort->getBm(130))/0.993351/10;
		double B_131_140 = (cohort->getBm(131)+cohort->getBm(132)+cohort->getBm(133)+cohort->getBm(134)+cohort->getBm(135)+cohort->getBm(136)+cohort->getBm(137)+cohort->getBm(138)+cohort->getBm(139)+cohort->getBm(140))/0.993351/10;
		double B_141_150 = (cohort->getBm(141)+cohort->getBm(142)+cohort->getBm(143)+cohort->getBm(144)+cohort->getBm(145)+cohort->getBm(146)+cohort->getBm(147)+cohort->getBm(148)+cohort->getBm(149)+cohort->getBm(150))/0.993351/10;
		double B_151_160 = (cohort->getBm(151)+cohort->getBm(152)+cohort->getBm(153)+cohort->getBm(154)+cohort->getBm(155)+cohort->getBm(156)+cohort->getBm(157)+cohort->getBm(158)+cohort->getBm(159)+cohort->getBm(160))/0.993351/10;
		double B_161_170 = (cohort->getBm(161)+cohort->getBm(162)+cohort->getBm(163)+cohort->getBm(164)+cohort->getBm(165)+cohort->getBm(166)+cohort->getBm(167)+cohort->getBm(168)+cohort->getBm(169)+cohort->getBm(170))/0.993351/10;
		double B_171_180= (cohort->getBm(171)+cohort->getBm(172)+cohort->getBm(173)+cohort->getBm(174)+cohort->getBm(175)+cohort->getBm(176)+cohort->getBm(177)+cohort->getBm(178)+cohort->getBm(179)+cohort->getBm(180))/0.993351/10;
		double B_181_190= (cohort->getBm(181)+cohort->getBm(182)+cohort->getBm(183)+cohort->getBm(184)+cohort->getBm(185)+cohort->getBm(186)+cohort->getBm(187)+cohort->getBm(188)+cohort->getBm(189)+cohort->getBm(190))/0.993351/10;
		double B_191_200= (cohort->getBm(191)+cohort->getBm(192)+cohort->getBm(193)+cohort->getBm(194)+cohort->getBm(195)+cohort->getBm(196)+cohort->getBm(197)+cohort->getBm(198)+cohort->getBm(199)+cohort->getBm(200))/0.993351/10;
		double B_201_210= (cohort->getBm(201)+cohort->getBm(202)+cohort->getBm(203)+cohort->getBm(204)+cohort->getBm(205)+cohort->getBm(206)+cohort->getBm(207)+cohort->getBm(208)+cohort->getBm(209)+cohort->getBm(210))/0.993351/10;
		double B_211_220= (cohort->getBm(211)+cohort->getBm(212)+cohort->getBm(213)+cohort->getBm(214)+cohort->getBm(215)+cohort->getBm(216)+cohort->getBm(217)+cohort->getBm(218)+cohort->getBm(219)+cohort->getBm(220))/0.993351/10;
			
*/
			  double B_221_230=(cohort->getBm(221)+cohort->getBm(222)+cohort->getBm(223)+cohort->getBm(224)+cohort->getBm(225)+cohort->getBm(226)+cohort->getBm(227)+cohort->getBm(228)+cohort->getBm(229)+cohort->getBm(230))/0.993351/10;
			  double B_231_240=(cohort->getBm(231)+cohort->getBm(232)+cohort->getBm(233)+cohort->getBm(234)+cohort->getBm(235)+cohort->getBm(236)+cohort->getBm(237)+cohort->getBm(238)+cohort->getBm(239)+cohort->getBm(240))/0.993351/10;
			  double B_241_250=(cohort->getBm(241)+cohort->getBm(242)+cohort->getBm(243)+cohort->getBm(244)+cohort->getBm(245)+cohort->getBm(246)+cohort->getBm(247)+cohort->getBm(248)+cohort->getBm(249)+cohort->getBm(250))/0.993351/10;
			  double B_251_260=(cohort->getBm(251)+cohort->getBm(252)+cohort->getBm(253)+cohort->getBm(254)+cohort->getBm(255)+cohort->getBm(256)+cohort->getBm(257)+cohort->getBm(258)+cohort->getBm(259)+cohort->getBm(260))/0.993351/10;
			  double B_261_270=(cohort->getBm(261)+cohort->getBm(262)+cohort->getBm(263)+cohort->getBm(264)+cohort->getBm(265)+cohort->getBm(266)+cohort->getBm(267)+cohort->getBm(268)+cohort->getBm(269)+cohort->getBm(270))/0.993351/10;		
		if (outfile.is_open()) {
			
			/*outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_0_10"<<"\t"<<CA_0_10<<"\t"<<B_0_10<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_11_20"<<"\t"<<CA_11_20<<"\t"<<B_11_20<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_21_30"<<"\t"<<CA_21_30<<"\t"<<B_21_30<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_31_40"<<"\t"<<CA_31_40<<"\t"<<B_31_40<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_41_50"<<"\t"<<CA_41_50<<"\t"<<B_41_50<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_51_60"<<"\t"<<CA_51_60<<"\t"<<B_51_60<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_61_70"<<"\t"<<CA_61_70<<"\t"<<B_61_70<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_71_80"<<"\t"<<CA_71_80<<"\t"<<B_71_80<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_81_90"<<"\t"<<CA_81_90<<"\t"<<B_81_90<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_91_100"<<"\t"<<CA_91_100<<"\t"<<B_91_100<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_101_110"<<"\t"<<CA_101_110<<"\t"<<B_101_110<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_111_120"<<"\t"<<CA_111_120<<"\t"<<B_111_120<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_121_130"<<"\t"<<CA_121_130<<"\t"<<B_121_130<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_131_140"<<"\t"<<CA_131_140<<"\t"<<B_131_140<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_141_150"<<"\t"<<CA_141_150<<"\t"<<B_141_150<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_151_160"<<"\t"<<CA_151_160<<"\t"<<B_151_160<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_161_170"<<"\t"<<CA_161_170<<"\t"<<B_161_170<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_171_180"<<"\t"<<CA_171_180<<"\t"<<B_171_180<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_181_190"<<"\t"<<CA_181_190<<"\t"<<B_181_190<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_191_200"<<"\t"<<CA_191_200<<"\t"<<B_191_200<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_201_210"<<"\t"<<CA_201_210<<"\t"<<B_201_210<<endl;
			outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tCA_211_220"<<"\t"<<CA_211_220<<"\t"<<B_211_220<<endl;			*/	
outfile<<asID<<"\t"<<X<<"\t"<<Y<<"\tBiom27"<<"\t"<<B_261_270<<endl;

		};	   //				***End of file
	  }     // End of biomass and erea extraction for UA
	  
	
//cout<< "countryIn = "<<Country<<"\tarea = "<<cohort->getArea()<<"\tbm = "<<cohort->getBm()<<"\tabC= "<<iter->CABOVEHA"][year]<<"\t stockingDegree = "<<stockingDegree<<"\t Th= "<<thinningForest.get(xi,yi)<<"\t oldestAgeGroup= "<<oldestAgeGroup<<"\t oag= "<<oag<<endl;
      cohort->setU(321);
      cohort->aging();
//cout<< "countryIn = "<<Country<<"\tarea = "<<cohort->getArea()<<"\tbm = "<<cohort->getBm()<<"\tabC= "<<iter->CABOVEHA"][year]<<"\t stockingDegree = "<<stockingDegree<<"\t Th= "<<thinningForest.get(xi,yi)<<"\t oldestAgeGroup= "<<oldestAgeGroup<<"\t oag= "<<oag<<endl;
    }
	//if ((Country=224)&&(OforestShare>0)) {
		//cout<< "country=" <<Country<<"		asId"<<asID<<"		OforestShare="<< OforestShare<<"	area="<<area<< endl; };  
    
//    cohort->setRotPeriod(species[int(iter->SPECIESTYPE[year])-1]->gUt(cohort->getBm()/cohort->getArea(), maiForest.get(xi,yi), 1));
    Rotation = species[int(iter->SPECIESTYPE[byear])-1]->gUSdTab(cohort->getBm()/cohort->getArea(), maiForest.get(xi,yi), thinningForest.get(xi,yi))+1;
    if (Rotation < MAIRot) Rotation = MAIRot;
    cohort->setU(Rotation);  
//    cohort->setRotPeriod(rotationForest.get(xi,yi));
//      rotationForest.set(xi,yi,Rotation);
//  } else {
  } else if (forFlag > 0) {
    if (Rotation < MAIRot) Rotation = MAIRot;
    cohort->createNormalForest(Rotation, forFlag, thinningForest.get(xi,yi));
    cohort->setStockingdegreeMin(thinningForest.get(xi,yi)*sdMinCoeff);         
    cohort->setStockingdegreeMax(thinningForest.get(xi,yi)*sdMaxCoeff);         
    cohort->setU(Rotation);
//    cohort->aging();    

//    cohort->setU(Rotation);    
  } else {cohort->createNormalForest(1, forFlag, 1);} // MG: create an existing forest with 0 area for consistency of the singleCell structure
        
  rotationForest.set(xi,yi,Rotation);      
        
/*        
        
        pair<g4m::ageStruct::v, g4m::ageStruct::v> res;    // MG: results vector for the set (old) forest
        res = cohort->aging();
//        res = cohort.aging();        
        sawnWpot = res.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
        restWpot = res.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
        sawnThWpot = res.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
        restThWpot = res.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
        bmHpot = res.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
        bmThPot = res.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
        harvResPot = (bmHpot+bmThPot - (sawnWpot + restWpot + sawnThWpot + restThWpot)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
        double potHarvest = (sawnWpot + restWpot + sawnThWpot + restThWpot + harvResPot) * iter->FTIMBER[byear] ;
*/


        double abBiomass0 = cohort->getBm(); // modelled biomass at time 0, tC/ha
 } //End for FOREST>0
//cout<<"biomass2=\t"<<abBiomass0<<endl;
// New forest (planted/afforested)
//        g4m::ageStruct *newCohort = new g4m::ageStruct(&fi, sws, hlv, hle, dbv, dbe, 0, 0, maiForest.get(xi,yi), Rotation, thinningForest.get(xi,yi),0, 0.75); 

//    ageStruct
//      (incrementTab *it    //Increment table which will be used, the time step width (simulation period legth) of *it will also be used in ageStruct
//       , ffipol<double> *sws //Sawnwood share of harvested wood depending on dbh
//       , ffipol<double> *hlv //1-harvesting losses thinning (Vornutzung) (depending on d) in relation to standing timber (Vorratsfestmeter)
//       , ffipol<double> *hle //1-harvesting losses final felling (Endnutzung) (depending on d) in relation to standing timber (Vorratsfestmeter)
//       , ffipolm<double> *cov //Thinning costs depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
//       , ffipolm<double> *coe //Harvesting costs depending on d and vol
//       , ffipolm<bool> *dov //Do thinning (depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
//       , ffipolm<bool> *doe //Do final felling (depending on d and stocking volume per hectare)
//       , double mai  //mean annual increment in tC stemmwood per hectar and year at increment optimal rotation time
//       , int objOfProd=3 //objective of production: 0..Rotation time in years, 1..ammount of wood which need to beharvested every year, 2..like 1 but the ammount will not be fulfilled if rotation time will be shorter than for Highest average increment; >2 ignore Value of u and claculate it instead  3 .. Highest average increment, 4 . .Maximum average Biomass, 5 .. Highest possible age, 6 .. Maximum harvest at final cut, 7 .. Average Maximum harvest at final cut
//       , double u=0.  //Rotation time if objOfProd 0
//       , double minSw=0. //if objOfProd 1,2 ammount of sawnwood to harvest
//       , double minRw=0. //if objOfProd 1,2 ammount of restwood to harvest
//       , double minHarv=0. //if objOfProd 1,2 ammount of total harvest
//       //Usage of stocking degree: 0..Keep all the time sdMax, 1..alternate between sdMin and sdMax, 3..alternate between sdMinH and sdMaxH
//       , int sdDef=0
//       //Stocking degree: if sd exceeds sdMax do thinning until sdMin. sd > 0 stockingDegree yield table, sd -1 to 0 natural stocking degree
//       //Maybe sdMin and sdMax can be made dependent on h/hmax and MAI
//       , double sdMax=1.
//       , double sdMin=1.
//       , unsigned int maiYears=30  //Years to calculate average mai
//       //Minimal rotation time in years or as share given in minRotRef which needs to be exceedes until final harvests are done
//       , double minRotVal=0.75
//       //meaning of minRotVal value
//       //0..use it as years, 1..minRotVal*u (u>0), 2..*uMaxAvgIncrement, 3..*uMaxAvgBiomass, 4..*uMaxAge, 5..*uMaxHarvest, 6..*uAvgMaxHarvest
//       , int minRotRef=2
//       //how fast should the stoking degree target be reached
//   //0..do only remove caused by stand density  to  1..do only typical removes
//       , double flexSd=0.
//       , ffipol<double> *sdMaxH=NULL//Stockindegree depending on max tree height
//       , ffipol<double> *sdMinH=NULL//Stockindegree depending on max tree height
//		 , unsigned int maxAge=300 // maximal age of forest considered. Number of cohorts = maxAge/timeStep. maxAge is altered by rotation in createNormalForest	
//       );

   g4m::ageStruct *newCohort;
   newCohort = new g4m::ageStruct(species[int(iter->SPECIESTYPE[byear])-1], &ffsws, &ffhlv, &ffhle, ffcov, ffcoe, ffdov, ffdoe, maiForest.get(xi,yi)
     , 0  //What stands the value of u for
     , Rotation  //Rotation time
     , 0, 0, 0
//     ,1 , thinningForest.get(xi,yi)*sdMaxCoef, thinningForest.get(xi,yi)*sdMinCoeff8, 30, MAIRot
//     ,0 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, MAIRot
       ,0 , thinningForest.get(xi,yi)*sdMaxCoeff, thinningForest.get(xi,yi)*sdMinCoeff, 30, 1.
	 , 1  //reference of minrot
     , 0. //Flexibility of stocking degree
     , &ffsdMaxH, &ffsdMinH,321);
   

 //  if((X==23.75)&&(Y==49.75)){cout<<"x= "<<X<<"\t y= "<<Y<< "\t CurrRotTime= "<<rotationTimeCurr<<endl;}

     newCohort->createNormalForest(Rotation,0.,thinningForest.get(xi,yi));

	 
//        pair<g4m::ageStruct::v, g4m::ageStruct::v> newRes; // MG: results vector for the new (planted/afforested) forest
 // if (Country==224) {	
		//if((X==25.75)&&(Y==49.75)){
	//  if(asID==1887){
        //    double siteIndex = maiForest.get(xi,yi);
			//cout << "checkpoint"<<endl;
		//	extractData(*cohort, dat_all[asID], siteIndex, X, Y);
			// cout << "checkpoint end"<<endl;			

 
	/*	double bio = abBiomass0;
	    double rtime = rotationForest.get(xi,yi);		
		cout<<"MAI  "<<siteIndex<<endl;
		double i;			
		//for (i=10.; i<=221.; i+=10){
		for (i=0.; i<=221.; i+=1){
			if (outfile3.is_open()) {
				//cohort->TEST(i, siteIndex, bio, rtime);
				//outfile3<<asID<<"\t"<<X<<"\t"<<Y<<"\t"<<i<<"\t"<<cohort->GS_t_ha<<"\t"<<cohort->R_ha<<"\t"<<cohort->SD_t<<"\t"<<cohort->GS_t_Nat_ha<<endl; 
				outfile3<<asID<<"\t"<<X<<"\t"<<Y<<"\t"<<i<<"\t"<<cohort->GS_t_ha<<"\t"<<cohort->IncTot_ha<<"\t"<<cohort->R_ha<<"\t"<<cohort->GS_t_Tab_ha<<endl; 
			};
		};*/
		       
 
//cout<<"start singleCell"<<endl;
        dat singleCell;
//************************************
//**** Element for preparing output for GUI 
//************************************
//        singleCell.simUnit = sMap.getSIMU(X,Y);        
        singleCell.simUnit = sMap.getSimu(X,Y);        
//cout<<"X=\t"<<X<<"\t Y=\t"<<Y<<"\t simUnit=\t"<<sMap.getSimu(X,Y)<<endl;
//        singleCell.Rotation = Rotation;
        singleCell.Rotation = rotationForest.get(xi,yi);
        singleCell.LandAreaHa = LandAreaHa;
//        singleCell.potHarvest = potHarvest;
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
//        singleCell.ObiomassPrev = iter->CABOVEHA"][year];
        singleCell.ObiomassPrev = abBiomass0;
        singleCell.Obiomass0 = abBiomass0;                   // Modelled biomass at time 0, tC/ha
//        singleCell.rotBiomass = rotationForest.get(xi,yi);
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
//cout<<"push_back singleCell"<<"\t asID"<<asID<<endl;
        dat_all.push_back(singleCell);
        iter->asID = asID;
//if (asID >=100) cout << "asID = " << asID << "\txi= "<<xi<<"\t yi= "<<yi<< "\tx= "<<X<<"\t y= "<<Y<< endl;
		asID++;
        numAgeStruct++;
//if (thinningForest.get(xi,yi)>0){
//     CountriesManagedForHa.inc(Country,iter->BYEAR"][year],OforestShare*singleCell.LandAreaHa);
//     CountriesManagedCount.inc(Country,iter->BYEAR"][year],1);
//     }
//          
////     CountriesWoodHarvestM3Year.inc(Country,iter->BYEAR"][year],(harvWood*OforestShare+harvWoodNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa);
//     CountriesWoodLoosCYear.inc(Country,iter->BYEAR"][year],(harvWoodLost*OforestShare+harvWoodLostNew*singleCell.AforestSharePrev)*singleCell.LandAreaHa);

//if (asID % 1000 == 0) cout << "asID = " << asID << endl;

    }     // End for IF unprotected 
   }      // countries
  }      // Test only some regions 
    iter++;
kk++;
//if (kk % 100 == 0)
//    cout << "kk = " << kk << endl;
  }       // End of WHILE cells loop
//for (int i=1;i<=209;i++){
//cout<<"Country= "<<i<<"\t harvWoodM3 = "<< woodHarvest[i-1]<<"\t woodPotM3 = "<< woodPot[i-1]<<"\t harvLostC = "<< woodLost[i-1]<<endl;
//}
//cout<< "End initLoop"<<endl;
 }
