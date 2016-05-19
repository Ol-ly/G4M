//******************************************************************************
// forestry NPV function
//******************************************************************************
#ifndef forNPV_cpp
#define forNPV_cpp
//#include <cmath>


// Value of Forestry during multiple rotation

double forNPV(g4m::dataStruct &plot, double forest, int year, int rotation, double maiV, 
              double priceW, double harvestedW, double harvestedArea, g4m::ipol<double,double> wpricecorr, double wpricecorr0, double Cprice)
 {

double cUptake = plot.NPP[2000] * 10. * coeff.FCuptake[2000]; //kg/m2 -> t/ha

double beta =                //Fraction of carbon costs during harvest
  //Depends on fraction of short and long term products
  (1. - coeff.decRateL[2000]/(coeff.decRateL[2000]+plot.R[2000])*plot.FRACLONGPROD[2000]
	 - coeff.decRateS[2000]/(coeff.decRateS[2000]+plot.R[2000])
	 *(1.-plot.FRACLONGPROD[2000]));


double cBenefit = //Carbon benefit (Eq. 3)
//	cout << "Cprice.g(year)"<<"\t"<<Cprice.g(year)<<endl;
  (Cprice * cUptake * (1. - coeff.baseline[2000]) *
	 ( ((1. - pow(1.+plot.R[2000],-rotation) ) /plot.R[2000]) -
	   rotation * (1.-beta) * pow(1.+plot.R[2000], -rotation)));




    //MG: Timber price external
  //double c4 = (priceWMax0.g(year) - priceWMin0.g(year))/9.;
  //double c3 = priceWMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));

  double sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  double c4 = (coeff.PriceTimberMaxR[2000] - coeff.PriceTimberMinR[2000])/99.;
  double c3 = coeff.PriceTimberMinR[2000] - c4;
double priceWExt =  ((c3 + c4 * plot.SPOPDENS[2000] * sfor)
	 * plot.PRICEINDEX[2000]/coeff.PriceIndexE[2000]
	 * wpricecorr.g(year)/wpricecorr0);


double priceHarvest =     //Price to harvest the timber
   (priceWExt * .0);
//Beside harvesting costs also thinning costs, branc-removal,... can be
//considered


    //Costs to plant 1 ha of forest
  //return(plantingCosts0.g(year)*priceIndex.g(year)/priceIndex0.g(year));
  //Maybe these costs do not ocure on the second rotation intervall
  //because of *)natural regeneration *)coppice forests
  double plantrate = (maiV-3.)/6.;
  if(plantrate > 1.) {plantrate = 1.;}
  if(plantrate < 0.) {plantrate = 0.;}
  double plantingCosts = (plantrate*coeff.PlantingCostsR[2000]*harvestedArea*plot.PRICEINDEX[year]/
	 coeff.PriceIndexE[2000]);




//double forestValueOneExt =  //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
//  //return(-plantingCosts() +
//  // (priceW()
//  //  -priceHarvest())*harvestedW()*pow(1+r.g(year), -rotation())
//  // +cBenefit());
//  ((-plantingCosts +
////   	 (priceW() //MG:  changed for external price correction
//  	 (priceWExt
//    	  -priceHarvest)*harvestedW
//  	 +cBenefit)/rotation
//);

double forestValueOneExt =  //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
  ((-plantingCosts +
  	 (priceWExt - priceHarvest)*harvestedW +cBenefit/rotation)
);


//double forValExt() { //MG: Value of Forestry during multiple rotation with External Timber price // Changed to multiple years!!!!
  double currF = 0.;
  double npvSum = 0.;
  int j=0;
  do {
       currF = 1./(pow(1+plot.R[2000],j));
       npvSum += currF;
       j++;
//cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
     } while (currF > 0.0000001 && j < 100);
cout<< "forestValueOneExt = "<<forestValueOneExt<<"  npvSum = "<<npvSum<<"  j = "<<j<<"  harvestedArea = "<<harvestedArea<<endl; 
 return(npvSum * forestValueOneExt);
//  return(forestValueOneExt()/(1.-pow(1.+r.g(year), -rotation())));
//  return(forestValueOne()/(1.-pow(1.+r.g(year), -rotation()))+cBenefit());  

}
//*****************************************************************************
double forNPVfd(g4m::dataStruct &plot, double forest, int year, int rotation, double maiV, 
              g4m::ipol<double,double> wpricecorr, double wpricecorr0, double harvestedW, double harvestedArea)
 {


double priceHarvest =     //Price to harvest the timber
//   (priceWExt * .0);
//Beside harvesting costs also thinning costs, branc-removal,... can be
//considered
//                 plot.HARVCOST[year];
                   0.;
    //Costs to plant 1 ha of forest
  //return(plantingCosts0.g(year)*priceIndex.g(year)/priceIndex0.g(year));
  //Maybe these costs do not ocure on the second rotation intervall
  //because of *)natural regeneration *)coppice forests
  double plantrate = (maiV-3.)/6.;
  if(plantrate > 1.) {plantrate = 1.;}
  if(plantrate < 0.) {plantrate = 0.;}
  double plantingCosts = (plantrate*coeff.PlantingCostsR[2000]*harvestedArea*plot.PRICEINDEX[year]/
	 coeff.PriceIndexE[2000]);

    //MG: Timber price external
  //double c4 = (priceWMax0.g(year) - priceWMin0.g(year))/9.;
  //double c3 = priceWMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));

  double sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  double c4 = (coeff.PriceTimberMaxR[2000] - coeff.PriceTimberMinR[2000])/99.;
  double c3 = coeff.PriceTimberMinR[2000] - c4;
double priceWExt =  ((c3 + c4 * plot.SPOPDENS[2000] * sfor)
	 * plot.PRICEINDEX[2000]/coeff.PriceIndexE[2000]
	 * wpricecorr.g(year)/wpricecorr0);


//double forestValueOneExt =  //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
//  //return(-plantingCosts() +
//  // (priceW()
//  //  -priceHarvest())*harvestedW()*pow(1+r.g(year), -rotation())
//  // +cBenefit());
//  ((-plantingCosts +
////   	 (priceW() //MG:  changed for external price correction
//  	 (priceWExt
//    	  -priceHarvest)*harvestedW
//  	 +cBenefit)/rotation
//);

double forestValueOneExt =  //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
  (-plantingCosts +
  	 (priceWExt - priceHarvest)*harvestedW);

return(forestValueOneExt);

////double forValExt() { //MG: Value of Forestry during multiple rotation with External Timber price // Changed to multiple years!!!!
//  double currF = 0.;
//  double npvSum = 0.;
//  int j=0;
//  do {
//       currF = 1./(pow(1+plot.R[2000],j));
//       npvSum += currF;
//       j++;
////cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
//     } while (currF > 0.0000001 && j < 100);
//cout<< "forestValueOneExt = "<<forestValueOneExt<<"  npvSum = "<<npvSum<<"  j = "<<j<<"  harvestedArea = "<<harvestedArea<<endl; 
// return(npvSum * forestValueOneExt);
////  return(forestValueOneExt()/(1.-pow(1.+r.g(year), -rotation())));
////  return(forestValueOne()/(1.-pow(1.+r.g(year), -rotation()))+cBenefit());  

}

//*****************************************************************************

//*****************************************************************************
double forNPVfdc(g4m::dataStruct &plot, double forest, int year, int rotation, double maiV, 
              g4m::ipol<double,double> wpricecorr, double wpricecorr0, double harvestedW, double harvestedArea, double biomass_cur)
 {


double priceHarvest =     //Price to harvest the timber
//   (priceWExt * .0);
//Beside harvesting costs also thinning costs, branc-removal,... can be
//considered
//                 plot.HARVCOST[year];
                   0.;
    //Costs to plant 1 ha of forest
  //return(plantingCosts0.g(year)*priceIndex.g(year)/priceIndex0.g(year));
  //Maybe these costs do not ocure on the second rotation intervall
  //because of *)natural regeneration *)coppice forests
 double plantingCosts = 0.;
 if (harvestedArea>0){
  double plantrate = (maiV-3.)/6.;
  if(plantrate > 1.) {plantrate = 1.;}
  if(plantrate < 0.) {plantrate = 0.;}
  plantingCosts = (plantrate*coeff.PlantingCostsR[2000]*harvestedArea*plot.PRICEINDEX.g(year)/
	 coeff.PriceIndexE[2000]);
 }
    //MG: Timber price external
  //double c4 = (priceWMax0.g(year) - priceWMin0.g(year))/9.;
  //double c3 = priceWMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));
 double priceWExt = 0.;
 if (harvestedW>0){
  double sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  double c4 = (coeff.PriceTimberMaxR[2000] - coeff.PriceTimberMinR[2000])/99.;
  double c3 = coeff.PriceTimberMinR[2000] - c4;
  double priceWExt =  ((c3 + c4 * plot.SPOPDENS[2000] * sfor)
	 * plot.PRICEINDEX[2000]/coeff.PriceIndexE[2000]
	 * wpricecorr.g(year)/wpricecorr0);
 }

//cout<<"Cbenefit start";
   int asID = plot.asID;
//cout<<"asID = "<<asID<<endl;
//cout<<"yearID = "<<year-refYear-1<<endl;
   double Cbenefit = coeff.PriceC[2000] * ((biomass_cur - biomass_bau[year-refYear-1][asID])*1.2); // extention to total biomass

//   double Cbenefit = cprice[cscenario].g(year)*plot.CORRUPTION.g(byear)*
//                   ((biomass_cur-biomass_bau[year-refYear-1][asID])*1.2); // extention to total biomass

//cout<<"year=\t"<<year<<"\tprice=\t"<<cprice[cscenario].g(year)<<"\tcorr=\t"<<plot.CORRUPTION.g(byear)<<"\tCbenefit(end)=\t"<<Cbenefit<<endl;

double forestValueOneExt =  //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
  (-plantingCosts +
  	 (priceWExt - priceHarvest)*harvestedW)
        + Cbenefit;


//cout<<"wpricecorr.g(year)= "<<wpricecorr.g(year)<<endl;
//cout<<"wpricecorr0= "<< wpricecorr0<<endl;
//cout<<"forest = "<<forest<<endl;
//cout <<"plantingCosts = "<<plantingCosts<<endl;
//cout<<"priceHarvest = "<<priceHarvest<<endl;
//cout<<"priceWExt = "<<priceWExt<<endl;
//cout<<"harvestedW = "<<harvestedW<<endl;
//cout<< "forestValueOneExt = "<< forestValueOneExt<<endl;
return(forestValueOneExt);

}
//*****************************************************************************

//*****************************************************************************
double forNPVfdc(g4m::dataStruct &plot, double forest, int year, int rotation, double maiV, 
              g4m::ipol<float,float> wpricecorr, float wpricecorr0, double harvestedW, double harvestedArea, double biomass_cur)
 {


double priceHarvest =     //Price to harvest the timber
//   (priceWExt * .0);
//Beside harvesting costs also thinning costs, branc-removal,... can be
//considered
//                 plot.HARVCOST[year];
                   0.;
    //Costs to plant 1 ha of forest
  //return(plantingCosts0.g(year)*priceIndex.g(year)/priceIndex0.g(year));
  //Maybe these costs do not ocure on the second rotation intervall
  //because of *)natural regeneration *)coppice forests
 double plantingCosts = 0.;
 if (harvestedArea>0){
  double plantrate = (maiV-3.)/6.;
  if(plantrate > 1.) {plantrate = 1.;}
  if(plantrate < 0.) {plantrate = 0.;}
  plantingCosts = (plantrate*coeff.PlantingCostsR[2000]*harvestedArea*plot.PRICEINDEX.g(year)/
	 coeff.PriceIndexE[2000]);
 }
    //MG: Timber price external
  //double c4 = (priceWMax0.g(year) - priceWMin0.g(year))/9.;
  //double c3 = priceWMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));
 double priceWExt = 0.;
 if (harvestedW>0){
  double sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  double c4 = (coeff.PriceTimberMaxR[2000] - coeff.PriceTimberMinR[2000])/99.;
  double c3 = coeff.PriceTimberMinR[2000] - c4;
  double priceWExt =  ((c3 + c4 * plot.SPOPDENS[2000] * sfor)
	 * plot.PRICEINDEX[2000]/coeff.PriceIndexE[2000]
	 * wpricecorr.g(year)/wpricecorr0);
 }

//cout<<"Cbenefit start";
   int asID = plot.asID;
//cout<<"asID = "<<asID<<endl;
//cout<<"yearID = "<<year-refYear-1<<endl;
   double Cbenefit = coeff.PriceC[2000] * ((biomass_cur - biomass_bau[year-refYear-1][asID])*1.2); // extention to total biomass

//   double Cbenefit = cprice[cscenario].g(year)*plot.CORRUPTION.g(byear)*
//                   ((biomass_cur-biomass_bau[year-refYear-1][asID])*1.2); // extention to total biomass

//cout<<"year=\t"<<year<<"\tprice=\t"<<cprice[cscenario].g(year)<<"\tcorr=\t"<<plot.CORRUPTION.g(byear)<<"\tCbenefit(end)=\t"<<Cbenefit<<endl;

double forestValueOneExt =  //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
  (-plantingCosts +
  	 (priceWExt - priceHarvest)*harvestedW)
        + Cbenefit;


//cout<<"wpricecorr.g(year)= "<<wpricecorr.g(year)<<endl;
//cout<<"wpricecorr0= "<< wpricecorr0<<endl;
//cout<<"forest = "<<forest<<endl;
//cout <<"plantingCosts = "<<plantingCosts<<endl;
//cout<<"priceHarvest = "<<priceHarvest<<endl;
//cout<<"priceWExt = "<<priceWExt<<endl;
//cout<<"harvestedW = "<<harvestedW<<endl;
//cout<< "forestValueOneExt = "<< forestValueOneExt<<endl;
return(forestValueOneExt);

}
//*****************************************************************************


// Calculate NPV for the period from current year till end of the modeling
double npv_calc(g4m::dataStruct &iter, g4m::ageStruct &cohortTmp, double maiV, int year, int rotation, string regprice, string regprice0,double OforestShare, short int used)
  {
          double currF = 0.;
          double npvSum = 0.;
          int j=0;
          int N = biomass_bau.size();  
          double harvestTmp = 0.;
          double harvestedArea = 0.;          
          do {
             pair<g4m::ageStruct::v, g4m::ageStruct::v> resTmp = cohortTmp.aging();
             double bm = cohortTmp.getBm();
             if (used >0){
                 double sawnW = resTmp.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                 double restW = resTmp.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                 double sawnThW = resTmp.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                 double restThW = resTmp.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                 double bmH = resTmp.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                 double bmTh = resTmp.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                 double harvL = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) ; // MG: harvest residues for the set (old) forest tC/ha      
                 harvestTmp = (sawnW + restW + sawnThW + restThW + resUse*harvL)  * iter.FTIMBER.g(byear);
                 harvestedArea = cohortTmp.getArea(0);
             }


//cout<< "NPVone start "<<endl;        
             double NPVone = forNPVfdc(iter, OforestShare, year+j, rotation, maiV, wprice[regprice], wprice[regprice0].g(2000), harvestTmp, harvestedArea, bm);
//cout<<"NPVone (end) = "<<NPVone<<endl;                
               currF = NPVone/(pow(1+iter.R[2000],j));
               npvSum += currF;
               j++;
        //cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
//             } while (abs(currF) > 0.001 && year + j <= coeff.eYear);
//             } while (abs(currF) > 0.01 && (year + j) <= coeff.eYear && (year+j) <= (refYear+N));
             } while ((year + j) <= coeff.eYear && (year+j) <= (refYear+N));
        return(npvSum);
  }



//*****************************************************************************
//*****************************************************************************
// Calculate NPV for the 10 year period
double npv_calc10(g4m::dataStruct &iter, g4m::ageStruct &cohortTmp, double maiV, int year, int rotation, string regprice, string regprice0,double OforestShare, short int used)
  {
          double currF = 0.;
          double npvSum = 0.;
          int j=0;
          int N = biomass_bau.size();  
          double harvestTmp = 0.;
          double harvestedArea = 0.;          
          do {
             pair<g4m::ageStruct::v, g4m::ageStruct::v> resTmp = cohortTmp.aging();
             double bm = cohortTmp.getBm();
             if (used >0){
                 double sawnW = resTmp.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                 double restW = resTmp.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                 double sawnThW = resTmp.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                 double restThW = resTmp.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                 double bmH = resTmp.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                 double bmTh = resTmp.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                 double harvL = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) ; // MG: harvest residues for the set (old) forest tC/ha        
                 harvestTmp = (sawnW + restW + sawnThW + restThW + resUse*harvL)  * iter.FTIMBER.g(byear);
                 harvestedArea = cohortTmp.getArea(0);
             }


//cout<< "NPVone start "<<endl;        
             double NPVone = forNPVfdc(iter, OforestShare, year+j, rotation, maiV, wprice[regprice], wprice[regprice0].g(2000), harvestTmp, harvestedArea, bm);
//cout<<"NPVone (end) = "<<NPVone<<endl;                
               currF = NPVone/(pow(1+iter.R[2000],j));
               npvSum += currF;
               j++;
        //cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
//             } while (abs(currF) > 0.001 && year + j <= year+10 && year + j <= coeff.eYear);
//             } while (abs(currF) > 0.01 && year + j <= year+10 && (year+j) <= (refYear+N));
             } while (year + j <= year+10 && (year+j) <= (refYear+N));
        return(npvSum);
  }
//******************************************************************************
//*****************************************************************************
// Calculate NPV for the 50 year period
double npv_calc50(g4m::dataStruct &iter, g4m::ageStruct &cohortTmp, double maiV, int year, int rotation, string regprice, string regprice0,double OforestShare, short int used)
  {
          double currF = 0.;
          double npvSum = 0.;
          int j=0;
          double harvestTmp = 0.;
          double harvestedArea = 0.;   
          int N = biomass_bau.size();       
          do {
             pair<g4m::ageStruct::v, g4m::ageStruct::v> resTmp = cohortTmp.aging();
             double bm = cohortTmp.getBm();
             if (used >0){
                 double sawnW = resTmp.second.sw;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
                 double restW = resTmp.second.rw;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
                 double sawnThW = resTmp.first.sw;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
                 double restThW = resTmp.first.rw;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
                 double bmH = resTmp.second.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
                 double bmTh = resTmp.first.bm;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
                 double harvL = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) ; // MG: harvest residues for the set (old) forest tC/ha        
                 harvestTmp = (sawnW + restW + sawnThW + restThW + resUse*harvL)  * iter.FTIMBER.g(byear);
                 harvestedArea = cohortTmp.getArea(0);
                 }


//cout<< "NPVone start "<<endl;        
             double NPVone = forNPVfdc(iter, OforestShare, year+j, rotation, maiV, wprice[regprice], wprice[regprice0].g(2000), harvestTmp, harvestedArea, bm);
//cout<<"NPVone (end) = "<<NPVone<<endl;                
               currF = NPVone/(pow(1+iter.R[2000],j));
               npvSum += currF;
               j++;
        //cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
//             } while (currF > 0.001 && year + j <= year+50 && year + j <= coeff.eYear);
//             } while (abs(currF) > 0.0001 && j <= 50 && j <= N);             
//cout<<"year=\t"<<year<<"\tj=\t"<<j<<"\tcurrF=\t"<<currF<<endl;
//             } while (abs(currF) > 0.01 && j <= 50 && (year+j) <= (refYear+N));             
             } while (j <= 50 && (year+j) <= (refYear+N));             
        return(npvSum);
  }
//******************************************************************************


void forNPV_init(void)
//initialisation of forestry NPV arrays
{
for (int i=0;i<NumberOfCountries;i++){
// NPV curves
    minRotNPV[i]=0.;   // country average NPV if max MAI rotation is applied to all forests
    minMedNPV[i]=0.;   // NPV at rotation between min and medium rotation
    medRotNPV[i]=0.;   // NPV at rotation between min and medium rotation
    medMaxNPV[i]=0.;     // NPV at rotation between medium and max biomass rotation
    maxRotNPV[i]=0.;  // country average NPV if max biomass rotation is applied to all forests
    minRot[i]=0.;    // country average min rotation (max harvest)
    maxRot[i]=0.;    // country average max biomass rotation
//---------------------
}
    
     }

#endif
