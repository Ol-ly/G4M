//   Name:      g4m + interaction with GUI & GDataView
//   Author:    Andriy Bun, based on works of ...
//   Date:      ** December, 2009

// Modified: Mykola Gusti, 4 Jan 2010

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <ctime>
#include <algorithm>

#include "IntToStr.cpp"

// Elements for preparing output for GDataView:
#include "simUnitsMap.cpp"
#include "simUnitsData.h"
#include "tableData.cpp"

// Rest of includes:
#include "vector2d2.h"        
#include "readSettings_Europe.h"
#include "misc.h"                          // ipol class for interpolation
#include "increment.cpp"
#include "ageStruct_float.cpp"                 // max stocking degree is limited to about 1.2
#include "dataStruct_Europe_float.h"                 // data structures
#include "griddata3_MG.h"                  // GridData class                       (v1.1, v1.2 or v1.3)
#include "griddata2.h"                     // data of any type on x-y grid by Andriy Bun
#include "countryData_tw.cpp"                   // time window smoothing
#include "forest_GUI_Europe.h"                     // definitions
#include "dima_float.h"
#include "readInput_Europe.cpp" 
#include "defineSpecies.cpp"
#include "initManagedForest5_ageStruct_Europe_bef.cpp" // BEF is taken into account
#include "fm_cpol_v12_1_Europe.cpp"                    // works with "adjustManagedForestCountryGradual_fmp_v6_12_0_lulucf.cpp"!  NPV function is used instead of NPV50 ; rotTime of others is proportional to C price till 50$/tC
#include "adjustManagedForestCountryGradual_fmp_v6_12_1_Europe_bef.cpp"  // the same as v6_6_1 + NoFMcpol countries are detected + Initial FM sink in unmanaged forests is adjusted taking into account BEF
//------------------------------------------------------
#include "calc_Europa_2Aug2011_v021_bef_withFW.cpp" // The same as calc_Europa_2Aug2011_v021.cpp + affor, defor and FM emissions are estimated taking into account BEF + Area and Biomass of new forests over 20 y.o. are tracked
#include "MAI_country_mcpfe_maiMax30June11_ROHU_nc.cpp" // 19 July 2011 maie for Romania is returned to the previous value and for Hungary is decreased acording to demand of Giacomo Grassi
#include "hurdle_and_deforaffor_EUCLIMIT_03052013_nc_05_fi1_pt1.h"
#include "woodProductionIndexes_1July2010_nc_lulucf_cr.cpp" // countryRegMix is extended to include Croatia under No4
#include "forNPV_Europe.cpp"
#include "cprice.cpp"
#include "countryCodes_new_cr.cpp"
#include "listsToConsider_countryNew192_eu27_cr.cpp"    // Croatia is added to the EU countries
#include "harvestWood.cpp" //
#include "biomassExpansionFactors.cpp" // Biomass Expansion Factors (BEF) functions for tree species (1-8 as in defineSpecies) and forest types (9-11)
#include "extractData.cpp"


//******************************************************************************
//***********************************  MAIN  ***********************************
//******************************************************************************
int main (int argc, char* argv[])
 {
// we read two parameters from the command line - 'price path' and 'price value'
int PriceC = atoi(argv[1]);
cout<<PriceC<<endl;
  time_t start=time(NULL);
  // Read settings from GUI:
  readSettings();
  listsToConsider();


//*******************************
//** Reading coefficients file **
//*******************************
  readCoeff(coeff);
  byear = coeff.bYear;
  eyear = coeff.eYear;
//*****************************************************
//** Reading input files (with resolution 0.5 x 0.5) **
//*****************************************************
//*******************************************
//** Reading detailed input data from file **
//*******************************************
  woodProductionIndexes();
  dataDetStruct plots;                         // structure with data
  numRecords = readInputDet(plots);        // plots[<elNum>].<variable>[year]  


//*****************************************
//** Add other parameters to output here **
//*****************************************
  countryData countriesForestArea = countryData();                              // table output
//************************************************
    string suffix="_lulucf_ref_it1_t10_t_05_fi1_pt1_bef_12062013_3_1";

  string suffix0=suffix;  
double priceC = 0.;
   if (PriceC > 0){ suffix += "_Pco2_" + IntToStr(PriceC);     
                   priceC =  PriceC * exchangeRate;}
   if (PriceC==0) {suffix += "_Pco2_0";}   
   string fileOut = settings.outputPath +"\\" + suffix + ".txt";

  fff.open(fileOut.c_str(),ios::out);
  if (fff.is_open()) {   
fff<<"asID"<<"\tX"<<"\tY"<<"\tyear"<<"\tCountry"<<"\tOforShare"<<"\tLandAreaHa"<< "\tAforestShare"<<"\tcab"<<"\tabFM"<<"\tabFMnew"<<"\tMAI_M3Ha"<<"\trotMAI";
fff<<"\t harvestThM3HaOld"<<"\t harvestFctM3HaOld"<<"\t harvestThM3HaNew"<<"\t harvestFctM3HaNew"<<"\tharvTotCur"<<"\tharvMAI"<<"\trotBiomass"<<"\tforNPV_CH"<<"\tfVal"<<"\t aval";
fff<<"\t rotationTimeCurr"<<"\t rotationTimeCurrNew";
fff<<"\tTimberPrice"<<"\tthinningForest"<<"\tmanagedForest";
fff<<"\tEmissionsSOCAfforCur"<<"\t abBiomassAffor"<<"\t EmissionsLitterAfforCurHa";
fff<<"\t LitterAfforHa"<<"\t OforBm"<<endl;

}

  string suffix1="AgeClassStruct"; 
 // string suffix1="GrowingStockTab_tC";
	   string fileOut1 = settings.outputPath +"\\"+"details" + suffix1 + ".txt";
	   outfile.open(fileOut1.c_str(),ios::out);  
	   if (outfile.is_open()) {   
//outfile<<"asID"<<"\tX"<<"\tY"<<"\tAgeClass"<<"\tForestArea"<<"\tBiomass"<<endl;
		   outfile<<"asID"<<"\tX"<<"\tY"<<"\tBiomClass"<<"\tBiomass"<<endl;
	//	   outfile<<"asID"<<"\tX"<<"\tY"<<"\ti"<<"\tBiomassPerHa"<<"\tAgeClassArea"<<"\tValueYieldTable"<<endl;
	   };

string suffix2="DataHarvWood";
string fileOut2 = settings.outputPath + "\\" +suffix2 + ".txt";
outfile1.open(fileOut2.c_str(),ios::out);  
if (outfile1.is_open()) {   
		   outfile1<<"\tYear"<<"\tX"<<"\tY"<<"\tsawlogs"<<"\trestlogs" <<endl;
	   };

string suffix3="dataG4M_290316";
string fileOut3 = settings.outputPath + "\\" +suffix3 + ".txt";
outfile3.open(fileOut3.c_str(),ios::out);  
if (outfile3.is_open()) {   
		//   outfile3<<"asID"<<"\tX"<<"\tY"<<"\tAge"<<"\tGrowingStockPerHa_gBmt"<<"\tTotalIncrement_gGwlt"<<"\tRemovalsTh_ha"<<"\tSDBiomass_ha"<<endl;
	//outfile3<<"PARAMETER THINNING_HA(g,i) thinned wood [m3 per ha]"<<endl;
	//outfile3<<"/"<<endl;
	outfile3<<"X"<<"\tY"<<"\tremovals_period"<<"\tAgeClassArea"<<"\tRemovalsTh_ha"<<"\tTotalIncrement_gGwlt"<<endl;
	   };

string suffix4="_dataG4M_270715";
string fileOut4 = settings.outputPath + "\\" +suffix4 + ".txt";
outfile4.open(fileOut4.c_str(),ios::out);  
if (outfile4.is_open()) {   
		//   outfile3<<"asID"<<"\tX"<<"\tY"<<"\tAge"<<"\tGrowingStockPerHa_gBmt"<<"\tTotalIncrement_gGwlt"<<"\tRemovalsTh_ha"<<"\tSDBiomass_ha"<<endl;
	//outfile3<<"PARAMETER THINNING_HA(g,i) thinned wood [m3 per ha]"<<endl;
	//outfile3<<"/"<<endl;
	outfile4<<"X"<<"\tY"<<"\tAge"<<"\tGrowingStockTab_ha"<<"\tDiameter"<<endl;
	   };


//** Initializing forest cover array by gridcells **
//**************************************************
  griddata harvestGrid = griddata(ResLongitude,ResLatitude,0);
  griddata maiForest = griddata(ResLongitude,ResLatitude,0);
  griddata rotationForest = griddata(ResLongitude,ResLatitude,0);
  griddata rotationForestNew = griddata(ResLongitude,ResLatitude,0);
  griddata thinningForest = griddata(ResLongitude,ResLatitude,0);
  griddata thinningForestNew = griddata(ResLongitude,ResLatitude,0);
  griddata2<float> OforestShGrid = griddata2<float>(ResLongitude,ResLatitude,0); 

  griddata2<char> decisionGrid = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> managedForest = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> manageChForest = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> rotationType = griddata2<char>(ResLongitude,ResLatitude,0);
  griddata2<char> unmanaged = griddata2<char>(ResLongitude,ResLatitude,0);
  
 
  g4m::incrementTab *species[8];
  defineSpecies(species);


// Setup forest management parameters similar for all countries (cells)
  sws.clear(); 
  sws.insert(10, .0);
  sws.insert(30, .6);

  hlv.clear();
  hlv.insert(0, .0);
  hlv.insert(25, .7);

  hle.clear();
  hle.insert(0, .0);
  hle.insert(25, .75);

  sdMaxH.clear();
  sdMaxH.insert(0.,1.);
  sdMinH.clear();
  sdMinH.insert(0.,1.);

  {
    vector<double> idx;
    idx.push_back(0); //diameter
    idx.push_back(0); //stocking volume [tC stemwood/ha]
    idx.push_back(0); //Share of harvest (0 - 1)
    idx[0]=0; idx[1]=2; idx[2]=.3;
    cov.insert(idx, 4);
    idx[0]=40; idx[1]=30; idx[2]=.2;
    cov.insert(idx, 2);
   
    idx[0]=4; idx[1]=2; idx[2]=.1;
    dov.insert(idx,false);
    idx[0]=5; idx[1]=3; idx[2]=.2;
    dov.insert(idx,true);

    idx.pop_back();
    idx[0]=0; idx[1]=2;
    coe.insert(idx, 3);
    idx[0]=40; idx[1]=30;
    coe.insert(idx, 1);

    idx[0]=15; idx[1]=10;
    doe.insert(idx,false);
    idx[0]=16; idx[1]=11;
    doe.insert(idx,true);
  }

  g4m::ffipolm<double> ffcov(cov); //Thinning costs depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
  g4m::ffipolm<double> ffcoe(coe); //Harvesting costs depending on d and vol
  g4m::ffipolm<bool> ffdov(dov); //Do thinning (depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
  g4m::ffipolm<bool> ffdoe(doe); //Do final felling (depending on d and stocking volume per hectare)
// Georg's recommendation
//7: 0
//25: 1-HarvestingLosesCountry
//50: 1-(0.7*HarvestingLosesCountry)
//This is for hle, for the thinnings (hlv) multiply the values of
//harvestable biomass with 0.8.

   ffsws.overwrite(sws);
   ffhlv.overwrite(hlv);
   ffhle.overwrite(hle);
   ffsdMinH.overwrite(sdMinH);
   ffsdMaxH.overwrite(sdMaxH);
// Initializing some arrays
  MAI_country();
  hurdle_aff_deff();
//  woodHarvestStatCountry();
  woodProductionIndexes();
  forNPV_init();
  countryCodes();

  for (int i=0; i < NumberOfCountries; i++) {  
    EmissionsCurCountry[i]=0.; 
    EmissionsCurAfforCountry[i]=0.;
  }
  ageStructVector cohort_all;
  ageStructVector newCohort_all;
  cohort_all.reserve(numRecords);
  newCohort_all.reserve(numRecords);
  datGlobal dat_all;
if (GUIcontainers){
 if (settings.produceMaps){                    
  ASU = simUnitsData();                      // Initializing data for simulation units
  ASU.rename("G4M parameters");
  // Adding dimensions' parameters:
  if (settings.maps[0]) ASU.addDim("Scenario", IntToStr(PriceC));
  if (settings.maps[1]) ASU.addDim("Year", years);
  if (settings.maps[2]) ASU.addDim("Parameter", settings.parametersMap);
  set<string>::iterator it = settings.parametersMap.begin();
  cout << "parameters " << ASU.dimN()<< endl;
  while (it != settings.parametersMap.end())
  {
     cout << *it << endl;    
     it++;
  }
 }
}
 if (fmpol && PriceC>0) {biomass_bau.readFromFile("biomass_bau"+suffix0); // reading BAU files
                         NPVbau.readFromFile("NPVbau"+suffix0);}   

cout<<biomass_bau.size()<<"\t"<<biomass_bau[1].size()<<endl;
cout<<NPVbau.size()<<"\t"<<NPVbau[1].size()<<endl;

// ****************************************************************************

//******************************************************************************
//***************************** start calculations *****************************
//******************************************************************************


cout<< "Start initialising cohorts"<< endl;
  initManagedForest(plots, species, &ffcov,&ffcoe, &ffdov, &ffdoe, dat_all, cohort_all, newCohort_all, 
	                                  maiForest, thinningForest, rotationType, managedForest, rotationForest, harvestGrid, OforestShGrid);

//************************
//**** loop by years *****
//************************
float TH_10;  
  int year = byear;
    do {
      cout << "Processing year " << year << endl;
      int Age = year-byear;
      if (year > refYear && PriceC>0) {priceC = PriceC * exchangeRate;} 
      if (year > byear)
       {
cout << "Adjusting FM .."<< endl;
       adjustManagedForest(plots, species, cohort_all, 
              newCohort_all, dat_all, maiForest, 
              thinningForest, rotationForest, managedForest,
              rotationForestNew, thinningForestNew, manageChForest,
              rotationType, harvestGrid, year, unmanaged,priceC);
       }


//************************************
//** processing data from all plots **
//************************************
      vector<double> tmpBm;
      vector<double> tmpPr;
  
                   
cout << "starting cell calculations .."<< endl;
      dataDetStruct::iterator it = plots.begin();
      while (it != plots.end()) {
      if (regions.find(it->POLESREG[2000]) != regions.end()) { // Test only some regions
       if (countriesList.find(it->COUNTRY[2000]) != countriesList.end()) { // Test only some countries    
        if (it->PROTECT[2000]==0) {

            int asID = it->asID;
			

{
            calc(*it, species, *cohort_all[asID], *newCohort_all[asID], dat_all[asID], managedForest,
                 maiForest, rotationForest, rotationForestNew, thinningForest, thinningForestNew, 
                 harvestGrid,year, priceC, asID, OforestShGrid);


}                 
             
            if (fmpol&&bau&&PriceC==0&&year>refYear){ tmpBm.push_back(dat_all[asID].ObiomassPrev);      
                                                 tmpPr.push_back(profit);                                            
                                                                                          }                
          }                        // End if not protected ...
         }                         // End if country
        }                          // End if region
    
        it++;
      }                            // End loop by plots
            if (fmpol&&bau&&PriceC==0&&year>refYear){ biomass_bau.push(tmpBm);     
                                                 Profit_bau.push(tmpPr); }  
      
//******************************************
//** writing to containers for GUI output **
//******************************************

//******************************************
//** writing to containers for GUI output **
//******************************************
string strtmp;
vector<string> point;                      // Point to be considered
if (GUIcontainers){
    // filling vector of coordinates
    strtmp = IntToStr(PriceC);
    if (settings.maps[0]) point.push_back(strtmp);
}

if (GUIcontainers){
 if (settings.produceMaps){                    
cout<<"writing to GUI containers"<<endl;
      if (years.find(year) != years.end()) {
        strtmp = IntToStr(year);
        if (settings.maps[1]) point.push_back(strtmp);
        dataDetStruct::iterator it = plots.begin();
        while (it != plots.end()) {
//          if (it->PROTECT[2000]==0) {
          if ((it->PROTECT[2000]==0)&&((it->POTVEG[2000]<10 && it->FOREST[2000]>0 && it->MAIE[2000]>0)||(it->POTVEG[2000]<10 && (it->NPP[2000]>0)||it->MAIN[2000]>0))){                          
            if (regions.find(it->POLESREG[2000]) != regions.end()) { // Test only some regions
              int asID = it->asID;
              if (settings.parametersMap.find("stocking degree") != settings.parametersMap.end()) {
                point.push_back("stocking degree");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].SD, point);
                point.pop_back();
              }
              if (settings.parametersMap.find("em_fm_ab mtco2hayear") != settings.parametersMap.end()) {
                point.push_back("em_fm_ab mtco2hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].FMsink*(-1), point);
                point.pop_back();
              }              
              if (settings.parametersMap.find("em_fm_bm mtco2hayear") != settings.parametersMap.end()) {
                point.push_back("em_fm_bm mtco2hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].FMsink_Bm*(-1), point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("area_forest_old ha") != settings.parametersMap.end()) {
                point.push_back("area_forest_old ha");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].OforestShare*dat_all[asID].LandAreaHa, point);
                point.pop_back();
              }               
              if (settings.parametersMap.find("area_forest_new ha") != settings.parametersMap.end()) {
                point.push_back("area_forest_new ha");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].AforestShare*dat_all[asID].LandAreaHa, point);
                point.pop_back();
              }                      
              if (settings.parametersMap.find("cai m3ha") != settings.parametersMap.end()) {
                point.push_back("cai m3ha");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].CAI, point);
                point.pop_back();
              }            
              if (settings.parametersMap.find("area_df hayear") != settings.parametersMap.end()) {
                point.push_back("area_df hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].deforestHaYear, point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("area_af hayear") != settings.parametersMap.end()) {
                point.push_back("area_af hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].afforestHaYear, point);
                point.pop_back();
              }                          
              if (settings.parametersMap.find("harvest_total m3hayear") != settings.parametersMap.end()) {
                point.push_back("harvest_total m3hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].harvestTot, point);
                point.pop_back();
              }      
              if (settings.parametersMap.find("harvest_fc m3hayear") != settings.parametersMap.end()) {
                point.push_back("harvest_fc m3hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].harvestFcM3Ha, point);
                point.pop_back();
              }      
              if (settings.parametersMap.find("harvest_th m3hayear") != settings.parametersMap.end()) {
                point.push_back("harvest_th m3hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].harvestThM3Ha, point);
                point.pop_back();
              }                    
              if (settings.parametersMap.find("biom_fm tcha") != settings.parametersMap.end()) {
                point.push_back("biom_fm tcha");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].oforestBm, point);
                point.pop_back();
              }                    
              if (settings.parametersMap.find("biom_af tcha") != settings.parametersMap.end()) {
                point.push_back("biom_af tcha");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].aforestBm, point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("em_df_bm mtco2hayear") != settings.parametersMap.end()) {
                point.push_back("em_df_bm mtco2hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].emissionsD_Bm, point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("em_df_sl mtco2hayear") != settings.parametersMap.end()) {
                point.push_back("em_df_sl mtco2hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].emissionsD_S, point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("em_af_bm mtco2hayear") != settings.parametersMap.end()) {
                point.push_back("em_af_bm mtco2hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].emissionsA_Bm, point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("em_af_sl mtco2hayear") != settings.parametersMap.end()) {
                point.push_back("em_af_sl mtco2hayear");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].emissionsA_S, point);
                point.pop_back();
              }  
              if (settings.parametersMap.find("rotation year") != settings.parametersMap.end()) {
                point.push_back("rotation year");
                ASU.insert(dat_all[asID].simUnit, dat_all[asID].Rotation, point);
                point.pop_back();
              }                                             
// Add parameters you want to output here as for "Forest area":
// The parameters you may want to output must be defined in listOfParameters.txt
// file used by the GUI. The parameters for output must be defined in settings.ini 
// created by the GUI for g4m. 
            }                      // End if Country ...
          }                        // End if not protected
          it++;
        }
        if (settings.maps[1]) point.pop_back();            // clearing year from vector of coordinates                          // End loop by plots
      }                            // End loop by years
   }  // end  if (settings.produceMaps)
}
      year++;
    } while (year <= eyear);       // End loop by years

if (GUIcontainers){

set<string> countriesNameList;      // coutry Names to be considered in the table GUI
{
set<int>::iterator it;
set<string>::iterator itn;
for (it=countriesList.begin(); it!=countriesList.end(); it++){
string tmp = countryOrderName[countryCodeOrder[(*it) - 1]]; 				   	
countriesNameList.insert(tmp);
}   	
}
	   
;

    // Clearing vector of coordinates:
    
    vector<string> point;                      // Point to be considered

//*************************************
//** table data
//*************************************
  tableData obj;
  if (settings.produceTabs) {
cout<<"start writing to GUI table"<<endl;
    if (settings.tabs[0]) obj.addDim("Scenario", IntToStr(PriceC));
    if (settings.tabs[1]) obj.addDim("Year", years);
    if (settings.tabs[2]) obj.addDim("Country", countriesNameList);
    /*if (settings.tabs[3])*/ obj.addDim("Parameter", settings.parametersTable);
    point.clear();
    for (int prc = 0; prc < 1; prc++) {
      if (settings.tabs[0]) point.push_back(IntToStr(PriceC));
      for (int year = byear; year <= eyear; year++) {
        if (years.find(year) != years.end()) {
          if (settings.tabs[1]) point.push_back(IntToStr(year));
          for (int countryCode = 1; countryCode < 244; countryCode++) {
            if (countriesList.find(countryCode) != countriesList.end()) {
			  std::string country_code = countryOrderName[countryCodeOrder[countryCode-1]]; 
              if (settings.tabs[2]) point.push_back(country_code);
              if (settings.parametersTable.find("em_fm_ab mtco2year") != settings.parametersTable.end()) {
                point.push_back("em_fm_ab mtco2year");
                obj.insert(CountriesFM.get(countryCode, year)*(-1), point);
                point.pop_back();
              }
              if (settings.parametersTable.find("em_fm_bm mtco2year") != settings.parametersTable.end()) {
                point.push_back("em_fm_bm mtco2year");
                obj.insert(CountriesFMbm.get(countryCode, year)*(-1), point);
                point.pop_back();
              }       
              if (settings.parametersTable.find("area_forest_old ha") != settings.parametersTable.end()) {
                point.push_back("area_forest_old ha");
                obj.insert(CountriesOforCover.get(countryCode, year), point);
                point.pop_back();
              }              
              if (settings.parametersTable.find("area_forest_new ha") != settings.parametersTable.end()) {
                point.push_back("area_forest_new ha");
                obj.insert(CountriesNforCover.get(countryCode, year), point);
                point.pop_back();
              }                         
              if (settings.parametersTable.find("area_forest_used ha") != settings.parametersTable.end()) {
                point.push_back("area_forest_used ha");
                obj.insert(CountriesManagedForHa.get(countryCode, year), point);
                point.pop_back();
              }                                       
              if (settings.parametersTable.find("cai_old m3ha") != settings.parametersTable.end()) {
                point.push_back("cai_old m3ha");
                obj.insert(CountriesCAI.get(countryCode, year)/CountriesOforCover.get(countryCode, year), point);
                point.pop_back();
              }                         
              if (settings.parametersTable.find("cai_new m3ha") != settings.parametersTable.end()) {
                point.push_back("cai_new m3ha");
                obj.insert(CountriesCAI_new.get(countryCode, year)/CountriesNforCover.get(countryCode, year), point);
                point.pop_back();
              }                               
              if (settings.parametersTable.find("area_df hayear") != settings.parametersTable.end()) {
                point.push_back("area_df hayear");
                obj.insert(CountriesDeforHaYear.get(countryCode, year), point);
                point.pop_back();
              }       
              if (settings.parametersTable.find("area_af hayear") != settings.parametersTable.end()) {
                point.push_back("area_af hayear");
                obj.insert(CountriesAfforHaYear.get(countryCode, year), point);
                point.pop_back();
              }  
              if (settings.parametersTable.find("harvest_demand m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_demand m3year");
                obj.insert(CountryregWprod.get(countryCode, year), point);
                point.pop_back();
              }                            	                                 
              if (settings.parametersTable.find("harvest_total m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_total m3year");
                obj.insert(CountriesWoodHarvestM3Year.get(countryCode, year), point);
                point.pop_back();
              }                    
              if (settings.parametersTable.find("harvest_old_fc m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_old_fc m3year");
                obj.insert(CountriesWoodHarvestFc_oldM3Year.get(countryCode, year), point);
                point.pop_back();
              }                                           	
              if (settings.parametersTable.find("harvest_old_th m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_old_th m3year");
                obj.insert(CountriesWoodHarvestTh_oldM3Year.get(countryCode, year), point);
                point.pop_back();
              }                
              if (settings.parametersTable.find("harvest_new_fc m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_new_fc m3year");
                obj.insert(CountriesWoodHarvestFc_newM3Year.get(countryCode, year), point);
                point.pop_back();
              }                                           	
              if (settings.parametersTable.find("harvest_new_th m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_new_th m3year");
                obj.insert(CountriesWoodHarvestTh_newM3Year.get(countryCode, year), point);
                point.pop_back();
              }                                                 
              if (settings.parametersTable.find("harvest_fm m3hayear") != settings.parametersTable.end()) {
                point.push_back("harvest_fm m3hayear");
                obj.insert(CountriesWoodHarvestFmM3Year.get(countryCode, year)/CountriesManagedForHa.get(countryCode, year), point);
                point.pop_back();
              }                            	
               if (settings.parametersTable.find("harvest_df m3hayear") != settings.parametersTable.end()) {
                point.push_back("harvest_df m3hayear");
                obj.insert(CountriesWoodHarvestDfM3Year.get(countryCode, year)/CountriesDeforHaYear.get(countryCode, year), point);
                point.pop_back();
              }        
              if (settings.parametersTable.find("harvest_fm m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_fm m3year");
                obj.insert(CountryregWoodHarvestFmM3Year.get(countryCode, year), point);
                point.pop_back();
              }                            	
               if (settings.parametersTable.find("harvest_df m3year") != settings.parametersTable.end()) {
                point.push_back("harvest_df m3year");
                obj.insert(CountryregWoodHarvestDfM3Year.get(countryCode, year), point);
                point.pop_back();
              }                                                               
              if (settings.parametersTable.find("biom_fm tc") != settings.parametersTable.end()) {
                point.push_back("biom_fm tc");
                obj.insert(CountriesOforC_biom.get(countryCode, year), point);
                point.pop_back();
              }                            	
              if (settings.parametersTable.find("biom_stem_fm tc") != settings.parametersTable.end()) {
                point.push_back("biom_stem_fm tc");
                obj.insert(CountriesOfor_stem_C.get(countryCode, year), point);
                point.pop_back();
              }                            	
              if (settings.parametersTable.find("biom_ab_fm tc") != settings.parametersTable.end()) {
                point.push_back("biom_ab_fm tc");
                obj.insert(CountriesOfor_ab_C.get(countryCode, year), point);
                point.pop_back();
              }                
              if (settings.parametersTable.find("biom_af tc") != settings.parametersTable.end()) {
                point.push_back("biom_af tc");
                obj.insert(CountriesNforTotC.get(countryCode, year), point);
                point.pop_back();
              }       
              if (settings.parametersTable.find("biom_stem_af tc") != settings.parametersTable.end()) {
                point.push_back("biom_stem_af tc");
                obj.insert(CountriesNfor_stem_C.get(countryCode, year), point);
                point.pop_back();
              }       
              if (settings.parametersTable.find("biom_ab_af tc") != settings.parametersTable.end()) {
                point.push_back("biom_ab_af tc");
                obj.insert(CountriesNfor_ab_C.get(countryCode, year), point);
                point.pop_back();
              }       

              if (settings.parametersTable.find("em_df_bm mtco2year") != settings.parametersTable.end()) {
                point.push_back("em_df_bm mtco2year");
                obj.insert(CountriesDeforCYear_biom.get(countryCode, year), point);
                point.pop_back();
              }      
              if (settings.parametersTable.find("em_df_sl mtco2year") != settings.parametersTable.end()) {
                point.push_back("em_df_sl mtco2year");
                obj.insert(CountriesDeforCYear_soil.get(countryCode, year), point);
                point.pop_back();
              }      
              if (settings.parametersTable.find("em_af_bm mtco2year") != settings.parametersTable.end()) {
                point.push_back("em_af_bm mtco2year");
                obj.insert(CountriesAfforCYear_biom.get(countryCode, year)*(-1), point);
                point.pop_back();
              }                    
              if (settings.parametersTable.find("em_af_sl mtco2year") != settings.parametersTable.end()) {
                point.push_back("em_af_sl mtco2year");
                obj.insert(CountriesAfforCYear_soil.get(countryCode, year)*(-1), point);
                point.pop_back();
              }       
              /*if (settings.parametersTable.find("rotation_avg year") != settings.parametersTable.end()) {
                point.push_back("rotation_avg year");
                obj.insert(CountryRotation.getAvg(countryCode, year), point);
                point.pop_back();
              } */
			  if (settings.parametersTable.find("rotation_avg year") != settings.parametersTable.end()) {
                point.push_back("rotation_avg year");
                obj.insert(CountryRotationMng.getAvg(countryCode, year), point);
                point.pop_back();
              }
              if (settings.parametersTable.find("area_affor_acc ha") != settings.parametersTable.end()) {
                point.push_back("area_affor_acc ha");
                obj.insert(CountriesAfforAccumHa.get(countryCode, year), point);
                point.pop_back();
              }   

              if (settings.parametersTable.find("area_affor_acc_o20 ha") != settings.parametersTable.end()) {
                point.push_back("area_affor_acc_o20 ha");
                obj.insert(CountriesAfforCover20.get(countryCode, year), point);
                point.pop_back();
              }   
              if (settings.parametersTable.find("biom_af_o20 tc") != settings.parametersTable.end()) {
                point.push_back("biom_af_o20 tc");
                obj.insert(CountriesAfforTotC20.get(countryCode, year), point);
                point.pop_back();
              }       
// Add parameters you want to output here as for "Forest area":
// The parameters you may want to output must be defined in listOfParameters.txt
// file used by the GUI. The parameters for output must be defined in settings.ini 
// created by the GUI for g4m. 
            if (settings.tabs[2]) point.pop_back();
            }
          }
          if (settings.tabs[1]) point.pop_back();
        }
      }
      if (settings.tabs[0]) point.pop_back();
    }
  } 
cout << "****" << endl;
//*************************************
//** Saving data to files
//*************************************
  if (settings.produceTabs) {
    obj.rename("G4M results for countries");
    point.clear();
    if (settings.tabs[0]) point.push_back("Scenario");
    if (settings.tabs[1]) point.push_back("Year");
    if (settings.tabs[2]) point.push_back("Country");
    /*if (settings.tabs[3])*/ point.push_back("Parameter");
    obj.renameDims(point);
//cout << "parameters " << ASU.dimN()<< endl;
    obj.SaveToFile(settings.outputPath, "tabs_gui"+suffix);
  }  
//**********************************
  if (settings.produceMaps) {
    ASU.SaveToFile(settings.outputPath, "-EU"+suffix);
  }
//  system("pause");
 //******************************
cout << "written to GUI files"<<endl;
} // End if GUIcontainer
//  cout << "parameters " << ASU.dimN()<< endl;
//------------------------------------------------------------------------------

if (fmpol&&bau&&PriceC==0){
cout << "writing bau files"<<endl;
dataDetStruct::iterator it;
int year = refYear+1;
 do {
    it = plots.begin();
    vector<double> tmp_npv;    
    while (it != plots.end()) {
        if (it->PROTECT[2000]==0) {
              if (regions.find(it->POLESREG[2000]) != regions.end()) { // Test only some regions
                  int asID = it->asID;                  
                  double sum = 0.;
                  int j=year;
                  do {
                      sum += Profit_bau[j-refYear-1][asID]/(pow(1+it->R[2000],(j-year)));
                      j++;
//                      } while (j <= eyear);
                      } while ((j-year<=50)&&(j <= eyear));                      
                      tmp_npv.push_back(sum);
              }
        }
       it++;
   }
                  NPVbau.push(tmp_npv);
                  year++;
 } while (year <= eyear);
          

                biomass_bau.saveToFile("biomass_bau"+suffix0);
                NPVbau.saveToFile("NPVbau"+suffix0);
            
cout<< "finished writing bau files"<<endl;
 }
//------------------------------------------------------------------------------ 
//
if (countryOutput){
     CountriesNforCover.printToFile("newForestHa"+suffix,coeff.bYear,coeff.eYear,1);
     CountriesAfforHaYear.printToFile("afforHaYear"+suffix,coeff.bYear,coeff.eYear,1);          

  
     CountriesNforTotC.printToFile("afforSinkC" + suffix,coeff.bYear,coeff.eYear,1); 
     CountriesAfforCYear.printToFile("afforSinkGgCO2Year" + suffix,coeff.bYear,coeff.eYear,1); 
     CountriesAfforCYear_ab.printToFile("afforSinkGgCO2Year_ab" + suffix,coeff.bYear,coeff.eYear,1);  
     CountriesAfforCYear_bl.printToFile("afforSinkGgCO2Year_bl" + suffix,coeff.bYear,coeff.eYear,1);   
     CountriesAfforCYear_biom.printToFile("afforSinkGgCO2Year_biom" + suffix,coeff.bYear,coeff.eYear,1);
     CountriesAfforCYear_dom.printToFile("afforSinkGgCO2Year_dom" + suffix,coeff.bYear,coeff.eYear,1);    
     CountriesAfforCYear_soil.printToFile("afforSinkGgCO2Year_soil" + suffix,coeff.bYear,coeff.eYear,1);
    
// //--------- 
     CountriesOforCover.printToFile("OforestHa" + suffix,coeff.bYear,coeff.eYear,1);

     CountriesDeforHaYear.printToFile("deforHaYear" + suffix,coeff.bYear,coeff.eYear,1); 
     
     CountriesOfor_ab_C.printToFile("OforC_ab" + suffix,coeff.bYear,coeff.eYear,1);
     CountriesOforC_biom.printToFile("OforC_biom" + suffix,coeff.bYear,coeff.eYear,1);  
     CountriesDeforCYear.printToFile("deforEmGgCO2Year_tot" + suffix,coeff.bYear,coeff.eYear,1);   
     CountriesDeforCYear_bl.printToFile("deforEmGgCO2Year_bl" + suffix,coeff.bYear,coeff.eYear,1);   
     CountriesDeforCYear_ab.printToFile("deforEmGgCO2Year_ab" + suffix,coeff.bYear,coeff.eYear,1); 
     CountriesDeforCYear_biom.printToFile("deforEmGgCO2Year_biom" + suffix,coeff.bYear,coeff.eYear,1);     
     CountriesDeforCYear_dom.printToFile("deforEmGgCO2Year_dom" + suffix,coeff.bYear,coeff.eYear,1);     
     CountriesDeforCYear_soil.printToFile("deforEmGgCO2Year_soil" + suffix,coeff.bYear,coeff.eYear,1);      

////--------- 

     CountriesWoodHarvestM3Year.printToFile("harvest_m3_year" + suffix,coeff.bYear,coeff.eYear,1); 
     CountriesWoodHarvestPlusM3Year.printToFile("harvestPlus_m3_year" + suffix,coeff.bYear,coeff.eYear,1);      
     CountriesWoodLoosCYear.printToFile("harvest_lostGgCO2Year" + suffix,coeff.bYear,coeff.eYear,1);
     CountriesHarvLossesYear.printToFile("harvLossesYear" + suffix,coeff.bYear,coeff.eYear,1);
////---------     
     CountriesManagedForHa.printToFile("managedForHa" + suffix,coeff.bYear,coeff.eYear,1);
     CountriesManagedCount.printToFile("managedCount" + suffix,coeff.bYear,coeff.eYear,1); 
     
     CountriesMAI.printToFile("mai" + suffix,coeff.bYear,coeff.bYear,1,"AVG");
     CountriesCAI.printToFile("cai" + suffix,coeff.bYear,coeff.eYear,1);     // total current increment m3/year
     CountriesCAI_new.printToFile("caiNew" + suffix,coeff.bYear,coeff.eYear,1); // total current increment m3/year

     CountriesFM.printToFile("FM_GgCO2Year" + suffix,coeff.bYear,coeff.eYear,1);
     CountriesFMbm.printToFile("FMbm_GgCO2Year" + suffix,coeff.bYear,coeff.eYear,1);     
//------------------------------------------------------------------------------ 

     CountriesWprod.printToFile("wprod_test" + suffix,coeff.bYear,coeff.eYear,1);
//------------------------------------------------------------------------------     
     CountriesProfit.printToFile("profit" + suffix,coeff.bYear,coeff.eYear,1); 

////------------------------------------------------------------------------------
}
    fff.close(); 
//  } 

//-- Writing ASCII grid to file -------------------------------------------------
if (gridOutput){ 
oForShare_grid_00.PrintToFile("oForShare_0" + suffix);
oForShare_grid_50.PrintToFile("oForShare_50" + suffix);
oForShare_grid_100.PrintToFile("oForShare_100" + suffix);
}
//-----------------------------------------------------------------------------
                               // End loop by prices
  cout << "> Working time is " << difftime(time(NULL),start)/60. << " min." << endl;
 // system("pause");
 }
//******************************************************************************
// end main
//******************************************************************************
