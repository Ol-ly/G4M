#ifndef G4M_DATASTRUCT_H
#define G4M_DATASTRUCT_H

#include <string>
#include "misc.h"
//#include "ipol.h"
namespace g4m{

struct dataStruct
 {
 
int	x	;
int	y	;
int asID;  // index of corresponding data in ageStruct vector

//  std::string classes;
ipol<float,float>	COUNTRY;
ipol<float,float>	POTVEG;
ipol<float,float>	PROTECT;
ipol<float,float>	USED;
ipol<float,float>	LANDAREA	;
ipol<float,float>	NPP	;
ipol<float,float>	POPDENS	;
ipol<float,float>	SAGRSUIT	;
ipol<float,float>	AGRSUIT	;
ipol<float,float>	PRICEINDEX	;
ipol<float,float>	BIOMASS	;
ipol<float,float>	FOREST	;
ipol<float,float>	FORESTTYPE	;
ipol<float,float>	SPECIESTYPE	;
ipol<float,float>	R	;
ipol<float,float>	GDP	;
ipol<float,float>	BUILTUP	;
ipol<float,float>	CROP	;
ipol<float,float>	FRACLONGPROD	;
ipol<float,float>	CORRUPTION	;
ipol<float,float>	SLASHBURN	;
ipol<float,float>	SPOPDENS	; // MG: added because it's used in forest_calculations
ipol<float,float>	BIOMASSBL	;
ipol<float,float>	DECHERB	;
ipol<float,float>	DECWOOD	;
ipol<float,float>	DECSOC	;
ipol<float,float>	CABOVEHA	;
ipol<float,float>	CBELOWHA	;
ipol<float,float>	CDEADHA	;
ipol<float,float>	CLITTERHA	;
ipol<float,float>	SOCHA	;
ipol<float,float>	FTIMBER	;
ipol<float,float>	IIASA_REGION	;
ipol<float,float>	POLESREG	;
ipol<float,float>	MAIE	;
ipol<float,float>	MAIN	;
ipol<float,float>	MANAGEDSHARE	;
ipol<float,float>	MANAGEDFLAG	;
ipol<float,float>	HARVESTCOSTS	;
ipol<float,float>	SIMUID	;
ipol<float,float>	COUNTRYREGMIX	;
ipol<float,float>	ROAD	;
ipol<float,float>	FORLOSS	;
ipol<float,float>	GLOBIOM_RESERVED	;
ipol<float,float>	AFFORMAX	;
//ipol<float,float>	ACCESSABILITY	;

//ipol<float,float>	LAND	;

// Constructor: initializing default values
  dataStruct()
   {
    x=0;
    y=0;
    asID=0;
//    classes = "AAA";
   }
 };


struct coeffStruct
 {
// Starting Year of simmulation
  short int bYear;
// Ending Year of simmulation
  short int eYear;
// Interaction between cells
  short int cellsInteract;
// Consider afforestation
  short int inclAffor;
// No pay
  short int noPay;
// Belowground biomass
  short int uBiomass;
// Litter
  short int litter;
// Soil organic carbon
  int SOC;
//******************************************************************************
//** parameters
//******************************************************************************
// Priceindex of reference country
  ipol<float,float> PriceIndexE;
// Minimum Landprice [cash/ha]
  ipol<float,float> PriceLandMinR;
// Maximum Landprice [cash/ha]
  ipol<float,float> PriceLandMaxR;
// Factor Carbon uptake (DIMA-Model)
  ipol<float,float> FCuptake;
// Comercial timbervolume per ton of carbon [m3/tC]
//  ipol<float,float> FTimber;
// HarvestingLosses (Share of losses during harvest)
  ipol<float,float> HarvLoos;
// Carbon price [Cash/tC] (9/25)
  ipol<float,float> PriceC;
// Share of Long-living products [0-1]
  ipol<float,float> FracLongProd;
// Decrease rate of long living products
  ipol<float,float> decRateL;
// Decrease rate of short living products
  ipol<float,float> decRateS;
// Share of SlashBurn deforestation [0-1]
  ipol<float,float> SlashBurn;
// Frequency of aid payments (PRICECAID) [Years]
  ipol<float,float> FreqAid;
// Aid Carbon Price [Cash/tC/FREQAID] (6)
  ipol<float,float> PriceCAid;
// Maximum rotation time im Years
  ipol<float,float> MaxRotInter;
// Minimum rotation time im Years
  ipol<float,float> MinRotInter;
// Baseline
  ipol<float,float> baseline;
// Maximum Timberprice [cash/m3]
  ipol<float,float> PriceTimberMaxR;
// Minimum Timberprice [cash/m3]
  ipol<float,float> PriceTimberMinR;
// Planting costs in reference country [Cash/ha]
  ipol<float,float> PlantingCostsR;
// Standardised Populationdensity [1-10]
  ipol<float,float> sPopDens;
 
  // Constructor: initializing default values
  coeffStruct()
   {
    bYear = 2000;
    eYear = 2005;
   }
 };

struct lwprice {
  int REG_ID;     
  string REG;
  string SCENARIO;
  ipol<float,float> LP;
  ipol<float,float> WP;
  };

struct wooddemand {
//  int REGMIX_ID;     
//  string REGMIX;
  int REG_ID;     
  string REG;
  string SCENARIO;
  ipol<float,float> WD;
  };


}
#endif
