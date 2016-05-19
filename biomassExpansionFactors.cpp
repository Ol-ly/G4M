// Biomass Expansion Factor (BEF) functions
// Mykola Gusti 11 April 2013
//The functions are from
//Teobaldelli M, Somogyi Z., Migliavacca M., Usoltsev V. (2009) Generalized functions of biomass expansion factors for conifers and broadleaved by stand age, growing stock and site index
//  Forest Ecology and Management, N257, 1004-1013
//Brown S. (1997) Estimating Biomass and Biomass Change of Tropical Forests: a Primer. (FAO Forestry Paper - 134, FAO 1997)
#include <math.h>

using namespace std;

//float BEF(int forType, float growingStockC, float tC_m3 = 4, float woodDensity = 1)
float BEF(int forType, float growingStockC, float tC_m3 = 4)
{
//	float growingStockBiomass = growingStock * woodDensity;
  float growingStock = growingStockC * tC_m3;
  float bef = 1;
  if (growingStock>0)
	{
	forType +=1;
	switch(forType)
		{
//1 / 0 - fir  
	     case 1:
                 bef = 1.069 + 1.919 / pow(growingStock,float(0.524));
				 if (bef>3.5) bef=3.5;
				 if (bef<1.1) bef=1.1;
	             break;

//2 / 1 - spruce
		 case 2:
                 bef = 1.204 + 0.903 * exp(-0.009 * growingStock);
				 if (bef>4) bef=4;
				 if (bef<1.1) bef=1.1;
				 break;

//3 / 2 - pine
		 case 3:
                 bef = 0.949 + 3.791/pow(growingStock,float(0.501));
				 if (bef>6) bef=6;
				 if (bef<1.1) bef=1.1;	
				 break;

//4 / 3 -  pinus Halepensis - pine
		 case 4:
                 bef = 0.949 + 3.791/pow(growingStock,float(0.501));
				 if (bef>6) bef=6;
				 if (bef<1.1) bef=1.1;	
				 break;

//5 / 4 - birch // alder // alnus incana
		 case 5:
	             bef = 1.105 + 9.793 / growingStock;
				 if (bef>1.6) bef=1.6;
				 if (bef<1.1) bef=1.1;
				 break;

//6 / 5 - beech
		 case 6:
	             bef = 1.197 + 0.386 * exp(-0.009 * growingStock);
				 if (bef>3.5) bef=3.5;
				 if (bef<1.1) bef=1.1;
				 break;
//7 / 6 - oak
		 case 7:
	             bef = 1.202 + 0.422 * exp(-0.013 * growingStock);
				 if (bef>3.5) bef=3.5;
				 if (bef<1.1) bef=1.1;
				 break;

//8 / 7 - larch
		 case 8:
	             bef = 1.023 + 2.058 / pow(growingStock,float(0.508));
				 if (bef>3.5) bef=3.5;
				 if (bef<1.1) bef=1.1;
				 break;
		}
	}     
	return bef;
}

float BEF(int forType, float growingStockC, float tC_m3 = 4, float woodDensity = 1)
{
//	float growingStockBiomass = growingStock * woodDensity;
  float growingStock = growingStockC * tC_m3;
  float bef = 1;
  if (growingStock>0)
	{
	forType +=1;
	switch(forType)
		{

// Tropics Broadleaf. 
// Brown S. (1997) Estimating Biomass and Biomass Change of Tropical Forests: a Primer. (FAO Forestry Paper - 134, FAO 1997)
//http://www.fao.org/docrep/w4095e/w4095e06.htm#3.1.2%20volume%20weighted%20average%20wood%20density%20%28wd%29
//
		 case 9:
                  if (growingStock * woodDensity < 190) {
					    bef = exp(3.213 - 0.506 * log(growingStock * woodDensity)); 
	              }else{bef = 1.74;}
				  break;
				  
// Broadleaf non-tropics Teobaldelli et al. 2009
		 case 10:
			      bef = 1.175 + 0.5 * exp(-0.018 * growingStock);
				  break;

// Conifers ALL, Teobaldelli et al. 2009
		 case 11:
	              bef = 1.172 + 0.739 * exp(-0.014 * growingStock);
				  break;
		}
	}     
	return bef;
}
