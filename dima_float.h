#ifndef DIMA_H_
#define DIMA_H_

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;
//////////
// DIMA //
//////////
class dima {
public:
  float forVal();      //Value of Forestry during multiple rotation
  float forValNC();      //Value of Forestry during multiple rotation
  float agrVal();      //Net present Value of Agriculture
  float amenVal();      //Net present Value of Agriculture
  float minPriceC();      //Get the minimum Carbonprice where Forest=Argic
  int setYear(int);
  float setForest(float);
  float forValExt();   //MG: Value of Forestry using External Timber price
  float forValComb();   //MG: Value of Forestry using combined G4M & External Timber prices
  float forValNCExt();      //MG: Value of Forestry during multiple rotation with 0 C price and using External wood price
  float agrVal2000();   //MG: Net present Value of Agriculture in the year 2000
  float plantingCosts();    //Costs to plant x ha of forest  
  dima(int ayear
       , g4m::ipol<float,float> &anpp          //npp [kg-C/m2/year]
       , g4m::ipol<float,float> &asPopDens   //Standardised (1-10) population density
       , g4m::ipol<float,float> &asAgrSuit   //Standardised (1-10) agricultural suitability
       , g4m::ipol<float,float> &apriceIndex  //Priceindex
       , g4m::ipol<float,float> &apriceIndex0 //Priceindex ref Country
       , g4m::ipol<float,float> &ar          //discount rate
       , g4m::ipol<float,float> &apriceC      //carbon price [$/tC]
       , g4m::ipol<float,float> &aplantingCosts0 //Planting costs ref Country [$/ha]
       , g4m::ipol<float,float> &apriceLandMin0   //Minimum Landprice in ref country [$/ha]
       , g4m::ipol<float,float> &apriceLandMax0  //Maximum Landprice in ref country [$/ha]
       , g4m::ipol<float,float> &amaxRotInter   //Maximal rotation intervall [years]
       , g4m::ipol<float,float> &aminRotInter       //Minimal rotation intervall [years]
       , g4m::ipol<float,float> &adecLongProd  //Decay rate for long lived prod
       , g4m::ipol<float,float> &adecShortProd  //Decay rate for short lived pro
       //Fraction of carbon stored in longterm products
       , g4m::ipol<float,float> &afracLongProd
       //Fraction of carbon substracted due to baseline considerations
       , g4m::ipol<float,float> &abaseline
       //Comercial timbervolume per ton of carbon [m3/tC]
//       , g4m::ipol<float,float> aftimber
       , g4m::ipol<float,float> &aftimber
//         , int aftimber
       //Maximal timber price in reference country [$/m3]
       , g4m::ipol<float,float> &apriceTimberMax0
       //Minimal timber price in reference country {$/m3]
       , g4m::ipol<float,float> &apriceTimberMin0
       , g4m::ipol<float,float> &afcuptake  //Factor of carbon uptake from npp
       , g4m::ipol<float,float> &agdp       //Gross domestic production
       , g4m::ipol<float,float> &aharvLoos  //Harvesting losses
       , float aforest      //Share of landarea used by forest [0-1]
       , g4m::ipol<float,float> awpricecorr //MG: Added for wood price correction
       , float awpricecorr0 //MG: Added for wood price correction in year 2000 when price of carbon=0
       , float arotInterM  //MG: Added rotation interval estimated from Georg's Forest Management Tool
       , float aharvWood   //MG: Added harvestable wood from Georg's Forest Management Tool
       ) : year(ayear)
	   , npp(anpp)
	   , sPopDens(asPopDens)
	   , sAgrSuit(asAgrSuit)
	   , priceIndex(apriceIndex)
	   , priceIndex0(apriceIndex0)
	   , r(ar)
	   , priceC(apriceC)
	   , plantingCosts0(aplantingCosts0)
	   , priceLandMin0(apriceLandMin0)
	   , priceLandMax0(apriceLandMax0)
	   , maxRotInter(amaxRotInter)
	   , minRotInter(aminRotInter)
	   , decLongProd(adecLongProd)
	   , decShortProd(adecShortProd)
	   , fracLongProd(afracLongProd)
	   , baseline(abaseline)
	   , ftimber(aftimber)
	   , priceTimberMax0(apriceTimberMax0)
	   , priceTimberMin0(apriceTimberMin0)
	   , fcuptake(afcuptake)
	   , gdp(agdp)
	   , harvLoos(aharvLoos)
	   , forest(aforest)
	   , wpricecorr(awpricecorr)
	   , wpricecorr0(awpricecorr0)
	   , rotInterM(arotInterM)
	   , harvWood(aharvWood)
  {
    ;
  }
  float priceTimber();      //Timber price
  float priceTimberExt();  //MG: Timber price External
  float priceTimberComb();  //MG: Timber price Combined (G4M + External)
  //private:
  int year;
  float rotInter();         //Rotation interval of a Forest in Years
  float woodHarvestVol();   //Harvest volume of the timber
  float cUptake();          //mean anual carbon uptake
  float vIncr();            //Harvestable wood-volum increment (m3/ha/year)
private:
  float beta();             //Fraction of carbon costs during harvest
  float cBenefit();         //Carbon benefit
  float priceHarvest();     //Price to harvest the timber
  float forestValueOne();   //Value of Forestry during one rotation
  float forestValueOneExt();//MG: Value of Forestry during one rotation  with External Timber price
  float forestValueOneComb();//MG: Value of Forestry during one rotation  with combined G4M & External Timber prices
  float forestValueOneNC();   //Value of Forestry during one rotation
  float forestValueOneNCExt();   //MG: Value of Forestry during one rotation with External Timber price
//  float plantingCosts();    //Costs to plant x ha of forest

  g4m::ipol<float,float> plantingCosts0;     //Costs to Plant 1ha of forests in ref country
  g4m::ipol<float,float> r;                  //Discount rate
  g4m::ipol<float,float> sAgrSuit;           //Standardised Agricultural Suitability
  g4m::ipol<float,float> sPopDens;           //Standardised Population Density (1-10)
  g4m::ipol<float,float> npp;                //neto primary production
  g4m::ipol<float,float> priceIndex;         //Price Index
  g4m::ipol<float,float> priceIndex0;        //Price Index of reference country
  g4m::ipol<float,float> priceC;             //Carbon Price
  g4m::ipol<float,float> priceLandMin0;      //Minimum Landprice in ref country
  g4m::ipol<float,float> priceLandMax0;      //Maximum Landprice in ref country
  g4m::ipol<float,float> maxRotInter;        //Maximal rotation intervall in years
  g4m::ipol<float,float> minRotInter;        //Minimal rotation intervall
  g4m::ipol<float,float> decLongProd;        //Decay rate for long lived products
  g4m::ipol<float,float> decShortProd;       //Decay rate for short lived products
  g4m::ipol<float,float> fracLongProd;       //Fraction of carbon stored in longterm products
  //Fraction of carbon substracted due to baseline considerations
  g4m::ipol<float,float> baseline;

//  g4m::ipol<float,float> ftimber;          //Comercial timbervolume per ton of carbon (m3/tC)
    g4m::ipol<float,float> ftimber;
//      int ftimber;

  g4m::ipol<float,float> priceTimberMax0;    //Maximal timber price in reference country
  g4m::ipol<float,float> priceTimberMin0;    //Minimal timber price in reference country
  g4m::ipol<float,float> fcuptake;            //Factor of carbon uptake from npp
  g4m::ipol<float,float> gdp;                 //Gross domestic production
  g4m::ipol<float,float> harvLoos;
  float forest;                //Share of landarea used by forest [0-1]
  g4m::ipol<float,float> wpricecorr;          //MG: Added for wood price correction
  float wpricecorr0;          //MG: Added for wood price correction in year 2000
  float rotInterM;            //MG: Added rotation interval estimated from Georg's Forest Management Tool
  float harvWood;             //MG: Added harvestable wood estimated from Georg's Forest Management Tool
};

float dima::setForest(float x) {
  forest = x;
  if(forest < 0.) {forest = 0.;}
  if(forest > 1.) {forest = 1.;}
  return(forest);
}

int dima::setYear(int i) {
  year = i;
  return(year);
}

float dima::cUptake() {     //mean anual carbon uptake (t-C/ha/year)
  //This value should depend on the NPP which is reduced by a factor
  //This factor should not be the same all over the world
  //Tropical takes not so much NPP because of:
  //   *)high rotation of leaves and litle brances
  //   *)Insects eat the leaves

  //  return(1.02*.13*(npp+.1)/10.);
//cout << "npp.g(year)= "<<npp.g(year)<<endl;
  return(npp.g(year) * 10. * fcuptake.g(year)); //kg/m2 -> t/ha

}

float dima::vIncr() {  //Harvestable wood-volume increment (m3/ha/year)
//cout << "ftimber.g(year)= "<<ftimber.g(year)<<endl;
  return(cUptake() * ftimber.g(year));
//  return(cUptake() * ftimber);

}

// MG: Use Georg's Optimal Rotation Time
float dima::rotInter() {          //Rotation interval of a Forest in Years

//  float t = 100;
//  float harvestVolume = 600. - abs(vIncr() -6.) * 50.;
//  if(cUptake() > 0.) {t = harvestVolume/vIncr();}
//  if(t < minRotInter.g(year)) {t = minRotInter.g(year);}
//  if(t > maxRotInter.g(year)) {t = maxRotInter.g(year);}

float t=rotInterM;  // MG: Use Georg's Optimal Rotation Time
//cout<<"t= "<<t<<endl;
return(t);

}

float dima::beta() {              //Fraction of carbon costs during harvest
  //Depends on fraction of short and long term products
  return(1. - decLongProd.g(year)/(decLongProd.g(year)+r.g(year))*fracLongProd.g(year)
	 - decShortProd.g(year)/(decShortProd.g(year)+r.g(year))
	 *(1.-fracLongProd.g(year)));
}

float dima::cBenefit() {       //Carbon benefit (Eq. 3)
//	cout << "priceC.g(year)"<<"\t"<<priceC.g(year)<<endl;
  return(priceC.g(year) * cUptake() * (1 - baseline.g(year)) *
	 ( ((1. - pow(1+r.g(year),-rotInter()) ) /r.g(year)) -
	   rotInter() * (1-beta()) * pow(1+r.g(year), -rotInter())));
}

float dima::woodHarvestVol() {  //Harvest volume of the timber during 1 rotation period
//  return(vIncr() * rotInter() * (1. - harvLoos.g(year)));
//cout << harvWood * rotInter() <<"   "<<vIncr() * rotInter() * (1. - harvLoos.g(year))<<endl;
return(harvWood * rotInter());  // MG: Georg's harvestable wood

  //Disturbances can also be mentioned
}

float dima::priceTimber() {     //Timber price internal
  //float c4 = (priceTimberMax0.g(year) - priceTimberMin0.g(year))/9.;
  //float c3 = priceTimberMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));

  float sfor = (1. - forest) * 9. + 1.;
  float c4 = (priceTimberMax0.g(year) - priceTimberMin0.g(year))/99.; 
  float c3 = priceTimberMin0.g(year) - c4;
  return((c3 + c4 * sPopDens.g(year) * sfor)
	 * priceIndex.g(year)/priceIndex0.g(year));

}

float dima::priceTimberExt() {     //MG: Timber price external
  //float c4 = (priceTimberMax0.g(year) - priceTimberMin0.g(year))/9.;
  //float c3 = priceTimberMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));

  float sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  float c4 = (priceTimberMax0.g(2000) - priceTimberMin0.g(2000))/99.;
  float c3 = priceTimberMin0.g(2000) - c4;
  return((c3 + c4 * sPopDens.g(2000) * sfor)
	 * priceIndex.g(2000)/priceIndex0.g(2000)
	 * wpricecorr.g(year)/wpricecorr0);
}

float dima::priceTimberComb() {     //MG: Combined timber price G4M+external
  //float c4 = (priceTimberMax0.g(year) - priceTimberMin0.g(year))/9.;
  //float c3 = priceTimberMin0.g(year) - c4;
  //return((c3 + c4 * sPopDens.g(year))
  //	 * priceIndex.g(year)/priceIndex0.g(year));

  float sfor = (1. - forest) * 9. + 1.;
// MG: use internal G4M wood price
//MG: Changed to external SawnLogsPrice
  float c4 = (priceTimberMax0.g(year) - priceTimberMin0.g(year))/99.;
  float c3 = priceTimberMin0.g(year) - c4;
  return((c3 + c4 * sPopDens.g(year) * sfor)
	 * priceIndex.g(year)/priceIndex0.g(year)
	 * wpricecorr.g(year)/wpricecorr0);
}


  //float shareFuelwood = 1./(1 + exp(-1.09525 + 0.30313*gdp.g(year)
  //			     + 1.50520*forest));
  //float harvestAmountYear = woodHarvestVol()/rotInter();
  //float harvFuel = harvestAmountYear * shareFuelwood;
  //if(harvFuel > 10.) {harvFuel = 10.;}
  //float harvInd = harvestAmountYear - harvFuel;
  //if(harvInd > 4.) {harvInd = 4.;}
  //float indWoodPrice = priceTimberMin0.g(year) + gdp.g(year)
  //  * exp(3.7619 - 1.098*forest - 0.5928 * harvInd);
  //float fuelWoodPrice = priceTimberMin0.g(year) + gdp.g(year)
  //  * exp(1.38493 - 0.97638*forest + 0.05744*harvFuel);
  //float price = fuelWoodPrice * shareFuelwood
  //  + indWoodPrice * (1. - shareFuelwood);
  //return(price);


float dima::priceHarvest() {     //Price to harvest the timber
  return(priceTimber() * .0);
//Beside harvesting costs also thinning costs, branc-removal,... can be
//considered
}

float dima::plantingCosts() {    //Costs to plant 1 ha of forest
  //return(plantingCosts0.g(year)*priceIndex.g(year)/priceIndex0.g(year));
  //Maybe these costs do not ocure on the second rotation intervall
  //because of *)natural regeneration *)coppice forests
  float plantrate = (vIncr()-3.)/6.;
  if(plantrate > 1.) {plantrate = 1.;}
  if(plantrate < 0.) {plantrate = 0.;}
  return(plantrate*plantingCosts0.g(year)*priceIndex.g(year)/
	 priceIndex0.g(year));
}


float dima::forestValueOneNC() {//Value of Forestry one rotation NoCarbonPrice // Changed to one year!!!
  //return(-plantingCosts() +
  // (priceTimber()
  // -priceHarvest())*woodHarvestVol()*pow(1+r.g(year), -rotInter()));
  return((-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimber()
    -priceHarvest())*woodHarvestVol())/rotInter());  // MG: I deleted *pow(1+r.g(year), -rotInter()) to make it similar to forestValueOne
}

float dima::forestValueOneNCExt() {//MG: Value of Forestry one rotation NoCarbonPrice using External wood price // Changed to one year!!!
  //return(-plantingCosts() +
  // (priceTimber()
  // -priceHarvest())*woodHarvestVol()*pow(1+r.g(year), -rotInter()));

  return((-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimberExt()
    -priceHarvest())*woodHarvestVol())/rotInter());  // MG: I deleted *pow(1+r.g(year), -rotInter()) to make it similar to forestValueOne
}

float dima::forValNC() { //Value of Forestry multiple rotation No Carbon Price // Changed to multiple years!!!!
  float currF = 0.;
  float npvSum = 0.;
  int j=0;
  do {
       currF = 1./(pow(1+r.g(year),j));
       npvSum += currF;
       j += modTimeStep;
//cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
     } while (currF > 0.00001 && j < 400);
 return(npvSum * forestValueOneNC() * modTimeStep);
//  return(forestValueOneNC()/(1.-pow(1.+r.g(year), -rotInter()))); // Georg's definition (like in Kindermann et al. 2007)
}

float dima::forValNCExt() { //MG: Value of Forestry multiple rotation No Carbon Price using External wood price// Changed to multiple years!!!!
  float currF = 0.;
  float npvSum = 0.;
  int j=0;
  do {
       currF = 1./(pow(1+r.g(year),j));
       npvSum += currF;
       j += modTimeStep;
//cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
     } while (currF > 0.00001 && j < 400);
 return(npvSum * forestValueOneNCExt() * modTimeStep);
//  return(forestValueOneNCExt()/(1.-pow(1.+r.g(year), -rotInter())));// Georg's definition (like in Kindermann et al. 2007)
}

float dima::forestValueOne() { //Value of Forestry during one rotation (Eq.1) // Changed to 1 year!!!!
  //return(-plantingCosts() +
  // (priceTimber()
  //  -priceHarvest())*woodHarvestVol()*pow(1+r.g(year), -rotInter())
  // +cBenefit());
  return((-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimber()
    	  -priceHarvest())*woodHarvestVol()
  	 +cBenefit())/rotInter()
);
}

float dima::forestValueOneExt() { //MG: Value of Forestry during one rotation External// Changed to 1 year!!!!
  //return(-plantingCosts() +
  // (priceTimber()
  //  -priceHarvest())*woodHarvestVol()*pow(1+r.g(year), -rotInter())
  // +cBenefit());


	return((-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimberExt()
    	  -priceHarvest())*woodHarvestVol()
  	 +cBenefit())/rotInter()
);
}

float dima::forestValueOneComb() { //MG: Value of Forestry during one rotation, combination of G4M + External// Changed to 1 year!!!!
  //return(-plantingCosts() +
  // (priceTimber()
  //  -priceHarvest())*woodHarvestVol()*pow(1+r.g(year), -rotInter())
  // +cBenefit());


	return((-plantingCosts() +
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimberComb()
    	  -priceHarvest())*woodHarvestVol()
  	 +cBenefit())/rotInter()
);


}

float dima::forVal() { //Value of Forestry during multiple rotation (Eq.4)// Changed to multiple years!!!!
  float currF = 0.;
  float npvSum = 0.;
  int j=0;
  do {
       currF = 1./(pow(1+r.g(year),j));
       npvSum += currF;
       j += modTimeStep;
//cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
     } while (currF > 0.00001 && j < 400);
 return(npvSum * forestValueOne() * modTimeStep);
//  return(forestValueOne()/(1.-pow(1.+r.g(year), -rotInter())));
//  return(forestValueOne()/(1.-pow(1.+r.g(year), -rotInter()))+cBenefit());  
}




float dima::forValExt() { //MG: Value of Forestry during multiple rotation with External Timber price // Changed to multiple years!!!!
  float currF = 0.;
  float npvSum = 0.;
  int j=0;
  do {
       currF = 1./(pow(1+r.g(year),j));
       npvSum += currF;
       j += modTimeStep;
//cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
     } while (currF > 0.00001 && j < 400);
 return(npvSum * forestValueOneExt() * modTimeStep);
//  return(forestValueOneExt()/(1.-pow(1.+r.g(year), -rotInter())));
//  return(forestValueOne()/(1.-pow(1.+r.g(year), -rotInter()))+cBenefit());  
}

float dima::forValComb() { //MG: Value of Forestry during multiple rotation with combined (G4M+External) Timber price // Changed to multiple years!!!!
  float currF = 0.;
  float npvSum = 0.;
  int j=0;
  do {
       currF = 1./(pow(1+r.g(year),j));
       npvSum += currF;
       j += modTimeStep;
//cout<<"j= "<<j<<"   currF= "<<currF<<"   npvSum= "<<npvSum<<endl;
     } while (currF > 0.00001 && j < 400);
 return(npvSum * forestValueOneComb() * modTimeStep);
//  return(forestValueOneExt()/(1.-pow(1.+r.g(year), -rotInter())));
//  return(forestValueOne()/(1.-pow(1.+r.g(year), -rotInter()))+cBenefit());  
}

float dima::agrVal() {   //Net present Value of Agriculture (Eq.5)
  float priceLevel = priceLandMin0.g(year)
    * priceIndex.g(year)/priceIndex0.g(year);
  //Importance of Population density
  float popImp = (log(priceLandMax0.g(year))
		   - log(priceLandMin0.g(year)))/(2. * log(10.));
  //Importance of the Suitable for Agriculture
  float agrImp = popImp;
//cout << "priceLandMin0.g(year)= "<<priceLandMin0.g(year)<<"\t priceIndex.g(year)= "<<priceIndex.g(year)<<"\t priceIndex0.g(year)= "<<priceIndex0.g(year);
//cout << "\t priceLevel= "<<priceLevel<<"\t sPopDens.g(year)= "<<sPopDens.g(year)<<"\t popImp= "<<popImp<<"\t sAgrSuit.g(year)= "<<sAgrSuit.g(year); 
//cout << "\t agrVal()= "<<priceLevel * pow(sPopDens.g(year),popImp) * pow(sAgrSuit.g(year),agrImp)<<endl;
  return(priceLevel * pow(sPopDens.g(year),popImp)
	 * pow(sAgrSuit.g(year),agrImp));

}

////MG: Attantion: agrVal changed! Only 2000 value is estimated here
float dima::agrVal2000() {   //MG: Net present Value of Agriculture in the year 2000
  float priceLevel = priceLandMin0.g(2000)
    * priceIndex.g(2000)/priceIndex0.g(2000);
  //Importance of Population density
  float popImp = (log(priceLandMax0.g(2000))
		   - log(priceLandMin0.g(2000)))/(2. * log(10.));
  //Importance of the Suitable for Agriculture
  float agrImp = popImp;
  return(priceLevel * pow(sPopDens.g(2000),popImp)
	 * pow(sAgrSuit.g(year),agrImp));

}

float dima::amenVal() {   //Value of amenity
  float priceLevel = priceLandMin0.g(year)
    * priceIndex.g(year)/priceIndex0.g(year);
  //Importance of Population density
  float popImp = (log(priceLandMax0.g(year))
		   - log(priceLandMin0.g(year)))/(2. * log(10.));
  //Importance of GDP
  float gdpImp = popImp;
  return(priceLevel * pow(sPopDens.g(year),popImp)
	 * pow(gdp.g(year),gdpImp));
}

float dima::minPriceC() {
  //Get the minimum Carbonprice where Forest=Argic (Eq.6)
  return((agrVal() * (1-pow(1+r.g(year),(-rotInter()))) +
	  plantingCosts() -
//   	 (priceTimber() //MG:  changed for external price correction
  	 (priceTimber()	  
	  
	   -priceHarvest())*woodHarvestVol()*pow(1+r.g(year), (-rotInter()))
	  )
	 /
	 (cUptake() * (1. - baseline.g(year)) *
	  ( ((1. - pow(1+r.g(year),(-rotInter())) ) /r.g(year)) -
	    rotInter() * (1.-beta()) * pow(1+r.g(year), (-rotInter())))));
}
////////////// END DIMA //////////////////////////////////////

#endif
