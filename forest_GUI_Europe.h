using namespace std;

//#include "interpol.h"
#include "dat.h"
//#include "initialisation.cpp"
//******************************************************************************
// types
//******************************************************************************
//typedef vector< g4m::ipol<double,double> > cellCol;
//typedef vector<cellCol> dataArray;
typedef vector<g4m::dataStruct> dataDetStruct;
typedef vector<g4m::lwprice> lwpricesStruct;
typedef vector<g4m::wooddemand> wooddemandStruct;

g4m::coeffStruct coeff;
map<string, g4m::ipol<double,double> > lprice; //datamap for land price corrections for current price  scenario (GLOBIOM)
//map<string, g4m::ipol<double,double> > wprice; //datamap for wood price corrections for current price  scenario (GLOBIOM)
map<string, g4m::ipol<float,float> > wprice; //datamap for wood price corrections for current price  scenario (GLOBIOM)
map<string, g4m::ipol<double,double> > wprod;  //datamap for wood production corrections for POLES regions	(GLOBIOM)
map<string, g4m::ipol<double,double> > cprice; // datamap for carbon price
typedef vector<g4m::ageStruct *> ageStructVector;
typedef vector<dat> datGlobal;
//******************************************************************************
// containers of data
//******************************************************************************
bool GUIcontainers = true;
//simUnitsMap sMap = simUnitsMap("D:\MGusti\CurrentWork\georgPrgs\dima\DeforAforCCurves_growth\G4M_EU\g4m_Europe\Debug\simu.bin");
//simUnitsMap sMap = simUnitsMap("D:/MGusti/CurrentWork/georgPrgs/dima/DeforAforCCurves_growth/G4M_EU/g4m_Europe/Debug/simu.bin");
//simUnitsMap sMap = simUnitsMap("E:/POLITECH/ASPIRANTURA/_G4M_versions/g4m_WoodTypes_UA_1/WT_1/Release/simu.bin");
simUnitsMap sMap = simUnitsMap("simu.bin");
simUnitsData ASU;

set<int> regions;
set<int> years;             // select years for results output
set<int> toAdjust;          // coutry where FM to be adjusted
set<int> doneList;          // countries already adjusted
set<int> countriesList;      // coutry to be considered
set<int> countryregList;      // country and region mixture to be considered
set<int> countriesFmcpol;     // List of Annex-1 countries for FMcpol
set<int> countriesNoFmcpol;   // List of countries where it's impossible to match demanded wood production in current year




//******************************************************************************
// constants and variables
//******************************************************************************
#ifdef unix
string homeDir = "./data/";
#else
string homeDir = "data\\";
#endif
//int ResLatitude, ResLongitude;    // resolutions of model
int eyear, byear;
const double GridStepLat = 0.5;   // step by latitude
const double GridStepLon = 0.5;   // step by longitude
// resolution of model  
const int ResLatitude = int(floor(180/GridStepLat));
const int ResLongitude = int(floor(360/GridStepLon));

const int nYears = 61;
const int NumberOfCountryregmix = 50; // number of POLES regions mixed with EU27 countries
const int NumberOfCountries = 244; // new country codes
const double unitConv = 1e-6; // unit conversion: 1e-6 - Mt
const int modTimeStep = 1;

//  double defIncome = 0.;
    //Get optimal rotation time
    //0 .. Highest average increment
    //1 . .Maximum average Biomass
    //2 .. Highest possible age
    //3 .. Maximum harvest at final cut
    //4 .. Average Maximum harvest at final cut
          int optimMAI = 0;
          int optimMaxBm = 1;
//          int optimMaxBmTh = 3;
          int optimHarvFin = 3;
          int optimHarvAve = 4;

		  double sdMinCoeff = 1; // minimum stocking degree multiplier: sdMin = SD * sdMinCoeff
		  double sdMaxCoeff = 1; // maximum stocking degree multiplier: sdMax = SD * sdMaxCoeff

double MAI_CountryUprotect[NumberOfCountries+1];
double MAI_CountryAll[NumberOfCountries+1];
double MAI_countryregmix_up_avg[NumberOfCountryregmix+1];
double woodHarvestStat[NumberOfCountries+1];
double Hurdle_opt[NumberOfCountries+1];
double afforRate_opt[NumberOfCountries+1];
double deforRate_opt[NumberOfCountries+1];
double EmissionsCurCountry[NumberOfCountries+1];
double EmissionsCurAfforCountry[NumberOfCountries+1];
short int countryNwp[NumberOfCountries+1];
short int yearNwp[25];							//OT: original value 25
double countryLosses[NumberOfCountries+1];
double FM_sink_stat[NumberOfCountries+1];
double FMs[NumberOfCountries+1];
short int countryRegion[NumberOfCountries+1];
short int countryRegMix[NumberOfCountries+1];

int numRecords = 0; // number of records in the input data file
//int country_asId[210][2]; //country - asId start and end  
//int xy_asID[3000][2]; // xi and yi correspondence to asID
int numAgeStruct = 0;
//int usedNumAgeStruct = 0;

short int countryCodeOrder[NumberOfCountries+1];
char counrtyOrderISO[NumberOfCountries+1][4]; 
string countryOrderName[NumberOfCountries+1];
string countryRegName[NumberOfCountryregmix+1];



//g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm // share of harvestable sawnwood per m3 (diameter, share)
//g4m::ipol<double,double> hlv;   //1-Ernteverluste Vornutzung // loosses after first prefinal cut (diameter, share of harvesting loses) ?
//g4m::ipol<double,double> hle;   //1-Ernteverluste Endnutzung // losses at final cut (diameter, share of harvesting loses)?
//g4m::ipol<double,double> dbv;  //Dekungsbeitrag vornutzung   // income per m3 for thinning (diameter,income)
//g4m::ipol<double,double> dbe;  //Dekungsbeitrag endnutzung   //  income per m3 for final harvest (diameter,income)
  g4m::ipol<double,double> sws;  //Schnittholzanteil an Vfm
  g4m::ipol<double,double> hlv;   //1-Ernteverluste Vornutzung
  g4m::ipol<double,double> hle;   //1-Ernteverluste Endnutzung
  g4m::ipol<vector<double>,double> cov;  //costs vornutzung
  g4m::ipol<vector<double>,double> coe;  //Costs endnutzung
  g4m::ipol<vector<double>,bool> dov;  //Do vornutzung
  g4m::ipol<vector<double>,bool> doe;  //Do endnutzung
  g4m::ipol<double,double> sdMaxH;  //sdMaxH
  g4m::ipol<double,double> sdMinH;  //sdMinH

//  initialisation(&sdMaxH,&sdMinH,&cov,&dov,&coe);
  
  g4m::ffipol<double> ffsws(0); //Sawnwood share of harvested wood depending on dbh
  g4m::ffipol<double> ffhlv(0); //1-harvesting losses thinning (Vornutzung) (depending on d) in relation to standing timber (Vorratsfestmeter)
  g4m::ffipol<double> ffhle(0); //1-harvesting losses final felling (Endnutzung) (depending on d) in relation to standing timber (Vorratsfestmeter)
  //g4m::ffipolm<double> ffcov(cov); //Thinning costs depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
  //g4m::ffipolm<double> ffcoe(coe); //Harvesting costs depending on d and vol
  //g4m::ffipolm<bool> ffdov(dov); //Do thinning (depending on d and removed volume per hectare) in relation to standing timber (Vorratsfestmeter)
  //g4m::ffipolm<bool> ffdoe(doe); //Do final felling (depending on d and stocking volume per hectare)
  g4m::ffipol<double> ffsdMaxH(0); //Stockindegree depending on max tree height
  g4m::ffipol<double> ffsdMinH(0); //Stockindegree depending on max (min?) tree height


/*
double LinPrice2050[51];
float CubPrice2050[51];
float MixPrice2050[51];
float LinPrice2030[51];
float CubPrice2030[51] ;
float MixPrice2030[51];
float LinPrice2020[51] ;
float CubPrice2020[51] ;
float MixPrice2020[51];
float LinPrice2015[51] ;
float CubPrice2015[51] ;
float MixPrice2015[51];
float LinPrice2010[51] ;
float CubPrice2010[51] ;
float MixPrice2010[51];
*/
// NPV curves
double profit = 0.;

float minRotNPV[NumberOfCountries+1];   // country average NPV if max MAI rotation is applied to all forests
float minMedNPV[NumberOfCountries+1];   // NPV at rotation between min and medium rotation
float medRotNPV[NumberOfCountries+1];   // NPV at rotation between min and medium rotation
float medMaxNPV[NumberOfCountries+1];     // NPV at rotation between medium and max biomass rotation
float maxRotNPV[NumberOfCountries+1];  // country average NPV if max biomass rotation is applied to all forests
float minRot[NumberOfCountries+1];    // country average min rotation (max harvest)
float maxRot[NumberOfCountries+1];    // country average max biomass rotation
bool forNPVcurves = false;
bool forNPVcuvvesDyn = false;
bool fmpol = false; // For testing FM response to C price incentive; requires bin files with BAU biomass and NPV
bool bau = false; // Write bin files with BAU biomass and NPV

//if (fmpol){
 vector2d NPVbau = vector2d(46000);
 vector2d biomass_bau = vector2d(46000);
 vector2d Profit_bau = vector2d(46000);
//} 
 int maxDiffCountry = 0;
 double harvDiff[NumberOfCountries+1];
 griddata NPVcGrid = griddata(ResLongitude,ResLatitude,0);  // grid to store current NPV of forest if management is adjusted

int refYear = 2015; // policies start the next year
//---------------------

// Adjusting (increasing) FM sink by disturbance management in "unmanaged" forests 
bool adjustFMsink = true;
//---------------------
double exchangeRate = 1.47; // Euro -> USD exchange rate average for 2008 (all prices in the model are in USD)
string cscenario="";
//int PriceCiS[12] = {0,10,20,30,50,70,100,200,300,500,1000,0};
//int PriceCiS[26] = {0,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,250,300,1000};
int PriceCiS[14] = {0,5,10,15,20,25,30,35,40,45,50, 70, 100, 150};
double deflator = 0.8807; // Deflator to convert 2000 USD prices to 1995 USD prices
double resUse = 0.; // share of harvest residuals that can be used for bioenergy
double tolerance = 0.03;//0.05;//0.1; // tolerance of forest management adjustment to match domestic wood demand
bool cellInteract = true;

//*****************************************************************************
// country outputs
//*****************************************************************************
  bool countryOutput = false;                                                                               

  countryData CountriesNforCover = countryData();
  countryData CountriesNforTotC = countryData();
  countryData CountriesNfor_stem_C = countryData(); 
  countryData CountriesNfor_ab_C = countryData(); 
  countryData CountriesAfforHaYear = countryData();
  countryData CountriesAfforAccumHa = countryData();

  countryData CountriesAfforCYear = countryData();  
  countryData CountriesAfforCYear_ab = countryData();  
  countryData CountriesAfforCYear_bl = countryData();  
  countryData CountriesAfforCYear_biom = countryData();  
  countryData CountriesAfforCYear_dom = countryData();    
  countryData CountriesAfforCYear_soil = countryData();  

  countryData CountriesAfforCover20 = countryData();
  countryData CountriesAfforTotC20 = countryData();
//---------  
  countryData CountriesOforCover = countryData();
  countryData CountriesDeforHaYear = countryData();  
    
  countryData CountriesOfor_ab_C = countryData();
  countryData CountriesOfor_stem_C = countryData();
  countryData CountriesOforC_biom = countryData(); 
  countryData CountriesDeforCYear = countryData();  
  countryData CountriesDeforCYear_bl = countryData();   
  countryData CountriesDeforCYear_ab = countryData(); 
  countryData CountriesDeforCYear_biom = countryData();
  countryData CountriesDeforCYear_dom = countryData();   
  countryData CountriesDeforCYear_soil = countryData();     
//---------  
  countryData CountriesWoodHarvestM3Year = countryData();    
  countryData CountriesWoodHarvestPlusM3Year = countryData(); 
  countryData CountriesWoodHarvestFmM3Year = countryData();
  countryData CountriesWoodHarvestDfM3Year = countryData();
  countryData CountriesWoodLoosCYear = countryData();   
  countryData CountriesHarvLossesYear = countryData();

  countryData CountriesWoodHarvestFc_oldM3Year = countryData();
  countryData CountriesWoodHarvestTh_oldM3Year = countryData();
  countryData CountriesWoodHarvestFc_newM3Year = countryData();
  countryData CountriesWoodHarvestTh_newM3Year = countryData();
  countryData CountriesWoodHarvestSawn = countryData ();
  countryData CountriesWoodHarvestRest = countryData ();
//---------
  countryData CountriesManagedForHa = countryData();     
  countryData CountriesManagedCount = countryData();  
  
  countryData CountriesMAI = countryData();    
  countryData CountriesCAI = countryData();    
  countryData CountriesCAI_new = countryData();      
  countryData CountriesFM = countryData();   
  countryData CountriesFMbm = countryData();     
//---------------------------------------------------------------------------  
  countryData CountryRotation =  countryData(); 
  countryData CountryRotationMng = countryData();
//----------
  countryData CountriesWprod =  countryData();  // test wood production input data
//-----------
  countryData CountriesProfit =  countryData();  // profit due to selling  harvested wood
//---------------------------------------------------------------------------  
  countryData CountryregWoodHarvestM3Year =  countryData(); 
  countryData CountryregWoodHarvestFmM3Year = countryData();
  countryData CountryregWoodHarvestDfM3Year = countryData();  
  countryData CountryregWprod =  countryData();
  //---------------------------------------------------------------------------  
  countryData CountryregRotation =  countryData(); 

  countryData CountryregMaxHarvest =  countryData(); // maximal possible sustainable harvest
//**************************
//Output file  
//**************************
  ofstream fff;  
  ofstream outfile;
  ofstream outfile1;
  ofstream outfile3;
  ofstream outfile4;
  
  
//******************************************************************************
// ASCII grid output
bool gridOutput = false; // to output maps (ascii grids)

//    griddata oForShare_grid_1990 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata oForShare_grid_00 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata oForShare_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata oForShare_grid_100 = griddata(360/GridStepLon,180/GridStepLat,-9999);            
	    
    griddata harvestm3_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999);     
    griddata harvestm3_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999); 
    griddata harvestm3_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);     
        
    griddata bmabtC_of_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999);         
    griddata bmabtC_of_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999);   
    griddata bmabtC_of_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);             

    griddata bmabtC_nf_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999);         
    griddata bmabtC_nf_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999);   
    griddata bmabtC_nf_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);

    griddata fmgGco2_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999); 
    griddata fmgGco2_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999);     
    griddata fmgGco2_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);         

    griddata mai_m3ha_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata mai_m3ha_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata mai_m3ha_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);        

    griddata SD_grid_1990 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata SD_grid_00 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata SD_grid_10 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata SD_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata SD_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata SD_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999);    

    griddata RL_grid_20 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata RL_grid_30 = griddata(360/GridStepLon,180/GridStepLat,-9999);        
    griddata RL_grid_50 = griddata(360/GridStepLon,180/GridStepLat,-9999); 


//******************************************************************************
// functions
//******************************************************************************
int readInputDet(dataDetStruct &);
vector<double> procPlots(g4m::dataStruct &,double,int,double, double,double, double, double);
vector<double> calcPlots(double [],int);
vector<double> plotsData(g4m::dataStruct &,int);
//void initManagedForest(dataDetStruct &, g4m::incrementTab &, datGlobal &,
//                       ageStructVector &, ageStructVector &,              
//                       griddata &, griddata2<char> &, griddata2<char> &, 
//                       griddata &, griddata &,griddata &);

void initManagedForest(dataDetStruct &, g4m::incrementTab *species[8],  
	                   g4m::ffipolm<double> *, 
                       g4m::ffipolm<double> *, 
                       g4m::ffipolm<bool> *, 
                       g4m::ffipolm<bool> *, 	                   
	                   datGlobal &,
                       ageStructVector &, ageStructVector &,              
                       griddata &, griddata2<char> &, griddata2<char> &, 
                       griddata &, griddata &,griddata &);


// void adjustManagedForest(dataDetStruct &data_all, g4m::incrementTab &fi, ageStructVector &cohort_all, 
//              ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
//              griddata &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
//              griddata &rotationForestNew, griddata &thinningForestNew, griddata2<char> &manageChForest,
//              griddata2<char> &rotationType, griddata &harvestGrid, int year); 

// void adjustManagedForest(dataDetStruct &data_all, g4m::incrementTab &fi, ageStructVector &cohort_all, 
//              ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
//              griddata &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
//              griddata &rotationForestNew, griddata &thinningForestNew, griddata2<char> &manageChForest,
//              griddata2<char> &rotationType, griddata &harvestGrid, int year, griddata2<char> &unmanaged);


void adjustManagedForest(dataDetStruct &data_all, g4m::incrementTab &, ageStructVector &cohort_all, 
              ageStructVector &newCohort_all, datGlobal &dat_all, griddata &maiForest, 
              griddata &thinningForest, griddata &rotationForest, griddata2<char> &managedForest,
              griddata &rotationForestNew, griddata &thinningForestNew, griddata2<char> &manageChForest,
              griddata2<char> &rotationType, griddata &harvestGrid, int year, griddata2<char> &unmanaged, double);

//void fm_cpol(dataDetStruct &, g4m::incrementTab &, ageStructVector &, 
//              ageStructVector &, datGlobal &, griddata &, 
//              griddata &, griddata &, griddata2<char> &,
//              griddata &, griddata &, griddata2<char> &,
//              griddata2<char> &, griddata &, int , griddata2<char> &, double,double &,double);
              
void fm_cpol(dataDetStruct &, g4m::incrementTab &, ageStructVector &, 
              ageStructVector &, datGlobal &, griddata &, 
              griddata &, griddata &, griddata2<char> &,
              griddata &, griddata &, griddata2<char> &,
              griddata2<char> &, griddata &, int , griddata2<char> &, double,double &,double,
              set<int> &, griddata &,griddata2<int> &,griddata2<int> &);
           
//void initLoop(int, dataDetStruct &, g4m::incrementTab *species[8], ageStructVector &, 
//              ageStructVector &, datGlobal &, griddata &, griddata &, griddata &);

void initLoop(int, dataDetStruct &, g4m::incrementTab *species[8], 
	          g4m::ffipolm<double> *, 
              g4m::ffipolm<double> *, 
              g4m::ffipolm<bool> *, 
              g4m::ffipolm<bool> *, 	          
     	      ageStructVector &, 
              ageStructVector &, datGlobal &, griddata &, griddata &, griddata &);

//void initLoop(int, dataDetStruct &, g4m::incrementTab &, 
//              datGlobal &, griddata &, griddata &, griddata &);

//void calc(g4m::dataStruct &, g4m::incrementTab &, g4m::ageStruct &, g4m::ageStruct &, 
//          dat &, griddata2<char> &, griddata &, griddata &, int, int, int);
//void calc(g4m::dataStruct &, g4m::incrementTab &, g4m::ageStruct &, g4m::ageStruct &,
//          dat &, griddata2<char> &, griddata &, griddata &,
//          griddata &, griddata2<char> &, griddata &, 
//          griddata &, int , int , int );  
void calc(g4m::dataStruct &, g4m::incrementTab &, g4m::ageStruct &, g4m::ageStruct &,
          dat &, griddata2<char> &, griddata &, griddata &,
          griddata &, griddata2<char> &, griddata &, 
          griddata &, int , double , int );           
          
//void initCohorts(dataDetStruct &, g4m::incrementTab &, ageStructVector &, 
//                 griddata &,griddata &, griddata &);           
void initCohorts(dataDetStruct &, g4m::incrementTab &, ageStructVector &, ageStructVector &,
                 griddata &,griddata &, griddata &);   
                  
void MAI_country(void);
void woodHarvestStatCountry(void);
void hurdle_aff_deff(void);
void cPrices(void);
void int2str(int i, char tmp[]);   // defined in misc.h

void woodProductionIndexes(void); // country code and year for wood production file

double forNPV(g4m::dataStruct &, double , int , int , double ,      // estimation of forestry NPV
              double , double , double, g4m::ipol<double,double> , double, double );
double forNPVfd(g4m::dataStruct &, double , int , int , double , 
                g4m::ipol<double,double> , double , double , double );
double forNPVfdc(g4m::dataStruct &, double , int , int , double , 
                g4m::ipol<double,double> , double , double , double,double);       
                        
double npv_calc(g4m::dataStruct &iter, g4m::ageStruct &cohortTmp, double maiV, int year, int rotation, string regprice, string regprice0,double OforestShare, short int used);

double npv_calc10(g4m::dataStruct &iter, g4m::ageStruct &cohortTmp, double maiV, int year, int rotation, string regprice, string regprice0,double OforestShare, short int used);

double npv_calc50(g4m::dataStruct &iter, g4m::ageStruct &cohortTmp, double maiV, int year, int rotation, string regprice, string regprice0,double OforestShare, short int used);

void forNPV_init(void); // initialisation of arrays for estimation of country average forestry NPV

void carbonPrice(void); // definition of carbon price time functions
//g4m::incrementCurves *species[8];
void defineSpecies(g4m::incrementTab &);

double cohortRes(double , pair<g4m::ageStruct::v, g4m::ageStruct::v> &);

float BEF(int , float , float);

//int extractData (dataDetStruct &, ageStructVector &, double);g4m::ageStruct &
int extractData (g4m::ageStruct &,dat &, double, double, double);