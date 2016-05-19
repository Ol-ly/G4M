#ifndef READINPUT_CPP
#define READINPUT_CPP
#include <algorithm>
#include <map>
#include "misc.h" 

using namespace std;

pair<double, string> getNumber(string str);

vector< pair<double, string> > header;

//******************************************************************************
// reads detailed data for plots from input file
//******************************************************************************
int readInputDet(dataDetStruct &data_all)
 {
//Opening a file  
  ifstream fp;

  string dataInputDir = settings.inputPath + "\\";
  string dataInputDir2 = settings.inputPath + "\\";

	  string FileName = dataInputDir +"b2_cramer_europe_euclimit_bas_ua_ONLY.txt";

  fp.open(FileName.c_str(), ios::in);
  if (!fp.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  
cout << "> Reading the rest of input data..." << FileName << endl;
  string line;
  getline(fp,line);
  stringstream ss(line);
  string buf;
  ss >> buf; ss >> buf;                      // skip x,y
  while (ss >> buf) {
    header.push_back(getNumber(buf));
    transform(header.back().second.begin(),header.back().second.end(),
    header.back().second.begin(), (int(*)(int)) toupper);
  }
  int LineNum = 0;   
  while (!fp.eof()) {
    getline(fp,line);
    if(line.size()>0 && line[0] != '#') {
      g4m::dataStruct d;
      stringstream ss(line);
      double val,x,y;
      ss >> x;
      ss >> y;
      d.x = int((x-GridStepLon/2+180)/GridStepLon);
      d.y = int((y-GridStepLat/2+90)/GridStepLat);
      int ElemNum = 0;
      while (ss >> val) {

         {
          if (!strcmp(header[ElemNum].second.c_str(),"COUNTRY"))
            d.COUNTRY.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"POTVEG"))
            d.POTVEG.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"PROTECT"))
            d.PROTECT.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"USED"))
            d.USED.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"LANDAREA"))
            d.LANDAREA.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"NPP"))
            d.NPP.insert(header[ElemNum].first,val);

          else if (!strcmp(header[ElemNum].second.c_str(),"POPDENS"))
            d.POPDENS.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"SPOPDENS")) //MG: Added
            d.SPOPDENS.insert(header[ElemNum].first,val);        
            
          else if (!strcmp(header[ElemNum].second.c_str(),"SAGRSUIT"))
            d.SAGRSUIT.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"AGRSUIT"))
            d.AGRSUIT.insert(header[ElemNum].first,val);

          else if (!strcmp(header[ElemNum].second.c_str(),"PRICEINDEX"))
            d.PRICEINDEX.insert(header[ElemNum].first,val);

          else if (!strcmp(header[ElemNum].second.c_str(),"BIOMASS"))
            d.BIOMASS.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"BIOMASSBL"))
            d.BIOMASSBL.insert(header[ElemNum].first,val);            
          else if (!strcmp(header[ElemNum].second.c_str(),"CABOVEHA"))
            d.CABOVEHA.insert(header[ElemNum].first,val);    
          else if (!strcmp(header[ElemNum].second.c_str(),"CBELOWHA"))
            d.CBELOWHA.insert(header[ElemNum].first,val);    
          else if (!strcmp(header[ElemNum].second.c_str(),"CDEADHA"))
            d.CDEADHA.insert(header[ElemNum].first,val);                                               
          else if (!strcmp(header[ElemNum].second.c_str(),"CLITTERHA"))
            d.CLITTERHA.insert(header[ElemNum].first,val);     
          else if (!strcmp(header[ElemNum].second.c_str(),"SOCHA"))
            d.SOCHA.insert(header[ElemNum].first,val);  

          else if (!strcmp(header[ElemNum].second.c_str(),"DECHERB"))
            d.DECHERB.insert(header[ElemNum].first,val);                                               
          else if (!strcmp(header[ElemNum].second.c_str(),"DECWOOD"))
            d.DECWOOD.insert(header[ElemNum].first,val);     
          else if (!strcmp(header[ElemNum].second.c_str(),"DECSOC"))
            d.DECSOC.insert(header[ElemNum].first,val);              
                          
          else if (!strcmp(header[ElemNum].second.c_str(),"FOREST"))
            d.FOREST.insert(header[ElemNum].first,val);

          else if (!strcmp(header[ElemNum].second.c_str(),"SPECIESTYPE"))
            d.SPECIESTYPE.insert(header[ElemNum].first,val);

          else if (!strcmp(header[ElemNum].second.c_str(),"R"))
            d.R.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"GDP"))
            d.GDP.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"BUILTUP"))
            d.BUILTUP.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"CROP"))
            d.CROP.insert(header[ElemNum].first,val); 

          else if (!strcmp(header[ElemNum].second.c_str(),"FRACLONGPROD"))
            d.FRACLONGPROD.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"CORRUPTION"))
            d.CORRUPTION.insert(header[ElemNum].first,val);
          else if (!strcmp(header[ElemNum].second.c_str(),"SLASHBURN"))
            d.SLASHBURN.insert(header[ElemNum].first,val);

          else if (!strcmp(header[ElemNum].second.c_str(),"IIASA_REGION"))
            d.IIASA_REGION.insert(header[ElemNum].first,val);       

          else if (!strcmp(header[ElemNum].second.c_str(),"FTIMBER"))
            d.FTIMBER.insert(header[ElemNum].first,val);                   

          else if (!strcmp(header[ElemNum].second.c_str(),"POLESREG"))
            d.POLESREG.insert(header[ElemNum].first,val);                        

          else if (!strcmp(header[ElemNum].second.c_str(),"MAIE"))
            d.MAIE.insert(header[ElemNum].first,val);        

          else if (!strcmp(header[ElemNum].second.c_str(),"MAIN"))
            d.MAIN.insert(header[ElemNum].first,val);                    

          else if (!strcmp(header[ElemNum].second.c_str(),"MANAGEDSHARE"))
            d.MANAGEDSHARE.insert(header[ElemNum].first,val);                    

          else if (!strcmp(header[ElemNum].second.c_str(),"MANAGEDFLAG"))
            d.MANAGEDFLAG.insert(header[ElemNum].first,val);                                

          else if (!strcmp(header[ElemNum].second.c_str(),"COUNTRYREGMIX"))
            d.COUNTRYREGMIX.insert(header[ElemNum].first,val);       

          else if (!strcmp(header[ElemNum].second.c_str(),"ROAD"))
            d.ROAD.insert(header[ElemNum].first,val);      

          else if (!strcmp(header[ElemNum].second.c_str(),"FORLOSS"))
            d.FORLOSS.insert(header[ElemNum].first,val);      

          else if (!strcmp(header[ElemNum].second.c_str(),"GLOBIOM_RESERVED"))
            d.GLOBIOM_RESERVED.insert(header[ElemNum].first,val);      

          else if (!strcmp(header[ElemNum].second.c_str(),"SIMUID"))
            d.SIMUID.insert(header[ElemNum].first,val);      

          else if (!strcmp(header[ElemNum].second.c_str(),"AFFORMAX"))
            d.AFFORMAX.insert(header[ElemNum].first,val);   
          // ...
         }
        ElemNum++;
      }
      data_all.push_back(d);
      LineNum++;
    }
  }
  fp.close();


cout << "Successfully read " << LineNum << " lines." << endl;

int ret = LineNum;

{
  ifstream fplp;

    FileName = dataInputDir2+"landPrice_bas_ua.txt"; 

cout<<"Reading file: "<<FileName<<endl;
  fplp.open(FileName.c_str(), ios::in);
  if (!fplp.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  getline(fplp,line); // we don't need 1st line
  int lineNum=0;
  while (!fplp.eof()) {
    getline(fplp,line);
    if(line.size()>0 && line[0] != '#') {
      stringstream ss(line);
      string regID;
      double val;
      ss >> regID;
      int ElemNum = 0;
      string colname;
      char regstr[3]; 
      int2str(lineNum+1,regstr);
      colname = "re" + string(regstr) + "price0";
      while (ss >> val) {
        lprice[colname].insert(2000+(ElemNum)*10,val);
        ElemNum++;
      }
      lineNum++;

	}
  }
  fplp.close();  
}
{
  ifstream fpwp;

    FileName = dataInputDir2+"woodPrice_bas_ua.txt";  

cout<<"Reading file: "<<FileName<<endl;
  fpwp.open(FileName.c_str(), ios::in);
  if (!fpwp.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  getline(fpwp,line);  // we don't need 1st line
  int LineNum=0;
  while (!fpwp.eof()) {
    getline(fpwp,line);
    if(line.size()>0 && line[0] != '#') {
      stringstream ss(line);
      string regID;
      double val;
      ss >> regID;
      int ElemNum = 0;
      string colname;
      char regstr[3]; 
      int2str(LineNum+1,regstr);
      colname = "re" + string(regstr) + "price0";
      while (ss >> val) {
        wprice[colname].insert(2000+(ElemNum)*10,val);
        ElemNum++;
      }
      LineNum++;
    }
  }

  fpwp.close();
}
{
  ifstream fppr;
	FileName = dataInputDir2+"SawnWoodProduction_UA_EU.txt";
	
cout<<"Reading file: "<<FileName<<endl;
  fppr.open(FileName.c_str(), ios::in);
  if (!fppr.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  getline(fppr,line);  // we don't need 1st line
  int LineNum=0;
  while (!fppr.eof()) {
    getline(fppr,line);
    if(line.size()>0 && line[0] != '#') 
    {
      stringstream ss(line);
      string regID;
      double val;
      ss >> regID;
      int ElemNum = 0;
      string colname;
      char regstr[4]; 
      int2str(countryNwp[LineNum],regstr);
      colname = "re" + string(regstr) + "price0";
      while (ss >> val) {
        wprod_sawnwood[colname].insert(yearNwp[ElemNum],val);  
//		cout << "colname	" << colname << endl;
//		cout << "yearNwp[ElemNum]	" << yearNwp[ElemNum] << endl;
//		cout << "wpor_sawnwood	" << val << endl;
        ElemNum++;
      }
      LineNum++;
    }
  }
  fppr.close();
}

{
  ifstream fppr;
	FileName = dataInputDir2+"RestWoodProduction_UA_EU.txt";
	
cout<<"Reading file: "<<FileName<<endl;
  fppr.open(FileName.c_str(), ios::in);
  if (!fppr.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  getline(fppr,line);  // we don't need 1st line
  int LineNum=0;
  while (!fppr.eof()) {
    getline(fppr,line);
    if(line.size()>0 && line[0] != '#') 
    {
      stringstream ss(line);
      string regID;
      double val;
      ss >> regID;
      int ElemNum = 0;
      string colname;
      char regstr[4]; 
      int2str(countryNwp[LineNum],regstr);
      colname = "re" + string(regstr) + "price0";
      while (ss >> val) {
        wprod_restwood[colname].insert(yearNwp[ElemNum],val);  
        ElemNum++;
      }
      LineNum++;
    }
  }
  fppr.close();
}

  return(ret);
// finished reading file
 }



//******************************************************************************
// reading coefficients
//******************************************************************************
void readCoeff(g4m::coeffStruct &coeff)
 {
//Opening a file  

  ifstream fp;

  string FileName = settings.coeffPath;
  fp.open(FileName.c_str(), ios::in);
  if (!fp.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  cout << "> Reading coefficients..." << endl;
  while (!fp.eof())
    {
     string line;
     getline(fp, line);
     if(line.size()>0 && line[0] != '#')                        //Jump over lines starting with #
       {
        transform(line.begin(),line.end(),line.begin(), (int(*)(int)) toupper);
        stringstream ss(line);
        string buf;
        double val;
        ss >> buf;
        pair<double,string> tmp=getNumber(buf);
        ss >> val;
        if (!strcmp(tmp.second.c_str(),"BYEAR"))
          coeff.bYear = int(val);
        else if (!strcmp(tmp.second.c_str(),"EYEAR"))
          coeff.eYear = int(val);
        else if (!strcmp(tmp.second.c_str(),"CELLSINTERACT"))
          coeff.cellsInteract = int(val);
        else if (!strcmp(tmp.second.c_str(),"INCLAFFOR"))
          coeff.inclAffor = int(val);
        else if (!strcmp(tmp.second.c_str(),"NOPAY"))
          coeff.noPay = int(val);
        else if (!strcmp(tmp.second.c_str(),"UBIOMASS"))
          coeff.uBiomass = int(val);
        else if (!strcmp(tmp.second.c_str(),"LITTER"))
          coeff.litter = int(val);
        else if (!strcmp(tmp.second.c_str(),"SOC"))
          coeff.SOC = int(val);
        else if (!strcmp(tmp.second.c_str(),"PRICELANDMINR"))
          coeff.PriceLandMinR.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PRICELANDMAXR"))
          coeff.PriceLandMaxR.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"FCUPTAKE"))
          coeff.FCuptake.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"HARVLOOS"))
          coeff.HarvLoos.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PRICEC"))
          coeff.PriceC.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"FRACLONGPROD"))
          coeff.FracLongProd.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"DECRATEL"))
          coeff.decRateL.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"DECRATES"))
          coeff.decRateS.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"SLASHBURN"))
          coeff.SlashBurn.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"FREQAID"))
          coeff.FreqAid.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PRICECAID"))
          coeff.PriceCAid.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"MAXROTINTER"))
          coeff.MaxRotInter.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"MINROTINTER"))
          coeff.MinRotInter.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"BASELINE"))
          coeff.baseline.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PRICETIMBERMAXR"))
          coeff.PriceTimberMaxR.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PRICETIMBERMINR"))
          coeff.PriceTimberMinR.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PRICEINDEXE"))
          coeff.PriceIndexE.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"PLANTINGCOSTSR"))
          coeff.PlantingCostsR.insert(tmp.first,val);
        else if (!strcmp(tmp.second.c_str(),"SPOPDENS"))
          coeff.sPopDens.insert(tmp.first,val);
       }
    }
  fp.close();

 }

//******************************************************************************
// splits string to pair <number,string>
//******************************************************************************
pair<double, string> getNumber(string str)
{ // Returns a pair of double, string;
  string text;
  string number;
  for(unsigned int i=0; i < str.size(); ++i) {
    if((str[i] >= '0' && str[i] <= '9') || str[i] == '.')
       number += str[i];
    else if((str[i] != '[') && (str[i] != ']'))
      text += str[i];
  }
  double x = 0;
  stringstream ss(number);
  ss >> x;
  return(pair<double, string>(x, text));
}

#endif
