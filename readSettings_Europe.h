#ifndef READSETTINGS_H_
#define READSETTINGS_H_

#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

struct {
  string coeffPath, inputPath, outputPath;
  set<string> parametersTable,parametersTableReg, parametersMap;
  bool produceTabs, produceMaps, tabs[3], maps[3];
} settings;

void readSettings()
 {
  ifstream f;
  //"..\..\..\..\g4mUAwoodtypes_v_1_0"+
  string FileName = "settings_Europe.ini";
  f.open(FileName.c_str(), ios::in);
  if (!f.is_open()) {
    cout << "Cannot read " << FileName << endl;
    system("pause");
    exit(0);
  }
  string line;
  int lineNum=0;
  while (!f.eof()) {
    getline(f,line);
    if(line.size()>0 && line[0] != '#') {
      stringstream ss(line);
      switch (lineNum) {
        case 0: settings.coeffPath = line;
                break;
        case 1: settings.inputPath = line;
                break;
        case 2: settings.outputPath = line;
                break;
        case 3: {int num;
                ss >> num;
cout<<"numMaps=\t"<<num<<endl;                
                for (int i = 0; i < num; i++) {
                  getline(f,line);
                  if(line.size()>0 && line[0] != '#' ) {
                    transform(line.begin(), line.end(), line.begin(), ::tolower);
                    settings.parametersMap.insert(line);
                  } else {
                    i--;
                  }
                }
               }
               break; 
        case 4: {int num;
                ss >> num;
cout<<"numTabs=\t"<<num<<endl;
                for (int i = 0; i < num; i++) {
                  getline(f,line);
                  if(line.size()>0 && line[0] != '#') {
                    transform(line.begin(), line.end(), line.begin(), ::tolower);
                    settings.parametersTable.insert(line);
                  } else {
                    i--;
                  }
                }
               }                
                break;
        case 5: {int num;
                ss >> num;
cout<<"numTabsReg=\t"<<num<<endl;
       
                for (int i = 0; i < num; i++) {
                  getline(f,line);
                  if(line.size()>0 && line[0] != '#') {
                    transform(line.begin(), line.end(), line.begin(), ::tolower);
                    settings.parametersTableReg.insert(line);
                  } else {
                    i--;
                  }
                }

              }                
                break;                
        case 6: settings.produceTabs = (line == "1");
cout<<"produceTabs=\t"<< settings.produceTabs<<endl;       
                break;
        case 7: settings.tabs[0] = (line == "1");
cout<<"Tabs[0]=\t"<< settings.tabs[0]<<endl;         
                break;
        case 8: settings.tabs[1] = (line == "1");
cout<<"Tabs[1]=\t"<< settings.tabs[1]<<endl;         
                break;
        case 9: settings.tabs[2] = (line == "1");
cout<<"Tabs[2]=\t"<< settings.tabs[2]<<endl;         
                break;
        case 10: settings.produceMaps = (line == "1");
cout<<"produceMaps=\t"<< settings.produceMaps<<endl;       
                break;
        case 11: settings.maps[0] = (line == "1");
cout<<"Maps[0]=\t"<< settings.maps[0]<<endl;           
                break;
        case 12: settings.maps[1] = (line == "1");
cout<<"Maps[1]=\t"<< settings.maps[1]<<endl;         
                break;
        case 13: settings.maps[2] = (line == "1");
cout<<"Maps[2]=\t"<< settings.maps[2]<<endl;         
                break;
      }
      lineNum++;
    }
  }
  f.close();
 }

#endif
