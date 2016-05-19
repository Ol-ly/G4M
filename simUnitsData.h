//  Name:         simUnitsData class
//  Author:       Andriy Bun
//  Date:         05.11.09
//  Modified:     11.01.10 (Andriy Bun)
//  Description:  Class for storing and producing data for viewing with GDataView
//                Software. Is a friend class for MDT class.

#ifndef SIMUNITSDATA_H_
#define SIMUNITSDATA_H_

#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <cstring>

//#include "endianness.h"
#include "endianness.cpp"
//#include "MDT.h"
#include "MDT.cpp"
#ifndef PATHSEPARATOR_
#define PATHSEPARATOR_

#ifdef unix
string pathSeparator = "/";
#else
string pathSeparator = "\\";
#endif

#endif

using namespace std;

class simUnitsData {
private:
  typedef vector<string> strVector;
  typedef vector<float> floatVector;
  map<int, floatVector> data;
  MDT descr;
  int N;                           // Number of data records per simulation unit
public:
  simUnitsData();
  simUnitsData(string fileName_MDT);
  ~simUnitsData();
  // Inserts a value "val" corresponding to an active simulation unit SIMU and
  // vector of coordinates "point" into the list.
  void insert(int SIMU, float val, strVector point);
  // Rename dataset:
  void rename(string name);
  // Rename dimensions of the dataset:
  void renameDims(strVector vec);
  // Add new dimension:
  void addDim(string dimName, set<string> elements);
  void addDim(string dimName, set<int> elements);
  void addDim(string dimName, string element);
  // Clearing dataset:
  void clear();
  // Writing dataset to files *.MSU and *.MDC files
  bool SaveToFile(string outDir, string fileName);
  
  int dimN();  
};

// default constructor
simUnitsData::simUnitsData()
 {
  descr = MDT();
  N = descr.getN();
 }

// constructor
simUnitsData::simUnitsData(string fileName_MDT)
 {
  descr = MDT(fileName_MDT);
  N = descr.getN();
 }

// destructor
simUnitsData::~simUnitsData()
 {
/*  map<int, floatVector>::iterator it = data.begin();
  while (it != data.end()) {
    delete it->second;
    it++;
  }*/
 }

//================
// Class methods:
//================

void simUnitsData::insert(int SIMU, float val, strVector point)
 {
  if (data.find(SIMU) == data.end()) {
    floatVector tmp;// = new floatVector;
    data.insert(make_pair(SIMU, tmp));
  }
  int card = descr.getHash(point);
  if (data[SIMU].size() < (card + 1)) data[SIMU].resize(card + 1);
  data[SIMU][card] = val;
 }

void simUnitsData::rename(string name)
 {
  descr.paramName = name;
 }

void simUnitsData::renameDims(strVector vec)
 {
  descr.dimNames = vec;
 }

void simUnitsData::addDim(string dimName, set<string> elements)
 {
  strVector tmp;
  set<string>::iterator it = elements.begin();
  while(it != elements.end()) {
    tmp.push_back(*it);
    it++;
  }
  descr.addDim(dimName,tmp);
 }

void simUnitsData::addDim(string dimName, set<int> elements)
 {
  strVector tmp;
  set<int>::iterator it = elements.begin();
  while(it != elements.end()) {
    tmp.push_back(IntToStr(*it));
    it++;
  }
  descr.addDim(dimName,tmp);
 }

void simUnitsData::addDim(string dimName, string element)
 {
  strVector tmp;
  tmp.push_back(element);
  descr.addDim(dimName,tmp);
 }

void simUnitsData::clear()
 {
  map<int, floatVector>::iterator it = data.begin();
  while (it != data.end()) {
    it->second.clear();
    it++;
  }
  data.clear();
  descr.clear();
 }

bool simUnitsData::SaveToFile(string outDir, string fileName)
 {
//****************************
//**** Writimg *.MSU file ****
//****************************
cout << "filename: " + outDir + pathSeparator + "main" << endl;
//descr.SaveToFile(outDir + pathSeparator + "main", "MAP");
  map<int, floatVector>::iterator it = data.begin();
  if (it == data.end()) {
    cout << "ERROR! Attempt to write empty data - msu - file." << endl;
    return false;
  }
  ofstream f;
  fileName = outDir + pathSeparator + fileName;
  string fileNameTmp = fileName + ".msu";
  f.open(fileNameTmp.c_str(), ios::out | ios::binary);
  if (f.is_open()) {
    while (it != data.end()) {
      int tmp = it->first;
      if (tmp < 0) {
        data.erase(tmp);
        it = data.begin();
      }
      else
        it++;
    }
    int numSimU = data.size();
    INV_BYTE_ORDER(numSimU);
    f.write(reinterpret_cast<char *>(&numSimU), sizeof(int));
// writing simulation units to file
    it = data.begin();
    while (it != data.end()) {
      int tmp = it->first;
        INV_BYTE_ORDER(tmp);
        f.write(reinterpret_cast<char *>(&tmp), sizeof(int));
      it++;
    }
    f.close();
    cout << "Successfully written to binary file: " << fileNameTmp << endl;
  } else {
    cout << "Unable to save to file!  "<< fileNameTmp  << endl;
    return false;
  }
//****************************
//**** Writimg *.MDC file ****
//****************************
  N = descr.getN();
  if (N==0) {cout << "ERROR! Attempt to write empty data (mdc - file)." << endl;}
  fileNameTmp = fileName + ".mdc";
  f.open(fileNameTmp.c_str(), ios::out | ios::binary);
  if (f.is_open()) {
    it = data.begin();
    while (it != data.end()) {
      if (it->first >= 0)
		if (N>it->second.size()) N=it->second.size(); //MG
        for (int i = 0; i < N; i++) {
          INV_BYTE_ORDER(it->second[i]);
          f.write(reinterpret_cast<char *>(&it->second[i]), sizeof(float));
        }
      it++;
    }
    f.close();
    cout << "Successfully written to binary file: " << fileNameTmp << endl;
  } else {
    cout << "Unable to save to file!  "  << fileNameTmp << endl;
    return false;
  }
  descr.SaveToFile(outDir + pathSeparator + "main", "MAP");
  return true;
 }
 
 int simUnitsData::dimN()
 {int dimn=descr.getN();
  return dimn;
 }
#endif
