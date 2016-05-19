#include "tableData.h"

// default constructor
tableData::tableData()
 {
  descr = MDT();
  N = descr.getN();
 }

// constructor
tableData::tableData(string fileName_GDT)
 {
  descr = MDT(fileName_GDT);
  N = descr.getN();
 }

// destructor
tableData::~tableData()
 {
/*  map<int, float_vector_t>::iterator it = data.begin();
  while (it != data.end()) {
    delete it->second;
    it++;
  }*/
 }

//================
// Class methods:
//================

void tableData::insert(float val, strVector point)
 {
  if (val != 0)
    data.insert(make_pair(descr.getHash(point), val));
 }

void tableData::update(float val, strVector point)
 {
  if (val != 0) {
    long long card = descr.getHash(point);
    data.erase(card);
    data.insert(make_pair(card, val));
  }
 }

void tableData::rename(string name)
 {
  descr.paramName = name;
 }

void tableData::renameDims(strVector vec)
 {
  descr.dimNames = vec;
 }

void tableData::addDim(string dimName, set<string> elements)
 {
  strVector tmp;
  set<string>::iterator it = elements.begin();
  while(it != elements.end()) {
    tmp.push_back(*it);
    it++;
  }
  descr.addDim(dimName,tmp);
 }

void tableData::addDim(string dimName, set<int> elements)
 {
  strVector tmp;
  set<int>::iterator it = elements.begin();
  while(it != elements.end()) {
    tmp.push_back(IntToStr(*it));
    it++;
  }
  descr.addDim(dimName,tmp);
 }

void tableData::addDim(string dimName, string element)
 {
  strVector tmp;
  tmp.push_back(element);
  descr.addDim(dimName,tmp);
 }

void tableData::clear()
 {
  data.clear();
  descr.clear();
 }

bool tableData::SaveToFile(string outDir, string fileName)
 {
//****************************
//**** Writimg *.GDC file ****
//****************************
  if (data.size() == 0) {
    cout << "ERROR! Attempt to write empty data. " << fileName << endl;
    return false;
  }
  fileName = outDir + pathSeparator + fileName;
  string fileNameTmp = fileName + ".gdc";
  ofstream f;
  f.open(fileNameTmp.c_str(), ios::out | ios::binary);
  if (f.is_open()) {
    // Writing info about dimensions:
    int tmp_int = descr.getNDims();
    INV_BYTE_ORDER(tmp_int);
    f.write(reinterpret_cast<char *>(&tmp_int), sizeof(int));
    for (int i = 0; i < descr.nDims; i++) {
      tmp_int = descr.dimCardinals[i];
      INV_BYTE_ORDER(tmp_int);
      f.write(reinterpret_cast<char *>(&tmp_int), sizeof(int));
    }
    // Writing hash data:
    tmp_int = data.size();
    int numVal = tmp_int;
    long long * vl = new long long [numVal];
    float * vf = new float [numVal];
    INV_BYTE_ORDER(tmp_int);
    f.write(reinterpret_cast<char *>(&tmp_int), sizeof(int));
    int i = 0;
    map<long long, float>::iterator it = data.begin();
    while (it != data.end()) {
      long long tmp_long = it->first;
      INV_BYTE_ORDER(tmp_long);
      vl[i] = tmp_long;
      float tmp_float = it->second;
      INV_BYTE_ORDER(tmp_float);
      vf[i] = tmp_float;
      i++;
      it++;
    }
    f.write(reinterpret_cast<char *>(vl), sizeof(long long) * numVal);
    f.write(reinterpret_cast<char *>(vf), sizeof(float) * numVal);
    delete []vl;
    delete []vf;
    f.close();
    cout << "Successfully written to binary file: " << fileNameTmp << endl;
  } else {
    cout << "Unable to save to file!  " << fileNameTmp << endl;
    return false;
  }
//****************************
//**** Writimg *.GDT file ****
//****************************
  descr.SaveToFile(fileName + ".gdt", "TABLE");
  return true;
 }