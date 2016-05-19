#include "MDT.h"

// Default constructor
MDT::MDT()
 {
  nDims = 0;
 }

// Constructor
MDT::MDT(string fileName)
 {
  ifstream f;
  f.open(fileName.c_str(), ios::in);
  if (f.is_open()) {
    string line;
    // Reading number of dimensions and parameter description:
    getline(f,line);
    stringstream s(line);
    s >> nDims;
    dimCardinals.resize(nDims);
    dimNames.resize(nDims);
    dimElements.resize(nDims);
    while (s >> line) paramName += line + " ";
    paramName.erase(paramName.length()-1);
    // Reading dimensions' descriptions:
    for (int i = 0; i < nDims; i++) {
      getline(f, line);
      stringstream ss(line);
      int tmp;
      ss >> tmp;
      dimCardinals[i] = tmp;
      while (ss >> line) dimNames[i] += line + " ";
      dimNames[i].erase(dimNames[i].length()-1);
      for (int j = 0; j < tmp; j++) {
        getline(f, line);
        dimElements[i].push_back(line);
      }
    }
//cout << paramName << " has " << nDims << " dimensions." << endl;
//for (int i = 0; i < nDims; i++) {
//  cout << dimNames[i] << " has " << dimCardinals[i] << " elements:" << endl;
//  for (int j = 0; j < dimCardinals[i]; j++) {
//    cout << "\t" << dimElements[i][j] << endl;
//  }
//}
    f.close();
    cout << "Successfully read file: " << fileName << endl;
  } else {
    cout << "Unable to open file! " << fileName << endl;
  }
 }

// Destructor
MDT::~MDT()
 {

 }

//================
// Class methods:
//================

void MDT::clear()
 {
  paramName = "";
  nDims = 0;
  dimCardinals.clear();
  dimNames.clear();
  dimElements.clear();
 }

void MDT::rename(string name)
 {
  paramName = name;
 }

bool MDT::addDim(string dimName, strVector elements)
 {
  for (int i = 0; i < nDims; i++)
    if (dimNames[i] == dimName)      // verifying if a dimension with the specified name exists;
      return false;                  // if exists, return false;
  nDims++;
  dimNames.push_back(dimName);
  dimCardinals.push_back(elements.size());
  dimElements.push_back(elements);
  return true;
 }

bool MDT::addDimEl(string dimName, string element)
 {
  bool dimExists = false;
  for (int i = 0; i < nDims; i++) {
    if (dimNames[i] == dimName) {            // verifying if a dimension with the specified name exists;
      dimExists = true;
      bool elementExists = false;
      for (int j = 0; j < dimElements[i].size(); j++) {
        if (dimElements[i][j] == element) elementExists = true;
      }                                      // verifying if an element with the specified name exists within dimension;
      if (elementExists) return false;       // if exists, don't add new element and return false;
      else {                                 // else - add new element
        dimCardinals[i]++;
        dimElements[i].push_back(element);
        return true;
      }
    }
  }
  if (!dimExists) {        // if a dimension with the specified name doesn't exist, add new dimension
    strVector tmp;
    tmp.push_back(element);
    addDim(dimName, tmp);
    return true;
  }
  return false;
 }

long long MDT::getN()
 {
  long long res = 1;
  for (int i = 0; i < nDims; i++) {
    res *= dimCardinals[i];
  }
  return res;
 }
 
int MDT::getNDims()
 {
  return nDims;
 }

long long MDT::getHash(strVector elements)
 {
  if (nDims < elements.size()) {
    for (int i = nDims; i < elements.size(); i++) {
      addDimEl("Dimension " + IntToStr(i + 1), elements[i]);
    }
  }
  int N = 1;
  int tmp = getCoordinate(0, elements[0]);
  if (tmp == -1) return -1;
  long long res = 0;
  res += N * tmp;
  for (int i = 1; i < nDims; i++) {
    N *= dimCardinals[i-1];
    tmp = getCoordinate(i, elements[i]);
//cout << i << "\t" << dimCardinals[i-1] << "\t" << elements[i] << "\t" << tmp << "\t" << N << endl;
    if (tmp == -1) return -1;
    res += N * tmp;
  }
//system("pause");
  return res;
 }

int MDT::getCoordinate(int dim, string element)
 {
  for (int i = 0; i < dimElements[dim].size(); i++) {
    if (dimElements[dim][i] == element) return i;
  }
  dimElements[dim].push_back(element);
  dimCardinals[dim]++;
  return (dimElements[dim].size() - 1);
 }

void MDT::SaveToFile(string fileName, string mode)
 {
  ofstream f;
  if (mode == "MAP") fileName += ".mdt";
  f.open(fileName.c_str(), ios::out);
  if (f.is_open()) {
    if (mode == "MAP") {
        f << nDims << " " << paramName << endl;
        cout<<"dimCardinalsMAP.size=\t"<<dimCardinals.size()<<endl;
        for (int i = nDims - 1; i >= 0; i--) {
          f << dimCardinals[i] << " " << dimNames[i] << endl;
          for (int j = 0; j < dimCardinals[i]; j++) {
            f << dimElements[i][j] << endl;
          }
        }
        cout << "Successfully written to file: " << fileName << endl;
    } else if (mode == "TABLE") {
        f << paramName << endl;
        cout<<"dimCardinalsTable.size=\t"<<dimCardinals.size()<<endl;        
        for (int i = 0; i < nDims; i++) {
          f << dimNames[i] << endl;
          for (int j = 0; j < dimCardinals[i]; j++) {
            f << dimElements[i][j] << endl;
          }
        }
        cout << "Successfully written to file: " << fileName << endl;
    } else {
        cout << "Error writing to file! Wrong mode! " << fileName << endl;
    }
    f.close();
    
  } else {
    cout << "Unable to save to file: " << fileName << endl;
  }
 }
