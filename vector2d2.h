#ifndef VECTOR2d
#define VECTOR2d
//  Name:           
//  Author:         Andriy Bun
//  Date:           01.04.10
//  Description:    

#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

class vector2d {
private:
  int dimensionsNum;
  vector<double *> dat_all;
  double *tmpVector;
public:
  vector2d();
  vector2d(int dims);
  ~vector2d();
  void push(vector<double> input);
  vector<double> pop();
  bool saveToFile(string fileName);
  bool readFromFile(string fileName);
  vector<double> operator [](int index);
  int size();
};

// Default constructor
vector2d::vector2d()
 {
  dimensionsNum = 1;
 }

// Class constructor
vector2d::vector2d(int dims)
 {
  dimensionsNum = dims;
 }

// Destructor
vector2d::~vector2d()
 {
  for(int i = dat_all.size() - 1; i >= 0; i--) {
    delete dat_all[i];
    dat_all.pop_back();
  }
 }

// Other methods:

void vector2d::push(vector<double> input)
 {
  tmpVector = new double[dimensionsNum];
  for (int i = 0; i < dimensionsNum; i++) {
    if (i >= input.size())
      tmpVector[i] = 0;
    else
      tmpVector[i] = input[i];  
  }
  dat_all.push_back(tmpVector);
 }

vector<double> vector2d::pop()
 {
  vector<double> res;
  int sz = dat_all.size() - 1;
  if (sz >= 0) {
    for (int j = 0; j < dimensionsNum; j++) {
      res.push_back(dat_all[sz][j]);
    }
    delete dat_all[sz];
    dat_all.pop_back();
  }
  return res;
 }

bool vector2d::saveToFile(string fileName)
 {
  if (dat_all.size() == 0) {
    cout << "ERROR! Attempt to write empty data." << endl;
    return false;
  }
  string fileNameTmp = fileName + ".bin";
  ofstream f;
  f.open(fileNameTmp.c_str(), ios::out | ios::binary);
  if (f.is_open()) {
    int tmp = dat_all.size();
    f.write(reinterpret_cast<char *>(&tmp), sizeof(int));
    f.write(reinterpret_cast<char *>(&dimensionsNum), sizeof(int));
    for (int i = 0; i < dat_all.size(); i++) {
      f.write(reinterpret_cast<char *>(dat_all[i]), sizeof(double)*dimensionsNum);
    }
  }
  f.close();
  return true;
 }

bool vector2d::readFromFile(string fileName)
 {
  for(int i = dat_all.size() - 1; i >= 0; i--) {
    delete dat_all[i];
    dat_all.pop_back();
  }
  string fileNameTmp = fileName + ".bin";
  ifstream f;
  f.open(fileNameTmp.c_str(), ios::in | ios::binary);
  if (f.is_open()) {
    int size;
    f.read(reinterpret_cast<char *>(&size), sizeof(int));
    f.read(reinterpret_cast<char *>(&dimensionsNum), sizeof(int));
    for (int i = 0; i < size; i++) {
      tmpVector = new double[dimensionsNum];
      f.read(reinterpret_cast<char *>(tmpVector), sizeof(double)*dimensionsNum);
      dat_all.push_back(tmpVector);
    }
    f.close();
  }
  else {
    cout << "ERROR! Unable to open the file: " << fileName << endl;
  }
//  for (int i = 0; i < dat_all.size(); i++) {
//    cout << i << " |\t";
//    for (int j = 0; j < dimensionsNum; j++) {
//      cout << dat_all[i][j] << "\t";
//    }
//    cout << endl;
//  }
  return true;
 }

vector<double> vector2d::operator [](int index)
 {
  vector<double> res;
  if (index < dat_all.size()) {
    for (int j = 0; j < dimensionsNum; j++) {
      res.push_back(dat_all[index][j]);
    }
  }
  return res;
 }

int vector2d::size()
 {
  return dat_all.size();
 }

#endif