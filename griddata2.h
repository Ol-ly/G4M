//   Name:          GridData class 2 + any type
//   Author:        Andriy Bun; Mykola Gusti
//   Date:          26.08.2009

#ifndef griddata2_h_
#define griddata2_h_

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <typeinfo>


using namespace std;

template <class TP>
class griddata2
 {
  private:
    int HorResolution;
    int VerResolution;
    int HorNeigh;
    int VerNeigh;
//    TP GridRows[];
    TP *GridRows;
    TP *grid, *gridPrev;  
  public:
      griddata2(int HR, int VR, TP val);
      griddata2();
      griddata2(const griddata2& g);
      ~griddata2();
      void ShowArray();                       // prints array
      void PrintToFile(string fileName, string rastrType);      // print array to file
      void PrintToFilePrev(string fileName, string rastrType);  // print previous year array to file
      void ShowArrayPrev();                   // prints array for the previous year
      void update();                          // updates values for previous year to current values
      void set(int x, int y, TP val);     // assigns value val to cell [x][y]
      void setPrev(int x, int y, TP val); // assigns previous year value val to cell [x][y]
      void inc(int x, int y, TP val);     // adds value val to the existing value in cell [x][y]
      TP get(int x, int y);               // returns value stored in cell [x][y]
      TP getPrev(int x, int y);           // returns value for the previous year stored in cell [x][y]
      void SetNeighNum(int n, int m);     // sets number of neighbour cells to be considered
      TP GetMax(int x, int y);            // returns maximum value of all neighbours for the previous year
      TP GetMin(int x, int y);            // returns minimum value of all neighbours for the previous year
      TP GetAvg(int x, int y);            // returns average value for the previous year
 };

  // using class constructor:
  // griddata2<DATA_TYPE> OBJECT_NAME = griddata2<DATA_TYPE>(LONGITUDE_RES,LATTITUDE_RES,DEFAULT_VALUE);
  //         EXAMPLE: griddata2<char> OBJ = griddata2<char>(36,18,-99.5);
  // using default class constructor:
  // griddata2<DATA_TYPE> OBJECT_NAME;
  // default parameters are:
  //         LONGITUDE_RES = 720
  //         LATTITUDE_RES = 360
  //         DEFAULT_VALUE = 0
  //         EXAMPLE: griddata2<double> OBJ_DEF;


template <class TP>
void griddata2<TP>::ShowArray()
 {
  TP val;
  int isChar = 0;
  char t1;
  unsigned char t2;
  if ((typeid(val).name() == typeid(t1).name()) ||
      (typeid(val).name() == typeid(t2).name())) isChar = 1;
  for (int j = 0; j < VerResolution; j++)
   {
    cout << j << "|\t";
    for (int i = 0; i < HorResolution; i++)
     {
      if (isChar)
        cout << (int)grid[j*HorResolution+i] << "\t";
      else
        cout << grid[j*HorResolution+i] << "\t";
     }
    cout << endl;  
   }   
 }

template <class TP>
void griddata2<TP>::ShowArrayPrev()
 {
  TP val;
  int isChar = 0;
  char t1;
  unsigned char t2;
  if ((typeid(val).name() == typeid(t1).name()) ||
      (typeid(val).name() == typeid(t2).name())) isChar = 1;
  for (int j = 0; j < VerResolution; j++)
   {
    cout << j << "|\t";
    for (int i = 0; i < HorResolution; i++)
     {
      if (isChar)
        cout << (int)gridPrev[j*HorResolution+i] << "\t";
      else
        cout << gridPrev[j*HorResolution+i] << "\t";
     }
    cout << endl;  
   }   
 }

template <class TP>
void griddata2<TP>::PrintToFile(string fileName, string rastrType = "ESRI")
 {
  TP val;
  int isChar = 0;
  char t1;
  unsigned char t2;
  if ((typeid(val).name() == typeid(t1).name()) ||
      (typeid(val).name() == typeid(t2).name())) isChar = 1;
  ofstream f;
  fileName = /*"output\\" +*/ fileName + ".asc";
  f.open(fileName.c_str(),ios::out);
  if (f.is_open()) {
/////////////////////// Select grid type:
if (rastrType == "GRASS"){
    f << "cols: " << HorResolution << endl;
    f << "rows: " << VerResolution << endl;
    f << "west: -180" << endl;
    f << "south: -90" << endl;
    f << "north: 90" << endl;
    f << "east: 180" << endl;
}else if (rastrType == "ESRI"){
//----------------------------    
// ESRI ascii Grid
    f << "NCOLS " << HorResolution << endl;
    f << "NROWS " << VerResolution << endl;
    f << "XLLCORNER -180" << endl;
    f << "YLLCORNER -90" << endl;
    f << "CELLSIZE " << 360./HorResolution  << endl;
    f << "NODATA_VALUE -9999" << endl;
}else {cout<<"griddata23 error message: Specify correct rastrType: ESRI or GRASS"<<endl;}
//-----------------------------        
    for (int j = 0; j < VerResolution; j++)
     {
      for (int i = 0; i < HorResolution; i++)
       {
        if (isChar)
          f << (int)grid[(VerResolution-j-1)*HorResolution+i] << " ";
        else
          f << grid[(VerResolution-j-1)*HorResolution+i] << " ";
       }
      f << endl;  
     }
    f.close(); 
  } else {
    cout << "Unable to save to file!" << endl;
  }
 }

template <class TP>
void griddata2<TP>::PrintToFilePrev(string fileName, string rastrType = "ESRI")
 {
  TP val;
  int isChar = 0;
  char t1;
  unsigned char t2;
  if ((typeid(val).name() == typeid(t1).name()) ||
      (typeid(val).name() == typeid(t2).name())) isChar = 1;
  ofstream f;
  fileName = "output\\" + fileName + ".asc";
  f.open(fileName.c_str(),ios::out);
  if (f.is_open()) {
/////////////////////// Select grid type:
//GRASS ascii Grid
if (rastrType == "GRASS"){
    f << "cols: " << HorResolution << endl;
    f << "rows: " << VerResolution << endl;
    f << "west: -180" << endl;
    f << "south: -90" << endl;
    f << "north: 90" << endl;
    f << "east: 180" << endl;
}else if (rastrType == "ESRI"){
//----------------------------    
// ESRI ascii Grid
    f << "NCOLS " << HorResolution << endl;
    f << "NROWS " << VerResolution << endl;
    f << "XLLCORNER -180" << endl;
    f << "YLLCORNER -90" << endl;
    f << "CELLSIZE " << 360./HorResolution  << endl;
    f << "NODATA_VALUE -9999" << endl;
}else {cout<<"griddata23 error message: Specify correct rastrType: ESRI or GRASS"<<endl;}
//-----------------------------        
    for (int j = 0; j < VerResolution; j++)
     {
      for (int i = 0; i < HorResolution; i++)
       {
        if (isChar)
          f << (int)gridPrev[(VerResolution-j-1)*HorResolution+i] << " ";
        else
          f << gridPrev[(VerResolution-j-1)*HorResolution+i] << " ";
       }
      f << endl;  
     }
    f.close(); 
  } else {
    cout << "Unable to save to file!" << endl;
  }
 } 

template <class TP>
void griddata2<TP>::set(int x, int y, TP val)
 {
  grid[y*HorResolution+x] = val;
 }

template <class TP>
void griddata2<TP>::setPrev(int x, int y, TP val)
 {
  gridPrev[y*HorResolution+x] = val;
 }

template <class TP>
void griddata2<TP>::inc(int x, int y, TP val)
 {
  grid[y*HorResolution+x] += val;
 }

template <class TP>
TP griddata2<TP>::get(int x, int y)
 {
  return (grid[y*HorResolution+x]);
 }

template <class TP>
TP griddata2<TP>::getPrev(int x, int y)
 {
  return (gridPrev[y*HorResolution+x]);
 }

template <class TP>
void griddata2<TP>::update()
 {
  memcpy(gridPrev,grid,VerResolution*HorResolution*sizeof(TP));
 }

template <class TP>
void griddata2<TP>::SetNeighNum(int n, int m)
 {
  HorNeigh = n;
  VerNeigh = m;
 }

template <class TP>
TP griddata2<TP>::GetMax(int x, int y)
 {
  int tmpx = x - HorNeigh;
  if (tmpx < 0) tmpx = HorResolution + tmpx;
  int tmpy = y - VerNeigh;
  if (tmpy < 0) tmpy = 0;
  TP maxv = gridPrev[tmpy*HorResolution+tmpx];
  for (int j = tmpy; j <= y + VerNeigh; j++)
   {
    if (j >= VerResolution) break;
    for (int i = -HorNeigh; i <= HorNeigh; i++)
     {
      int ii = x + i;
      if (ii >= HorResolution) ii -= HorResolution;
      else if (ii < 0) ii += HorResolution;
      if ((gridPrev[j*HorResolution+ii] > maxv) && !((ii == x) && (j == y)))
        maxv = gridPrev[j*HorResolution+ii];
     }
   }  
  return(maxv);
 }

template <class TP>
TP griddata2<TP>::GetMin(int x, int y)
 {
  int tmpx = x - HorNeigh;
  if (tmpx < 0) tmpx = HorResolution + tmpx;
  int tmpy = y - VerNeigh;
  if (tmpy < 0) tmpy = 0;
  TP minv = gridPrev[tmpy*HorResolution+tmpx];
  for (int j = tmpy; j <= y + VerNeigh; j++)
   {
    if (j >= VerResolution) break;
    for (int i = -HorNeigh; i <= HorNeigh; i++)
     {
      int ii = x + i;
      if (ii >= HorResolution) ii -= HorResolution;
      else if (ii < 0) ii += HorResolution;
      if ((gridPrev[j*HorResolution+ii] < minv) && !((ii == x) && (j == y)))
        minv = gridPrev[j*HorResolution+ii];
     }
   }  
  return(minv);
 }

template <class TP>
TP griddata2<TP>::GetAvg(int x, int y)
 {
  int count = 0;
  int tmpx = x - HorNeigh;
  if (tmpx < 0) tmpx = HorResolution + tmpx;
  int tmpy = y - VerNeigh;
  if (tmpy < 0) tmpy = 0;
  TP sumv = 0;
  for (int j = tmpy; j <= y + VerNeigh; j++)
   {
    if (j >= VerResolution) break;
    for (int i = -HorNeigh; i <= HorNeigh; i++)
     {
      int ii = x + i;
      if (ii >= HorResolution) ii -= HorResolution;
      else if (ii < 0) ii += HorResolution;
      count++;
      sumv += gridPrev[j*HorResolution+ii];
     }
   }  
  return(sumv/count);
 }

// Class constructor
template <class TP>
griddata2<TP>::griddata2(int HR, int VR, TP val)
 {
  HorResolution = HR;
  VerResolution = VR;
  HorNeigh = 1;
  VerNeigh = 1;
  grid = new TP[HorResolution*VerResolution];
  gridPrev = new TP[HorResolution*VerResolution];
  for (int j = 0; j < VerResolution; j++)
    for (int i = 0; i < HorResolution; i++)
     {
      grid[j*HorResolution+i] = val;
      gridPrev[j*HorResolution+i] = val;
     }
 }

// Default constructor
template <class TP>
griddata2<TP>::griddata2()
 {
  HorResolution = 720;
  VerResolution = 360;
  HorNeigh = 1;
  VerNeigh = 1;
  TP val = 0;
  grid = new TP[HorResolution*VerResolution];
  gridPrev = new TP[HorResolution*VerResolution];
  for (int j = 0; j < VerResolution; j++)
    for (int i = 0; i < HorResolution; i++)
     {
      grid[j*HorResolution+i] = val;
      gridPrev[j*HorResolution+i] = val;
     }
 }

// Copy constructor
template <class TP>
griddata2<TP>::griddata2(const griddata2& g)
 {
  HorResolution = g.HorResolution;
  VerResolution = g.VerResolution;
  HorNeigh = g.HorNeigh;
  VerNeigh = g.VerNeigh;
  grid = new TP[HorResolution*VerResolution];
  gridPrev = new TP[HorResolution*VerResolution];
  memcpy(grid,g.grid,HorResolution*VerResolution*sizeof(TP));
  memcpy(gridPrev,g.gridPrev,HorResolution*VerResolution*sizeof(TP));
 }

// Destructor
template <class TP>
griddata2<TP>::~griddata2()
 {
  delete []grid;
  delete []gridPrev;
 }

#endif
