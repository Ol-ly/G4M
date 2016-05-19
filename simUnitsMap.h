//  Name:           simUintsMap class
//  Author:         Andriy Bun
//  Date:           04.11.09
//  Description:    class for (pre-)processing binary maps of simulation units in 
//                  GLOBIOM GUI format
//  Modified by Mykola Gusti 22 June 2012 - correction of 0.5x0.5 SIMU shift by one 5-min cell row up; correction of one 5-min cell to the right.
//                                          correction of ASCII grid left corner coordinates


#ifndef SIMUNITSMAP_H_
#define SIMUNITSMAP_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include <string.h>
#include <cmath>

#define X_RES				4320		// Longitude resolution 
#define Y_RES				1674		// Latitude resolution
#define X_RES_BIG			720
#define Y_RES_BIG			360
#define RESOLUTION_RATIO	6			// proportion of resolutions (bigger cells to smaller)
#define NSIMU				212707		// Total number of simulation units

using namespace std;

class simUnitsMap
{
private:
	typedef vector<int> int_vector_t;

	// Some constants:
	static const float xMin;		// Minimum longitude value
	static const float yMin;		// Minimum latitude value

	// Data containers:
	int * ptr;									// ptr[i] indicates the number of points (cells) that belong to i-th simulation unit;
	int * simUMapInt;							// Internal simulation units map of resolution 0.5x0.5 for use within model
	int * simUMap;								// Map of simulation units as in simu.bin

	vector<int_vector_t> xPoints, yPoints;
	int round(float val);
public:
	typedef struct
	{
		int simu;
		char numCells;
		double simuFraction;
	} simu_info_struct_t;

	simUnitsMap();
	simUnitsMap(string fileName);
	simUnitsMap(const simUnitsMap & sMap);
	simUnitsMap& operator=(const simUnitsMap & sMap);
	~simUnitsMap();
	int getSimu(double x, double y);          //
	bool getSimuInfoByXY(double x, double y, vector<simu_info_struct_t> & result);
	int simuPerCell(double x, double y);    //
	void saveToFile();
	void saveToFile_ESRIGrid();
};

#endif
