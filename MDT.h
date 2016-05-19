//  Name:           MDT class (Map Data Textual description)
//  Author:         Andriy Bun
//  Date:           06.11.09
//  Modified:       11.01.10 (Andriy Bun)
//  Description:    Class for managing textual description of datasets

#ifndef MDT_H_
#define MDT_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "common.h"
#include "IntToStr.h"
//#include "simUnitsData.h"
//#include "tableData.h"

using namespace std;

class MDT {
	friend class simUnitsData;
	friend class tableData;
private:
	typedef vector<string> strVector;
	typedef vector<int> int_vector_t;
	// Vectors of dimensions' description:
	int_vector_t dimCardinals;
	strVector dimNames;
	vector<strVector> dimElements;
	string paramName;
	int nDims;
public:
	MDT();
	MDT(string fileName);
	~MDT();
	void clear();
	// Rename parameter name "paramName":
	void rename(string name);
	// Adding a new dimension "dimName" to the space, with vector of "elements".
	// Returns true if successfully added, and false, if dimension already exists
	bool addDim(string dimName, strVector elements);
	// Adding a new element "element" to a dimension "dimName". If dimension doesn't
	// exist, a new dimension with "dimName" is created with one "element". Returns
	// true if successfully added
	bool addDimEl(string dimName, string element);
	// Returns maximum hash value
	long long getN();
	// Returns number of dimensions
	int getNDims();
	// Returns HASH value of a point with vector of coordinates "elements". If
	// vector "elements" has more dimensions that are currently available, new 
	// dimensions are created
	long long getHash(strVector elements);
	// Returns coordinate of "element" in dimension "dim". If "element" doesn't
	// exist within "dim", adds it to the dimension
	int getCoordinate(int dim, string element);
	// Saves description to *.MDT file (if mode == "MAP") or *.GDT file (if mode == "TABLE")
	void SaveToFile(string fileName, string mode);
};

#endif
