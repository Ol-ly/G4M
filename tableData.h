//  Name:         tableData class
//  Author:       Andriy Bun
//  Date:         13.11.09
//  Modified:     11.01.10 (Andriy Bun)
//  Description:  Class for storing and producing data for viewing with
//                GDagaView Software. Is a friend class for MDT class.

#ifndef TABLEDATA_H_
#define TABLEDATA_H_

#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <string>

#include "common.h"
#include "endianness.h"
#include "IntToStr.h"
#include "MDT.h"
//#include "MDT.cpp"

using namespace std;

class tableData
{
private:
	typedef vector<string> strVector;
	map<long long, float> data;
	MDT descr;
	long long N;                           // Number of data records per simulation unit
public:
	tableData();
	tableData(string fileName_GDT);
	~tableData();
	// Inserts a value "val" corresponding to a vector of coordinates "point"
	// into the list.
	void insert(float val, strVector point);
	// Updates a value "val" corresponding to a vector of coordinates "point"
	// into the list.
	void update(float val, strVector point);
	// Rename dataset:
	void rename(string name);
	// Rename dimensions of the dataset:
	void renameDims(strVector vec);
	// Add new dimension:
	void addDim(string dimName, set<string> elements);
	void addDim(string dimName, set<int> elements);
	void addDim(string dimName, string element);
	// Clearing object
	void clear();
	// Writing object to files *.GDT and *.GDC files
	bool SaveToFile(string outDir, string fileName);
};

#endif
