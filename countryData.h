//  Name: CountryData class with statistics by regions
//  Author: Andriy Bun
//  Date: 29.10.09
//  Description: 
//
//	Modifications:
//      29.10.2009: added selecting countries to print; if no countries specified - prints all countries.
//		21.02.2011: added getTimeAvg, getSmoothAvg methods; renamed PrintToFile method to printToFile, now it takes full file name,
//					doesn't add extension; added copy constructor and assignment operator.


#ifndef COUNTRYDATA_H_
#define COUNTRYDATA_H_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>

//#include "ipol.h"

#define xmin(a,b) (((a)<(b))?(a):(b))
#define xmax(a,b) (((a)<(b))?(b):(a))

const int countriesNum = 245;

using namespace std;

class countryData
{
private:
	typedef set<int> set_t;
	vector< g4m::ipol<int, double> > values;
	vector< g4m::ipol<int, int> > count;
	vector< unsigned char > regions;
	set_t countriesToPrint;
	bool printAllCountries;
	void setRegions();
	double computeAvg(double *values, int *count, int timeIdx, int timePeriodWidth, int numYears);
public:
	countryData();
	countryData(const countryData & g);
	countryData & operator = (const countryData & g);
	~countryData();
	void reset(void);
	void setVal(int, int, double);
	void inc(int, int, double);
	double get(int countryIdx, int year);
	double getAvg(int countryIdx, int year);
	countryData getTimeAvg(int timePeriodWidth);
	countryData getSmoothAvg(int timePeriodWidth, int timeStep = 1);
	double getRegionSum(unsigned char regIdx, int year);
	void insertCountryToPrint(int countryIdx);
	set_t getListOfCountries();
	void setListOfCountries(set_t S);
	void printToFile(string fileName, int firstYear, int lastYear, int step, string statType = "VAL");  // print array to file
};

#endif
