#include "countryData.h"

void countryData::reset(void)
{
	int ind = 0;
	do
	{
		values[ind].clear();
		count[ind].clear();
		ind++;
	}
	while (ind <countriesNum);
}


void countryData::setVal(int ind, int year, double val)
{
	values[ind].insert(year,val);
	count[ind].insert(year,1);
}

void countryData::inc(int ind, int year, double val)
{
	values[ind].inc(year,val);
	count[ind].inc(year,1);
}

double countryData::get(int countryIdx, int year)
{
	return(values[countryIdx][year]);
}

double countryData::getAvg(int countryIdx, int year)
{
	return(values[countryIdx][year]/count[countryIdx][year]);
}

countryData countryData::getTimeAvg(int timePeriodWidth)
{
	countryData result = getSmoothAvg(timePeriodWidth, timePeriodWidth);
	return result;
}

countryData countryData::getSmoothAvg(int timePeriodWidth, int timeStep)
{
	countryData result;

	if (timePeriodWidth % 2 != 1) 
	{
		cout << "Error! Please provide odd number for time period!" << endl;
		return result;
	}

	map<int, double>::iterator valuesIter, valuesIterLastEl;
	map<int, int>::iterator countIter, countIterLastEl;

	int sumValues = 0;
	int sumCount = 0;

	for (int countryIdx = 0; countryIdx < values.size(); countryIdx++)
	{
		if (values[countryIdx].aMap.size() > 0)
		{
			valuesIter = values[countryIdx].aMap.begin();
			valuesIterLastEl = --values[countryIdx].aMap.end();
			countIter = count[countryIdx].aMap.begin();
			countIterLastEl = --count[countryIdx].aMap.end();
			int firstYear = valuesIter->first;
			int lastYear = valuesIterLastEl->first;
			int numYears = lastYear - firstYear + 1;
			double *valuesByYears = new double[numYears];
			int *countByYears = new int[numYears];
			memset(countByYears, 0, numYears * sizeof(int));

			do
			{
				valuesByYears[valuesIter->first - firstYear] = valuesIter->second;
				countByYears[countIter->first - firstYear] = countIter->second;
			}
			while ((valuesIter++ != valuesIterLastEl) && (countIter++ != countIterLastEl));

			for (int yearIdx = 0; yearIdx < numYears; yearIdx += timeStep)
			{
				result.setVal(countryIdx, firstYear + yearIdx, 
					computeAvg(valuesByYears, countByYears, yearIdx, timePeriodWidth, numYears));
			}

			delete [] valuesByYears;
			delete [] countByYears;
		}
	}
	return result;
}

double countryData::computeAvg(double *values, int *count, int yearIdx, int timePeriodWidth, int numYears)
{
	double sum = (double)0;
	int num = 0;
	int timeHalfPeriodWidth = (timePeriodWidth - 1) / 2;
	for(int idx = xmax(0, yearIdx - timeHalfPeriodWidth); idx <= xmin(numYears-1, yearIdx+timeHalfPeriodWidth); idx++)
	{
		num += count[idx];
		sum += (count[idx] == 0 ? 0 : values[idx]);
	}
	return sum / num;
}

double countryData::getRegionSum(unsigned char regIdx, int year)
{
	double sum = 0;
	int ind = 0;
	do
	{
		if (regions[ind] == regIdx)
		{
			sum += get(ind, year);
		}
		ind++;
	} 
	while (ind <countriesNum);  
	return sum;
}

void countryData::insertCountryToPrint(int countryIdx)
{
	countriesToPrint.insert(countryIdx);
	printAllCountries = false;
}

set<int> countryData::getListOfCountries()
{
	return countriesToPrint;
}

void countryData::setListOfCountries(set_t S)
{
	countriesToPrint = S;
	if (countriesToPrint.size() > 0) printAllCountries = false;
}

void countryData::printToFile(string fileName, int firstYear, int lastYear, int step, string statType)
{
	ofstream f;
	f.open(fileName.c_str(),ios::out);
	if (f.is_open())
	{
		f << "CountryN/Year"; 
		for (int j = firstYear; j <= lastYear; j += step)
		{
			f <<"\t"<<j;
		}
		f << endl;
		//-----------------------------        
		for (int i = 0; i < countriesNum; i++)
		{
			if ((printAllCountries) || (countriesToPrint.find(i) != countriesToPrint.end()))
			{
				f << i ;
				for (int j = firstYear; j <= lastYear; j += step)
				{
					if (statType == "VAL")
						f<<"\t"<<values[i][j];
					else if (statType == "AVG")
						f << "\t" << getAvg(i,j);
				}
				f << endl;  
			}
		}
		f.close(); 
		cout << "Successfully written to:" << fileName << endl;
	}
	else
	{
		cout << "Unable to save to file!" << endl;
	}
}

// default constructor
countryData::countryData()
{
	printAllCountries = true;
	values.resize(countriesNum);
	regions.resize(countriesNum);
	count.resize(countriesNum);
	setRegions();
}

// copy constructor
countryData::countryData(const countryData & g)
{
	values = g.values;
	count = g.count;
	regions = g.regions;
	countriesToPrint = g.countriesToPrint;
	printAllCountries = g.printAllCountries;
}

// assignment operator
countryData & countryData::operator = (const countryData & g)
{
	if (this != &g)
	{
		values = g.values;
		count = g.count;
		regions = g.regions;
		countriesToPrint = g.countriesToPrint;
		printAllCountries = g.printAllCountries;
	}
	return *this;
}

// destructor
countryData::~countryData()
{

}

void countryData::setRegions()
{
	regions.assign(countriesNum, 0);
	regions[	1	]=	0	;
	regions[	2	]=	20	;
	regions[	3	]=	25	;
	regions[	4	]=	0	;
	regions[	5	]=	0	;
	regions[	6	]=	17	;
	regions[	7	]=	0	;
	regions[	8	]=	16	;
	regions[	9	]=	14	;
	regions[	10	]=	19	;
	regions[	11	]=	10	;
	regions[	12	]=	0	;
	regions[	13	]=	0	;
	regions[	14	]=	0	;
	regions[	15	]=	0	;
	regions[	16	]=	1	;
	regions[	17	]=	7	;
	regions[	18	]=	10	;
	regions[	19	]=	25	;
	regions[	20	]=	7	;
	regions[	21	]=	25	;
	regions[	22	]=	25	;
	regions[	23	]=	20	;
	regions[	24	]=	6	;
	regions[	25	]=	14	;
	regions[	26	]=	16	;
	regions[	27	]=	17	;
	regions[	28	]=	10	;
	regions[	29	]=	16	;
	regions[	30	]=	16	;
	regions[	31	]=	19	;
	regions[	32	]=	2	;
	regions[	33	]=	16	;
	regions[	34	]=	21	;
	regions[	35	]=	20	;
	regions[	36	]=	0	;
	regions[	37	]=	25	;
	regions[	38	]=	25	;
	regions[	39	]=	3	;
	regions[	40	]=	0	;
	regions[	41	]=	18	;
	regions[	42	]=	19	;
	regions[	43	]=	4	;
	regions[	44	]=	25	;
	regions[	45	]=	25	;
	regions[	46	]=	25	;
	regions[	47	]=	25	;
	regions[	48	]=	0	;
	regions[	49	]=	19	;
	regions[	50	]=	25	;
	regions[	51	]=	25	;
	regions[	52	]=	16	;
	regions[	53	]=	16	;
	regions[	54	]=	0	;
	regions[	55	]=	0	;
	regions[	56	]=	9	;
	regions[	57	]=	6	;
	regions[	58	]=	7	;
	regions[	59	]=	25	;
	regions[	60	]=	16	;
	regions[	61	]=	8	;
	regions[	62	]=	16	;
	regions[	63	]=	14	;
	regions[	64	]=	19	;
	regions[	65	]=	14	;
	regions[	66	]=	25	;
	regions[	67	]=	0	;
	regions[	68	]=	9	;
	regions[	69	]=	5	;
	regions[	70	]=	25	;
	regions[	71	]=	8	;
	regions[	72	]=	15	;
	regions[	73	]=	0	;
	regions[	74	]=	7	;
	regions[	75	]=	0	;
	regions[	76	]=	0	;
	regions[	77	]=	25	;
	regions[	78	]=	8	;
	regions[	79	]=	10	;
	regions[	80	]=	25	;
	regions[	81	]=	18	;
	regions[	82	]=	25	;
	regions[	83	]=	0	;
	regions[	84	]=	25	;
	regions[	85	]=	25	;
	regions[	86	]=	0	;
	regions[	87	]=	9	;
	regions[	88	]=	16	;
	regions[	89	]=	0	;
	regions[	90	]=	16	;
	regions[	91	]=	0	;
	regions[	92	]=	0	;
	regions[	93	]=	19	;
	regions[	94	]=	0	;
	regions[	95	]=	0	;
	regions[	96	]=	16	;
	regions[	97	]=	17	;
	regions[	98	]=	16	;
	regions[	99	]=	6	;
	regions[	100	]=	21	;
	regions[	101	]=	11	;
	regions[	102	]=	0	;
	regions[	103	]=	8	;
	regions[	104	]=	14	;
	regions[	105	]=	14	;
	regions[	106	]=	18	;
	regions[	107	]=	14	;
	regions[	108	]=	9	;
	regions[	109	]=	16	;
	regions[	110	]=	14	;
	regions[	111	]=	12	;
	regions[	112	]=	10	;
	regions[	113	]=	25	;
	regions[	114	]=	10	;
	regions[	115	]=	22	;
	regions[	116	]=	15	;
	regions[	117	]=	0	;
	regions[	118	]=	24	;
	regions[	119	]=	14	;
	regions[	120	]=	22	;
	regions[	121	]=	14	;
	regions[	122	]=	25	;
	regions[	123	]=	14	;
	regions[	124	]=	16	;
	regions[	125	]=	0	;
	regions[	126	]=	20	;
	regions[	127	]=	25	;
	regions[	128	]=	5	;
	regions[	129	]=	7	;
	regions[	130	]=	5	;
	regions[	131	]=	0	;
	regions[	132	]=	14	;
	regions[	133	]=	0	;
	regions[	134	]=	10	;
	regions[	135	]=	25	;
	regions[	136	]=	20	;
	regions[	137	]=	13	;
	regions[	138	]=	0	;
	regions[	139	]=	17	;
	regions[	140	]=	25	;
	regions[	141	]=	9	;
	regions[	142	]=	21	;
	regions[	143	]=	17	;
	regions[	144	]=	22	;
	regions[	145	]=	0	;
	regions[	146	]=	25	;
	regions[	147	]=	25	;
	regions[	148	]=	0	;
	regions[	149	]=	25	;
	regions[	150	]=	0	;
	regions[	151	]=	25	;
	regions[	152	]=	21	;
	regions[	153	]=	0	;
	regions[	154	]=	25	;
	regions[	155	]=	0	;
	regions[	156	]=	25	;
	regions[	157	]=	0	;
	regions[	158	]=	25	;
	regions[	159	]=	16	;
	regions[	160	]=	0	;
	regions[	161	]=	7	;
	regions[	162	]=	18	;
	regions[	163	]=	20	;
	regions[	164	]=	0	;
	regions[	165	]=	1	;
	regions[	166	]=	14	;
	regions[	167	]=	20	;
	regions[	168	]=	16	;
	regions[	169	]=	0	;
	regions[	170	]=	19	;
	regions[	171	]=	21	;
	regions[	172	]=	0	;
	regions[	173	]=	15	;
	regions[	174	]=	6	;
	regions[	175	]=	0	;
	regions[	176	]=	22	;
	regions[	177	]=	9	;
	regions[	178	]=	19	;
	regions[	179	]=	0	;
	regions[	180	]=	0	;
	regions[	181	]=	14	;
	regions[	182	]=	0	;
	regions[	183	]=	6	;
	regions[	184	]=	10	;
	regions[	185	]=	25	;
	regions[	186	]=	14	;
	regions[	187	]=	25	;
	regions[	188	]=	25	;
	regions[	189	]=	21	;
	regions[	190	]=	0	;
	regions[	191	]=	0	;
	regions[	192	]=	0	;
	regions[	193	]=	15	;
	regions[	194	]=	25	;
	regions[	195	]=	25	;
	regions[	196	]=	0	;
	regions[	197	]=	25	;
	regions[	198	]=	0	;
	regions[	199	]=	17	;
	regions[	200	]=	25	;
	regions[	201	]=	19	;
	regions[	202	]=	6	;
	regions[	203	]=	6	;
	regions[	204	]=	8	;
	regions[	205	]=	25	;
	regions[	206	]=	25	;
	regions[	207	]=	14	;
	regions[	208	]=	0	;
	regions[	209	]=	25	;
	regions[	210	]=	25	;
	regions[	211	]=	21	;
	regions[	212	]=	10	;
	regions[	213	]=	0	;
	regions[	214	]=	10	;
	regions[	215	]=	0	;
	regions[	216	]=	0	;
	regions[	217	]=	16	;
	regions[	218	]=	14	;
	regions[	219	]=	26	;
	regions[	220	]=	0	;
	regions[	221	]=	4	;
	regions[	222	]=	25	;
	regions[	223	]=	25	;
	regions[	224	]=	10	;				
	regions[	225	]=	0	;
	regions[	226	]=	19	;
	regions[	227	]=	27	;
	regions[	228	]=	10	;
	regions[	229	]=	0	;
	regions[	230	]=	16	;
	regions[	231	]=	19	;
	regions[	232	]=	0	;
	regions[	233	]=	0	;
	regions[	234	]=	22	;
	regions[	235	]=	15	;
	regions[	236	]=	0	;
	regions[	237	]=	0	;
	regions[	238	]=	14	;
	regions[	239	]=	23	;
	regions[	240	]=	25	;
	regions[	241	]=	25	;
	regions[	242	]=	0	;
	regions[	243	]=	0	;
	regions[	244	]=	0	;

}
