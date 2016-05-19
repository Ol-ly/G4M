#include "simUnitsMap.h"

// Some constants:
const float simUnitsMap::xMin = -179.9583332f;   // Minimum longitude value
const float simUnitsMap::yMin = -55.87501355f;   // Minimum latitude value
//const float simUnitsMap::yMin = -55.958333333f;   // Minimum latitude value

// Default constructor
simUnitsMap::simUnitsMap()
{
	// Creating internal coordinate system:
	ptr = new int [NSIMU+1];
	simUMap = new int[X_RES * Y_RES];
	simUMapInt = new int[X_RES_BIG * Y_RES_BIG];
	for (int yy = 0; yy < Y_RES; yy++)
	{
		for (int xx = 0; xx < X_RES; xx++)
		{
			// internal simulation units:
//			int x = int((float)xx/RESOLUTION_RATIO);
//			int y = int((float)yy/RESOLUTION_RATIO);
			int x = int(float(xx)/float(RESOLUTION_RATIO));
			int y = int(float(yy)/float(RESOLUTION_RATIO));			
			int intSimUnit = y * X_RES_BIG + x;
			simUMapInt[intSimUnit] = intSimUnit;
			// model simulation units:
			simUMap[yy * X_RES + xx] = intSimUnit;
			ptr[intSimUnit+1]++;
		}
	}
	ptr[0] = 0;
	for (int i = 1; i <= NSIMU; i++)
	{
		ptr[i] = ptr[i-1] + ptr[i];
	}
}

// Constructor
simUnitsMap::simUnitsMap(string fileName)
{
	ptr = new int [NSIMU+1];
	simUMap = new int[X_RES * Y_RES];
	simUMapInt = new int[X_RES_BIG * Y_RES_BIG];
	for (int i = 0; i < X_RES * Y_RES; i++) simUMap[i] = -1;
	ifstream f;
	f.open(fileName.c_str(), ios::in | ios::binary);
	if (f.is_open())
	{
		f.read(reinterpret_cast<char *>(ptr), sizeof(int) * (NSIMU + 1));
		int NPTS = ptr[NSIMU];
		int SIMU = 0;
		for (int i = 0; i < NPTS; i++)
		{
			if (i > ptr[SIMU]) SIMU++;
			int xx, yy;
			if (f.eof())
			{
				cout << "End of file!" << endl;
				break;
			}
			f.read(reinterpret_cast<char *>(&xx), sizeof(int));
			if (f.eof())
			{
				cout << "End of file!" << endl;
				break;
			}
			f.read(reinterpret_cast<char *>(&yy), sizeof(int));
			// model simulation units:
			simUMap[yy * X_RES + xx] = SIMU;
			// internal simulation units:
                       
			int x = int(float(xx)/float(RESOLUTION_RATIO));
//			int y = int(float(yy)/float(RESOLUTION_RATIO)); MG: corrected to avoid 5 min shift on lat
            int y=0;
            if (yy==0) {y = 0;} else {y = int(float(yy+5.)/float(RESOLUTION_RATIO));} //MG: corrected to avoid one-pixel tail in the 0.5x0.5 cells

			int intSimUnit = y * X_RES_BIG + x;
			simUMapInt[intSimUnit] = intSimUnit;
			simUMap[yy * X_RES + xx] = intSimUnit;
//			ptr[intSimUnit+1]++;
		}

//MG: Corrected (commented) to allow creation of G4M simu.bin file from the original GUI simu.bin
/*
//For creation of G4M simu.bin from original GLOBIOM simu.bin uncomment the lines below.
// // // Read the original GLOBIOM simu.bin: simUnitsMap sMap = simUnitsMap("simu_GLOBIOM.bin");
// // // Save G4M simu.bin:                  sMap.saveToFile();
  delete [] ptr;
  ptr = new int[NSIMU+1];
  for (int j = 0; j < Y_RES; j++) {
    for (int i = 0; i < X_RES; i++) {
         int tmp = simUMap[j * X_RES + i];
         if (tmp > 0)
          ptr[tmp+1]++;
//cout<<"SIMU[\t"<<j * xRes + i<<"\t]=\t"<<tmp<<"\tptr1=\t"<<ptr[tmp+1]<<endl;
    }
  }
//cout<<"SIMU[\t"<<(yRes-1) * xRes + (xRes-1)<<"\t]=\t"<<simUMap[(yRes-1) * xRes + (xRes-1)]<<"\tptr1=\t"<<ptr[(yRes-1) * xRes + (xRes-1)+1]<<endl;

    ptr[0] = 0;
    for (int i = 1; i <= NSIMU; i++) {
      ptr[i] = ptr[i-1] + ptr[i];
    }
*/
//----------------------------------------


		f.close();
		cout << "Successfully read from binary file: " << fileName << endl;
	}
	else
	{
		cout << "Unable to open file!" << endl;
	}
}

// Copy constructor
simUnitsMap::simUnitsMap(const simUnitsMap & sMap)
{
	ptr = new int [NSIMU+1];
	simUMap = new int[X_RES * Y_RES];
	simUMapInt = new int[X_RES_BIG * Y_RES_BIG];
	memcpy(ptr, sMap.ptr, (NSIMU+1) * sizeof(int));
	memcpy(simUMap, sMap.simUMap, (X_RES * Y_RES) * sizeof(int));
	memcpy(simUMapInt, sMap.simUMapInt, (X_RES_BIG * Y_RES_BIG) * sizeof(int));
}

// Assignment operator
simUnitsMap& simUnitsMap::operator=(const simUnitsMap & sMap)
{
	if (this != &sMap)
	{
		memcpy(this->ptr, sMap.ptr, (NSIMU+1) * sizeof(int));
		memcpy(this->simUMap, sMap.simUMap, (X_RES * Y_RES) * sizeof(int));
		memcpy(this->simUMapInt, sMap.simUMapInt, (X_RES_BIG * Y_RES_BIG) * sizeof(int));
	}
	return *this;
}

// Destructor
simUnitsMap::~simUnitsMap()
{
	delete []ptr;
	delete []simUMap;
	delete []simUMapInt;
}

//================
// Class methods:
//================

int simUnitsMap::round(float val)
{
	if ((val + 0.5) >= (int(val) + 1))
		return int(val)+1;
	else
		return int(val);
}

int simUnitsMap::getSimu(double x, double y)
{
//	int xID = round(2. * (x - xMin)); //MG:
	int xID = round(12. * (x - xMin));
//	int yID = Y_RES_BIG - 1 - round(2. * (y - yMin)) ; //MG:
	int yID = Y_RES - 1 - round(12. * (y - yMin)) ;	
//	if ((xID < 0) || (xID >= X_RES_BIG) || (yID < 0) || (yID >= Y_RES_BIG)) return -2;//MG:
	if ((xID < 0) || (xID >= X_RES) || (yID < 0) || (yID >= Y_RES)) return -2;	
//	return simUMap[yID * X_RES_BIG + xID];//MG:
	return simUMap[yID * X_RES + xID];	
}


bool simUnitsMap::getSimuInfoByXY(double x, double y, vector<simUnitsMap::simu_info_struct_t> & result)
{
	int numActiveCells = 0;
	int xID = round(2. * RESOLUTION_RATIO * (x - xMin));
	int yID = Y_RES - 1 - round(2. * RESOLUTION_RATIO * (y - yMin));
	if ((xID < 0) || (xID >= X_RES) || (yID < 0) || (yID >= Y_RES)) return false;

	for (int i = RESOLUTION_RATIO - 1; (i >= 0) && ((xID + i) < X_RES); i--)
	{
		for (int j = RESOLUTION_RATIO - 1; (j >= 0) && ((yID + j) < Y_RES); j--)
		{
			int simu = simUMap[(yID + j) * X_RES + (xID + i)];
			if (simu >= 0)
			{
				numActiveCells++;
				bool simuIsInBlock = false;
				for (int it = 0; it < result.size(); it++)
				{
					if (simu == result[it].simu)
					{
						result[it].numCells++;
						simuIsInBlock = true;
					}
				}
				if (!simuIsInBlock)
				{
					simu_info_struct_t simuInfo;
					simuInfo.simu = simu;
					simuInfo.numCells = 1;
					result.push_back(simuInfo);
				}
			}
		}
	}
	for (int i = 0; i < result.size(); i++)
	{
		result[i].simuFraction = double(result[i].numCells) / numActiveCells;
	}
	return true;
}

int simUnitsMap::simuPerCell(double x, double y)
{
	set<int> res;
	int xID = round(2. * RESOLUTION_RATIO * (x - xMin));
	int yID = Y_RES - 1 - round(2. * RESOLUTION_RATIO * (y - yMin)) ;
	if ((xID < 0) || (xID >= X_RES) || (yID < 0) || (yID >= Y_RES)) ;
	else
	{
		for (int i = RESOLUTION_RATIO - 1; i >= 0; i--)
		{
			if ((xID + i) < X_RES)
			{
				for (int j = RESOLUTION_RATIO - 1; j >= 0; j--)
				{
					if ((yID + j) < Y_RES)
					{
						if (simUMap[(yID + j) * X_RES + (xID + i)] >= 0)
						{
							res.insert(simUMap[(yID + j) * X_RES + (xID + i)]);
						}
					}
				}
			}
		}
	}
	int result = res.size();
	return result;
}

void simUnitsMap::saveToFile()
{
	xPoints.resize(NSIMU);
	yPoints.resize(NSIMU);
	for (int j = 0; j < Y_RES; j++)
	{
		for (int i = 0; i < X_RES; i++)
		{
			int tmp = simUMap[j * X_RES + i];
			if (tmp >=0)
			{
				xPoints[tmp].push_back(i);
				yPoints[tmp].push_back(j);
			}
		}
	}
	int vSize = 2 * ptr[NSIMU];
	// temporary array for output of geographic indeces:
	int * v = new int[vSize];
	int id = 0;
	for (int i = 0; i < NSIMU; i++)
	{
		for (int j = 0; j < xPoints[i].size(); j++)
		{
			v[2 * id]     = xPoints[i][j];
			v[2 * id + 1] = yPoints[i][j];
			id++;
		}
	}
	// Writing to file:
	ofstream f;
	string fileName = "simu.bin";
	f.open(fileName.c_str(), ios::out | ios::binary);
	if (f.is_open())
	{
		f.write(reinterpret_cast<char *>(ptr), (NSIMU + 1) * sizeof(int));
		f.write(reinterpret_cast<char *>(v), vSize * sizeof(int));
		f.close();
		cout << "Successfully written to binary file: " << fileName << endl;
	}
	else
	{
		cout << "Unable to save to file!" << endl;
	}
	delete []v;
}

void simUnitsMap::saveToFile_ESRIGrid()
{
	ofstream f;
	string fileName = "SIMU.asc";
	f.open(fileName.c_str(),ios::out);
	if (f.is_open())
	{
		f << "NCOLS " << X_RES << endl;
		f << "NROWS " << Y_RES << endl;
		f << "XLLCORNER " << xMin-0.5/(2*RESOLUTION_RATIO) << endl;
		f << "YLLCORNER " << yMin-0.5/(2*RESOLUTION_RATIO) << endl;
		f << "CELLSIZE " << 0.5/RESOLUTION_RATIO << endl;
//		f << "NODATA_VALUE -9999" << endl;
		f << "NODATA_VALUE -1" << endl;
		for (int j = 0; j < Y_RES; j++)
		{
			for (int i = 0; i < X_RES; i++)
			{
				int tmp = simUMap[j * X_RES + i];
				f << tmp << " ";
			}
			f << endl;
		}
		f.close();     
		cout << "Successfully written to ASCII file: " << fileName << endl;
	}
	else
	{
		cout << "Unable to save to ASCII file!" << endl;
	}
}
