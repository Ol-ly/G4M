#pragma once

#include <map>
#include <vector>
#include <deque>
#include <utility>
#include <cmath>
#include <limits>

#include "ageStruct_float.h"
#include "dataStruct_Europe_float.h"
#include "misc.h"
#include "forest_GUI_Europe_wt.h"
#include "griddata3_MG.h"


class HarvestData
{
public:
	float sawnW;
	float restW;
	float sawnThW;
	float restThW;
	float bmH;
	float bmTh;
	float harvL;

	void InitializeData (pair<g4m::ageStruct::v, g4m::ageStruct::v>&, double, double, int);
	void CleanData();
	HarvestData(void);
	~HarvestData(void);
};

