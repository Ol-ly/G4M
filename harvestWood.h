#pragma once

#include <utility>

#include "misc.h"


using namespace std;
class harvestWood
{
private:
	float harvAreaO;       //	
	float sawnW ;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
	float restW;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
	float sawnThW;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
	float restThW;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
	float bmH;
	float bmTh;

	void ClearData ();

	//double bmH;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
	//double bmTh;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning 		
	//double harvRes; // MG: usable harvest residues for the set (old) forest tC/ha
	//double harvestSO;       // harvested sawnwood, tC/ha
	//double harvestRO;       //harvested restwood, tC/ha	
	//double harvSawnWood; //harvested sawnwood, m3/ha
	//double harvRestWood; //harvested restwood, m3/ha 
	//double cohortRes_SW(double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int, double, int);		// return sawnwood in m3/ha
	//double cohortRes_RW(double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int, double, int);		// return restwood in m3/ha
	//double UsableResidues(double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int);
	//double cohortRes_SW_adj (double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int);
	//double cohortRes_RW_adj (double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int);

public:
	float harvestWS;
	float harvestWR;
	float sumSW;	// total sawnwood in the cell, tC
	float sumRW;	// total restwood in the cell, tC
	float harvL;
	double thinning;
	double fcut;

	harvestWood(void);	//Constructor
	~harvestWood(void);		//Destructor
	void calcHarvM3 (double, std::pair <g4m::ageStruct::v, g4m::ageStruct::v>&, int, double, int);	 // calculate harvested wood in m3/ha
	void calculateHarvest (double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int);	 // calculate harvested wood per 1 ha
	void InitData (pair<g4m::ageStruct::v, g4m::ageStruct::v>&, double, double, int);
	//void tmpWoodHarv (double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int);
	//	void Inicialize(double, pair<g4m::ageStruct::v, g4m::ageStruct::v> &, int);
};


