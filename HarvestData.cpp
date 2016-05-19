#include "HarvestData.h"


HarvestData::HarvestData(void)
{
	CleanData();
}
HarvestData::~HarvestData(void)
{
}
void HarvestData::CleanData(void)
{
	sawnW=0.;
	restW=0.;
	sawnThW=0.;
	restThW=0.;
	bmH=0.;
	bmTh=0.;
	harvL=0.;
}
void HarvestData::InitData (pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, double realArea, double harvArea, int modTimeStep)
{
	sawnW = res.second.sw * harvArea/realArea/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
	restW = res.second.rw * harvArea/realArea/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
	sawnThW = res.first.sw/realArea/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
	restThW = res.first.rw/realArea/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
	bmH = res.second.bm * harvArea/realArea/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
	bmTh = res.first.bm/realArea/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
	harvL = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) ; // MG: harvest residues for the set (old) forest tC/ha
}
