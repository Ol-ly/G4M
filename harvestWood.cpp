#include "harvestWood.h"

#include <utility>


#include "misc.h"
harvestWood::harvestWood(void)
{
	
}
void harvestWood::calculateHarvest(double realArea, pair<g4m::ageStruct::v, g4m::ageStruct::v> & res, int modTimeStep)
{
	ClearData ();
	harvAreaO=res.second.area;
	if (realArea > 0) {
		sawnW = res.second.sw * harvAreaO/realArea;      			
		sawnThW = res.first.sw/realArea;    
		restW = res.second.rw * harvAreaO/realArea;       
		restThW = res.first.rw/realArea;   
		harvestWS = (sawnW+sawnThW)/modTimeStep;
		harvestWR = (restW+restThW)/modTimeStep;

	}

}

void harvestWood::calcHarvM3 (double realArea, std::pair<g4m::ageStruct::v, g4m::ageStruct::v> & res, int modTimeStep, double units, int biomasRot)
{
	ClearData ();
	harvAreaO=res.second.area;
	if (realArea > 0) {
		sawnW = res.second.sw * harvAreaO/realArea;      			
		sawnThW = res.first.sw/realArea;    
		restW = res.second.rw * harvAreaO/realArea;       
		restThW = res.first.rw/realArea;   
		harvestWS = (sawnW+sawnThW)/modTimeStep * units;
		harvestWR = (restW+restThW)/modTimeStep * units;
	}
	if ((biomasRot ==0)||(harvestWS < 0)) harvestWS = 0.;
	if ((biomasRot ==0)||(harvestWR < 0)) harvestWR = 0.;
}

void harvestWood::InitData (std::pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, double harvArea, double realArea, int modTimeStep)
{
	ClearData ();
//	cout << "CLASS: harvArea	" << harvArea << endl;
//	cout << "CLASS: realArea	" << realArea << endl;
	if (realArea > 0) {
	sawnW = res.second.sw * harvArea/realArea/modTimeStep;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
	restW = res.second.rw * harvArea/realArea/modTimeStep;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
	sawnThW = res.first.sw/realArea/modTimeStep;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
	restThW = res.first.rw/realArea/modTimeStep;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
	bmH = res.second.bm * harvArea/realArea/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
	bmTh = res.first.bm/realArea/modTimeStep;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning 
	sumSW = sawnW + sawnThW;
	sumRW = restW + restThW;
	//cout << "CLASS: sawnW	" << sawnW << endl;
	//cout << "CLASS: sawnThW	" << sawnThW << endl;
	
	harvL = (bmH+bmTh - (sumSW + sumRW)) ; // MG: harvest residues for the set (old) forest tC/ha
	thinning = sawnThW + restThW;
	fcut = sawnW + restW;
	}
}

void harvestWood::ClearData ()
{
	sawnW = 0.;
	restW = 0.;
	sawnThW = 0.;
	restThW = 0.;
	harvestWS = 0.;
	harvestWR = 0.;
	harvAreaO = 0.;
	sumSW = 0.;
	sumRW = 0.;
	harvL = 0.;
	thinning = 0.;
	fcut = 0.;
}
harvestWood::~ harvestWood(void)
{ 
}


	//void harvestWood::tmpWoodHarv (double realArea, pair<g4m::ageStruct::v, g4m::ageStruct::v> & res, int modTimeStep)
	//{
	//	calculateHarvest (realAreaO, resTmp, modTimeStep);
	//}
	//void harvestWood::Inicialize(double area, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, int modTimeStep){
	//	//TODO add clear fields;
	//	harvAreaO = res.second.area;
	//	harvestRO = 0.; 
	//	harvestSO = 0.;
	//	if (area>0) {
	//		sawnW = res.second.sw * harvAreaO/area;      			
	//		sawnThW = res.first.sw/area;    						
	//		harvestSO = (sawnW + sawnThW)/modTimeStep;	
	//
	//		restW = res.second.rw * harvAreaO/area;       
	//		restThW = res.first.rw/area;    
	//		harvestRO = (restW + restThW)/modTimeStep;
	//	}
	//	so = harvestSO * units; 
	//	if ((biomasRot ==0)||(so < 0)) 
	//		so = 0.;
	//	ro = harvestRO * units; 
	//	if ((biomasRot ==0)||(ro < 0)) 
	//		ro = 0.;
	//}

	//double harvestWood::cohortRes_SW (double realAreaO, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, int modTimeStep, double units, int biomasRot)
	//{
	//	harvAreaO = res.second.area;
	//	harvestSO = 0.;
	//	if (realAreaO>0) {
	//		sawnW = res.second.sw * harvAreaO/realAreaO;      			
	//		sawnThW = res.first.sw/realAreaO;    						
	//		harvestSO = (sawnW + sawnThW)/modTimeStep;	
	//	}
	//	harvSawnWood = harvestSO * units; 
	//	if ((biomasRot ==0)||(harvSawnWood < 0)) harvSawnWood = 0.;		
	//	return(harvSawnWood) ;
	//}
	//
	//double harvestWood::cohortRes_RW (double realAreaO, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, int modTimeStep, double units, int biomasRot)
	//{	
	//	harvAreaO = res.second.area;
	//	harvestRO = 0.; 
	//	if (realAreaO>0) {
	//		restW = res.second.rw * harvAreaO/realAreaO;       
	//		restThW = res.first.rw/realAreaO;    
	//		harvestRO = (restW + restThW)/modTimeStep;
	//	}
	//	harvRestWood = harvestRO * units; 
	//	if ((biomasRot ==0)||(harvRestWood < 0)) harvRestWood = 0.;		
	//	return(harvRestWood) ;
	//
	//}
	//
	//double harvestWood::UsableResidues (double realAreaO, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, int modTimeStep)
	//{
	//	harvAreaO = res.second.area;
	//	if (realAreaO>0) {
	//		bmH = res.second.bm * harvAreaO/realAreaO;       
	//		bmTh = res.first.bm/realAreaO;                 
	//		harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; 		  
	//	}
	//	return(harvRes);
	//}
	//
	//double harvestWood::cohortRes_SW_adj (double Area, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, int modTimeStep)
	//{
	//	harvAreaO = res.second.area;
	//	harvestSO = 0.;
	//	if (Area>0) {
	//		sawnW = res.second.sw * harvAreaO/Area;      			
	//		sawnThW = res.first.sw/Area;    						
	//		harvestSO = (sawnW + sawnThW)/modTimeStep;	
	//	}			
	//	return(harvestSO) ;
	//}
	//
	//double harvestWood::cohortRes_RW_adj (double Area, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res, int modTimeStep)
	//{	
	//	harvAreaO = res.second.area;
	//	harvestRO = 0.; 
	//	if (Area>0) {
	//		restW = res.second.rw * harvAreaO/Area;       
	//		restThW = res.first.rw/Area;    
	//		harvestRO = (restW + restThW)/modTimeStep;
	//	}		
	//	return(harvestRO) ;
	//}
	//

