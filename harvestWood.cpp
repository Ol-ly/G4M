double cohortRes(double realAreaO, pair<g4m::ageStruct::v, g4m::ageStruct::v> &res)
{
		double harvestW = 0.;		
		double harvAreaO = res.second.area;
        if (realAreaO>0) {
			double sawnW = res.second.sw * harvAreaO/realAreaO;      // MG: get harvestable sawnwood for the set (old) forest tC/ha for final cut. 
			double restW = res.second.rw * harvAreaO/realAreaO;      // MG: get harvestable restwood for the set (old) forest tC/ha for final cut. 
			double sawnThW = res.first.sw/realAreaO;    // MG: get harvestable sawnwood for the set (old) forest tC/ha for thinning. 
			double restThW = res.first.rw/realAreaO;    // MG: get harvestable restwood for the set (old) forest tC/ha for thinning.
			double bmH = res.second.bm * harvAreaO/realAreaO;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor final cut
			double bmTh = res.first.bm/realAreaO;       // MG: get total harvestable biomass including harvest losses for the set (old) forest tC/hafor thinning          
			//harvL = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) ; // MG: harvest residues for the set (old) forest tC/ha
		    double harvRes = (bmH+bmTh - (sawnW + restW + sawnThW + restThW)) * resUse; // MG: usable harvest residues for the set (old) forest tC/ha
		    harvestW = (sawnW + restW + sawnThW + restThW + harvRes)/modTimeStep;
		}

		return(harvestW);
}