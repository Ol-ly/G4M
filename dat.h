#ifndef DAT_H_
#define DAT_H_

struct dat {
  int simUnit;
  int Rotation;
  double LandAreaHa;
  double potHarvest;
  double forestShare;
  double OforestShare;
  double AforestShare;
  double prevOForShare;
  double prevOForShareRP;
  double prevAForShareRP;
  double AforestSharePrev;
  double savedCarbonPrev;
  double gainedCarbonPrev;
  double EmissionsTotPrev;
  double EmissionsAfforPrev;
  double prevPlantPhytHaBmGr;
  double prevPlantPhytHaBmGrBef;
  double prevPlantPhytHaBmGrBef20;
  double prevPlantPhytHaBlGr;
  double deforestHaTot;
  double afforestHaTot;
  double EmissionsProduct;  
  double EmissionsLitter;  
  double EmissionsSOC;      
  double EmissionsSlashBurn;
  double EmissionsDeadBurn;
  double EmissionsCRootBurn;    
  double EmissionsTot;     
  double EmLitterAffor;
  double EmSOCAffor; 
  double EmissionsAffor;
  double forestAgeShare[110];
  double BDeadA[110];
  double LitterA[110];
  double SOCA[110];
  double ProdLongA[110];
  double ProdShortA[110];
  double deforestA[110];  
  double FineRootA[110];
  double LitterAffor[110];
  double SOCaffor[110];
  double prevReportYear;
  double ObiomassPrev;
  int rotBiomass;  
  int ireportYear;
  double SD;  
  double usedForest;
  double FMsink;
  double FMsink_Bm;
  double Obiomass0; // modelled biomass at time 0
  pair <double, double> DeforWoodM3;	// total deforested sawn and rest wood
  double deforWoodTotM3; // stem wood obtained from deforestation averaged for last 5 years
  double deforWoodArray[5]; //stem wood obtained from deforestation for last 5 years
  double deforAreaArray[5]; // deforested area for the last 5 years
  double deforestShare;
  double afforestShare;  
  double CAI;
  double afforestHaYear;
  double deforestHaYear;
  double harvSawnWood;
  double harvRestWood;
  double harvestTot;
  double harvestFcM3Ha;
  double harvestThM3Ha;
  double oforestBm;
  double aforestBm;
  double emissionsD_Bm;  
  double emissionsD_S;
  double emissionsA_Bm;  
  double emissionsA_S;
  double R0; // discount rate at time 0
  double kR; // a coefficient for modification the discount rate
  double deforPrev;
  double road;
  double deforRateCoeffCell;
  };

#endif
