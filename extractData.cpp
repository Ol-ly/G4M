#include <iostream>
using namespace std;


int extractData (g4m::ageStruct &cohort,dat &singleCell, double siteIndex, double X, double Y)
{

	const int size1=23;
	double thinning [size1]={0.};
	double forestarea[size1]={0.};
	double increment[size1]={0.};
	double diameter [size1]={0.};
	double biomass [size1]={0.};
    double growingStockTab[size1] = {0.};
	const int size= 221;
	double age [size]={0.};
	
	double currentAge =0.;
	int period = 0;


	for (int iter=0; iter<=220; iter+=1){
		age[iter]=age[iter]+iter;		
		currentAge=age[iter];	
		cohort.TEST(currentAge, siteIndex);
		forestarea[period] = forestarea[period] + cohort.ClassArea;			// ADD *(singleCell.LandAreaHa*OforestShare)/0.993351;
		thinning[period] = thinning[period] + cohort.R_ha;
		increment[period] = increment[period] + cohort.IncTot_ha;		 
		biomass[iter]=cohort.GS_t_Tab_ha;

		if((X==25.75)&&(Y==49.75)){
			cout<<iter<<"	biomass "<<biomass[iter]<<endl;}
		if ( (currentAge+1)==((period+1)*10) ){
			forestarea[period]=forestarea[period]*(singleCell.LandAreaHa*singleCell.OforestShare)/0.993351;;
			period++;
			if (period==23){
				forestarea[23] = forestarea[22];
				thinning[23] = thinning[22];
				increment[23] = increment[22];
				}
			if (period==24){				
				break;}
		}

		currentAge=0.;
	}

	for(int iter=0;iter<=220;iter+=10){

		cohort.TEST(iter, siteIndex);
		growingStockTab[iter] = cohort.GS_t_Tab_ha;
		//diameter[iter] = cohort.D_t_ha;
		//cout << iter<<"\t"<<diameter[iter]<<endl;
	//outfile4<<X<<"\t"<<Y<<"\t"<<iter<<"\t"<<growingStockTab[iter]<<endl;
	//"\t"<<diameter[iter]<<
	if (iter==220){break;}

	}

	
	for (int k=0;k<=22;k++){
		 
	outfile3<<X<<"\t"<<Y<<"\t"<<k<<"\t"<<forestarea[k]<<"\t"<<thinning[k]<<"\t"<<increment[k]<<endl;
	//"\t"<<growingStockTab[k*10]<<
	//	if (k==22){break;}
	}

	return 0;
}