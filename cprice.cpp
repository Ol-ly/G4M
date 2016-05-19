// Carbon values ETS, Prices in EURO 2008.

#ifndef CPRICE_CPP
#define CPRICE_CPP

#include <map>
#include "misc.h"
//#include "ipol.h"

void carbonPrice(void)
{
//bau - BAU scenario
string colname = "bau";
cprice[colname].insert(2010,0.);
cprice[colname].insert(2015,0.);
cprice[colname].insert(2020,0.);
cprice[colname].insert(2025,0.);
cprice[colname].insert(2030,0.);
cprice[colname].insert(2035,0.);
cprice[colname].insert(2040,0.);
cprice[colname].insert(2045,0.);
cprice[colname].insert(2050,0.);

//Ref - Reference scenario
colname = "ref";
cprice[colname].insert(2010,11.1*  exchangeRate);
cprice[colname].insert(2015,13.6*  exchangeRate);
cprice[colname].insert(2020,16.5*  exchangeRate);
cprice[colname].insert(2025,20.0*  exchangeRate);
cprice[colname].insert(2030,36.0*  exchangeRate);
cprice[colname].insert(2035,50.0*  exchangeRate);
cprice[colname].insert(2040,51.5*  exchangeRate);
cprice[colname].insert(2045,50.5*  exchangeRate);
cprice[colname].insert(2050,50.0*  exchangeRate);

//Ref-SH - Ref + Oil shock
colname = "refsh";
cprice[colname].insert(2010,11.1* exchangeRate);
cprice[colname].insert(2015,13.6* exchangeRate);
cprice[colname].insert(2020,16.5*  exchangeRate);
cprice[colname].insert(2025,20.0*  exchangeRate);
cprice[colname].insert(2030,36.0*  exchangeRate);
cprice[colname].insert(2035,43.0*  exchangeRate);
cprice[colname].insert(2040,45.0*  exchangeRate);
cprice[colname].insert(2045,45.0*  exchangeRate);
cprice[colname].insert(2050,45.3*  exchangeRate);

//Decarb-ETGCA (effective technologies)
colname = "etgca";
cprice[colname].insert(2010,11.1*  exchangeRate);
cprice[colname].insert(2015,15.0*  exchangeRate);
cprice[colname].insert(2020,25.0*  exchangeRate);
cprice[colname].insert(2025,50.0*  exchangeRate);
cprice[colname].insert(2030,64.0*  exchangeRate);
cprice[colname].insert(2035,127.0*  exchangeRate);
cprice[colname].insert(2040,170.0*  exchangeRate);
cprice[colname].insert(2045,195.0*  exchangeRate);
cprice[colname].insert(2050,240.0*  exchangeRate);


//Decarb-DGCA (delayed climate action)
colname = "dgca";
cprice[colname].insert(2010,11.1*  exchangeRate);
cprice[colname].insert(2015,13.6*  exchangeRate);
cprice[colname].insert(2020,16.5*  exchangeRate);
cprice[colname].insert(2025,20.0*  exchangeRate);
cprice[colname].insert(2030,36.0*  exchangeRate);
cprice[colname].insert(2035,200.0*  exchangeRate);
cprice[colname].insert(2040,406.4*  exchangeRate);
cprice[colname].insert(2045,515.7*  exchangeRate);
cprice[colname].insert(2050,540.0*  exchangeRate);
}
#endif
