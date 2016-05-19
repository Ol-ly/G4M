#ifndef G4M_AGESTRUCT_CPP
#define G4M_AGESTRUCT_CPP

#include "ageStruct_float.h"

#include <iostream>

namespace g4m {
  ageStruct::ageStruct
  (incrementTab *ait
   , ffipol<double> *asws
   , ffipol<double> *ahlv, ffipol<double> *ahle
   , ffipolm<double> *acov, ffipolm<double> *acoe
   , ffipolm<bool> *adov, ffipolm<bool> *adoe
   , double amai
   , int aobjOfProd, double au
   , double aminSw, double aminRw, double aminHarv
   , int asdDef, double asdMax, double asdMin
   , unsigned int amaiYears
   , double aminRotVal, int aminRotRef
   , double aflexSd
   , ffipol<double> *sdMaxH, ffipol<double> *sdMinH
   , unsigned int amaxAge)
  {
    it = ait;
    sws = asws;
    hlv = ahlv; hle = ahle;
    cov = acov; coe = acoe; dov = adov; doe = adoe;
    mai = amai;
    if(!qMai.empty()) {qMai.clear();}
    qMai.resize(amaiYears, mai);
    avgMai = mai;
    objOfProd = aobjOfProd;
    uRef = au;
    u = setRotationTime();
    minSw = aminSw; minRw = aminRw; minHarv = aminHarv;
    sdDef = asdDef;
    sdMin = asdMin; sdMax = asdMax;
    minRotVal = aminRotVal;
    minRotRef = aminRotRef;
    setMinRot();
    flexSd=aflexSd;
    timeStep = it->gtimeframe();  //Check if this should be tStep or timeframe
    area = 0.;
    maxNumberOfAgeClasses = ceil(amaxAge/timeStep);
  }

  ageStruct::~ageStruct() {
    it = NULL;
    sws = NULL;
    hlv = NULL;
    hle = NULL;
    cov = NULL;
    coe = NULL;
    dov = NULL;
    doe = NULL;
  }

  double ageStruct::calcAvgMai() {
    if(qMai.size() > 0) {
      double product = 1.;
      double sWeight = 0.;
      double weight = 0.;
      double step = 2./qMai.size();
      std::deque<double>::iterator iter = qMai.begin();
      while(iter != qMai.end()) {
        weight += step;
        if(*iter <= 0.) {   //The years before 0 are unimportant
          product = 1.;
          sWeight = 0.;
        } else {
          product *= std::pow(*iter, weight);
          sWeight += weight;
        }
        ++iter;
      }
      if(sWeight > 0.) {avgMai = std::pow(product, 1./sWeight);}
      else {avgMai = 0.;}
    } else {
      avgMai = mai;
    }
    return(avgMai);
  }
 
  double ageStruct::setRotationTime() {
    if(objOfProd == 0) {u=uRef;}
    else if(objOfProd == 1 || objOfProd == 2) {u=-1.;}
    else if(objOfProd > 2 && objOfProd < 8) {
      const int type = objOfProd - 3;
      double sd = (sdMin + sdMax)/2.;
      if(sd > 0) {
        if(sd == 1.) {u = it->gToptt(avgMai, type);}
        else {u = it->gToptSdTab(avgMai, sd, type);}
      } else {
        if(sd == -1.) {u = it->gTopt(avgMai, type);}
        else {u = it->gToptSdNat(avgMai, -1. * sd , type);}
      }
    }
    else {u=0.;}
    return(u);
  }

  double ageStruct::setMinRot() {
    if(minRotRef == 0) {minRot = minRotVal;}
    else if(minRotRef == 1 && u > 0.) {minRot = minRotVal * u;}
    else if(minRotRef < 7) {
      const int type = minRotRef - 2;
      double sd = (sdMin + sdMax)/2.;
      if(sd > 0) {
        if(sd == 1.) {minRot = minRotVal * it->gToptt(avgMai, type);}
        else {minRot = minRotVal * it->gToptSdTab(avgMai, sd, type);}
      } else {
        if(sd == -1.) {minRot = minRotVal * it->gTopt(avgMai, type);}
        else {minRot = minRotVal * it->gToptSdNat(avgMai, -sd , type);}
      }
    } else {minRot = -1.;}
    return(minRot);
  }

  double ageStruct::createNormalForest
  (double rotationPeriod, double aarea, double sd) {
    if(aarea < 0.) {aarea = 0.;}
    area = aarea;
    if(rotationPeriod < 1.) {rotationPeriod = 1;}
    if(rotationPeriod >= it->gTmax()) {rotationPeriod = it->gTmax()-1;}
    unsigned int ageClasses = std::ceil(0.5 + rotationPeriod/timeStep);
    if(dat.size() < ageClasses) {dat.resize(ageClasses);}
    if(maxNumberOfAgeClasses < ageClasses) {maxNumberOfAgeClasses = ageClasses;}
    aarea /= rotationPeriod;
    cohort tmp;
    tmp.area = aarea*timeStep;
    tmp.bm = 0.; tmp.d = 0.; tmp.h = 0.;
    for(unsigned int i=0; i<ageClasses; ++i) {
      double age = i*timeStep;
      if(sd > 0.) {
        if(sd == 1.) {
          tmp.bm = it->gBmt(age, avgMai);
          tmp.d = it->gDbht(age, avgMai);
        }
        else {
          tmp.bm = it->gBmSdTab(age, avgMai, sd);
          tmp.d = it->gDbhSdTab(age, avgMai, sd);
        }
      } else {
        if(sd == -1.) {
          tmp.bm = it->gBm(age, avgMai);
          tmp.d = it->gDbh(age, avgMai);
        }
        else {
          tmp.bm = it->gBmSdNat(age, avgMai, -sd);
          tmp.d = it->gDbhSdNat(age, avgMai, -sd);
        }
      }
      tmp.h =  it->gHeight(age, avgMai);
      dat[i] = tmp;
    }
    dat[0].area /= 2.;
    double share = 0.5 + rotationPeriod/timeStep - static_cast<int>(0.5 + rotationPeriod/timeStep);
    if(share > 0.) {dat[ageClasses-1].area = aarea * timeStep * share;}
    for(unsigned int i=ageClasses; i<dat.size(); ++i) {
      dat[i].area = 0.;
    }
    return(area);
  }

  double ageStruct::getD(double age) {
    age /= timeStep;
    unsigned int ageH = ceil(age);
    if(ageH >= dat.size()) {return(dat[dat.size()-1].d);}
    if(ageH < 1) {return(dat[0].d);}
    --ageH;
    return(dat[ageH].d + (dat[ageH+1].d - dat[ageH].d) * (age - ageH));
  }

  double ageStruct::getH(double age) {
    age /= timeStep;
    unsigned int ageH = ceil(age);
    if(ageH >= dat.size()) {return(dat[dat.size()-1].h);}
    if(ageH < 1) {return(dat[0].h);}
    --ageH;
    return(dat[ageH].h + (dat[ageH+1].h - dat[ageH].h) * (age - ageH));
  }

  double ageStruct::getBm(double age) {
    age /= timeStep;
    unsigned int ageH = ceil(age);
    if(ageH >= dat.size()) {return(dat[dat.size()-1].bm);}
    if(ageH < 1) {return(dat[0].bm);}
    --ageH;
    return(dat[ageH].bm + (dat[ageH+1].bm - dat[ageH].bm) * (age - ageH));
  }

  double ageStruct::getBm() {
    double area=0.;
    double bm=0.;
    for(unsigned int i=0; i<dat.size(); ++i) {
      if(dat[i].area==dat[i].area && dat[i].bm==dat[i].bm) { //Sould be Improved
        area += dat[i].area;
        bm += dat[i].area * dat[i].bm;
      }
    }
    if(area > 0.) {return(bm/area);}
    return(0.);
  }

  double ageStruct::getArea(int age) {
    if(age >= static_cast<int>(dat.size())) {return(0.);}
    return(dat[age].area);
  }

  double ageStruct::getArea(double age) {
    //return(dat[0].area);
    //Here seems to be an error - Above just to allow some colculations
    age /= timeStep;
    unsigned int ageH = static_cast<int>(age+0.5);
    if(ageH >= dat.size()) {return(0.);}
    if(ageH <= 0) {return(dat[0].area/timeStep);}
    return(dat[ageH].area/timeStep);
  }

  double ageStruct::getArea() {
    return(area);
  }

  double ageStruct::calcArea() { //Calculates the area
    area = 0.;
    for(unsigned int i=0; i<dat.size(); ++i) {
      area += dat[i].area;
    }
    return(area);
  }

  double ageStruct::setArea(unsigned int ageClass, double aarea) {
    //area += aarea - dat[0].area;
    //dat[0].area = aarea;
    //return(dat[0].area);
    //Here seems to be an error - Above just to allow some colculations
    if(aarea < 0.) {return(-1.);}
    if(ageClass >= static_cast<unsigned int>(it->gTmax()/timeStep)) {
      ageClass = static_cast<unsigned int>(it->gTmax()/timeStep)-1;
    }
    if(ageClass < dat.size()) {
      area += aarea - dat[ageClass].area;
      dat[ageClass].area = aarea;
      return(dat[ageClass].area);
    }
    area += aarea;
    unsigned int oldSize = dat.size();
    dat.resize(ageClass+1);
    initCohort(oldSize, ageClass);
    dat[ageClass].area = aarea;
    return(dat[ageClass].area);
  }

  double ageStruct::setBm(unsigned int ageClass, double biomass) {
    if(biomass < 0.) {return(0.);}
    if(ageClass >= static_cast<unsigned int>(it->gTmax()/timeStep)) {
      ageClass = static_cast<unsigned int>(it->gTmax()/timeStep)-1;
    }
    double maxBm = it->gBm(ageClass*timeStep, avgMai);
    if(biomass > maxBm) {biomass = maxBm;}
    if(ageClass < dat.size()) {
      dat[ageClass].bm = biomass;
      return(dat[ageClass].bm);
    }
    unsigned int oldSize = dat.size();
    dat.resize(ageClass+1);
    initCohort(oldSize, ageClass);
    dat[ageClass].bm = biomass;
    return(dat[ageClass].bm);
  }

  double ageStruct::setD(unsigned int ageClass, double dbh) {
    if(ageClass >= static_cast<unsigned int>(it->gTmax()/timeStep)) {
      ageClass = static_cast<unsigned int>(it->gTmax()/timeStep)-1;
    }
    double age = ageClass*timeStep;
    double minD = it->gDbh(age, avgMai);
    double maxD = it->gDbhSdNat(age, avgMai, 0.);
    if(dbh < minD) {dbh = minD;}
    if(dbh > maxD) {dbh = maxD;}
    if(ageClass < dat.size()) {
      dat[ageClass].d = dbh;
      return(dat[ageClass].d);
    }
    unsigned int oldSize = dat.size();
    dat.resize(ageClass+1);
    initCohort(oldSize, ageClass);
    dat[ageClass].d = dbh;
    return(dat[ageClass].d);
  }

  unsigned int ageStruct::initCohort(unsigned int ageClassL, unsigned int ageClassH) {
    cohort tmp;
    tmp.area = 0.;
    double sd = (sdMin + sdMax)/2.;
    for(unsigned int i = ageClassL; i < ageClassH; ++i) {
      double age = i*timeStep;
      if(sd > 0.) {
        if(sd == 1.) {
          tmp.bm = it->gBmt(age, avgMai);
          tmp.d = it->gDbh(age, avgMai);
        }
        else {
          tmp.bm = it->gBmSdTab(age, avgMai, sd);
          tmp.d = it->gDbhSdTab(age, avgMai, sd);
        }
      } else {
        if(sd == -1.) {
          tmp.bm = it->gBm(age, avgMai);
          tmp.d = it->gDbh(age, avgMai);
        }
        else {
          tmp.bm = it->gBmSdNat(age, avgMai, -sd);
          tmp.d = it->gDbhSdNat(age, avgMai, -sd);
        }
      }
      tmp.h =  it->gHeight(age, avgMai);
      dat[i] = tmp;
    }
    return(ageClassH - ageClassL);
  }

  double ageStruct::setMai(double aMai) {
    mai = aMai;
    return(mai);
  }

  unsigned int ageStruct::setMaiYears(unsigned int maiYears) {
    while(maiYears > qMai.size()) {qMai.push_front(avgMai);}
    while(maiYears < qMai.size()) {qMai.pop_front();}
    calcAvgMai();
    setRotationTime();
    setMinRot();
    return(qMai.size());
  }

  double ageStruct::setAvgMai(double aavgMai) {
    avgMai = aavgMai;
    if(avgMai < 0.) {avgMai = 0.;}
    std::deque<double>::iterator iter = qMai.begin();
    while(iter != qMai.end()) {
      *iter = avgMai;
      ++iter;
    }
    setRotationTime();
    setMinRot();
    return(avgMai);
  }

  double ageStruct::getMai() {
    return(mai);
  }

  double ageStruct::getAvgMai() {
    return(avgMai);
  }

  int ageStruct::setObjOfProd(int aobjOfProd) {
    objOfProd = aobjOfProd;
    return(objOfProd);
  }

  double ageStruct::setU(double au) {
    uRef = au;
    setRotationTime();
    setMinRot();
    return(u);
  }

  double ageStruct::setStockingdegreeMin(double sd) {
    sdMin=sd;
    return(sdMin);
  }

  double ageStruct::setStockingdegreeMax(double sd) {
    sdMax=sd;
    return(sdMax);
  }

  int ageStruct::setMinRotRef(int aminRotRef) {
    minRotRef = aminRotRef;
    return(minRotRef);
  }

  double ageStruct::setMinRotVal(double aminRotVal) {
    minRotVal = aminRotVal;
    setMinRot();
    return(minRotVal);
  }

  double ageStruct::setFlexSd(double aflexSd) {
    flexSd=aflexSd;
    return(flexSd);
  }

  double ageStruct::afforest(double aarea) {
//    if(area > 0) { //MG:
    if(aarea > 0) {
      dat[0].area += aarea;
      area += aarea;
    }
    return(dat[0].area);
  }
/*
  double ageStruct::reforest(double aarea) {
    if(dat.size() < 2) {
      int oldSize=dat.size();
      dat.resize(2);
      initCohort(oldSize, 2);
    }
    if(area > 0) {
      dat[0].area += aarea/2.;
      dat[1].area += aarea/2.;
      area += aarea;
    }
    return(dat[0].area + dat[1].area);
  }
  */

    double ageStruct::reforest(double aarea) {
	double reforestArea=0;
    if(dat.size() < 2) {
      int oldSize=dat.size();
      dat.resize(2);
      initCohort(oldSize, 2);
    }
    if(aarea > 0 && dat[1].bm>0) {
      dat[0].area += aarea/2.;
      dat[1].area += aarea/2.;
	  reforestArea=dat[0].area + dat[1].area;
      area += aarea;
    }else{
		dat[0].area += aarea;
		area += aarea;
		reforestArea=dat[0].area;
	}

    return(reforestArea);
  }

  ageStruct::v ageStruct::deforest(double aarea, int type) {
    v ret = {0., 0., 0., 0., 0.};
    if(type == 0) { //Take from all age classes
      if(area < aarea) {aarea = area;}
      if(area > 0.) {
        double mul = 1. - aarea/area;
        static std::vector<double> dbhBm(2,0); //Key of dbh and biomass
        for(unsigned int i=0; i<dat.size(); ++i) {
          if(dat[i].area > 0.) {
            double sdNat = it->gSdNat(i*timeStep, avgMai, dat[i].bm);
            double dbm = it->gIncBmSdNat(i*timeStep, avgMai, sdNat)/2.;
            double id = it->gIncDbhSdNat(i*timeStep, avgMai, sdNat)/2.;
//if (dbm<0 || id<0 || dbm>100 || id > 100) {cout<<"agestruct:"<<"\tcellid="<<cellid<<"\tdbm/id 1:"<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<"\tavgMai="<<avgMai<<"\ti*timeStep="<<i*timeStep<<"\tdat.size()=\t"<<dat.size()<<endl;} //MG: test
if (dbm>100 || id > 100) {cout<<"agestruct:"<<"\tcellid="<<cellid<<"\tdbm/id 1:"<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<"\tavgMai="<<avgMai<<"\ti*timeStep="<<i*timeStep<<"\tdat.size()=\t"<<dat.size()<<endl;} //MG: test
//if (dbm<0) {dbm=0;} if (dbm>dat[i].bm*(1./sqrt(i+1.))) {dbm=dat[i].bm*(1./sqrt(i+1.));}// MG:: temporal fix. Must correct the tables!
//if (id<0) {id=0;} if (id>dat[i].d*(1./sqrt(i+1.))) {id=dat[i].d*(1./sqrt(i+1.));}// MG:: temporal fix. Must correct the tables!
if (dbm<0) {dbm=0;} //MG: we don't allow negative increments (suggestion of Geogr 25 March 2013) 
if (id<0) {id=0;} //MG: we don't allow negative increments (suggestion of Geogr 25 March 2013) 
if (dat[i].bm<0) dat[i].bm=0; // MG: we don't allow negative biomass in any age group

            dbhBm[0] = dat[i].d+id; dbhBm[1] = dat[i].bm+dbm;
            double totalWood = dat[i].area * (1.-mul) * dbhBm[1];
            ret.bm += totalWood;
//            ret.bm += totalWood * 10000; MG: test
            double harvestedWood = totalWood * hle->g(dbhBm[0]);
            double sawnWood = harvestedWood * sws->g(dbhBm[0]);
            ret.sw += sawnWood;
            ret.rw += harvestedWood - sawnWood;
            ret.co += totalWood * coe->g(dbhBm);
            area -= dat[i].area * (1. - mul);
            dat[i].area *= mul;
          }
        }
        ret.area=aarea;
      }
      if(ret.area > 0.) { //Values per hectare
        ret.sw /= area; ret.rw /= area; ret.co /= area; ret.bm /= area;}
    } else { //Take from the old age classes
      ret =  finalCut(aarea, false);
    }
    return(ret);
  }

  ageStruct::v ageStruct::finalCut(bool eco, double aarea, double minSw, double minRw, double minHarv, bool sustainable) {
    v ret = {0., 0., 0., 0., 0.};
    int endYear = 0;
    if(eco || sustainable) {endYear = minRot/timeStep;}
    static std::vector<double> dbhBm(2,0); //Key to ask if harvest is economic
    for(int i=dat.size()-1; i>=endYear && (ret.area<aarea || (ret.sw<minSw || ret.rw<minRw || (ret.sw+ret.rw)<minHarv)); --i) {
      if(dat[i].area > 0.) {
        //The Stands get half increment of the next growing period
        double sdNat = it->gSdNat(i*timeStep, avgMai, dat[i].bm);
        double dbm = it->gIncBmSdNat(i*timeStep, avgMai, sdNat)/2.;
        double id = it->gIncDbhSdNat(i*timeStep, avgMai, sdNat)/2.;
if (dat[i].bm<0 || dat[i].bm>400 || dbm>30) {cout<<"agestruct:"<<"\tcellid="<<cellid<<"\tdbm/id 2:"<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<"\tavgMai="<<avgMai<<"\ti*timeStep="<<i*timeStep<<"\tdat.size()=\t"<<dat.size()<<endl;} //MG: test
//if (dbm<0 || id<0 || dbm>100 || id > 100) {cout<<"agestruct:"<<"\tcellid="<<cellid<<"\tdbm/id 3:"<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<"\tavgMai="<<avgMai<<"\ti*timeStep="<<i*timeStep<<"\tdat.size()=\t"<<dat.size()<<endl;} //MG: test
//if (dbm>100 || id > 100) {cout<<"agestruct:"<<"\tcellid="<<cellid<<"\tdbm/id 3:"<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<"\tavgMai="<<avgMai<<"\ti*timeStep="<<i*timeStep<<"\tdat.size()=\t"<<dat.size()<<endl;} //MG: test
//if (dbm<0) {dbm=0;} if (dbm>dat[i].bm*(1./sqrt(i+1.))) {dbm=dat[i].bm*(1./sqrt(i+1.));}// MG:: temporal fix. Must correct the tables!
//if (id<0) {id=0;} if (id>dat[i].d*(1./sqrt(i+1.))) {id=dat[i].d*(1./sqrt(i+1.));}// MG:: temporal fix. Must correct the tables!
if (dbm<0) {dbm=0;} //MG: we don't allow negative increments (suggestion of Geogr 25 March 2013) 
if (id<0) {id=0;} //MG: we don't allow negative increments (suggestion of Geogr 25 March 2013) 
if (dat[i].bm<0) dat[i].bm=0; // MG: we don't allow negative biomass in any age group
        dbhBm[0] = dat[i].d+id; dbhBm[1] = dat[i].bm+dbm;
        if(doe->g(dbhBm) || eco==false) { //do harvest if it is economic
          if(aarea >=0) { //Given area to harvest
            double totalWood = 0.;
            if(ret.area+dat[i].area < aarea) { //Harvest all of this age class
              totalWood = dat[i].area * dbhBm[1];
			   
              ret.area += dat[i].area;
			  //double harvArea[dat.size()-1]=dat[i].area;
              area -= dat[i].area;
              dat[i].area = 0.;
            } else {
              totalWood = (aarea - ret.area) * dbhBm[1];
              area -= aarea - ret.area;
              dat[i].area -= aarea - ret.area;
              ret.area = aarea;
            }
            ret.bm += totalWood;
            double harvestedWood = totalWood * hle->g(dbhBm[0]);
            double sawnWood = harvestedWood * sws->g(dbhBm[0]);
//cout<<"agestruct: sawnWood=\t"<<sawnWood<<"\tharvestedWood=\t"<<harvestedWood<<"\ttotalWood=\t"<<totalWood<<"\thle=\t"<<hle->g(dbhBm[0])<<"\tsws=\t"<<sws->g(dbhBm[0])<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<endl;
if (harvestedWood<0 || harvestedWood>200 || sawnWood<0 || sawnWood>200 || dat[i].bm<0 || dat[i].bm>400) cout<<"agestruct:"<<"\tcellid="<<cellid<<"\tsawnWood=\t"<<sawnWood<<"\tharvestedWood=\t"<<harvestedWood<<"\ttotalWood=\t"<<totalWood<<"\thle=\t"<<hle->g(dbhBm[0])<<"\tsws=\t"<<sws->g(dbhBm[0])<<"\tsdNat=\t"<<sdNat<<"\tdbm=\t"<<dbm<<"\tid=\t"<<id<<"\tdat[i].bm=\t"<<dat[i].bm<<"\tharea=\t"<<aarea - ret.area<<endl;
            ret.sw += sawnWood;
            ret.rw += harvestedWood - sawnWood;
            ret.co += totalWood * coe->g(dbhBm);
          } else { //given amount of harvest
            double harvestShare = 1.;
            double totalWood = dat[i].area * dbhBm[1];
            double harvestedWood = totalWood * hle->g(dbhBm[0]);
            double sawnWood = harvestedWood * sws->g(dbhBm[0]);
            if((ret.sw+sawnWood) < minSw || (ret.rw+harvestedWood-sawnWood) < minRw || (ret.sw+ret.rw+harvestedWood) < minHarv) {//Harvest all
              ret.area += dat[i].area;
              area -= dat[i].area;
              dat[i].area = 0.;
            } else {//Harvest part of ageclass
              harvestShare = 0.;
              double tmp = minSw - ret.sw;
              if(tmp > 0. && sawnWood>0.) {tmp /= sawnWood;
                if(harvestShare < tmp) {harvestShare = tmp;}}
              tmp = minRw - ret.rw;
              if(tmp > 0. && (harvestedWood-sawnWood)>0.) {tmp /= harvestedWood-sawnWood;
                if(harvestShare < tmp) {harvestShare = tmp;}}
              tmp = minHarv - (ret.sw+ret.rw);
              if(tmp > 0. && harvestedWood>0.) {tmp /= harvestedWood;
                if(harvestShare < tmp) {harvestShare = tmp;}}
              if(harvestShare < 0.) {harvestShare = 0.;}
              if(harvestShare > 1.) {harvestShare = 1.;}
              ret.area += dat[i].area * harvestShare;
              area -= dat[i].area * harvestShare;
              dat[i].area *= 1. - harvestShare;
            }
            ret.bm += totalWood * harvestShare;
            ret.sw += sawnWood * harvestShare;
            ret.rw += (harvestedWood - sawnWood) * harvestShare;
            ret.co += (totalWood * coe->g(dbhBm)) * harvestShare;
          }
        }
      }
    }
    if(ret.area > 0.) { //Values per hectare
      ret.sw /= ret.area; ret.rw /= ret.area; ret.co /= ret.area; ret.bm /= ret.area;}
    return(ret);
  }

  ageStruct::v ageStruct::finalCut(double aarea, bool eco) {
    return(finalCut(eco, aarea, -1., -1., -1.));
  }

  ageStruct::v ageStruct::finalCut(double minSw, double minRw, double minHarv
                                   , bool eco, bool sustainable) {
    return(finalCut(eco, -1., minSw, minRw, minHarv,sustainable));
  }

  std::pair<ageStruct::v, ageStruct::v> ageStruct::aging() {
    return(aging(mai));
  }

  std::pair<ageStruct::v, ageStruct::v> ageStruct::aging(double amai) {
    v retThin = {0., 0., 0., 0., 0.};
    v retHarvest = {0., 0., 0., 0., 0.};
    mai = amai;
    qMai.push_back(mai);
    qMai.pop_front();
    calcAvgMai();
    setRotationTime();
    setMinRot();
    if(objOfProd == 1 || objOfProd == 2) { //Fulfill an amount of harvest
      retThin = thinAndGrow();
      if(objOfProd == 1) {retHarvest = finalCut(minSw-retThin.sw, minRw-retThin.rw, minHarv-(retThin.sw+retThin.rw), true, false);
      } else {retHarvest = finalCut(minSw-retThin.sw, minRw-retThin.rw, minHarv-(retThin.sw+retThin.rw), true, true);}
    } else { //We have a rotation time to fulfill
      retHarvest = finalCut(area*timeStep/u, true); //do final cut
      retThin = thinAndGrow();
      //retThin = thinAndGrowOLD();
    }
    //Make reforestations on final harvested area
    reforest(retHarvest.area); //Maybe include also reforestation costs
    return(std::make_pair(retThin, retHarvest));
  }

  ageStruct::v ageStruct::thinAndGrow() {
    if(sdDef==0) { //Constant stocking degre
      return(thinAndGrowStatic());
    } else if(sdDef==1) { //Alternate between constant stocking degrees
      return(thinAndGrowOLD());
    }
    //Alternate between varying stocking degrees
    return(thinAndGrowOLD());
  }

  ageStruct::v ageStruct::thinAndGrowStatic() {
    v ret = {0., 0., 0., 0., 0.};
    for(int i = static_cast<int>(dat.size())-2; i>-1; --i) {
      if(dat[i].area > 0.) {
        double sd, iGwl, bmT, id;
        incStatic(i, sd, iGwl, bmT, id);
        if(flexSd > 0.) { //The typical amount of harvest
          double bmTCom =  dat[i].bm + incCommon(i, sd, iGwl);
          bmT = bmT * (1. - flexSd) + bmTCom * flexSd;
        }
        double totalWood = dat[i].area * (iGwl - (bmT-dat[i].bm));
        if(totalWood < 0.) {totalWood = 0.;}
        static std::vector<double> dbhBmSh(3,0);
        dbhBmSh[0] = dat[i].d+id/2.;
        dbhBmSh[1] = dat[i].bm+iGwl/2.;
        dbhBmSh[2] = totalWood/dbhBmSh[1];
        if(dov->g(dbhBmSh)) { //Do Thinning if it is economic
          double harvestedWood = totalWood * hlv->g(dbhBmSh[0]);
          double sawnWood = harvestedWood * sws->g(dbhBmSh[0]);
          ret.area += dat[i].area;
          ret.bm += totalWood;
          ret.sw += sawnWood;
          ret.rw += harvestedWood - sawnWood;
          ret.co += totalWood * cov->g(dbhBmSh);
          dat[i].bm += iGwl - totalWood/dat[i].area;
        } else {  //No thinning
          dat[i].bm += iGwl;
          double bmMax = it->gBm((i+1)*timeStep, avgMai);
          if(dat[i].bm > bmMax) {dat[i].bm = bmMax;}
        }
        dat[i].d += id;
        dat[i].h += it->gIncHeight(i*timeStep, avgMai);
      }
    }
    cohortShift();
    return(ret);
  }

  int ageStruct::incStatic(int i, double &sd, double &iGwl, double &bm, double &id) {
    if(sdMax > 0) {  //Yield table stocking degree
      if(sdMax == 1.) {bm = it->gBmt((i+1)*timeStep, mai);}
      else {bm = it->gBmSdTab((i+1)*timeStep, mai, sdMax);}
      if(i) {sd = it->gSdTab(i*timeStep, avgMai, dat[i].bm);}
      else {sd = 1.;} //New plantations have a stocking degree of 1
      if(sd == 1.) {
        iGwl = it->gIncGwlt(i*timeStep, mai);
        id = it->gIncDbht(i*timeStep, mai);
      } else {
        iGwl = it->gIncGwlSdTab(i*timeStep, mai, sd);
        id = it->gIncDbhSdTab(i*timeStep, mai, sd);
      }
    } else {  //Natural stocking degree
      if(sdMax == -1.) {bm = it->gBm((i+1)*timeStep, mai);}
      else {bm = it->gBmSdNat((i+1)*timeStep, mai, -sdMax);}
      if(i) {sd = it->gSdNat(i*timeStep, avgMai, dat[i].bm);}
      else {sd = 1.;} //New plantations have a stocking degree of 1
      if(sd == 1.) {
        iGwl = it->gIncGwl(i*timeStep, mai);
        id = it->gIncDbh(i*timeStep, mai);
      } else {
        iGwl = it->gIncGwlSdNat(i*timeStep, mai, sd);
        id = it->gIncDbhSdNat(i*timeStep, mai, sd);
      }
    }
    return(0);
  }

  double ageStruct::incCommon(const int i, const double &sd, const double &iGwl) {
    double sdTarget = fabs(sdMax);
    double bmInc;
    if(sd > 0.) {
      if(sdMax > 0) {  //Yield table stocking degree
        if(sd == 1.) {bmInc = it->gIncBmt((i+1)*timeStep, mai);}
        else {bmInc = it->gIncBmSdTab((i+1)*timeStep, mai, sd);}
      } else {  //Natural stocking degree
        if(sd == 1.) {bmInc = it->gIncBm((i+1)*timeStep, mai);}
        else {bmInc = it->gIncBmSdNat((i+1)*timeStep, mai, sd);}
      }
      if(bmInc > iGwl) {bmInc = iGwl;}
      if(sd > sdTarget) {
        bmInc *= sdTarget/sd;
      } else {
        bmInc += (iGwl - bmInc) * (1. - 1./(1. + sdTarget-sd));
      }
    } else { //Current stocking degree is 0
      if(sdTarget > 0.) {
        bmInc = iGwl;
      } else {
        bmInc = 0.;
      }
    }
    return(bmInc);
  }

  int ageStruct::cohortShift() {
    if(static_cast<int>(dat.size()) > 1) {
      int i = static_cast<int>(dat.size())-2;
      if(dat[i+1].area > 0.) {
        if(static_cast<unsigned int>(i+2) < maxNumberOfAgeClasses) { //Add one more age class
          unsigned int oldSize = dat.size();
          dat.resize(i+3);
          initCohort(oldSize, i+2);
          dat[i+1] = dat[i];
        } else { //It would be good to allow one more age class
          dat[i+1].d = (dat[i+1].d * dat[i+1].area + dat[i].d * dat[i].area) / (dat[i+1].area + dat[i].area);
          dat[i+1].h = (dat[i+1].h * dat[i+1].area + dat[i].h * dat[i].area) / (dat[i+1].area + dat[i].area);
          dat[i+1].bm = (dat[i+1].bm * dat[i+1].area + dat[i].bm * dat[i].area) / (dat[i+1].area + dat[i].area);
          dat[i+1].area += dat[i].area;
        }
      } else {
        dat[i+1] = dat[i];
      }
      for(i = static_cast<int>(dat.size())-3; i>-1; --i) {
        dat[i+1] = dat[i];
      }
      dat[0].area=0; dat[0].bm=0; dat[0].d=0; dat[0].h=0;
    } else if(maxNumberOfAgeClasses > 1 && static_cast<int>(dat.size()) == 1) {
      dat.resize(2);
      initCohort(0, 1);
      dat[1] = dat[0];
      dat[0].area=0; dat[0].bm=0; dat[0].d=0; dat[0].h=0;
    }
  return(0);
  }

  unsigned int ageStruct::getMaxAge() {return(maxNumberOfAgeClasses*timeStep);}

  unsigned int ageStruct::setMaxAge(unsigned int maxAge) {
    maxNumberOfAgeClasses = ceil(maxAge/timeStep);
    return(maxNumberOfAgeClasses*timeStep);
  }

  ageStruct::v ageStruct::thinAndGrowOLD() {
    v ret = {0., 0., 0., 0., 0.};
    bool constSd = (sdDef==0) ? true : false;
//    for(int i = static_cast<int>(dat.size())-1; i>-1; --i) { //MG: GK: Georg's quick solution to the "out of range" problem
    for(int i = static_cast<int>(dat.size())-2; i>-1; --i) {
      if(dat[i].area > 0.) {
        double sdNat;
        if(i > 0) {sdNat = it->gSdNat(i*timeStep, avgMai, dat[i].bm);}
        else {  //Afforestations have a stocking degree of sdMin
          if(sdMin > 0.) {sdNat = sdMin * it->gSdNat(0, avgMai);
          } else {sdNat = -sdMin;}
        }
        double gwl = it->gIncGwlSdNat(i*timeStep, mai, sdNat);
        double sd=0.;  //Stocking degree after growing period
        if(sdMax > 0.) {sd = (gwl+dat[i].bm) / it->gBmt((i+1)*timeStep,avgMai);}
        else {sd = (gwl+dat[i].bm) / it->gBm((i+1)*timeStep, avgMai);}
        double id = it->gIncDbhSdNat(i*timeStep, avgMai, sdNat);
        //The oldest age class has no increment
        if(i == static_cast<int>(dat.size())-1) {
          gwl = 0.;
          id = 0.;
        } else {dat[i].h += it->gIncHeight(i*timeStep, avgMai);}
        bool thinningWasDone = false;
        //Key to ask if harvest is economic
        static std::vector<double> dbhBmSh(3,0);
 //Stocking degree to high or typical thinnigs are forced -> make thinning
        if(sd > fabs(sdMax) || flexSd > 0.) {
          //Thinning caused by stand density
          double reduce = 0.;
          if(constSd) {
            reduce = fabs(sdMax) / sd;
          } else {
            if(sd > fabs(sdMax)) {
              if(sdMin > 0) {
                reduce = sdMin / ((gwl/2.+dat[i].bm) / it->gBmt((i+1)*timeStep,avgMai));
              } else {
                reduce = -sdMin / ((gwl/2.+dat[i].bm) / it->gBm((i+1)*timeStep, avgMai));
              }
            }
          }
          if(reduce < 0.) {reduce = 0.;}
          if(reduce > 1.) {reduce = 1.;}
          //Thinning caused by typical thinnigs
          double reduceB = 0.;
          if(flexSd > 0.) {
            double cutVol = gwl - it->gIncBmSdNat(i*timeStep, mai, sdNat);
            if(cutVol > 0.) {
              //A weighting by sd / sdMax is done to come sometime to sdMax
              double weight;
              if(sdMax > 0) {weight = it->gSdTab(i*timeStep, avgMai, dat[i].bm) / sdMax;}
              else {weight = sdNat / (-sdMax);}
              reduceB = 1. - (cutVol / (gwl/+dat[i].bm)) * weight;
            }
          }
          if(reduceB < 0.) {reduceB = 0.;}
          if(reduceB > 1.) {reduceB = 1.;}
          //Bring in the thinning fluctuation softener
          //flexSd > 0. && flexSd <= 1.
          reduce = reduce * (1. - flexSd) + reduceB * flexSd;
          dbhBmSh[0] = dat[i].d+id/2.;
          dbhBmSh[1] = dat[i].bm+gwl/2.;
          dbhBmSh[2] = (1. - reduce);
          if(reduce > 0. && dov->g(dbhBmSh)) { //Do Thinning if it is economic
            thinningWasDone = true;
            if(constSd) {
              dat[i].bm += gwl;
              dat[i].d += id;
            } else {
              //The Stands get half increment at high stand density
              dat[i].bm += gwl/2.;
              dat[i].d += id/2.;
            }
            double totalWood = dat[i].area * dat[i].bm * (1. - reduce);
            dat[i].bm *= reduce;
            double harvestedWood = totalWood * hlv->g(dat[i].d);
            double sawnWood = harvestedWood * sws->g(dat[i].d);
            ret.area += dat[i].area;
            ret.bm += totalWood;
            ret.sw += sawnWood;
            ret.rw += harvestedWood - sawnWood;
            ret.co += totalWood * cov->g(dbhBmSh);
            if(!constSd) {
              //The Stands get half increment at low stand density
              sdNat = it->gSdNat((i+.5)*timeStep, avgMai, dat[i].bm);
              dat[i].bm += it->gIncBmSdNat((i+.5)*timeStep, avgMai, sdNat)/2.;
              dat[i].d += it->gIncDbhSdNat((i+.5)*timeStep, avgMai, sdNat)/2.;
            }
          } else {  //No thinning
            dat[i].bm += gwl;
            dat[i].d += id;
          }
        } else {  //No thinning
          dat[i].bm += gwl;
          dat[i].d += id;
        }
        //Check if Bm is not higher than maximum possible at end of period
        double bmMax = it->gBm((i+1)*timeStep, avgMai);
        if(dat[i].bm > bmMax) {
          if(constSd) {
            dat[i].bm = bmMax;
          } else {
            if(thinningWasDone) {
              double totalWood = dat[i].area * (dat[i].bm - bmMax);
              dat[i].bm = bmMax;
              double harvestedWood = totalWood * hlv->g(dat[i].d);
              double sawnWood = harvestedWood * sws->g(dat[i].d);
              ret.bm += totalWood;
              ret.sw += sawnWood;
              ret.rw += harvestedWood - sawnWood;
              ret.co += totalWood * cov->g(dbhBmSh);
            } else {
              dat[i].bm = bmMax;
            }
          }
        }
        //Bring the Data to the next age class
//        if(i < static_cast<int>(dat.size())-1) {//MG: GK: Georg's quick solution to the "out of range" problem
        if(i < static_cast<int>(dat.size())-2) {
          dat[i+1] = dat[i];
        } else {
          dat[i+1].d = (dat[i+1].d * dat[i+1].area + dat[i].d * dat[i].area) / (dat[i+1].area + dat[i].area);
          dat[i+1].h = (dat[i+1].h * dat[i+1].area + dat[i].h * dat[i].area) / (dat[i+1].area + dat[i].area);
          dat[i+1].bm = (dat[i+1].bm * dat[i+1].area + dat[i].bm * dat[i].area) / (dat[i+1].area + dat[i].area);
          dat[i+1].area += dat[i].area;
        }
        dat[i].area = 0.; dat[i].bm = 0.; dat[i].d = 0.; dat[i].h = 0.;
      }
    }
    if(ret.area > 0.) { //Values per hectare
      ret.sw /= ret.area; ret.rw /= ret.area; ret.co /= ret.area; ret.bm /= ret.area;}
    return(ret);
  }
//MG: added
  int ageStruct::getActiveAge()
  { // MG: Find "active Age" - the oldest cohort with area>0
      int i = static_cast<int>(dat.size())-1;
      while(dat[i].area == 0 && i>0) {--i;}
     while(dat[i].area != 0 && i<static_cast<int>(dat.size())-1) {++i;}
    //int activeAge = int((static_cast<int>(dat.size())-1) * timeStep);
	int activeAge = (i-1) * timeStep;
    return(activeAge);

  }

  double ageStruct::BiomSdTab(int i){
	  return (it->gBmSdNat(i, mai, 1.));
  }

  double ageStruct::Biomass(int i)
  {
	  cout <<"mai "<<mai<<endl;
	  return (it->gBmt(i, mai));
	  
	 // double mai_avg=calcAvgMai();
	  //cout << "Avarage MAI  " << mai_avg << endl;
	  //double increment = it->gGwl(10., mai_avg);
	//  double biomStock = it->gBm(10., mai);


	  /*  double refStock_0_10 = (it->gBmt(0., mai)+it->gBmt(1., mai)+it->gBmt(2., mai)+it->gBmt(3., mai)+it->gBmt(4., mai)+it->gBmt(5., mai)+it->gBmt(6., mai)+it->gBmt(7., mai)+it->gBmt(8., mai)+it->gBmt(9., mai)+it->gBmt(10., mai))/0.993351/11;
	    double refStock_11_20 = (it->gBmt(11., mai)+it->gBmt(12., mai)+it->gBmt(13., mai)+it->gBmt(14., mai)+it->gBmt(15., mai)+it->gBmt(16., mai)+it->gBmt(17., mai)+it->gBmt(18., mai)+it->gBmt(19., mai)+it->gBmt(20., mai))/0.993351/10;
		double refStock_21_30 = (it->gBmt(21., mai)+it->gBmt(22., mai)+it->gBmt(23., mai)+it->gBmt(24., mai)+it->gBmt(25., mai)+it->gBmt(26., mai)+it->gBmt(27., mai)+it->gBmt(28., mai)+it->gBmt(29., mai)+it->gBmt(30., mai))/0.993351/10;
		double refStock_31_40 = (it->gBmt(31., mai)+it->gBmt(32., mai)+it->gBmt(33., mai)+it->gBmt(34., mai)+it->gBmt(35., mai)+it->gBmt(36., mai)+it->gBmt(37., mai)+it->gBmt(38., mai)+it->gBmt(39., mai)+it->gBmt(40., mai))/0.993351/10;
		double refStock_41_50 = (it->gBmt(41., mai)+it->gBmt(42., mai)+it->gBmt(43., mai)+it->gBmt(44., mai)+it->gBmt(45., mai)+it->gBmt(46., mai)+it->gBmt(47., mai)+it->gBmt(48., mai)+it->gBmt(49., mai)+it->gBmt(50., mai))/0.993351/10;
		double refStock_51_60 = (it->gBmt(51., mai)+it->gBmt(52., mai)+it->gBmt(53., mai)+it->gBmt(54., mai)+it->gBmt(55., mai)+it->gBmt(56., mai)+it->gBmt(57., mai)+it->gBmt(58., mai)+it->gBmt(59., mai)+it->gBmt(60., mai))/0.993351/10;
		double refStock_61_70 = (it->gBmt(61., mai)+it->gBmt(62., mai)+it->gBmt(63., mai)+it->gBmt(64., mai)+it->gBmt(65., mai)+it->gBmt(66., mai)+it->gBmt(67., mai)+it->gBmt(68., mai)+it->gBmt(69., mai)+it->gBmt(70., mai))/0.993351/10;
		double refStock_71_80 = (it->gBmt(71., mai)+it->gBmt(72., mai)+it->gBmt(73., mai)+it->gBmt(74., mai)+it->gBmt(75., mai)+it->gBmt(76., mai)+it->gBmt(77., mai)+it->gBmt(78., mai)+it->gBmt(79., mai)+it->gBmt(80., mai))/0.993351/10;
		double refStock_81_90 = (it->gBmt(81., mai)+it->gBmt(82., mai)+it->gBmt(83., mai)+it->gBmt(84., mai)+it->gBmt(85., mai)+it->gBmt(86., mai)+it->gBmt(87., mai)+it->gBmt(88., mai)+it->gBmt(89., mai)+it->gBmt(90., mai))/0.993351/10;
		double refStock_91_100 = (it->gBmt(91., mai)+it->gBmt(92., mai)+it->gBmt(93., mai)+it->gBmt(94., mai)+it->gBmt(95., mai)+it->gBmt(96., mai)+it->gBmt(97., mai)+it->gBmt(98., mai)+it->gBmt(99., mai)+it->gBmt(100., mai))/0.993351/10;
		double refStock_101_110 = (it->gBmt(101., mai)+it->gBmt(102., mai)+it->gBmt(103., mai)+it->gBmt(104., mai)+it->gBmt(105., mai)+it->gBmt(106., mai)+it->gBmt(107., mai)+it->gBmt(108., mai)+it->gBmt(109., mai)+it->gBmt(110., mai))/0.993351/10;
		double refStock_111_120 = (it->gBmt(111., mai)+it->gBmt(112., mai)+it->gBmt(113., mai)+it->gBmt(114., mai)+it->gBmt(115., mai)+it->gBmt(116., mai)+it->gBmt(117., mai)+it->gBmt(118., mai)+it->gBmt(119., mai)+it->gBmt(120., mai))/0.993351/10;
		double refStock_121_130 = (it->gBmt(121., mai)+it->gBmt(122., mai)+it->gBmt(123., mai)+it->gBmt(124., mai)+it->gBmt(125., mai)+it->gBmt(126., mai)+it->gBmt(127., mai)+it->gBmt(128., mai)+it->gBmt(129., mai)+it->gBmt(130., mai))/0.993351/10;
		double refStock_131_140 = (it->gBmt(131., mai)+it->gBmt(132., mai)+it->gBmt(133., mai)+it->gBmt(134., mai)+it->gBmt(135., mai)+it->gBmt(136., mai)+it->gBmt(137., mai)+it->gBmt(138., mai)+it->gBmt(139., mai)+it->gBmt(140., mai))/0.993351/10;
		double refStock_141_150 = (it->gBmt(141., mai)+it->gBmt(142., mai)+it->gBmt(143., mai)+it->gBmt(144., mai)+it->gBmt(145., mai)+it->gBmt(146., mai)+it->gBmt(147., mai)+it->gBmt(148., mai)+it->gBmt(149., mai)+it->gBmt(150., mai))/0.993351/10;
		double refStock_151_160 = (it->gBmt(151., mai)+it->gBmt(152., mai)+it->gBmt(153., mai)+it->gBmt(154., mai)+it->gBmt(155., mai)+it->gBmt(156., mai)+it->gBmt(157., mai)+it->gBmt(158., mai)+it->gBmt(159., mai)+it->gBmt(160., mai))/0.993351/10;
		double refStock_161_170 = (it->gBmt(161., mai)+it->gBmt(162., mai)+it->gBmt(163., mai)+it->gBmt(164., mai)+it->gBmt(165., mai)+it->gBmt(166., mai)+it->gBmt(167., mai)+it->gBmt(168., mai)+it->gBmt(169., mai)+it->gBmt(170., mai))/0.993351/10;
		double refStock_171_180 = (it->gBmt(171., mai)+it->gBmt(172., mai)+it->gBmt(173., mai)+it->gBmt(174., mai)+it->gBmt(175., mai)+it->gBmt(176., mai)+it->gBmt(177., mai)+it->gBmt(178., mai)+it->gBmt(179., mai)+it->gBmt(180., mai))/0.993351/10;
		double refStock_181_190 = (it->gBmt(181., mai)+it->gBmt(182., mai)+it->gBmt(183., mai)+it->gBmt(184., mai)+it->gBmt(185., mai)+it->gBmt(186., mai)+it->gBmt(187., mai)+it->gBmt(188., mai)+it->gBmt(189., mai)+it->gBmt(190., mai))/0.993351/10;
		double refStock_191_200 = (it->gBmt(191., mai)+it->gBmt(192., mai)+it->gBmt(193., mai)+it->gBmt(194., mai)+it->gBmt(195., mai)+it->gBmt(196., mai)+it->gBmt(197., mai)+it->gBmt(198., mai)+it->gBmt(199., mai)+it->gBmt(200., mai))/0.993351/10;
		double refStock_201_210 = (it->gBmt(201., mai)+it->gBmt(202., mai)+it->gBmt(203., mai)+it->gBmt(204., mai)+it->gBmt(205., mai)+it->gBmt(206., mai)+it->gBmt(207., mai)+it->gBmt(208., mai)+it->gBmt(209., mai)+it->gBmt(210., mai))/0.993351/10;
		double refStock_211_220 = (it->gBmt(211., mai)+it->gBmt(212., mai)+it->gBmt(213., mai)+it->gBmt(214., mai)+it->gBmt(215., mai)+it->gBmt(216., mai)+it->gBmt(217., mai)+it->gBmt(218., mai)+it->gBmt(219., mai)+it->gBmt(220., mai))/0.993351/10;

		
//
	  
*/
	  
	  
 }
  // int ageStruct::TEST (double age_test, double mai_test, double bio_test, double rtime_test)
  int ageStruct::TEST (double age_test, double mai_test)
  {
	 // cout<<"checkpoint test"<<endl;
	 // return (it->gBmSdTab(age_test, mai_test, sd_test));
	  //return(it->gBmt(age_test, mai_test));
	SD_t=1.;
	//GS_t_ha = it->gBmt(age_test, mai_test);	// growing stock with thinning per ha
	ClassArea = getArea(age_test);
	IncTot_ha = it->gIncGwlt(age_test, mai_test);// total increment with thinning per ha 
	  R_ha =  it->gRemBmt (age_test, mai_test); // removals for thinning or mortality per ha
	GS_t_Tab_ha = it->gBmSdTab(age_test, mai_test,SD_t); // growing stock from yield table per ha
	D_t_ha =it -> gDbht(age_test, mai_test);//gHeight(age_test, mai_test); 
	   //SD_t = it->gAvgSdTab(rtime_test, mai_test,bio_test); // stocking degree for selected growing stock and MAI
	   
	  return (0);
  }
}

#endif
