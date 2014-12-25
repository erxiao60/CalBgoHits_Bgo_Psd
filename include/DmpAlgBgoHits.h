/* for BGO Hits                           */
/* yong zhzhy@mail.ustc.edu.cn 08-10-2014 */
#ifndef DmpAlgBgoHits_H
#define DmpAlgBgoHits_H

#include "DmpVAlg.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtNudRaw.h"
#include "DmpEvtBgoDyCoe.h"
#include "DmpEvtBgoMips.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtNudHits.h"
#include "TVector3.h"
class DmpEvtBgoRaw;
class DmpEvtBgoDyCoe;
class DmpEvtBgoMips;
class DmpEvtBgoHits;
class DmpEvtHeader;

class DmpAlgBgoHits : public DmpVAlg{
/*
 *  DmpAlgBgoHits
 *
 */
public:
  DmpAlgBgoHits();
  ~DmpAlgBgoHits();

  //void Set(const std::string &type,const std::string &value);
  // if you need to set some options for your algorithm at run time. Overload Set()
  bool Initialize();
  bool GetDyCoePar();
  bool GetMipsPar();
  bool GetAttPar();
  bool ProcessThisEvent();    // only for algorithm
  bool ProcessPsdHits();    // only for algorithm
  bool ProcessNudHits();    // only for algorithm
  bool Finalize();
  bool Reset(); 
  bool GetPsdPed();
  bool GetPsdMips();
  bool GetPsdDyCoe(); 
  void SaveHeader(){fSaveEvtHeader = true;}
private:
DmpEvtHeader *fEvtHeader;
  bool fSaveEvtHeader;
DmpEvtBgoRaw *fBgoRaw;
DmpEvtBgoDyCoe *fBgoDyCoe;
DmpEvtBgoMips  *fBgoMips;
DmpEvtBgoHits *fBgoHits;

DmpEvtBgoRaw *fPsdRaw;
DmpEvtBgoHits *fPsdHits;

DmpEvtNudRaw *fNudRaw;
DmpEvtNudHits *fNudHits;

double DyCoePar_58[14][22][2][2];//layer,bar,side, 0:Slope and 1:Intercept
double DyCoePar_25[14][22][2][2];//layer,bar,side, 0:Slope and 1:Intercept
double MipsPar[14][22][3][3];//layer,bar,side :2 mean combined value,0:MPV 1:Gsigma 2:Lwidth
double AttPar[14][22][2];//layer,bar,0:Slope 1:Intercept

double PsdPedMean[2][41][2][2];//layer,bar,side, dy
double PsdPedSigma[2][41][2][2];//
double PsdMips[2][41][2][3];//layer,bar,side, 0:MPV 1:Gsigma 2:Lwidth

//time cut
  int timecut;

//dynode choise
  double adc_dy2[14][24][2];
  double adc_dy5[14][24][2];
  double adc_dy8[14][24][2];
//Mipsenergy
  double Mipsenergy;

};

#endif
