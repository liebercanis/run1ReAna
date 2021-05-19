/**
** MG, March 2020 
**/
#ifndef TBACONEVENT_DEFINED
#define TBACONEVENT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TPulse.hxx"

using namespace std;

// class to store info for the event

class TBaconEvent : public TNamed
{
public:
  TBaconEvent();
  //~TBaconEvent();

  void clear();
  // data elements
  Int_t run;
  Int_t event;
  Int_t npulse;
  Int_t npmt;
  Int_t nspe;
  Double_t totQ;
  Double_t spe;
  Double_t F40;
  Double_t muVmax;
  //
  Double_t T0;
  Double_t totalCharge;
  Double_t tMax;
  Double_t vMax;
  Double_t cMax;
  Double_t baseline;
  Double_t sDev;

  std::vector<TPulse> hits;
  ClassDef(TBaconEvent, 4)
};
#endif
