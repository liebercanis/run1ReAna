#include "TBaconEvent.hxx"
ClassImp(TBaconEvent)

TBaconEvent::TBaconEvent(): TNamed("TBaconEvent","TBaconEvent")
{
  clear();
}

//TBaconEvent::~TBaconEvent(){}

void TBaconEvent::clear()
{
  run=0;
  event=0;
  npmt=0;
  npulse=0;
  nspe=0;
  totQ=0;
  spe=0;
  F40=0;
  muVmax=0;
  //
  T0=0;
  totalCharge=0;
  tMax=0;
  vMax=0;
  cMax=0;
  baseline=0;
  sDev=0; 

  hits.clear();
}

