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
  npulse=0;
  nspe=0;
  wsum=0;
  qsum=0;
  q900=0;
  hits.clear();
}

