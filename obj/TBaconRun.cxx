#include "TBaconRun.hxx"
ClassImp(TBaconRun)

TBaconRun::TBaconRun(): TNamed("TBaconRun","TBaconRun")
{
  clear();
}

//TBaconRun::~TBaconRun(){}

void TBaconRun::clear()
{
  run=0;
  totTrigger=0;
  goodTrigger=0;
  //
  nevn=0;
  mpvn=0;
  mpvnErr=0;
  widn=0;
  widnErr=0;
  spenSum=0;
  //
  nevc=0;
  maxc=0;
  maxcErr=0;
  mpvc=0;
  mpvcErr=0;
  widc=0;
  widcErr=0;
  specSum=0;
  //
  nevq=0;
  mpvq=0;
  mpvqErr=0;
  widq=0;
  widqErr=0;
  speqSum=0;
  //
  aveSpe=0;
  aveSpeErr=0;
}

