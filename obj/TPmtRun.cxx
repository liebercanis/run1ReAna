#include "TPmtRun.hxx"
ClassImp(TPmtRun)

TPmtRun::TPmtRun(): TNamed("TPmtRun","TPmtRun")
{
  clear();
}

//TPmtRun::~TPmtRun(){}

void TPmtRun::clear()
{
  //Run Info
  run=0;
  gpsYear=0;
  gpsDay=0;
  gpsSec=0;
  gpsNs=0;  
  event.clear();

  //Event Info
  entry = 0;
  deltaT = 0;
  nPMTs = 0;
  //Pulse Info
  charge.clear();
  startTimes.clear();
  peakWidths.clear();
  peakHeights.clear();
  peakTimes.clear();
  //Waveform info
  T0.clear();
  totalCharge.clear();
  tMax.clear();
  vMax.clear();
  cMax.clear();
  baseline.clear();
  sDev.clear();
}

void TPmtRun::resize()
{

  //Pulse Info
  charge.resize(nPMTs);
  startTimes.resize(nPMTs);
  peakWidths.resize(nPMTs);
  peakHeights.resize(nPMTs);
  peakTimes.resize(nPMTs);
  //Waveform Info
  T0.resize(nPMTs);
  totalCharge.resize(nPMTs);
  tMax.resize(nPMTs);
  vMax.resize(nPMTs);
  cMax.resize(nPMTs);
  baseline.resize(nPMTs);
  sDev.resize(nPMTs);
}

