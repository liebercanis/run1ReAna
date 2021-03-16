#include "TPmtHit.hxx"
ClassImp(TPmtHit)

TPmtHit::TPmtHit(): TNamed("TPmtHit","TPmtHit")
{
  clear();
}

//TPmtHit::~TPmtHit(){}

void TPmtHit::clear()
{
  firstBin=0;
  lastBin=0;
  startTime=0;
  peakWidth=0;
  qpeak=0;
  peakt=0;
  peakMaxTime=0;
  peakBin=0;
  qsum=0;
  qerr=0;
  good=0;
  kind=-1; // unassigned
}

