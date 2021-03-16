#include "TPmtRun.hxx"
ClassImp(TPmtRun)

TPmtRun::TPmtRun(): TNamed("TPmtRun","TPmtRun")
{
  clear();
}

//TPmtRun::~TPmtRun(){}

void TPmtRun::clear()
{
  run=0;
  gpsYear=0;
  gpsDay=0;
  gpsSec=0;
  gpsNs=0;  
  event.clear();	 
 }

