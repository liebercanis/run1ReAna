#include "TPulse.hxx"
ClassImp(TPulse)

    TPulse::TPulse() : TNamed("TPulse", "TPulse")
{
  clear();
}

//TPulse::~TPulse(){}

void TPulse::clear()
{
  istart = 0;
  iend = 0;
  good = 0;
  kind = 0;
  time = 0;
  tpeak = 0;
  peak = 0;
  pwidth = 0;
  q = 0;
  qerr = 0;
  baseline = 0;
  baselineErr = 0;
  slope = 0;
  slopeErr = 0;
  baseChi = 0;
  digi.clear();
}
