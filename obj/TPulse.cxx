#include "TPulse.hxx"
ClassImp(TPulse)

TPulse::TPulse(): TNamed("TPulse","TPulse")
{
  clear();
}

//TPulse::~TPulse(){}

void TPulse::clear()
{
    istart=0;
    nwidth=0;
    good=-1;
    kind=-1;
    time=0;
    thit=0;
    q=0;
    qerr=0;
}

