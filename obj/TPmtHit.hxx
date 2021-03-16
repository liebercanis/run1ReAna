/**
** MG, July 2018 
**/
#ifndef TPMTHIT_DEFINED
#define TPMTHIT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store pmt hit 

class TPmtHit: public TNamed {
  public:
    TPmtHit();
    ~TPmtHit(){;}

    void clear();
    std::vector<Double_t> getPulse(Int_t nwid, std::vector<Double_t> v ) {
      std::vector<Double_t> pulse;
      if(v.size()<peakBin+nwid) return pulse;
      for(Int_t i=peakBin-nwid; i<= peakBin+nwid; ++i) pulse.push_back(v[i]);
      return pulse;
    }
    // data elements
    Int_t firstBin;
    Int_t lastBin;
    Double_t startTime;
    Double_t peakWidth;
    Double_t qpeak;
    UInt_t peakt;
    Double_t peakMaxTime;
    Int_t peakBin;
    Double_t qsum;
    Double_t qerr;
    Int_t good;
    Int_t kind;

    ClassDef(TPmtHit,7)
};
#endif

