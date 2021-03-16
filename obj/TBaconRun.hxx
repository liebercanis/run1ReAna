/**
** MG, May 2020 
**/
#ifndef TBACONRUN_DEFINED
#define TBACONRUN_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>

using namespace std;

// class to store info for the run

class TBaconRun: public TNamed {
  public:
    TBaconRun();
    //~TBaconRun();

    void clear();
    Int_t run;
    Long64_t totTrigger;
    Long64_t goodTrigger;

    Int_t nevn;
    Double_t mpvn;
    Double_t mpvnErr;
    Double_t widn;
    Double_t widnErr;
    Double_t spenSum;
    //
    Int_t nevc;
    Double_t maxc;
    Double_t maxcErr;
    Double_t mpvc;
    Double_t mpvcErr;
    Double_t widc;
    Double_t widcErr;
    Double_t specSum;
    //
    Int_t nevq;
    Double_t mpvq;
    Double_t mpvqErr;
    Double_t widq;
    Double_t widqErr;
    Double_t speqSum;
    //
    Double_t aveSpe;
    Double_t aveSpeErr;

    ClassDef(TBaconRun,3)
};
#endif

