/**
** MG, March 2020 
**/
#ifndef TBACONEVENT_DEFINED
#define TBACONEVENT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TPulse.hxx"

using namespace std;

// class to store info for the event

class TBaconEvent: public TNamed {
	public:
		TBaconEvent();
    //~TBaconEvent();

		void clear();
		// data elements
    Int_t    run;
    Int_t    event;
    Int_t    npulse;
    Int_t    nspe;
    Double_t wsum;
    Double_t qsum;
    Double_t q900;
    std::vector<TPulse> hits;
		ClassDef(TBaconEvent,1)
};
#endif

