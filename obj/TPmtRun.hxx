/**
** MG, August 2017 
**/
#ifndef TPMTRUN_DEFINED
#define TPMTRUN_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TPmtEvent.hxx"

using namespace std;

// class to store info for the event

class TPmtRun: public TNamed {
	public:
		TPmtRun();
    //		~TPmtRun();

		void clear();
		// data elements
    Int_t    run;
    UShort_t  gpsYear;
    UShort_t  gpsDay;
    UInt_t    gpsSec;
    UInt_t   gpsNs;
		std::vector<TPmtEvent> event;		 
		ClassDef(TPmtRun,1)
};
#endif

