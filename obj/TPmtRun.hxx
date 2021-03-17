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
    void resize();
		// data elements
    //
    // Event Info
    Int_t    run;
    UShort_t  gpsYear;
    UShort_t  gpsDay;
    UInt_t    gpsSec;
    UInt_t   gpsNs;
		std::vector<TPmtEvent> event;
    // Event Info
    //
    Int_t entry;
    Double_t deltaT;
    Int_t nPMTs;

    // Pulse Info
    //
    std::vector<std::vector<Double_t>> charge;
    std::vector<std::vector<Double_t>> startTimes;
    std::vector<std::vector<Double_t>> peakWidths;
    std::vector<std::vector<Double_t>> peakHeights;
    std::vector<std::vector<Double_t>> peakTimes;
    // Waveform info
    //
    std::vector<Double_t> T0;
    std::vector<Double_t> totalCharge;
    std::vector<Double_t> tMax;
    std::vector<Double_t> vMax;
    std::vector<Double_t> cMax;
    std::vector<Double_t> baseline;
    std::vector<Double_t> sDev; 

		ClassDef(TPmtRun,1)
};
#endif

