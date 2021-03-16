/**
** MG, June 2018 
**/
#ifndef TPMTEVENT_DEFINED
#define TPMTEVENT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store info for the event

class TPmtEvent: public TNamed {
	public:
		TPmtEvent();
    //		~TPmtEvent();

		void clear();
		// data elements
    Int_t    event;
		std::vector<Double_t> time;		 
		std::vector<Double_t> volt1;		 
		std::vector<Double_t> volt2;		 
		ClassDef(TPmtEvent,2)
};
#endif

