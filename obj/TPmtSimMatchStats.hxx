/**
** MG December 2018 
**/
#ifndef TPMTSIMMATCHSTATS_DEFINED
#define TPMTSIMMATCHSTATS_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store info for the event

class TPmtSimMatchStats: public TNamed {
	public:
		TPmtSimMatchStats();

		void clear();
    void fill(unsigned nsim, unsigned nhit, unsigned nmatch, unsigned nnot, unsigned nmiss);
    void print();
		// data elements
    unsigned entries;
    unsigned tSim;
    unsigned tHit;
    unsigned tMatch;
    unsigned tNot;
    unsigned tMiss;
    Double_t eff;
    Double_t effError;
    Double_t over;
    Double_t overError;
		ClassDef(TPmtSimMatchStats,1)
};
#endif

