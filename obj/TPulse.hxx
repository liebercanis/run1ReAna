/**
** MG, March 2020 
**/
#ifndef TPULSE_DEFINED
#define TPULSE_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store pmt hit 

class TPulse: public TNamed {
  public:
    TPulse();
    ~TPulse(){;}

    void clear();
    // data elements
    //  int nev;
    int istart;
    int good;
    int kind;
    float time;
    float tpeak;
    float peak;
    float pwidth;
    float q;
    float qerr;

    ClassDef(TPulse,7)
};
#endif

