#include <iostream>
#include <fstream>
#include <numeric>
#include "compiled/BaconAnalysis.hh"


// time is in microseconds
void testNtuples()
{


  TFile* inFile = new TFile("BaconAnalysisHighCut-2.root", "READONLY");
  //TFile* fout = new TFile("testNtuples.root", "RECREATE");

  if (!inFile)
    return;

  inFile->ls();

  TString sumString = setSumNames();
  TNtuple *tSummary=NULL;
  inFile->GetObject("Summary", tSummary);
  printf(" sum %lld \n", tSummary->GetEntries());

  TString preString = setPreNames();
  TNtuple *tEvPre = NULL;
  inFile->GetObject("EvPre", tEvPre);
  printf(" evpre %lld \n", tEvPre->GetEntries());

  TString evString = setEvNames();
  TNtuple *tEvent = NULL;
  inFile->GetObject("Event", tEvent);
  printf(" ev %lld \n", tEvent->GetEntries());

  printf("sumvars string %s \n", sumString.Data());
  for(int iv=0; iv< SUMVARS; ++iv ) tSummary->SetBranchAddress(sumNames[iv],&sumVars[iv]);
  printf("preEvent string %s \n", preString.Data());
  for(int iv=0; iv< PREVARS; ++iv ) tEvPre->SetBranchAddress(preNames[iv],&preVars[iv]);
  printf("evEvent string %s \n", evString.Data());
  for(int iv=0; iv< EVVARS; ++iv ) tEvent->SetBranchAddress(evNames[iv],&evVars[iv]);

  //fout->Write();

} 
