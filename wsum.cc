/* MG sum waveforms for reanalysis Jan 13 2021*/
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <valarray>
//
#include <TROOT.h>
#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
//
#include "TPmtEvent.hxx"

class wsum
{
public:
  enum
  {
    NPMT = 2
  };
  double microSec = 1.0E6;
  wsum(TString tag = "XenonDoping1ppm_20038", Int_t maxEvents = 0);
  virtual ~wsum() { ; }

  // summed wave histograms
  TH1D *hPMTRaw[NPMT];
  TH1D *hPMTSum[NPMT];
};

wsum::wsum(TString tag, Int_t maxEvents)
{
  printf(" starting wsum tag %s \n", tag.Data());

  // open ouput file and make some histograms

  TString fileName;
  fileName.Form("rootData/DS2/%s.root", tag.Data());
  printf(" looking for file %s\n", fileName.Data());
  TFile *fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
  {
    printf(" couldnt open file %s\n", fileName.Data());
    return;
  }
  else
    printf("  found file %s \n", fileName.Data());
  // get pmtTree from file
  TTree *pmtTree = NULL;
  fin->GetObject("pmtTree", pmtTree);
  Long64_t nentries = pmtTree->GetEntries();

  // set up memory for reading
  TPmtEvent *pmtEvent = new TPmtEvent();
  Int_t iset = pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
  if (iset != 0)
  {
    printf(" couldnt find pmtEvent in file  %s\n", fileName.Data());
    return;
  }

  //switch to output file
  TString outfileName;
  outfileName.Form("/data1/bacon/waveSum/wsum_%s_%i.root", tag.Data(), maxEvents);
  printf(" opening output file %s \n", outfileName.Data());
  TFile *outfile = new TFile(outfileName, "recreate");
  outfile->cd();
  // loop over entries
  if (maxEvents > 0)
    nentries = maxEvents;
  printf(" STARTING RUN %s with  events  %lld of %lld  \n", tag.Data(), nentries, pmtTree->GetEntries());
  for (Long64_t ientry = 0; ientry < nentries; ientry++)
  {
    pmtTree->GetEntry(ientry);
    if (pmtEvent->time.size() == 0)
      continue;
    Int_t nSamples = pmtEvent->time.size();

    // check how many pmts we have
    int gotPMT = 0;
    if (pmtEvent->volt1.size() > 0)
      ++gotPMT;
    //if (pmtEvent->volt2.size() > 0)
    //  ++gotPMT;
    if (ientry == 0)
      printf(" .... events %lld samples %i PMT0 %zu PMT1 %zu \n", pmtTree->GetEntries(), nSamples, pmtEvent->volt1.size(), pmtEvent->volt2.size());
    if (gotPMT < 1)
      break;

    // define pmt signal histograms
    Double_t pmtXLow = pmtEvent->time[0] * microSec;
    Double_t pmtXHigh = pmtEvent->time[nSamples - 1] * microSec;
    // define histos
    if (ientry == 0)
    {
      for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
      {
        TString hname;
        hname.Form("PMTRaw%i-Ev%lld-%s", ipmt, ientry, tag.Data());
        hPMTRaw[ipmt] = new TH1D(hname, hname, nSamples, pmtXLow, pmtXHigh);
        hname.Form("PMTSum%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hPMTSum[ipmt] = new TH1D(hname, hname, nSamples, pmtXLow, pmtXHigh);
      }
    }

    // zero histos
    for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
      hPMTRaw[ipmt]->Reset();
    // loop over PMT 1
    for (unsigned isample = 0; isample < pmtEvent->volt1.size(); isample++)
    {
      Double_t volt0 = pmtEvent->volt1[isample];
      Double_t digi0 = -1.0 * (double(volt0));
      hPMTRaw[0]->SetBinContent(isample + 1, digi0);
      hPMTSum[0]->SetBinContent(isample + 1, digi0 + hPMTSum[0]->GetBinContent(isample + 1));
    }

    /* loop over PMT 2
    for (unsigned isample = 0; isample < pmtEvent->volt2.size(); isample++)
    {
      Double_t volt1 = pmtEvent->volt2[isample];
      Double_t digi1 = -1.0 * (double(volt1));
      hPMTRaw[1]->SetBinContent(isample + 1, digi1);
      hPMTSum[1]->SetBinContent(isample + 1, digi1 + hPMTSum[1]->GetBinContent(isample + 1));
    }
    */
  }

  printf(" END RUN %s  %lld of %lld \n", tag.Data(), nentries, pmtTree->GetEntries());

  outfile->ls();
  outfile->Write();
}