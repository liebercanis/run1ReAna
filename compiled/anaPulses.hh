#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm> // std::sort
#include "TSpectrum.h"
#include <math.h>
//
#include "TPmtEvent.hxx"

class anaPulses
{
public:
  anaPulses(TString tag, Int_t maxEvents);
  virtual ~anaPulses() { ; }
  enum
  {
    NPMT = 2
  };
  TTree *pmtTree;
  TPmtEvent *pmtEvent;
  void anaEntry(Long64_t ientry);
  std::vector<Double_t> Derivative(std::vector<Double_t>, Int_t NSteps);
  std::vector<Int_t> PeakFinding(std::vector<Double_t>, std::vector<Double_t>, Int_t pmtNum, Int_t entry);
  std::vector<Double_t> BubbleSort(std::vector<Double_t> A);
  Double_t CalculateMean(std::vector<Double_t> vec);
  Double_t CalculateMean(std::vector<Double_t> vec, Int_t startBin, Int_t stopBin);
  Double_t CalculateMean(std::vector<Int_t> vec);
  Double_t CalculateSdev(std::vector<Double_t> vec, Int_t startBin, Int_t stopBin, Double_t mean);
  std::vector<Double_t> TrigFilter(std::vector<Double_t> sig, Int_t N, Int_t pmtNum, Int_t ientry);
  std::pair<Double_t, Double_t> BaselineSubtraction(std::vector<Double_t> sig, Int_t pmtNum, Int_t entry);
  std::pair<Double_t, Double_t> BaselineSubtraction(std::vector<Double_t> sig, std::vector<Int_t> weight, Int_t pmtNum, Int_t entry);
  std::vector<Double_t> BaselineSubtraction(std::vector<Double_t> sig, Int_t weight, Int_t pmtNum, Int_t ientry);
  std::vector<Double_t> Integral(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry);
  std::vector<Double_t> TrapFilter(std::vector<Double_t> sig, Int_t ramp, Int_t flat, Int_t pmtNum, Int_t ientry);

  Int_t nSamples;
  Int_t NEvents;
  Int_t NHistograms = 1000;
  double microSec = 1.0E6;
  Int_t irun = 1;

  Double_t minPeakWidth = 2.5e-9;
  Double_t maxPeakWidth = 100e-9;
  std::vector<std::vector<Double_t>> signal;
  std::vector<std::vector<Double_t>> derivative;
  std::vector<std::vector<Double_t>> integral;

  // summed wave histograms
  TH1D *hPMTRaw[NPMT];
  TH1D *hPMTSum[NPMT];
  TH1D *hSignal[2];
  TH1D *hDerivative[2];
  TH1D *hIntegral[2];
  //    TH1D * hFilter[2];
  TH1D *hBaseline[2];
  TH1D *hSum[2];
  TH1D *hSumBaseline[2];
  TNtuple *ntuplePulse;
  TNtuple *ntupleEvent;
};
