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
#include <TDatime.h>
#include "TRandom2.h"

//
#include "TPmtEvent.hxx"
#include "TPmtSimulation.hxx"
#include "TPmtRun.hxx"

class anaSim
{
public:
  anaSim(Int_t irun, Int_t nDer, Int_t dsNum);
  anaSim();
  ~anaSim() { ; };

  std::vector<Double_t> DownSampler(std::vector<Double_t>, Int_t NSteps);
  std::vector<Double_t> Derivative(std::vector<Double_t>, Int_t NSteps);
  Double_t RMSCalculator(std::vector<Double_t> vec, Int_t ientry);
  std::vector<Double_t> BaselineWMA(std::vector<Double_t> sig, std::vector<Int_t> peaks, Int_t ientry, Int_t N, Int_t pmtNum, Double_t waveformFraction);
  std::vector<Int_t> PeakFinding(std::vector<Double_t>, std::vector<Double_t>, Int_t pmtNum, Int_t entry);
  std::vector<Double_t> BubbleSort(std::vector<Double_t> A);
  Double_t CalculateMean(std::vector<Double_t> vec);
  Double_t CalculateMean(std::vector<Int_t> vec);
  Double_t CalculateMean(std::vector<Double_t> vec, Int_t startBin, Int_t stopBin);
  Double_t CalculateSdev(std::vector<Double_t> vec, Int_t startBin, Int_t stopBin, Double_t mean);
  std::vector<Double_t> TrigFilter(std::vector<Double_t> sig, Int_t N, Int_t pmtNum, Int_t ientry);
  std::pair<Double_t, Double_t> BaselineSubtraction(std::vector<Double_t> sig, Int_t pmtNum, Int_t entry);
  std::pair<Double_t, Double_t> BaselineSubtraction(std::vector<Double_t> sig, std::vector<Int_t> weight, Int_t pmtNum, Int_t entry);
  std::vector<Double_t> BaselineSubtraction(std::vector<Double_t> sig, Int_t weight, Int_t pmtNum, Int_t ientry);
  std::vector<Double_t> Integral(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry);
  std::vector<Double_t> TrapFilter(std::vector<Double_t> sig, Int_t ramp, Int_t flat, Int_t pmtNum, Int_t ientry);
  std::vector<Double_t> LowPassFrequency(std::vector<Double_t> input, Double_t cutoff, Double_t sampleRate);
  std::vector<Double_t> NotchFilter(std::vector<Double_t> input, Double_t notch, Double_t BW, Double_t sampleRate);
  TPmtSimulation *simEvent;
  TPmtEvent *pmtEvent;
  TPmtEvent *processedEvent;
  TPmtRun *pmtRun;

  Int_t NEvents;
  Int_t NHistograms = 100;
  Double_t pi = 3.14159265359;
  Double_t deltaT = 0;
  bool peakFindingDebug = false;

  Double_t minPeakWidth = 3e-9;
  Double_t maxPeakWidth = 5000e-9;
  std::vector<std::vector<Double_t>> signal;
  std::vector<std::vector<Double_t>> derivative;
  std::vector<std::vector<Double_t>> integral;

  TH1D *hSignal[2];
  TH1D *hDerivative[2];
  TH1D *hBaseline[2];
  TH1D *hWindowSum[2];
  TH1D *hSumWaveform[2];
  TH1D *hSum[2];
  TH1D *hSumFresh[2];
  TH1D *hSumBaseline[2];
  TH1D *hIntegral[2];
  TH1D *hFlatBaseline[2];
  TH1D *hRawSig[2];
};
