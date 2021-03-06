/* uses TBaconEvent Class */
//////////////////////////////////////////////////////////
//  M.Gold May 2020
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>
#include <map>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm> // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"

#include "TBaconEvent.hxx"
#include "TBaconRun.hxx"

/**
** Double_t tMaxCut = 5e-6,tMinCut = 0.00e-6,vMaxEventCut = 10,vMinCut = 1e-3,peakWidthCut = 0,nHits = 10;//25e-9;
**  startTime > 9e-7
**        double peakTimeMean = 1.53e-6,peakTimeSigma = 100e-9;
          if(peakTime > peakTimeMean - peakTimeSigma && peakTime < peakTimeMean + peakTimeSigma ) continue;
          if(charge < meanSPE[j] - 2*sigmaSPE[j]) continue;
*/

class TBReader
{
public:
  TBReader(Int_t runstart = 3000, Int_t runstop = 29999);
  virtual ~TBReader() { ; }
  enum
  {
    MAXSETS = 6
  };
  TString runTag[MAXSETS];
  TString runRange[MAXSETS];
  int setEvents[MAXSETS];
  TFile *fout;
  TDirectory *trigDir;
  TDirectory *setDir;
  void newRun();
  void getFileList();
  void getMaxValue(TH1D *h, int &maxbin);
  void setMaxValue(TH1D *h, double vmax);
  void getRunSet();
  int irun;
  int iset;
  std::map<int, string> runList;
  Int_t currentRun;
  Int_t currentSet;
  Long64_t treeNumber;
  Long64_t nloop;
  TChain *tree;
  TBaconEvent *bEvent;
  TBaconRun *bRun;
  TTree *treeRun;
  TNtuple *ntTrig;

  TH1D *hLife;
  TH1D *hChargeSum;
  TH1D *hLifeRun;
  TH1D *hLifeSet[MAXSETS];
  TH1D *hLifeACutRun;
  TH1D *hChargeSumRun;
  TH1D *hChargeCutRun;
  TH1D *hLifeQCutRun;
  TH1D *hChargeQCutRun;
};
void TBReader::getRunSet()
{
  iset = -1;
  irun = -1;
  irun = bEvent->run;
  if (bEvent->run >= 3000 && bEvent->run <= 3020)
  {
    irun = irun - 3000 + 20000 - 20;
    iset = 0;
  }
  if (bEvent->run >= 20005 && bEvent->run < 20020)
    iset = 1;
  if (bEvent->run >= 20025 && bEvent->run <= 20040)
    iset = 2;
  if (bEvent->run >= 20045 && bEvent->run <= 20060)
    iset = 3;
  if (bEvent->run >= 20065 && bEvent->run <= 20080)
    iset = 4;
  if (bEvent->run >= 20220 && bEvent->run <= 20222)
    iset = 5;
  if (iset < 0)
  {
    printf(" ..... skipping run  %i  \n", bEvent->run);
  }
}

void TBReader::getFileList()
{
  TString dirname("/data2/mgold/run1RootData/");
  TString ext("_0.root");
  TString prefix("anaPulsesDecon");
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files)
  {
    TSystemFile *file;
    TString tfname;
    string fname;
    TIter next(files);
    std::map<int, string>::iterator mapIter = runList.begin();
    while ((file = (TSystemFile *)next()))
    {
      fname = string(file->GetName());
      tfname = TString(file->GetName());
      if (!file->IsDirectory() && tfname.EndsWith(ext) && tfname.BeginsWith(prefix))
      {
        vector<string> tokens;
        string intermediate;
        stringstream check1(fname);
        // Tokenizing w.r.t. space ' '
        while (std::getline(check1, intermediate, '_'))
          tokens.push_back(intermediate);
        int irun = atoi(tokens[tokens.size() - 2].c_str());
        //cout << fname << "  run " << irun << '\n';
        runList.insert(mapIter, std::pair<int, string>(irun, fname));
      }
    }
  }
  printf(" getFileList found %lu files: \n", runList.size());
  for (std::map<int, string>::iterator mapIter = runList.begin(); mapIter != runList.end(); ++mapIter)
    std::cout << mapIter->first << " => " << mapIter->second << '\n';
}

void TBReader::getMaxValue(TH1D *h, int &maxbin)
{
  maxbin = 0;
  double vmax = 0;
  for (int j = 0; j < h->GetNbinsX(); ++j)
  {
    if (h->GetBinContent(j) > vmax)
    {
      maxbin = j;
      vmax = h->GetBinContent(j);
    }
  }
  double fmin = h->GetBinLowEdge(maxbin - 15);
  double fmax = h->GetBinLowEdge(maxbin + 15);
  // gauss fit to region
  h->Fit("gaus", "Q0R", "", fmin, fmax);
  TF1 *fgaus = h->GetFunction("gaus");
  double fmaxcErr = h->GetBinWidth(maxbin) / sqrt(12.0);
  double fmaxc = h->GetBinLowEdge(maxbin) + h->GetBinWidth(maxbin) / 2.0;
  if (fgaus)
  {
    bRun->maxc = fgaus->GetParameter(1);
    bRun->maxcErr = fgaus->GetParError(1);
  }
  else
  {
    bRun->maxc = fmaxc;
    bRun->maxcErr = fmaxcErr;
  }
  printf(">> get max value %s maxbin %i (%.3f+/-%.3f) (%.3f+/-%.3f) \n ", h->GetName(), maxbin, fmaxc, fmaxcErr, bRun->maxc, bRun->maxcErr);
  return;
}

void TBReader::setMaxValue(TH1D *h, double vmax)
{

  int maxbin = h->FindBin(vmax);
  double fmin = h->GetBinLowEdge(TMath::Max(1, maxbin - 10));
  double fmax = h->GetBinLowEdge(TMath::Min(h->GetNbinsX(), maxbin + 10));
  bRun->maxc = h->GetBinLowEdge(maxbin) + h->GetBinWidth(maxbin) / 2.0;
  bRun->maxcErr = h->GetBinWidth(maxbin) / sqrt(12.0);
  // gauss fit to region
  h->Fit("gaus", "0R", "", fmin, fmax);
  TF1 *fgaus = h->GetFunction("gaus");
  double maxFit = 0;
  double maxFitErr = 0;
  if (fgaus)
  {
    maxFit = fgaus->GetParameter(1);
    maxFitErr = fgaus->GetParError(1);
  }
  if (fmax > fmin && fmax < fmax)
  {
    bRun->maxc = maxFit;
    bRun->maxcErr = maxFitErr;
  }
  printf(">> set max value %s bin %i  (%.3f to %.3f) max (%.3f+/-%.3f) lan (%.3f+/-%.3f) \n ", h->GetName(), maxbin, fmin, fmax, bRun->maxc, bRun->maxcErr, bRun->mpvc, bRun->mpvcErr);
  return;
}

void TBReader::newRun()
{

  double fmax = 0;
  double fmaxErr = 0;
  int maxbin = 0;
  if (hChargeCutRun)
  {
    getMaxValue(hChargeCutRun, maxbin);
    hChargeCutRun->Fit("landau");
    TF1 *lfit = hChargeCutRun->GetFunction("landau");
    //fill bRun branch
    if (lfit)
    {
      bRun->mpvc = lfit->GetParameter(1);
      bRun->mpvcErr = lfit->GetParError(1);
      bRun->widc = lfit->GetParameter(2);
      bRun->widcErr = lfit->GetParError(2);
      //setMaxValue(hChargeCutRun,bRun->mpvc);
    }
    hChargeSumRun->Fit("landau");
    TF1 *lfit2 = hChargeSumRun->GetFunction("landau");
    if (lfit2)
    {
      bRun->mpvn = lfit2->GetParameter(1);
      bRun->mpvnErr = lfit2->GetParError(1);
      bRun->widn = lfit2->GetParameter(2);
      bRun->widnErr = lfit2->GetParError(2);
    }
    hChargeQCutRun->Fit("landau");
    TF1 *lfit3 = hChargeQCutRun->GetFunction("landau");
    if (lfit3)
    {
      bRun->mpvq = lfit3->GetParameter(1);
      bRun->mpvqErr = lfit3->GetParError(1);
      bRun->widq = lfit3->GetParameter(2);
      bRun->widqErr = lfit3->GetParError(2);
    }
    bRun->nevn = hChargeSumRun->GetEntries();
    bRun->spenSum = hLifeRun->Integral();
    bRun->nevc = hChargeSumRun->GetEntries(); // ten hit cut
    bRun->specSum = hLifeACutRun->Integral(); // afterpulse cut
    bRun->nevq = hChargeQCutRun->GetEntries();
    bRun->speqSum = hLifeQCutRun->Integral();
    //
    if (bRun->nevq > 0)
    {
      bRun->aveSpe = bRun->specSum / bRun->nevq;
      // root n counting errors
      bRun->aveSpeErr = bRun->aveSpe * sqrt(1. / bRun->specSum + 1. / bRun->nevq);
    }
    treeRun->Fill();
    printf("   at entry %lld end of run %i tot(good) %lld(%lld) trigs  ", nloop, currentRun, bRun->totTrigger, bRun->goodTrigger);
    printf("   ending run %d size %lld events %d totSPE %E ", bRun->run, treeRun->GetEntries(), bRun->nevc, bRun->specSum);
    printf("   max  %.3f+/-%.3f mpv %.3f+/-%.3f \n", bRun->maxc, bRun->maxcErr, bRun->mpvc, bRun->mpvcErr);
  }

  bRun->clear();
  bRun->run = bEvent->run;
  currentRun = bEvent->run;
  hLifeRun = (TH1D *)hLife->Clone(Form("LifeRun%05i", int(currentRun)));
  hLifeACutRun = (TH1D *)hLife->Clone(Form("LifeQACutRun%05i", int(currentRun)));
  hLifeQCutRun = (TH1D *)hLife->Clone(Form("LifeQCutRun%05i", int(currentRun)));
  hChargeSumRun = (TH1D *)hChargeSum->Clone(Form("ChargeSumRun%05i", int(currentRun)));
  hChargeCutRun = (TH1D *)hChargeSum->Clone(Form("ChargeCutRun%05i", int(currentRun)));
  hChargeQCutRun = (TH1D *)hChargeSum->Clone(Form("ChargeQCutRun%05i", int(currentRun)));
  /*
  fout->Add(hLifeRun);
  fout->Add(hLifeACutRun);
  fout->Add(hLifeQCutRun);
  fout->Add(hChargeSumRun);
  fout->Add(hChargeCutRun);
  fout->Add(hChargeQCutRun);
  */
  getRunSet();
  if (iset != currentSet)
  {
    currentSet = iset;
    printf(" starting  set %d \n", iset);
    setDir->cd();
    hLifeSet[iset] = (TH1D *)hLife->Clone(Form("LifeSet-%s-%s", runRange[iset].Data(), runTag[iset].Data()));
    hLifeSet[iset]->SetTitle(Form("LifeSet %s %s", runRange[iset].Data(), runTag[iset].Data()));
    fout->cd();
  }
  printf(" starting run %d irun %d set %d \n", bRun->run, irun, iset);
}

TBReader::TBReader(Int_t runstart, Int_t runstop)
{
  for (int jset = 0; jset < MAXSETS; ++jset)
    setEvents[jset] = 0;

  runTag[0] = TString("00PPM-MU");
  runTag[1] = TString("01PPM-Ran");
  runTag[2] = TString("02PPM-Ran");
  runTag[3] = TString("05PPM-Ran");
  runTag[4] = TString("10PPM-Ran");
  runTag[5] = TString("10PPM-MU");
  runRange[0] = TString("03000-03020");
  runRange[1] = TString("20005-20020");
  runRange[2] = TString("20025-20040");
  runRange[3] = TString("20045-20060");
  runRange[4] = TString("20065-20080");
  runRange[5] = TString("20220-20222");

  printf("TBReader starting run %i to %i \n", runstart, runstop);
  getFileList();
  currentRun = -1;
  currentSet = -1;
  treeNumber - 1;
  hLifeRun = NULL;
  hChargeSumRun = NULL;
  hChargeCutRun = NULL;
  hLifeQCutRun = NULL;
  hLifeACutRun = NULL;
  hChargeQCutRun = NULL;

  double maxLife = 10.0;
  int lifeBins = int(1.0E4 / 8.0); // ns bins
  double QSumCutValue = 200;
  double startTimeMax = 9.0E-7;
  double startRatioCut = 0.05;

  //int lifeBins=400;

  //bool VmaxCut = false;

  tree = new TChain("TBacon");
  bEvent = new TBaconEvent();
  tree->SetBranchAddress("bevent", &bEvent);
  tree->GetListOfFiles()->ls();
  for (std::map<int, string>::iterator mapIter = runList.begin(); mapIter != runList.end(); ++mapIter)
  {
    TString addFile = TString("run1RootData/") + TString(mapIter->second.c_str());
    if (mapIter->first >= runstart && mapIter->first <= runstop)
    {
      std::cout << " add run " << mapIter->first << '\n';
      tree->Add(addFile);
    }
  }
  cout << " TBacon has " << tree->GetEntries() << endl;

  fout = new TFile(Form("TBReader-%i-%i-%i-bins-notrig.root", runstart, runstop, lifeBins), "RECREATE");
  cout << " output file is  " << fout->GetName() << endl;
  trigDir = fout->mkdir("trigDir");
  setDir = fout->mkdir("setDir");
  fout->cd();

  treeRun = new TTree("tRun", " by Run data ");
  bRun = new TBaconRun();
  treeRun->Branch("brun", &bRun);

  TH2D *hNHitTotQ = new TH2D("NHitByTotQ", " nhit by tot Q ", 50, 0, 50, 50, 0, 50);
  TH1D *hNHits = new TH1D("NHits", " hits in event  ", 350, 0, 350);
  TH1D *hTotQ = new TH1D("TotQ", " totQ in event  after 10 cut ", 350, 0, 350);
  TH2D *hStartCut = new TH2D("StartCut", " pre trigger hit SPE ", 110, 0., 1.1, 100, 0., 100.);
  TH1D *hStartSumCut = new TH1D("StartSumCut", " 700ns pre trigger sum SPE < 0.05 ", 1000, 0., 5.);

  // trigger
  std::map<double, double> mapDer;
  std::map<double, double> mapMax;
  ntTrig = new TNtuple("ntTrig", "", "event:sigTime:sigQ:derTime:derQ:maxTime:maxQ");

  TH1D *hTrigQhit = new TH1D("TrigQhit", " trig hit charge ", 1000, .9, 1.1);
  TH1D *hTrigQsum = new TH1D("TrigQsum", " trig sum charge ", 1000, .9, 1.1);
  TH2D *hTriqQsumvsTime = new TH2D("TrigQsumvsTime", " charge vs time  ", 100, .9, 1.1, 1000, 0, 200);

  TH1D *hTrigTimeQsig = new TH1D("TrigTimeQsig", " 3 Qspe Trigger time ", 240, .9, 1.5);
  TH1D *hTrigTimeQder = new TH1D("TrigTimeQder", " max der Trigger time ", 240, .9, 1.5);
  TH1D *hTrigTimeQmax = new TH1D("TrigTimeQmax", " max Q   Trigger time ", 240, .9, 1.5);

  // event
  TH1D *hMuVmax = new TH1D("MuVmax", " muVmax ", 1000, 0, 1);
  hMuVmax->GetXaxis()->SetTitle("muon Vmax");
  TNtuple *ntEvent = new TNtuple("ntEvent", " event ", "entry:maxBin:maxVal:run:ntrig:nhits:totQ:totHit");

  //hit
  TH1D *hHitResidual = new TH1D("HitResidual7microsec", " residual  ", 1000, 0, 5);

  // prototypes
  hLife = new TH1D("LifeCut", " lifetime PMT >10 hits ", lifeBins, 0, maxLife);
  hLife->GetXaxis()->SetTitle(" micro-seconds ");
  hLife->SetMarkerColor(kBlack);
  hLife->SetMarkerStyle(22);
  hLife->SetMarkerSize(.2);
  hLife->Sumw2();

  hChargeSum = new TH1D("ChargeSum", " good summed event charge > 10 pulses  ", 1000, 0, 10000);
  hChargeSum->GetXaxis()->SetTitle(" run sum dt-Q (x10^9) ");
  hChargeSum->Sumw2();

  nloop = tree->GetEntries();

  for (Long64_t entry = 0; entry < nloop; ++entry)
  {
    Long64_t itree = tree->LoadTree(entry);
    if (itree != treeNumber)
    {
      tree->SetBranchAddress("bevent", &bEvent);
      treeNumber = itree;
    }
    tree->GetEntry(entry);
    // find run range
    if (bEvent->run < runstart)
      continue;
    if (bEvent->run > runstop)
      break;

    // if new run
    if (bEvent->run != currentRun)
      newRun();
    if (bEvent->npmt != 0) // in anaPules, i set to one.. will fix
      continue;            // only use PMT zero

    setEvents[iset] += 1;
    ++bRun->totTrigger;
    if (entry % 1000 == 0)
      printf(" ... %lld run %d  bevent %lld nhits %lu \n", entry, bEvent->run, ntEvent->GetEntries(), bEvent->hits.size());

    hNHitTotQ->Fill(bEvent->totQ, bEvent->hits.size());
    hNHits->Fill(bEvent->hits.size());
    if (bEvent->hits.size() < 10)
      continue;
    hTotQ->Fill(bEvent->totQ);

    double qspe = 1E9 * bEvent->spe;
    // early trigger data
    TH1D *hTrigQhitEvent = NULL;
    TH1D *hTrigQsumEvent = NULL;
    if (bRun->totTrigger < 10)
    {
      trigDir->cd();
      hTrigQhitEvent = (TH1D *)hTrigQhit->Clone(Form("TrigQRun%iEvent%i", int(bEvent->run), int(entry)));
      hTrigQsumEvent = (TH1D *)hTrigQsum->Clone(Form("TrigQSumRun%iEvent%i", int(bEvent->run), int(entry)));
      fout->cd();
    }

    // find start time
    double lastQhit = 0;
    double qtrigSum = 0;
    double triggerTime = -1;
    double trigQ = 0;
    double trigSig = 3;
    mapDer.clear();
    mapMax.clear();

    // hit loop to find trigger time
    for (unsigned ip = 0; ip < bEvent->hits.size(); ++ip)
    {
      double hitTime = bEvent->hits[ip].time * 1E6;
      if (hitTime < 0.9)
        continue;
      if (hitTime > 1.5)
        continue;
      double hitq = 1.0E9 * bEvent->hits[ip].q / qspe; // spe
      mapDer.insert(pair<double, double>(hitq - lastQhit, hitTime));
      lastQhit = hitq;
      mapMax.insert(pair<double, double>(hitq, hitTime));
      qtrigSum += hitq;
      int ibin = hTrigQhit->FindBin(hitTime);
      hTrigQhit->SetBinContent(ibin, hitq);
      hTrigQsum->SetBinContent(ibin, qtrigSum);
      hTriqQsumvsTime->Fill(hitTime, qtrigSum);

      if (hitq > trigSig && triggerTime == -1)
      {
        triggerTime = hitTime;
        trigQ = hitq;
      }

      if (hTrigQhitEvent)
      {
        hTrigQhitEvent->SetBinContent(ibin, hitq);
        hTrigQsumEvent->SetBinContent(ibin, qtrigSum);
        //printf("\t xxxx event %i hit %i ibin %i time %f q %f sum %f  \n",entry, ip, ibin, hitTime, hitq, qtrigSum);
      }
    }

    /* start time ratio cut. already made
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      if(bEvent->hits[ip].time>startTimeMax) break;
      startCharge += bEvent->hits[ip].q;
    }
    hStartRatio->Fill(startCharge/bEvent->totalCharge);
    if(startCharge/bEvent->totalCharge > startRatioCut) continue;
    */

    // hit loop to find QSum and pre trigger charge
    double QSumCut = 0;
    double QSum = 0;
    double startCharge = 0;
    for (unsigned ip = 0; ip < bEvent->hits.size(); ++ip)
    {
      double hitTime = bEvent->hits[ip].time * 1E6;
      double hitq = 1.0E9 * bEvent->hits[ip].q / qspe; // spe
      if (hitq <= 0)
        continue;
      bool peakCut = (bEvent->hits[ip].peak > 0.795 && bEvent->hits[ip].peak < 0.804) || (bEvent->hits[ip].peak > 1.604 && bEvent->hits[ip].peak < 1.611);
      QSum += hitq;
      if (!peakCut)
        QSumCut += hitq;
      double relTime = hitTime - triggerTime + 1;
      if (relTime < 1.1)
        hStartCut->Fill(relTime, hitq);
      if (relTime < 1. && relTime > 0.3)
        startCharge += hitq;
    }
    hStartSumCut->Fill(startCharge);
    if (startCharge > 0.05)
      continue;

    hChargeSumRun->Fill(QSum);
    hChargeCutRun->Fill(QSumCut);
    if (QSumCut > QSumCutValue)
      hChargeQCutRun->Fill(QSumCut);

    std::map<double, double>::iterator derIter;
    std::map<double, double>::iterator maxIter;

    derIter = mapDer.begin();
    maxIter = mapMax.begin();

    hTrigTimeQsig->Fill(triggerTime);
    ntTrig->Fill(float(entry), triggerTime, trigQ, derIter->second, derIter->first, maxIter->second, maxIter->first);

    // throw out events with bad trigger
    //if(triggerTime<0.9||triggerTime>1.1) continue;

    ++bRun->goodTrigger;

    // muonVmax cut
    hMuVmax->Fill(bEvent->muVmax);
    //if(VmaxCut&&bEvent->muVmax  < 0.05) continue;
    // loop over pulses
    for (unsigned ip = 0; ip < bEvent->hits.size(); ++ip)
    {
      //cout << entry << "  " << hitq << endl;
      //if(theHit.order==0) triggerTime = theHit.time;
      double hitTime = 1.0E6 * bEvent->hits[ip].time - triggerTime + 1;
      double hitq = 1.0E9 * bEvent->hits[ip].q / qspe; // spe

      // throw out negative hits
      if (hitq <= 0)
        continue;

      // guess  at noise of tenth qspe
      //  ***** setting time unit 1.000000E-09 maxLife 10.000000 # digis 10000
      double width = (qspe * .5) * 1.0E9 * bEvent->hits[ip].pwidth;
      //double hitqerr= sqrt(2*(abs(hitq)+ width*width));
      //gain is 4.36e6 electrons per photon
      double hitqerr = sqrt(abs(hitq) + width * width);
      double hitres = hitq / hitqerr;

      if (hitTime > 7)
        hHitResidual->Fill(hitres);

      int hitBin = hLifeRun->FindBin(hitTime);
      double hitqerr2 = pow(hitqerr, 2);
      if (std::isnan(double(hitqerr2)))
        printf(" XXXXXXX %f %f \n", hitq, hitqerr);
      bool peakCut = (bEvent->hits[ip].peak > 0.795 && bEvent->hits[ip].peak < 0.804) || (bEvent->hits[ip].peak > 1.604 && bEvent->hits[ip].peak < 1.611);
      // cut bad hits
      //if(peakCut||afterCut) continue;
      hLifeRun->SetBinContent(hitBin, hLifeRun->GetBinContent(hitBin) + hitq);
      hLifeRun->SetBinError(hitBin, sqrt(pow(hLifeRun->GetBinError(hitBin), 2) + pow(hitqerr, 2)));
      hLifeSet[iset]->SetBinContent(hitBin, hLifeSet[iset]->GetBinContent(hitBin) + hitq);
      hLifeSet[iset]->SetBinError(hitBin, sqrt(pow(hLifeSet[iset]->GetBinError(hitBin), 2) + pow(hitqerr, 2)));
      if (QSumCut > QSumCutValue)
      {
        hLifeQCutRun->SetBinContent(hitBin, hLifeQCutRun->GetBinContent(hitBin) + hitq);
        hLifeQCutRun->SetBinError(hitBin, sqrt(pow(hLifeQCutRun->GetBinError(hitBin), 2) + pow(hitqerr, 2)));
        hLifeQCutRun->SetBinContent(hitBin, hLifeQCutRun->GetBinContent(hitBin) + hitq);
        hLifeQCutRun->SetBinError(hitBin, sqrt(pow(hLifeQCutRun->GetBinError(hitBin), 2) + pow(hitqerr, 2)));
      }
      if (!peakCut)
      { // for counting without afterpulses
        hLifeACutRun->SetBinContent(hitBin, hLifeACutRun->GetBinContent(hitBin) + hitq);
        hLifeACutRun->SetBinError(hitBin, sqrt(pow(hLifeACutRun->GetBinError(hitBin), 2) + pow(hitqerr, 2)));
        hLifeACutRun->SetBinContent(hitBin, hLifeACutRun->GetBinContent(hitBin) + hitq);
        hLifeACutRun->SetBinError(hitBin, sqrt(pow(hLifeACutRun->GetBinError(hitBin), 2) + pow(hitqerr, 2)));
      }
    }

    ntEvent->Fill(float(entry), bEvent->maxBin, bEvent->maxVal, float(bEvent->run), float(bRun->goodTrigger), float(bEvent->hits.size()), float(bEvent->totQ), float(bEvent->totHit));
  }
  // run landau fit for last run
  newRun();
  //fout->ls();
  // normalize
  double xbin;
  double ebin;
  for (int jset = 0; jset < MAXSETS; ++jset)
  {
    printf("set %i events %i \n", jset, setEvents[jset]);
    for (int ibin = 0; ibin < hLifeSet[jset]->GetNbinsX(); ++ibin)
    {
      xbin = hLifeSet[jset]->GetBinContent(ibin) / double(setEvents[jset]);
      //ebin = hLifeSet[jset]->GetBinError(ibin) / sqrt(double(setEvents[jset]));
      ebin = sqrt(xbin);
      hLifeSet[jset]->SetBinContent(ibin, xbin);
      hLifeSet[jset]->SetBinError(ibin, ebin);
    }
  }

  fout->Write();

  return;
}
