#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TProof.h"
#include "TFile.h"
#include "TParameter.h"
#include "TPmtEvent.cxx"
#include "BaconAnalysis.hh"

#include <iostream>
#include <string>
#include <vector>
#include <dirent.h>
#include <chrono>
#include <algorithm>
#include <bitset>


using namespace std;
mt19937 gen( chrono::system_clock::now().time_since_epoch().count() );

TH1D *h_v1 = NULL;
TH1D *h_v1Fix = NULL;
TH1D *h_v1histo = NULL;
TH1D *h_v1histobaseline = NULL;
TH1D *h_v1histobaselineend = NULL;
TH1D *h_runqualityB = NULL;
double baseline_mean = 0;
double baseline_rms = 0;
double baseline_min = 0;
double baselineend_mean = 0;
double baselineend_rms = 0;

enum
{
  NBITS =5
};
std::bitset<NBITS> rejectBits;
TH1D *hCode;
// event cut routine
int rejectEvent(){
  ///////////////////////////
  
  int event_code = 0;
  rejectBits.reset();
  int bin_a = h_v1histo->FindFirstBinAbove(0);
  int bin_b = h_v1histo->FindFirstBinAbove(3);
  if (bin_a == bin_b)
  {
    event_code += 1; //skip saturated events
    rejectBits.set(0);
    h_runqualityB->Fill(1);
  }
  ///////////////////////////
  baseline_mean = h_v1histobaseline->GetMean();
  baseline_rms = h_v1histobaseline->GetRMS();
  baseline_min = h_v1histobaseline->GetXaxis()->GetBinCenter(h_v1histobaseline->FindFirstBinAbove(0));
  if (baseline_mean - 5 * baseline_rms > baseline_min)
  {
    event_code += 2; //skip events with peaks on baseline
    rejectBits.set(1);
    //printf("... baseline  %f \n",baseline_min);
    h_runqualityB->Fill(2);
  }
  ///////////////////////////
  baselineend_mean = h_v1histobaselineend->GetMean();
  baselineend_rms = h_v1histobaselineend->GetRMS();
  if (baseline_rms > 2 * baselineend_rms)
  {
    event_code += 4; //slope in baseline + early events
    //printf("... baseline rms %f \n",baseline_rms);
    rejectBits.set(2);
    h_runqualityB->Fill(3);
  }
  ///////////////////////////
  double wf_counter = h_v1histo->Integral(0, h_v1histo->GetXaxis()->FindBin(baseline_mean - 5 * baseline_rms));
  if (wf_counter < 10)
  {
    event_code += 8; // events that are only noise (or tiny signals)
    rejectBits.set(3);
    //printf("... wf_counter %f \n",wf_counter);
    h_runqualityB->Fill(4);
  }
  ///////////////////////////
  if ((h_v1->GetMinimumBin() > 1100) || (h_v1->GetMinimumBin() < 900))
  {
    event_code += 16; //skip if singlet is not the highest peak (likely pileup, check with Michael)
    rejectBits.set(4);
    //printf("... min bin  %i \n", (int) h_v1->GetMinimumBin());
    h_runqualityB->Fill(5);
  }
  if(event_code==0) h_runqualityB->Fill(15);

  hCode->Fill(NBITS+1); // total events tested 
  for (int ib = 0; ib < NBITS; ++ib )
    if(rejectBits.test(ib))
      hCode->Fill(ib);
  return event_code;
}

enum
{
  NSETS = 6
};

int getRunSet(int irun)
{
  int iset = -1;
  if (irun >= 3000 && irun <= 3020)
  //if (irun >= 3000 && irun <= 5700)
  {
    //irun = irun - 3000 + 20000 - 20;
    iset = 0;
  }
  if (irun >= 20005 && irun < 20020)
    iset = 1;
  if (irun >= 20025 && irun <= 20040)
    iset = 2;
  if (irun >= 20045 && irun <= 20060)
    iset = 3;
  if (irun >= 20065 && irun <= 20080)
    iset = 4;
  if (irun >= 20220 && irun <= 20222)
    iset = 5;
  if (iset < 0)
  {
    printf(" ..... skipping run  %i  not in a run set \n", irun);
  }
  return iset;
}

double getMaxTime(TH1D* h, int istart, double& vmax) {
  int max = istart;
  int ilow = h->FindBin(afterLow*1.E3);
  int ihigh = h->FindBin(afterHigh * 1.E3);
  int ilow2 = h->FindBin(afterLow2 * 1.E3);
  int ihigh2 = h->FindBin(afterHigh2 * 1.E3);

  vmax = 0;
  for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin) {
    if ((ibin >= ilow && ibin <= ihigh) || (ibin >= ilow2 && ibin <= ihigh2))
      continue;
    if (h->GetBinContent(ibin) > vmax)
    {
      vmax = h->GetBinContent(ibin);
      max = ibin;
    }
  }
  return h->GetBinCenter(max);
}

double getSum(TH1D *h, double start, double end )
{
  int istart = h->FindBin(start);
  int iend = h->FindBin(end);
  iend = min(iend, h->GetNbinsX());
  double sum = 0;
  for (int ibin = istart; ibin < iend; ++ibin) {
    sum += h->GetBinContent(ibin);
  }
  return sum;
}

double afterPulsing(TH1D *h, int i1, int i2, int j)
{

  double y1 = 0;
  for (int ib = i1 - 9; ib <= i1; ++ib)
    y1 += h->GetBinContent(ib);
  y1 /= double(10);

  double y2 = 0;
  for (int ib = i2; ib <= i2 + 9; ++ib)
    y2 += h->GetBinContent(ib);
  y2 /= double(10);

  double x1 = h->GetBinCenter(i1);
  double x2 = h->GetBinCenter(i2);
  double xj = h->GetBinCenter(j);
  double m = (y2 - y1) / (x2 - x1);
  double b = y2 - m * x2;
  return m * xj + b;
}

void afterFix(TH1D *h, TH1D* hFix ) {
  TH1D *hTemp  = (TH1D *) h->Clone(Form("h%s-temp",h->GetName()));
<<<<<<< HEAD
  hTemp->SetDirectory(0);
=======
>>>>>>> ad97c1d945828312f4c52be34216249877ca6650
  // shift trigger times
  // align all the trigger times
  int maxBin = h->GetMaximumBin();
  int startBin = h->FindBin(1.0E3);
  for (int k = 0; k < h->GetNbinsX(); k++) {
    int jBin = k - maxBin + startBin;
<<<<<<< HEAD
    // baseline subtraction
=======
>>>>>>> ad97c1d945828312f4c52be34216249877ca6650
    hTemp->SetBinContent(jBin, h->GetBinContent(k));
    hTemp->SetBinError(jBin, h->GetBinError(k));
  }

  int ilow = hTemp->FindBin(afterLow * 1.E3);
  int ihigh = hTemp->FindBin(afterHigh * 1.E3);
  int ilow2 = hTemp->FindBin(afterLow2 * 1.E3);
  int ihigh2 = hTemp->FindBin(afterHigh2 * 1.E3);

  for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
  {
    double val = hTemp->GetBinContent(ibin);
    if (ibin >= ilow && ibin <= ihigh)
      val = afterPulsing(h, ilow, ihigh, ibin);
    if (ibin >= ilow2 && ibin <= ihigh2)
      val = afterPulsing(h, ilow2, ihigh2, ibin);
    hFix->SetBinContent(ibin, val );
    hFix->SetBinError(ibin,  h->GetBinError(ibin) );
  }
}

int currentSet = -1;
int currentRun = -1;
//////////////////////////////////////////////////////////////////////////
void BaconAnalysis(int maxFiles ){
  
  double microToNano = 1.0E3;
  auto start = chrono::steady_clock::now();

  double singleCut[NSETS] = {4.,4.,4.,4.,4.,4.};
  cout<< " max files " << maxFiles << endl;
  printf(" single cuts are :\n");
  for (int iset = 0; iset < NSETS; ++iset)
    printf(" \t set %i cut %.1f \n", iset, singleCut[iset]);

  TH1D *hWave[NSETS];
  for (int iset = 0; iset < NSETS; ++iset )
      hWave[iset] = NULL;

  const char *inDir = "/data1/bacon/rootData/DS2";

  TFile *fout = new TFile(Form("BaconAnalysisHighCut-%i.root",maxFiles),"recreate");
  TDirectory *qualDir = fout->mkdir("quality");
  TDirectory *waveDir = fout->mkdir("waveForms");
  TDirectory *rejectDir = fout->mkdir("rejectWaveForms");

  fout->cd();

  cout << inDir << endl;
  // /mnt/229fa628-4721-4c78-b916-eca13f19f0d8/296464/BaconData";
  char* dir = gSystem->ExpandPathName(inDir);
  void* dirp = gSystem->OpenDirectory(dir);
  const char* ext = ".root";


  
  const char* entry;
  TString str;
  
  vector<TString> files;
  while((entry = (char*)gSystem->GetDirEntry(dirp))) {
    str = entry;
    bool isRunFile = str.EndsWith(ext) && (str.Contains("Recirculate") || str.Contains("XenonDoping"));
    if (isRunFile)
    {
      //cout << gSystem->ConcatFileName(dir, entry) << endl;
      files.push_back(gSystem->ConcatFileName(dir, entry));
    }
  }
  
  cout << "---" << files.size() << " files " << endl;
  sort (files.begin(), files.end());
  //shuffle(files.begin(), files.end(), gen); //random runs for testing
<<<<<<< HEAD
  printf("--- from %s to %s  \n", files[0].Data(), files[files.size()-1].Data());
=======
  /*cout << " list of files to run " << endl;
  for (unsigned j = 0; j < files.size(); ++j)
    printf("\t %i %s \n", j, files[j].Data());
    */
>>>>>>> ad97c1d945828312f4c52be34216249877ca6650
  TString filename, f2;
  int filenumber;
  int nentries;
  int ntotal = 0;
  int ntotalgood = 0;
  int nrungood;
  TTree *t1 = 0;
  int m;
  double time_a;
  double time_b;
  
  double wf_min, wf_singlet, wf_width;
  double min_baseline;
  //vector<int> good_entries;
  bool skipfile;

  
  Float_t sumVars[SUMVARS];
  Float_t evVars[EVVARS];
  TH1D *h_runquality = NULL;
  TH1D *h_runqualityhisto = NULL;
  TH1D *h_runqualityS = NULL;
  TH1D *h_runqualityhistoS = NULL;
  TH1D *h_runqualityT = NULL;
  TH1D *h_runqualityhistoT = NULL;
  TH2D *h_runquality2D = NULL;
  TH1D *h_runqualityR = NULL;
  TH1D *h_minimum = NULL;
  TH1D *h_maximum = NULL;
  TH1D *h_mean = NULL;

  TH1D *h_v1b = NULL;
  TH1D *h_v1c = NULL;

  // these for cuts
  h_v1histo = new TH1D("h_v1histo", "h_v1histo", 5000, -5, 5);
  h_v1histobaseline = new TH1D("h_v1histobaseline","h_v1histobaseline",5000,-5,5);
  h_v1histobaselineend = new TH1D("h_v1histobaselineend","h_v1histobaselineend",5000,-5,5);

  hCode = new TH1D("rejectCode","rejectCode",NBITS+1,0,NBITS+1);
  TNtuple *ntSummary = new TNtuple("Summary", " summary","run:set:base:baseend:accept:total:singlet:dublet:triplet:ngood:over:minbase");
  TNtuple *ntEvPre = new TNtuple("EvPre", " event ","run:set:flag:sum:singlet:triplet:late:latetime");
  TNtuple *ntEvent = new TNtuple("Event", " event ","run:set:flag:sum:singlet:triplet:late:latetime:wfsinglet:wfmin");
  TH1D* h_acceptance = new TH1D("h_acceptance","h_acceptance",files.size(),0,files.size());
  TH1D* h_Triplet = new TH1D("h_Triplet","h_Triplet",files.size(),0,files.size());
  TH1D* h_Dublet = new TH1D("h_Dublet","h_Dublet",files.size(),0,files.size());
  TH1D* h_Singlet = new TH1D("h_Singlet","h_Singlet",files.size(),0,files.size());
  TH1D* h_Total = new TH1D("h_Total","h_Total",files.size(),0,files.size());
  TH2D* h_wf_integral = new TH2D("h_wf_integral","h_wf_integral",files.size(),0,files.size(),400,-100,300);
  TH2D* h_wf_max = new TH2D("h_wf_max","h_wf_max",files.size(),0,files.size(),100,0,5);

  TString canName;
  canName.Form("c1-files-%i",maxFiles);
  TCanvas *c1 = new TCanvas(canName, canName, 1200, 900);
  c1->Divide(3, 2);
  canName.Form("c2-files-%i",maxFiles);
  TCanvas *c2 = new TCanvas(canName, canName, 1200, 500);
  c2->Divide(4);

  /*
  h_v1histo->SetDirectory(0);
  h_v1histobaseline->SetDirectory(0);
  h_v1histobaselineend->SetDirectory(0);
  
  h_acceptance->SetDirectory(0);
  h_max->SetDirectory(0);
  h_wf_integral->SetDirectory(0);
  */
    
  
  /////////////////////////////////////////////////////////////////////////////
  //run each file seperately
  int nFilesRead = 0;
  int runSet = -1;
  int nDigi = 0;
  for (int i = 0; i < int(files.size()); i++)
  {
    int rejectSaved = 0;
    filename = files[i];
    f2 = filename(filename.First("_")+1,filename.Length());
    f2 = f2(0, f2.Length() - f2.First("_") - 2);
    f2 = f2(0,f2.First("."));
    filenumber = f2.Atoi();
    runSet = getRunSet(filenumber);

    cout << "=============================" << endl;
    cout << "----" << i << " : " << filename;
    cout << " -- " << filenumber << " set " << runSet << endl;
    
    skipfile = false;
    if(runSet<0)
      skipfile = true;

    //if((filenumber>=3000) && (filenumber<=3020)) skipfile = false;
    /*if((filenumber>=5580) && (filenumber<=5600)) skipfile = false;
    if((filenumber>=20005) && (filenumber<=20020)) skipfile = false;
    if((filenumber>=20025) && (filenumber<=20040)) skipfile = false;
    if((filenumber>=20045) && (filenumber<=20060)) skipfile = false;
    if((filenumber>=20065) && (filenumber<=20080)) skipfile = false;
    if((filenumber>=20220) && (filenumber<=20222)) skipfile = false;
    */
                
    if(skipfile){
      cout << " -- " << filenumber << " skipped !" << endl;
      continue; 
    }

    if (nFilesRead  >= maxFiles)
      break;
    ++nFilesRead;

    TFile *f = new TFile(filename);
    //f->ls();
    //skip empty files
    if (f->GetNkeys() == 0){
      cout << " -- file empty" << endl;
      continue;
    }
    t1 = (TTree*)f->Get("pmtTree");
    //sometimes header corrupted
    if (t1 == 0){
      cout << " -- header corrupted" << endl;
      continue;
    }
    
    TPmtEvent* oneEvent = new TPmtEvent();
    t1->SetBranchAddress("pmtEvent",&oneEvent);

    nentries = t1->GetEntries();
    t1->GetEntry(0);
    // collect some info 
    nDigi = oneEvent->time.size();
    time_a = oneEvent->time.at(0)*1.E6;
    time_b = oneEvent->time.at(nDigi- 1)*1.E6;
    cout << "\t *********  starting run " << filename << " set  " << runSet << " entries =  " << nentries 
      << " ndigi " << nDigi << " time_a " << time_a << " time_b " << time_b << endl;
    if (runSet != currentSet) {
      currentSet = runSet;
      fout->cd();
      hWave[runSet] = new TH1D(Form("Wave-%i", runSet), Form("wave %i", runSet), nDigi, time_a, time_b);
    }
    if(!h_minimum) {
      h_minimum = new TH1D("h_minimum", "h_minimum", nDigi, 0, nDigi); // the minimum
      h_maximum = new TH1D("h_maximum", "h_maximum", nDigi, 0, nDigi); // the minimum
      h_mean = new TH1D("h_mean", "h_mean", nDigi, 0, nDigi);          // the minimum
      h_minimum->SetDirectory(0);
      h_maximum->SetDirectory(0);
      h_mean->SetDirectory(0);
    }
    /////////////////////////////////////////////////////////////////////////////
    //event loop
    nrungood = 0;
    TString hname;
    if(currentRun!=filenumber) {
      currentRun = filenumber;
      qualDir->cd();
      hname.Form("runquality-%i", filenumber);
      h_runquality = new TH1D(hname, hname, nentries, 0, nentries);
      hname.Form("runqualityHisto-%i", filenumber);
      h_runqualityhisto = new TH1D(hname, hname, 500, -100, 400);
      hname.Form("runqualityS-%i", filenumber);
      h_runqualityS = new TH1D(hname, hname, nentries, 0, nentries);
      hname.Form("runqualityHistoS-%i", filenumber);
      h_runqualityhistoS = new TH1D(hname, hname, 500, -100, 400);
      hname.Form("runqualityT-%i", filenumber);
      h_runqualityT = new TH1D(hname, hname, nentries, 0, nentries);
      hname.Form("runqualityHistoT-%i", filenumber);
      h_runqualityhistoT = new TH1D(hname, hname, 500, -100, 400);
      hname.Form("runquality2D-%i", filenumber);
      h_runquality2D = new TH2D(hname, hname, 500, -100, 400, 500, -100, 400);
      hname.Form("runqualityR-%i", filenumber);
      h_runqualityR = new TH1D(hname, hname, 200, 0, 100);
      hname.Form("runqualityB-%i", filenumber);
      h_runqualityB = new TH1D(hname, hname, 20, 0, 20);
      TH1D *h_minimumRun = (TH1D *)h_minimum->Clone(Form("MinimumRun%i",currentRun) );
      TH1D *h_maximumRun = (TH1D *)h_maximum->Clone(Form("MaximumRun%i",currentRun) );
      waveDir->Append(h_minimumRun);
      waveDir->Append(h_maximumRun);
      h_minimum->Reset("ICESM");
      h_maximum->Reset("ICESM");
      h_mean->Reset("ICESM");
      for (int k = 0; k < nDigi; k++)
      {
        h_minimum->SetBinContent(k, 5);
        h_maximum->SetBinContent(k, -5);
      }
      waveDir->cd();
      hname.Form("va%i", filenumber);
      h_v1 = new TH1D(hname, hname, nDigi, 0, nDigi);
      h_v1->Sumw2();
      hname.Form("vbFix%i", filenumber);
      h_v1Fix = new TH1D(hname, hname, nDigi, 0, nDigi);
      h_v1Fix->Sumw2();
      hname.Form("vb%i", filenumber);
      h_v1b = new TH1D(hname, hname, nDigi, 0, nDigi);
      h_v1b->Sumw2();
      hname.Form("vc%i", filenumber);
      h_v1c = new TH1D(hname, hname, nDigi, 0, nDigi);
      h_v1c->Sumw2();
      fout->cd();
    }

    //good_entries.clear();
    cout << " loop over " << nentries << endl;
    for (int j = 0; j < nentries; j++)
    {
      //if(j<2650) continue;
     
      t1->GetEntry(j);
      m = oneEvent->time.size();
      ntotal++;
      // reset cut histos 
      h_v1histo->Reset("ICESM");
      h_v1histobaseline->Reset("ICESM");
      h_v1histobaselineend->Reset("ICESM");

      
      ///loop ver wf, set v1, fill the histo for baseline etc
      if(j%1000==0) cout << " first loop   " << j << endl;
      for(int k = 0;k<m;k++){
        h_v1->SetBinContent(k,oneEvent->volt1.at(k));
        h_v1histo->Fill(oneEvent->volt1.at(k));
        if(wf_min > oneEvent->volt1.at(k)){
          wf_min = oneEvent->volt1.at(k);
        }
        if(k<900){
          h_v1histobaseline->Fill(oneEvent->volt1.at(k));
        }
        if(k>m-900){
          h_v1histobaselineend->Fill(oneEvent->volt1.at(k));
        }
      }
      // start cuts
      int event_code = rejectEvent();
      std::string mystring;
      if (j % 1000 == 0)
      {
        printf("... %i hex code %X ", j, event_code);
        mystring = rejectBits.to_string<char, std::string::traits_type, std::string::allocator_type>();
        std::cout << " bits  " << mystring << '\n';
      }
      // save some reject examples
      if (event_code > 0 && rejectSaved< 1)
      {
        mystring = rejectBits.to_string<char, std::string::traits_type, std::string::allocator_type>();
        TH1D *hTemp = (TH1D *)h_v1->Clone(Form("run-%i-ev-%i-bits-%s",filenumber,j,mystring.c_str()));
        rejectDir->Append(hTemp);
        ++rejectSaved;
        printf(" clone %s total %i \n", hTemp->GetName(), rejectSaved);
      }

      float evPre[8];
      double vmaxPre = 0;
      double maxTimePre = getMaxTime(h_v1, singletEnd, vmaxPre);
      // fill evvars
      evPre[0] = filenumber;
      evPre[1] = runSet;
      evPre[2] = event_code;
      evPre[3] = getSum(h_v1, singletStart, double(nDigi));
      evPre[4] = getSum(h_v1, singletStart , singletEnd);
      evPre[5] = getSum(h_v1, singletEnd , double(nDigi));
      evPre[6] = vmaxPre;
      evPre[7] = maxTimePre;
      ntEvPre->Fill(evPre);

      if (event_code > 0)
        continue;

      nrungood++;
      wf_min = 0;
      //done with cleaning
      //good_entries.push_back(j);

      wf_min = -1 * (wf_min - baseline_mean); //minimum of singlet
      wf_singlet = -1 * (h_v1->Integral(900, 1100) - 301 * baseline_mean) / wf_min;
      wf_width = 5;

      //////////////////////////
      //second loop
      if (j % 1000 == 0)
        cout << " second loop   " << j << " set  " << currentSet << endl;
      for (int k = 0; k < m; k++)
      {
        double l = -1 * (h_v1->GetBinContent(k) - baseline_mean);

        h_v1b->SetBinContent(k, l);
        h_v1c->SetBinContent(k, l / wf_min);

        if ((wf_singlet > 15) && (k > wf_width) && (k < m - wf_width))
        {
          double local_min = h_minimum->GetBinContent(k);
          double local_max = h_maximum->GetBinContent(k);
          double new_value = -1 * (h_v1->Integral(k - wf_width, k + wf_width) / (2 * wf_width + 1) - baseline_mean) / wf_min;
          if (new_value < local_min)
          {
            h_minimum->SetBinContent(k, new_value);
          }
          else if (new_value > local_max)
          {
            h_maximum->SetBinContent(k, new_value);
          }
        }

      }

      // align wavforms to time = 1000 ns and fix after pulsing
      afterFix(h_v1b, h_v1Fix);
      double vmax = 0;
      double maxTime = getMaxTime(h_v1Fix, singletEnd, vmax);
      // fill evvars
      evVars[EVEVENT] = float(j);
      evVars[EVRUN] = filenumber;
      evVars[EVSET] = runSet;
      evVars[EVFLAG] = event_code;
      evVars[EVSUM] = getSum(h_v1Fix, singletStart*microToNano , double(nDigi));
      evVars[EVSINGLET] = getSum(h_v1Fix, singletStart*microToNano, singletEnd*microToNano);
      evVars[EVTRIPLET] = getSum(h_v1Fix, singletEnd*microToNano , double(nDigi));
      evVars[EVLATE] = vmax;
      evVars[EVLATETIME] = maxTime;
      evVars[EVWFSINGLET] = wf_singlet;
      evVars[EVWFMIN] = wf_min;

      //printf(" ev %i run %i singlet %f triplet %f \n", int(evVars[EVEVENT]),int(evVars[EVRUN]), evVars[EVSINGLET] , evVars[EVTRIPLET]);
      
      ntEvent->Fill(evVars);

      if (evVars[EVSINGLET]<singleCut[runSet])
        continue;

      // sum waveforms by set
      for (int jbin = 0; jbin < m; ++jbin)
      {
        hWave[currentSet]->SetBinContent(jbin, hWave[currentSet]->GetBinContent(jbin) + h_v1Fix->GetBinContent(jbin));
      }

      if (wf_singlet <= 15)
          h_runqualityB->Fill(6);
      else

        {
          h_wf_integral->Fill(i, h_v1c->Integral());
          h_wf_max->Fill(i, h_v1b->GetMaximum());
          h_mean->Add(h_v1c);
        }

        h_runquality->Fill(j, h_v1c->Integral());
        h_runqualityhisto->Fill(h_v1c->Integral());
        h_runqualityS->Fill(j, h_v1c->Integral(900, 1100));
        h_runqualityhistoS->Fill(h_v1c->Integral(900, 1100));
        h_runqualityT->Fill(j, h_v1c->Integral(1600, 9900));
        h_runqualityhistoT->Fill(h_v1c->Integral(1600, 9900));
        h_runquality2D->Fill(h_v1c->Integral(900, 1100), h_v1c->Integral(1600, 9900));
        h_runqualityR->Fill(h_v1c->Integral(900, 1100) - 0.15 * h_v1c->Integral(1600, 99000));

        if ((nrungood % 100 == 0))
        {
          c1->cd(1);
          //h_v1->Draw();
          h_v1c->GetXaxis()->SetRangeUser(0.5E-6, 9E-6);
          h_v1c->Draw();
          //h_v1b->Rebin(10);
          h_v1c->SetLineColor(kBlack);
          h_v1c->GetYaxis()->SetRangeUser(-0.5, 2);
          h_v1b->Draw("SAME");
          h_v1b->SetLineColor(kRed);
          h_minimum->Draw("SAME");
          h_minimum->SetLineColor(kBlue);
          h_maximum->Draw("SAME");
          h_maximum->SetLineColor(kCyan);

          c1->cd(2);
          h_runquality->Draw("histo p *");
          h_runquality->SetMarkerColor(kBlack);
          h_runqualityS->Draw("SAME histo p *");
          h_runqualityS->SetMarkerColor(kRed);

          c1->cd(3);
          h_runqualityhistoS->Draw("histo");
          h_runqualityhistoS->SetLineColor(kRed);
          h_runqualityhistoS->GetXaxis()->SetRangeUser(-50, 100);
          h_runqualityhistoT->Draw("SAME histo");
          h_runqualityhistoT->SetLineColor(kGreen - 2);
          h_runqualityhisto->Draw("SAME histo");
          h_runqualityhisto->SetLineColor(kBlack);
          gPad->SetLogy();

          c1->cd(4);
          h_runqualityR->Draw("histo");
          h_runqualityR->SetLineColor(kBlack);
          gPad->SetLogy();

          c1->cd(5);
          h_runquality2D->Draw("COLZ");
          h_runquality2D->GetXaxis()->SetRangeUser(0, 100);
          h_runquality2D->GetYaxis()->SetRangeUser(-100, 100);

          c1->cd(6);
          h_runqualityB->Draw("histo ");
          h_runqualityB->SetLineColor(kBlack);
          gPad->SetLogy();

          c1->Update();

          cout << " -- entry " << j << " file " << filenumber << " ngood " << nrungood << " int " << h_v1c->Integral(900, 1100) << endl;
          //cout << " " << h_runqualityhistoS->GetMean() <<" " << h_runqualityhistoS->GetRMS() << endl;
          //cin.get();
        }

    } // end of event loop

    min_baseline=0;
    if(nrungood){
      h_mean->Scale(1.0/nrungood);  
      min_baseline=h_minimum->Integral(wf_width+1,900);
      min_baseline = min_baseline/(900-wf_width-1);
      for(int k = wf_width; k<h_mean->GetNbinsX();k++){
        h_minimum->SetBinContent(k, h_minimum->GetBinContent(k)-min_baseline);
      }
    }
   
    c1->cd(1);
    h_mean->Draw("histo");
    h_mean->SetLineColor(kGreen+2);
    h_mean->GetXaxis()->SetRangeUser(0.5E-6,3E-6);
    h_mean->GetYaxis()->SetRangeUser(-0.05,0.2);
    h_minimum->Draw("SAME");
    h_minimum->SetLineColor(kBlue);
    h_maximum->Draw("SAME");
    h_maximum->SetLineColor(kCyan);
    c1->Update();
    ////////////////////////////////////////////////////////////
    ntotalgood+=nrungood;    
    h_acceptance->Fill(i,nrungood*1.0/nentries);
    h_Triplet->Fill(i,h_mean->Integral(1600,9900) - h_minimum->Integral(1600,9900));
    h_Dublet->Fill(i,(h_mean->Integral(1100,1400) - h_minimum->Integral(1100,1400))*500/300);
    h_Singlet->Fill(i,h_mean->Integral(950,1100));
    //1.25 to cover the left out zone
    h_Total->Fill(i,h_mean->Integral(950,1100)
                    +h_mean->Integral(1600,9900) - h_minimum->Integral(1600,9900)
                    +(h_mean->Integral(1100,1400) - h_minimum->Integral(1100,1400))*500/300 
                    );
   
    c2->cd(1);
    h_wf_integral->Draw("COLZ"); 
    c2->cd(2);
    h_wf_max->Draw("COLZ"); 
    c2->cd(3);
    h_acceptance->Draw("hist p");
    h_acceptance->SetMarkerStyle(3);
    c2->cd(4);
    h_Total->Draw("hist p *");
    h_Total->SetMarkerColor(kBlack);
    h_Triplet->Draw("SAME hist p *");
    h_Triplet->SetMarkerColor(kRed);
    h_Dublet->Draw("SAME hist p *");
    h_Dublet->SetMarkerColor(kGreen+2);
    h_Singlet->Draw("SAME hist p *");
    h_Singlet->SetMarkerColor(kBlue);
    c2->Update();    
    
    cout << " -- triplet light " << h_mean->Integral(1600,9900)<<endl;
    cout << " -- triplet base " << h_minimum->Integral(1600,9900) << endl;
    cout << " -- good wfs: " << nrungood*1.0/nentries << endl;
    

    // fill summary
    sumVars[ERUN] = filenumber;
    sumVars[ESET] = runSet;
    sumVars[EACCEPT] = h_acceptance->GetBinContent(i);
    sumVars[EBASE] = h_v1histobaseline->GetMean();
    sumVars[EBASEEND] = h_v1histobaselineend->GetMean();
    sumVars[ESINGLET] = h_Singlet->GetBinContent(i);
    sumVars[EDOUBLET] = h_Dublet->GetBinContent(i);
    sumVars[ETRIPLET] = h_Triplet->GetBinContent(i);
    sumVars[ETOTAL] = h_Total->GetBinContent(i);
    sumVars[NTOTAL] = ntotal;
    sumVars[NGOOD] = nrungood;
    sumVars[OVERSHOOT] = h_minimum->Integral(1100, 10000);
    sumVars[EMINBASE] = min_baseline;
    ntSummary->Fill(sumVars);
    f->Close();
    //
    auto end = chrono::steady_clock::now();
    cout << " \t finished " << filenumber << " after " 
      << chrono::duration_cast<chrono::seconds>(end - start).count() 
      << " sec events " << nentries << " good " 
      << nrungood << endl;

    printf("\t  singlet %f doublet %f triplet %f \n ", sumVars[ESINGLET], sumVars[EDOUBLET], sumVars[ETRIPLET]);
    TH1D *h_minimumRun = (TH1D *)h_minimum->Clone(Form("MinimumRun%i", currentRun));
    TH1D *h_maximumRun = (TH1D *)h_maximum->Clone(Form("MaximumRun%i", currentRun));
    waveDir->Append(h_minimumRun);
    waveDir->Append(h_maximumRun);

    //if(i>10) break;
  } // end file loop

  cout << " Total Files " << nFilesRead << endl;
  cout << " Total Events " << ntotal << endl;
  cout << " Total Events accepted " << ntotalgood << endl;

  fout->Append(c1);
  fout->Append(c2);
  fout->ls();
  fout->Write();
  fout->Close();
  }
