#include <iostream>
#include <fstream>
#include <numeric>
#include "compiled/BaconAnalysis.hh"
// time is in microseconds
using namespace TMath;
double microToNano = 1.0E3;
enum { NSETS = 6};
int nDigi = 10000;
TH1D* hWave[NSETS];
TH1D* hWaveS[NSETS];
TH1D* hWaveB[NSETS];

int setTotal[NSETS];
std::vector<int> runGood;
std::vector<int> fileNumber;
std::vector<TH1D*> waveList;
std::vector<TH1D*> waveNormed;
std::vector<TH1D*> baseList;
std::vector<TH1D*> baseScaled;
std::vector<int> runList;
vector<double> singletSum;
vector<double> tripletSum;
vector<double> setNumber;

std::vector<double> wEvents;
std::vector<double> wMinimum;
double allMinimum;
double aveEvents;


TNtuple* ntRun;
TDirectory *wdir;
TDirectory *bdir;
TFile* fout;


vector<TF1*> baseFits;


int getRunSet(int irun)
{
  int iset = -1;
  if (irun >= 3000 && irun <= 5700)
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


void getHistStats(int i) {

  TH1D* h = waveList[i];
  double wevents = double(h->GetEntries())/double(h->GetNbinsX());
  double min = 1.0E9;
  for(int ibin = 1; ibin<h->GetNbinsX(); ++ibin ) {
    if( h->GetBinContent(ibin) == 0) continue;
    if( h->GetBinContent(ibin) < min) min= h->GetBinContent(ibin);
  }

  wEvents.push_back(wevents);
  wMinimum.push_back(min);

}



void histNorm(int i)
{
  wdir->cd();
  TString hname;
  hname.Form("waveNorm-%i",runList[i]);
  TH1D *hNorm = (TH1D*) waveList[i]->Clone(hname);
  double wevents =  double(waveList[i]->GetEntries())/double(waveList[i]->GetNbinsX());
  double wnorm  = aveEvents/wevents;

  hNorm->Reset();
  hNorm->Setw2();
  for (int ibin = 0; ibin < waveList[i]->GetNbinsX(); ++ibin)
  {
    double val = waveList[i]->GetBinContent(ibin) - allMinimum;
    hNorm->SetBinContent(ibin, val*wnorm);
    hNorm->SetBinError(ibin, sqrt(abs(waveList[i]->GetBinContent(ibin)))*wnorm);
  }

  waveNormed.push_back(hNorm);

  double singlet = hNorm->Integral(900,1100);
  double triplet = hNorm->Integral(1100,10000);
  int set = getRunSet(fileNumber[i]);
  //printf(" normed run %i set %i events %i fact %.3f singlet %.3E triplet %.3E \n", i,set, runGood[i], wnorm, singlet, triplet  );
  setNumber.push_back(double(set));
  singletSum.push_back(singlet);
  tripletSum.push_back(triplet);
  ntRun->Fill(float(fileNumber[i]),float(set),float(wevents),float(singlet),float(triplet));
  fout->cd();
}


TF1 *baseFit(int ifile)
{
  double start=3000.; double end =10000.;
  TH1D* hwave =  waveNormed[ifile];
  double binwidth = hwave->GetBinWidth(1);
  printf(" fitting to %s bw %f from %f to %f \n",hwave->GetName(), binwidth,start,end);
  TString fname; fname.Form("fbaseRun-%i", ifile);
  TF1 *fbase=NULL;
  fbase = new TF1(fname,"[0]*exp( -[1]*x)+[2]",start,end);
  fbase->SetNpx(1000); // numb points for function
  fbase->SetParNames("const", "invtime","off");
  fbase->SetParameters(1,0,0);
  //fbase->Print();

  TFitResultPtr fptr=hwave->Fit(fbase,"RSO+","",start,end);
  printf(" \t\t !!!! fit to %s status %i !!! \n", hwave->GetName(),fbase->GetObjectStat() ); 

  fptr->Print();
  
  if(fbase->GetObjectStat()) {
    printf(" fit to %s failed \n",hwave->GetName() ); 
    return fbase;
  }

  new TCanvas(fname,fname);
  gStyle->SetOptFit(1111);
  hwave->Draw("");
  fbase->Draw("same");
  return fbase;
}



TH1D* reScale(int i)
{
  bdir->cd();
  TString hname;
  if(baseList.size()<1) return NULL;
  hname.Form("BaselineRun%i-Scaled",i);

  TH1D *hScale = (TH1D*) baseList[i]->Clone(hname);
  double bevents = double(baseList[i]->GetEntries());
  double wevents = double(waveList[i]->GetEntries())/double(waveList[i]->GetNbinsX());

  double normFact=  wevents/bevents;
  for (int ibin = 0; ibin < baseList[i]->GetNbinsX(); ++ibin)
  {
    double val = baseList[i]->GetBinContent(ibin);
    hScale->SetBinContent(ibin, val * normFact);
    //hScale->SetBinError(ibin, normFact * sqrt(val));
  }
  printf(" \t norm factor run %i file %i ngood %i for %s base int %.f norm %.2f scaled %.f  \n",
      i, fileNumber[i],  int(wevents) ,  hScale->GetName(),waveList[i]->Integral() ,normFact, hScale->Integral() );
  fout->cd();
  return hScale;
}




double peakInRange(TH1D* h, double xmin, double xmax){
  int ilow = h->FindBin(xmin);
  int ihigh = h->FindBin(xmax);
  double val=0;
  double xpeak=0;
  for(int i=ilow; i<= ihigh; ++i) if(h->GetBinContent(i)>val) {
    val = h->GetBinContent(i);
    xpeak = h->GetBinLowEdge(i);
  }
  return xpeak;
}

void weightedMean(  std::vector<double> val, std::vector<double> err, double &mean, double& emean) 
{

  mean =0;
  double w =0;
  for(unsigned i=0; i< val.size(); ++i ) {
    double wi = 1./pow(err[i] ,2.);
    mean += wi*val[i];
    w += wi;
  }
  mean /= w;
  emean = 1./sqrt(w);
}

void ratioE(double x, double xe, double y, double ye, double& r,double& re){
  r = x/y;
  re = r * sqrt( pow(xe/x,2.)+pow(ye/y,2.) );
}

void getStats(std::vector<double> v, double& mean, double& stdev)
{

  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  mean = sum / v.size();

  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(),
      std::bind2nd(std::minus<double>(), mean));
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  stdev = std::sqrt(sq_sum / v.size());
}

void baselineFix()
{
  
  TString cname;
  TCanvas *can;
 
    
  TFile *inFile = new TFile("BaconAnalysisHighCut-1000.root", "READONLY");
  if (!inFile)
    return;

  TString sumString = setSumNames();
  printf("sumvars string %s \n", sumString.Data());



  TTree* tSummary=NULL;
  for(int iset=0; iset<NSETS; ++iset ) {
   setTotal[iset]=0;
  }


  inFile->GetObject("Summary",tSummary);
  for(int iv=0; iv< SUMVARS; ++iv ) tSummary->SetBranchAddress(sumNames[iv],&sumVars[iv]);


  int setRuns[NSETS]={0,0,0,0,0,0};

  // all entries sum events
  for (Long64_t jent=0;jent<tSummary->GetEntries() ;jent ++) {
    tSummary->GetEntry(jent);
    int jset = int(sumVars[ESET]);
    //printf(" %lli %i %i \n", jent, jset ,int(sumVars[NGOOD]) );
    setTotal[jset] +=  int(sumVars[NGOOD]);  // fix me
    ++setRuns[jset];
    fileNumber.push_back(int(sumVars[ERUN]));
    runGood.push_back(int(sumVars[NGOOD]));
    //printf(" ... file %i  total runs %i good runs  %i set %i set total %i  \n ", 
    //    fileNumber[runGood.size()-1] , int(sumVars[NTOTAL]), runGood[runGood.size()-1], jset, setTotal[jset] );
  }

  printf(" set totals  \n");
  for(int iset=0; iset<NSETS; ++iset ) {
    printf(" set %i  runs in set %i good events  %i \n ", iset, setRuns[iset], setTotal[iset]);
  }

  // end of  get summary
  //


  fout = new TFile("baselineFix.root", "RECREATE");
  ntRun = new TNtuple("ntRun","run","run:set:nev:singlet:triplet");
  TDirectory* waveDir=NULL;
  wdir= fout->mkdir("waveDir");
  bdir= fout->mkdir("baseDir");
  inFile->GetObject("waveForms",waveDir);
  TList *list = waveDir->GetListOfKeys();

  TIter next(list);
  TKey *key;
  TString theName("blah");
  TH1D* hErr;

  while (TKey *key = (TKey *)next()) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *) key->ReadObj();
    TString hname(h->GetName());
    if(hname==theName) continue;
    theName=hname;
    if(hname.Contains("Fix")){
      wdir->Append(h);
      waveList.push_back(h);
      TString snum = (TString) hname(5,hname.Length());
      runList.push_back(atoi(snum.Data()));

    }
    /*
    else if(hname.Contains("MinimumRun")) {
      hErr=fillErrors(h);
      bdir->Append(hErr);
      baseList.push_back(hErr);
    }
    */
  }

  fout->cd();

  printf(" num waveList %lu num baseList %lu \n", waveList.size(), baseList.size());

  
  for(int i=0; i<waveList.size(); ++i) getHistStats(i);


  allMinimum = *std::min_element(wMinimum.begin(), wMinimum.end());


  aveEvents = std::accumulate(wEvents.begin(), wEvents.end(), 0)/double(wEvents.size());

  printf(" ave events %f all min %f \n", aveEvents,allMinimum);


  for(int i=0; i<waveList.size(); ++i)  histNorm(i);

  
  //baseList[0]->Print("all");
  //return;
  //
  TF1 *fit =NULL;

  //for(unsigned ifile=0; ifile< waveList.size(); ++ifile) {
  for(int ifile=0; ifile<1; ++ifile) {
    fit =  baseFit(ifile);
    if(fit->GetObjectStat()) break;
    fit->Print();
    for (int ii = 0; ii < 3; ++ii)
    {
      printf("\t  param %i %s %.3f +/- %.3f \n", ii, fit->GetParName(ii), fit->GetParameter(ii), fit->GetParError(ii));
    }
    baseFits.push_back(fit);
  }

  //unsigned filesInList = Min(waveList.size(),baseList.size());
  //for(unsigned ih=0; ih<waveList.size() ; ++ih) printf( " %s wave %E   \n",   
  //    waveList[ih]->GetName(), waveList[ih]->GetEntries()/1.E4 );


  for(int iset=0; iset<NSETS; ++iset) {
    hWave[iset] = new TH1D(Form("Wave%i", iset), Form("wave %i", iset), nDigi,0,10);
    hWave[iset]->Sumw2();
  }

  /*rescale baseline
  for(int i=0; i<baseList.size(); ++i) {
    baseScaled.push_back(reScale(i));
  }
  */
  // norm
  


  fout->Write();


  TMultiGraph *mgTriplet = new TMultiGraph();

   TGraph *gSinglet = new TGraphErrors(setNumber.size()-1, &setNumber[0], &singletSum[0]);
   TGraph *gTriplet = new TGraphErrors(setNumber.size()-1, &setNumber[0], &tripletSum[0]);
   TCanvas *cTriplet0 = new TCanvas("set-sums","set-sums");
   cTriplet0->SetGridx(); cTriplet0->SetGridy();
   gSinglet->SetName("singlet");
   gSinglet->SetTitle("singlet");
   gSinglet->SetMarkerColor(kRed);
   gSinglet->SetMarkerStyle(22);
   gSinglet->SetMarkerSize(1.4);
   gSinglet->GetHistogram()->GetXaxis()->SetTitle(" set");
   gSinglet->GetHistogram()->GetYaxis()->SetTitle(" singlet yield");

   gTriplet->SetName("triplet");
   gTriplet->SetTitle("triplet");
   gTriplet->SetMarkerColor(kBlack);
   gTriplet->SetMarkerStyle(22);
   gTriplet->SetMarkerSize(1.4);
   gTriplet->GetHistogram()->GetXaxis()->SetTitle(" set");
   gTriplet->GetHistogram()->GetYaxis()->SetTitle(" triplet yield");

   mgTriplet->Add(gTriplet,"p");
   mgTriplet->Add(gSinglet,"p");
   mgTriplet->SetTitle("Light Yield; set ; yield");
   mgTriplet->Draw("a");
   cTriplet0->Print(".png");
   cTriplet0->BuildLegend();



  int theSet = -1;
  vector<int> goodRunList;
  vector<int> goodRunNumber;
  vector<int> badRunList;
  double singletCut[NSETS]={0,0,0,0,0,0};
  double tripletCut[NSETS]={0,80.E3,120.E3,120.E3,150.E3,150.E3};
  bool first=false;
  can=NULL;
  for(unsigned i=0; i< waveList.size(); ++i) {
    int iset = getRunSet(runList[i]);
    if(iset!=theSet) {
      theSet=iset;
      if(can) can->Print(".png");
      cname.Form("WavesForSet-%i",iset);
      can = new TCanvas(cname,cname);
      first=true;
    }
    //printf(" %i file %i set %i singlet %E\n",i,fileNumber[i],iset,singletSum[i]);
    if(singletSum[i]<singletCut[iset]|| tripletSum[i]<tripletCut[iset]) {
      badRunList.push_back(runList[i]);
      continue;
    } else 
      goodRunList.push_back(runList[i]);
      goodRunNumber.push_back(i);

    if(first) {
      waveNormed[i]->Draw("");
      first=false;
    }
    else 
      waveNormed[i]->Draw("same");
  }
  can->Print(".png");

 

   printf(" RUN SUMMARY: good %lu bad %lu \n", goodRunList.size(), badRunList.size());

  
  double setTotal[NSETS]={0,0,0,0,0,0};
  for(int iset=0; iset<NSETS; ++iset) setRuns[iset]=0;



  for(unsigned igood=0; igood< goodRunList.size(); ++igood) {
    unsigned irun = goodRunList[igood];
    unsigned i = goodRunNumber[igood];

    // sum waveforms by set
    int iset = getRunSet(irun);
    if(iset<0) continue;
    
    double wevents =  double(waveList[i]->GetEntries())/double(waveList[i]->GetNbinsX());
    setTotal[iset]+= wevents;
    setRuns[iset]++;
    for (int jbin = 0; jbin < hWave[iset]->GetNbinsX() ; ++jbin) {
      hWave[iset]->SetBinContent(jbin, hWave[iset]->GetBinContent(jbin) + waveList[i]->GetBinContent(jbin) - allMinimum);
    }
  }



  // normalize  sets
  for(int iset=0; iset<NSETS; ++iset) {

    for (int jbin = 0; jbin < hWave[iset]->GetNbinsX() ; ++jbin) {
      hWave[iset]->SetBinContent(jbin, hWave[iset]->GetBinContent(jbin)/setTotal[iset]);
      hWave[iset]->SetBinError(jbin,sqrt(abs(hWave[iset]->GetBinContent(jbin)))/setTotal[iset]);
    }
  }


  for(int i=0; i<NSETS; ++i) printf(" set %i set runs %i total events %.0f \n",i,setRuns[i],setTotal[i]);


  cname.Form("SumWaveSets");
  can = new TCanvas(cname,cname);

  int setColor[NSETS];
  setColor[0] = kOrange + 4;
  setColor[1] = kBlue;
  setColor[2] = kRed;
  setColor[3] = kGreen;
  setColor[4] = kCyan;
  setColor[5] = kMagenta + 1;
  for(int iset=0; iset<NSETS; ++iset) {
    hWave[iset]->SetLineColor(setColor[iset]);
    hWave[iset]->SetMarkerColor(setColor[iset]);
    hWave[iset]->SetMarkerStyle(7);
  }


  hWave[0]->GetYaxis()->SetRangeUser(0.001,0.5);
  hWave[0]->Draw(0);
  for(int iset=1; iset<NSETS; ++iset) {
    hWave[iset]->Draw("same");
  }

  can->BuildLegend();
  can->Print(".png");



}

