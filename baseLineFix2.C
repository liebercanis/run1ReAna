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



TDirectory *wdir;
TDirectory *bdir;
TFile* fout;

TNtuple* ntRun;

vector<TF1*> baseFits;

TF1 *baseFit(int ifile)
{
  double start=4000.; double end =10000.;
  // fit exponential tail
  TH1D* hwave =  waveList[ifile];
  double binwidth = hwave->GetBinWidth(1);
  TString fname; fname.Form("fbaseRun-%i", ifile);
  TF1 *fbase = new TF1(fname,"[0]*exp( -[1]*x)+[2]",start,end);
  fbase->SetNpx(1000); // numb points for function
  fbase->SetParNames("const", "time","off");
  fbase->SetParameters(1,0,0);
  fbase->Print();
  printf(" fit to %s bw %f from %f to %f \n",hwave->GetName(), binwidth,start,end);
  hwave->Fit(fbase,"RLO+","",start,end);
  new TCanvas(fname,fname);
  gStyle->SetOptFit();
  hwave->Draw("");
  fbase->Draw("same");
  return fbase;
}



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


TH1D* reScale(int i)
{
  bdir->cd();
  TString hname;
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



void histNorm(int i)
{
  bdir->cd();
  TString hname;
  hname.Form("WaveRun%i-normed",i);
  double wevents = double(waveList[i]->GetEntries())/double(waveList[i]->GetNbinsX());
  TH1D *hNorm = (TH1D*) waveList[i]->Clone(hname);
  for (int ibin = 0; ibin < baseList[i]->GetNbinsX(); ++ibin)
  {
    double val = waveList[i]->GetBinContent(ibin);
    hNorm->SetBinContent(ibin, val/wevents );
  }
  waveNormed.push_back(hNorm);
  double singlet = hNorm->Integral(900,1100);
  double triplet = hNorm->Integral(1100,10000);
  int set = getRunSet(fileNumber[i]);
  printf(" normed run %i set %i events %i fact %.0f singlet %.3E triplet %.3E \n", i,set, runGood[i], wevents, singlet, triplet  );
  singletSum.push_back(singlet);
  tripletSum.push_back(triplet);
  ntRun->Fill(float(fileNumber[i]),float(set),float(wevents),float(singlet),float(triplet));
  fout->cd();
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


  // all entries sum events
  for (Long64_t jent=0;jent<tSummary->GetEntries() ;jent ++) {
    tSummary->GetEntry(jent);
    int jset = int(sumVars[ESET]);
    printf(" %lli %i %i \n", jent, jset ,int(sumVars[NGOOD]) );
    setTotal[jset] +=  int(sumVars[NGOOD]);  // fix me
    fileNumber.push_back(int(sumVars[ERUN]));
    runGood.push_back(int(sumVars[NGOOD]));
    printf(" ... file %i  total runs %i good runs  %i set %i set total %i  \n ", 
        fileNumber[runGood.size()-1] , int(sumVars[NTOTAL]), runGood[runGood.size()-1], jset, setTotal[jset] );
  }

  printf(" set totals  \n");
  for(int iset=0; iset<NSETS; ++iset ) {
    printf(" set %i  good events  %i \n ", iset, setTotal[iset]);
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
    else if(hname.Contains("MinimumRun")) {
      bdir->Append(h);
      baseList.push_back(h);
    }
  }

  fout->cd();

  printf(" num waveList %lu num baseList %lu \n", waveList.size(), baseList.size());

  unsigned filesInList = Min(waveList.size(),baseList.size());
  for(unsigned ih=0; ih<filesInList ; ++ih) printf( " %s wave %E %s base %E  \n",   
      waveList[ih]->GetName(), waveList[ih]->GetEntries()/1.E4 ,baseList[ih]->GetName(), baseList[ih]->GetEntries()  );


  for(int iset=0; iset<NSETS; ++iset) {
    hWave[iset] = new TH1D(Form("Wave-%i", iset), Form("wave %i", iset), nDigi,0,10);
    hWaveS[iset] = new TH1D(Form("WaveSubtracted-%i", iset), Form("wave %i", iset), nDigi,0,10);
    hWaveB[iset] = new TH1D(Form("NormalizedBaseline-%i", iset), Form("wave %i", iset), nDigi,0,10);
  }

  //baseList[0]->Print("all");
  //return;
  //
  TF1 *fit =NULL;

  for(unsigned ifile=0; ifile< waveList.size(); ++ifile) {
  //for(int ifile=0; ifile<2; ++ifile) {
    fit =  baseFit(ifile);
    if(!fit->GetObjectStat()) break;
    baseFits.push_back(fit);
  }

  return;


  //rescale baseline
  for(int i=0; i<baseList.size(); ++i) {
       baseScaled.push_back(reScale(i));
  }
  // norm
  for(int i=0; i<baseList.size(); ++i) {
       histNorm(i);
  }



  printf(" waves %lu base %lu baseScale %lu \n",waveList.size(),baseList.size(),baseScaled.size());

  for(unsigned i=0; i< waveList.size(); ++i) {
    TString hname = waveList[i]->GetName();
    //printf(" %i %s %i \n", i,waveList[i]->GetName() , runList[i] );
    //if(i<10) tailFit(i,waveList[i]);
    // sum set
    // sum waveforms by set
    int iset = getRunSet(runList[i]);
    if(iset>-1) 
      for (int jbin = 0; jbin < hWave[iset]->GetNbinsX() ; ++jbin)
      {
        hWave[iset]->SetBinContent(jbin, hWave[iset]->GetBinContent(jbin) + waveList[i]->GetBinContent(jbin));
        hWaveS[iset]->SetBinContent(jbin, hWaveS[iset]->GetBinContent(jbin) + waveList[i]->GetBinContent(jbin) - baseScaled[i]->GetBinContent(jbin));
        hWaveB[iset]->SetBinContent(jbin, hWaveB[iset]->GetBinContent(jbin) + baseScaled[i]->GetBinContent(jbin));
      }
  }

  fout->Write();




  TCanvas *can;
  TString cname;
  for(int iset=0; iset<NSETS; ++iset) {
    cname.Form("Set-%i",iset);
    can = new TCanvas(cname,cname);
    //can->Divide(1,2);
    hWave[iset]->SetLineColor(kBlack);
    hWaveS[iset]->SetLineColor(kGreen);
    hWaveB[iset]->SetLineColor(kRed);
    //can->cd(1);
    hWave[iset]->Draw();
    hWaveS[iset]->Draw("same");
    //can->cd(2);
    hWaveB[iset]->Draw("same");
  }


  
  int theSet = -1;
  vector<int> goodRunList;
  vector<int> badRunList;
  double singletCut[NSETS]={8,10,10,10,10,10};
  double tripletCut[NSETS]={50,50,50,50,55,55};
  bool first=false;
  for(unsigned i=0; i< waveList.size(); ++i) {
    int iset = getRunSet(runList[i]);
    if(iset!=theSet) {
      theSet=iset;
      cname.Form("WavesForSet-%i",iset);
      can = new TCanvas(cname,cname);
      first=true;
    }
    printf(" %i file %i set %i singlet %E\n",i,fileNumber[i],iset,singletSum[i]);
    if(singletSum[i]<singletCut[iset]|| tripletSum[i]<tripletCut[iset]) {
      badRunList.push_back(runList[i]);
      continue;
    } else 
      goodRunList.push_back(runList[i]);

    if(first) {
      waveNormed[i]->Draw("");
      first=false;
    }
    else 
      waveNormed[i]->Draw("same");
  }

  printf(" RUN SUMMARY: good %lu bad %lu \n", goodRunList.size(), badRunList.size());

}

