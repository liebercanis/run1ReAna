#include <iostream>
#include <fstream>
#include <numeric>
#include "compiled/BaconAnalysis.hh"


// time is in microseconds
using namespace TMath;
double microToNano = 1.0E3;
enum { NSETS = 6};
enum { NPARS = 7};

int ngood[NSETS];
int nall[NSETS];
TH1D* hSingleSet[NSETS];
TH1D* hTripleSet[NSETS];
TH1D* hLate[NSETS];
TH1D* hLateTime[NSETS];
TH1D* hSingleSetN[NSETS];
TH1D* hTripleSetN[NSETS];


TH1D* hSingletSetPass[NSETS];
TH1D* hTripletSetPass[NSETS];



int setColor[NSETS];
int setStyle[NSETS];
std::vector<std::vector<int>> runsForSet;

TGraph *gModel[3];
TGraphErrors *gModelNorm[3];
TString modelGraphName[3];
TDirectory *wdir;
TFile *inFile;
TFile* fout;
double singletSetCut[NSETS];
double tripletSetCut[NSETS];

double singletPass[NSETS];
double tripletPass[NSETS];
double singletPassErr[NSETS];
double tripletPassErr[NSETS];


enum
{
  NBITS =5
};

std::bitset<NBITS> rejectBits;

std::vector<int> runList;
std::vector<TH1D*> waveList;
std::vector<TH1D*> waveNormed;
std::vector<double> wEvents;
std::vector<double> wMinimum;
double allMinimum;
double aveEvents;
std::vector<int> fileNumber;
vector<double> runNumber;
vector<double> singletRun;
vector<double> tripletRun;

vector<double> setNumber;
vector<double> singletSum;
vector<double> tripletSum;


TString cname;
TCanvas *can;
TNtuple *ntRun;

TH1D* rescale(TH1D* h,int nevents ,TString hname)
{
  TH1D *hScale = (TH1D*) h->Clone(hname);
  double normFact = 1.0 / double(nevents);
  for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin)
  {
    double val = h->GetBinContent(ibin);
    hScale->SetBinContent(ibin, val * normFact);
    //hScale->SetBinError(ibin, normFact * sqrt(val));
  }
  printf(" \t norm factor set for %s set n=%i integral %f \n", hScale->GetName(), nevents ,hScale->Integral());
  return hScale;
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


void teventAna()
{
   
  int pcolor[NSETS];
  pcolor[5] = kBlack;
  pcolor[4] = kBlue;
  pcolor[3] = kYellow - 2;
  pcolor[2] = kGreen;
  pcolor[1] = kMagenta + 2;
  pcolor[0] = kRed;


  double PPM[NSETS]={0,1,2,5,10,10};

  double singleCut[NSETS] = {4.,4.,4.,4.,4.,4.};
  for (int i = 0; i < NSETS; ++i)
  {
    setColor[i] = i + 1;
    setStyle[i] = 20 + i;
  }
  setStyle[4] = 3;
  setStyle[5] = 29;

  setColor[0] = kOrange + 4;
  setColor[1] = kBlue;
  setColor[2] = kRed;
  setColor[3] = kGreen - 2;
  setColor[4] = kCyan - 2;
  setColor[5] = kMagenta + 1;

  printf(" single cuts are :\n");
  for (int iset = 0; iset < NSETS; ++iset)
    printf(" \t set %i cut %.1f \n", iset, singleCut[iset]);


  inFile = new TFile("BaconAnalysisHighCut-1000.root", "READONLY");
  if (!inFile)
    return;

  //inFile->ls();
  

  // get ntuples 

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


  fout = new TFile("teventAna.root", "RECREATE");
  
  TString hname;
  TString htitle;
  for(int iset=0; iset<NSETS; ++iset) {
    hSingleSetN[iset] = NULL;
    hTripleSetN[iset] = NULL;
    hLate[iset] = NULL;
    hLateTime[iset] =  NULL;
    hSingleSet[iset] =  NULL;
    hTripleSet[iset] =  NULL;
    
    hname.Form("Singlet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Singlet yield set %i, %i PPM",iset,int(PPM[iset]));
    hSingleSet[iset] = new TH1D(hname,htitle,1000,0.,30000.);
    hSingleSet[iset]->SetLineColor(pcolor[iset]);

    hname.Form("Triplet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Triplet yield set %i, %i PPM",iset,int(PPM[iset]));
    hTripleSet[iset] = new TH1D(hname,htitle,1000,0.,250000.);
    hTripleSet[iset]->SetLineColor(pcolor[iset]);
  }

  vector<vector<double>> single;
  vector<vector<double>> triple;
  vector<double> tripleRun;
  vector<double> aveTripleRun;
  vector<double> errTripleRun;
  vector<double> singleRun;
  vector<double> aveSingleRun;
  vector<double> errSingleRun;

  vector<double> vRun;
  vector<double> vRunErr;
  vector<int> vSet;



  single.resize(NSETS);
  triple.resize(NSETS);

  vector<vector<double>> singleSet;
  vector<vector<double>> singleSetErr;
  vector<vector<double>> tripleSet;
  vector<vector<double>> tripleSetErr;

  singleSet.resize(NSETS);
  singleSetErr.resize(NSETS);
  tripleSet.resize(NSETS);
  tripleSetErr.resize(NSETS);

  for (Long64_t jent = 0; jent < tEvent->GetEntries(); jent++)
  {
    tEvent->GetEntry(jent);
    int iset = evVars[EVSET];
    if (evVars[EVFLAG] == 0 )
    {
      ++nall[iset];
      hSingleSet[iset]->Fill(evVars[EVSINGLET]);
    } 
  }
  int thisRun = -1;
  for (Long64_t jent = 0; jent < tEvent->GetEntries(); jent++)
  {
    tEvent->GetEntry(jent);
    int iset = evVars[EVSET];
    if (evVars[EVRUN]<2500)
      continue;
    int run = evVars[EVRUN];
    if(run<4000) run = run - 3000;
    else if(run<6000) run = run - 5000;
    else if(run< 20100) run = run - 20000 + 700;
      else break;
      if (run < 0)
         run = evVars[EVRUN] - 3000;
      if (run != thisRun)
      {
        if (thisRun > 0)
        {
          vRun.push_back(double(thisRun));
          vRunErr.push_back(0.0);
          vSet.push_back(iset);
          double ave, err;
          getStats(tripleRun, ave, err);
          aveTripleRun.push_back(ave);
          errTripleRun.push_back(err);
          tripleSet[iset].push_back(ave);
          tripleSetErr[iset].push_back(err);
          tripleRun.clear();
          //printf("run %i size %lu triplet ave %f  err %f   ",int(evVars[EVRUN]), tripleRun.size(), ave , err );
          getStats(singleRun, ave, err);
          aveSingleRun.push_back(ave);
          errSingleRun.push_back(err);
          singleSet[iset].push_back(ave);
          singleSetErr[iset].push_back(err);
          //printf("single  ave %f  err %f  \n ",  ave , err );
          singleRun.clear();
        }
        thisRun = run;
      }
      bool pass = evVars[EVSINGLET] > singleCut[iset];
      if (evVars[EVFLAG] == 0 && pass)
      {
        ++ngood[iset];
        single[iset].push_back(evVars[EVSINGLET]);
        triple[iset].push_back(evVars[EVTRIPLET]);
        tripleRun.push_back(evVars[EVTRIPLET]);
        singleRun.push_back(evVars[EVSINGLET]);
        hSingleSet[iset]->Fill(evVars[EVSINGLET]);
        hTripleSet[iset]->Fill(evVars[EVTRIPLET]);
        //hRunTriplet->Fill(run, evVars[ETRIPLET]);
      } 
    }



    // make plots
  for (int iset = 0; iset < NSETS; ++iset) {
    if (!hSingleSet[iset])
      continue;
    TString newName;
    newName.Form("%s-norm", hSingleSet[iset]->GetName());
    printf(" iset %i %s nall %i \n", iset, newName.Data(),  nall[iset]);
    hSingleSetN[iset] = rescale(hSingleSet[iset], nall[iset], newName);
    hSingleSetN[iset]->GetYaxis()->SetTitle(" frequency/event");
    hSingleSetN[iset]->GetXaxis()->SetTitle("sum singlet light yield in event");
    if(iset==NSETS-1) {
      hSingleSetN[iset]->SetFillColor(kBlack);
      hSingleSetN[iset]->SetFillStyle(3001);
    }
  }


  TCanvas *canSingleSet = new TCanvas("canSingleSetNormalized", "canSingleSetNormalize");
  hSingleSetN[0]->Draw("");
  for (int iset = 1; iset < NSETS; ++iset){
    if(!hSingleSetN[iset]) continue;
    hSingleSetN[iset]->Draw("sames");
  }
  canSingleSet->BuildLegend();
  canSingleSet->Print(".png");



  for (int iset = 0; iset < NSETS; ++iset) {
    if(!hTripleSet[iset])
      continue;
    TString newName; newName.Form("%s-norm",hTripleSet[iset]->GetName());
    hTripleSetN[iset] = rescale( hTripleSet[iset],ngood[iset],newName);
    hTripleSetN[iset]->GetYaxis()->SetTitle(" frequency/event");
    hTripleSetN[iset]->GetXaxis()->SetTitle("sum triplet light yield");
    if(iset==NSETS-1) {
      hTripleSetN[iset]->SetFillStyle(3003);
      hTripleSetN[iset]->SetFillColor(kBlack);
    }
  }


  TCanvas *canTripleSet = new TCanvas("canTripleSetNormalized", "canTripleSetNormalized");
  hTripleSetN[0]->Draw("");
  for (int iset = 1; iset < NSETS; ++iset) {
    if(!hTripleSetN[iset]) continue;
    hTripleSetN[iset]->Draw("sames");
  }
  canTripleSet->BuildLegend();
  canTripleSet->Print(".png");




  fout->Write();
}
