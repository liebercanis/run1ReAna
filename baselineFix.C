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

TDirectory *wdir;
TFile* fout;


TF1 *tailFit(int irun, TH1D *hTail)
{
  double start=4000.; double end =10000.;
  // fit exponential tail
  double binwidth = hTail->GetBinWidth(1);
  TString fname; fname.Form("ftailEv-%i", irun);
  TF1 *ftail = new TF1(fname,"[0]+x*[1]",start,end);
  ftail->SetNpx(1000); // numb points for function
  ftail->SetParNames("const", "slope");
  ftail->SetParameters(0,0);
  ftail->Print();
  printf(" fit to %s bw %f from %f to %f \n",hTail->GetName(), binwidth,start,end);
  hTail->Fit(ftail,"RLO+","",start,end);
  new TCanvas(fname,fname);
  gStyle->SetOptFit();
  hTail->Draw("");
  ftail->Draw("same");
  return ftail;
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


TH1D* rescale(TH1D* h,int nevents ,TString hname)
{
  TH1D *hScale = (TH1D*) h->Clone(hname);
  double normFact = 1.0 / double(nevents);
  for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
  {
    double val = h->GetBinContent(ibin);
    hScale->SetBinContent(ibin, val * normFact);
    //hScale->SetBinError(ibin, normFact * sqrt(val));
  }
  printf(" \t norm factor set for %s set n=%i integral %f \n", hScale->GetName(), nevents ,hScale->Integral());
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
   
  std::vector<TH1D*> hlist;
  std::vector<int> runList;
  
  TFile *inFile = new TFile("BaconAnalysisHighCutAll-10000.root", "READONLY");
  if (!inFile)
    return;
 

  fout = new TFile("baselineFix.root", "RECREATE");
  TDirectory* waveDir=NULL;
  wdir= fout->mkdir("wavedir");
  inFile->GetObject("waveForms",waveDir);
  TList *list = waveDir->GetListOfKeys();

  TIter next(list);
  TKey *key;
  while (TKey *key = (TKey *)next()) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *) key->ReadObj();
    TString hname(h->GetName());
    if(hname.Contains("Fix")) {
      ///h->Rebin(5);
      wdir->Append(h);
      hlist.push_back(h);
      TString snum = (TString) hname(5,hname.Length());
      runList.push_back(atoi(snum.Data()));
    }
  }

  fout->cd();

  
  for(int iset=0; iset<NSETS; ++iset) {
    hWave[iset] = new TH1D(Form("Wave-%i", iset), Form("wave %i", iset), nDigi,0,10);
  }


  for(unsigned i=0; i< hlist.size(); ++i) {
    TString hname = hlist[i]->GetName();
    //printf(" %i %s %i \n", i,hlist[i]->GetName() , runList[i] );
    if(i<10) tailFit(i,hlist[i]);
    // sum set
    // sum waveforms by set
    int iset = getRunSet(runList[i]);
    if(iset>-1) 
      for (int jbin = 0; jbin < hWave[iset]->GetNbinsX() ; ++jbin)
      {
        hWave[iset]->SetBinContent(jbin, hWave[iset]->GetBinContent(jbin) + hlist[i]->GetBinContent(jbin));
      }
  }

  /*

     float_t sumVars[SUMVARS];
     int setTotal[NSETS];
     int setRuns[NSETS];

     TDirectory *waveForms;
  inFile->GetObject(waveForms, hInWave[i]);



  for (int i = 0; i < NSETS; ++i)
  {
    hInWave[i] = NULL;
    inFile->GetObject(Form("Wave-%i", i), hInWave[i]);
    if (hInWave[i])
    {
      cout << " wave hist " << hInWave[i]->GetName();
      fout->Append(hInWave[i]);
    }
  }
*/
  
    fout->Write();
  }

