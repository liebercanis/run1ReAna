#include <iostream>
#include <fstream>
#include <numeric>
#include "compiled/BaconAnalysis.hh"
// time is in microseconds
using namespace TMath;
double microToNano = 1.0E3;
enum { NSETS = 6};
int ngood[NSETS];
TH1D *hInWave[NSETS];
TH1D* hWave[NSETS];
int setColor[NSETS];
int setStyle[NSETS];

TH1D* rescale(TH1D* h,int iset)
{
  TH1D *hScale = (TH1D*) h->Clone(Form("NormWave-%i",iset));
  printf(" \t\t norm factor set %s %i \n", h->GetName(), ngood[iset]);
  double normFact = 1.0 / double(ngood[iset]);
  for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
  {
    double val = h->GetBinContent(ibin);
    hScale->SetBinContent(ibin, val * normFact);
    //hScale->SetBinError(ibin, normFact * sqrt(val));
  }
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

void eventAna()
{
  int pcolor[NSETS];
  pcolor[5] = kBlack;
  pcolor[4] = kBlue;
  pcolor[3] = kYellow - 2;
  pcolor[2] = kGreen;
  pcolor[1] = kMagenta + 2;
  pcolor[0] = kRed;

  double PPM[NSETS]={0,1,2,5,10,10};

  double singleCut[NSETS] = {1.,4.,4.,4.,4.,4.};
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

  TFile *inFile = new TFile("BaconAnalysisHighCut-1000.root", "READONLY");
  //TFile *inFile = new TFile("BaconAnalysis-646.root", "READONLY");
  if (!inFile)
    return;
  TFile *fout = new TFile("eventAna.root", "RECREATE");
  //
  TH2D* hRunTriplet= new TH2D("RunTriplet","run triplet",200,0,200,300,-100,200);
  //
  TH1D* hLate[NSETS];
  TH1D* hLateTime[NSETS];
  TH1D* hSingleSet[NSETS];
  TH1D* hTripleSet[NSETS];
  TString hname;
  TString htitle;
  for(int iset=0; iset<NSETS; ++iset) {
    hname.Form("Late%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Late set %i, %i PPM",iset,int(PPM[iset]));
    hLate[iset] = new TH1D(hname,htitle,40,0.,.2);
    hLate[iset]->SetLineColor(pcolor[iset]);
    hname.Form("LateTime%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("LateTime set %i, %i PPM",iset,int(PPM[iset]));
    hLateTime[iset] = new TH1D(hname,htitle,200,0.,10000.);
    hLateTime[iset]->SetLineColor(pcolor[iset]);

    hname.Form("Singlet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Singlet yield set %i, %i PPM",iset,int(PPM[iset]));
    hSingleSet[iset] = new TH1D(hname,htitle,310,-1.,31.);
    hSingleSet[iset]->SetLineColor(pcolor[iset]);

    hname.Form("Triplet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Triplet yield set %i, %i PPM",iset,int(PPM[iset]));
    hTripleSet[iset] = new TH1D(hname,htitle,400,-100.,300.);
    hTripleSet[iset]->SetLineColor(pcolor[iset]);

  }

  Float_t sumVars[SUMVARS];
  int setTotal[NSETS];
  int setRuns[NSETS];
  TTree* tSummary=NULL;


  for(int i=0; i<NSETS; ++i) {
    hInWave[i]=NULL;
    inFile->GetObject(Form("Wave-%i",i),hInWave[i]);
    if(hInWave[i]) {
      cout << " wave hist " << hInWave[i]->GetName();
      fout->Append(hInWave[i]);
    }
  }


  inFile->GetObject("Summary", tSummary);
  if (tSummary)
    {
      //TNtuple *ntSummary = new TNtuple("Summary", " Summary ", "run:set:base:baseend:accept:total:singlet:dublet:triplet:ngood");
      TString sumNames[SUMVARS];
      sumNames[0] = TString("run");
      sumNames[1] = TString("set");
      sumNames[2] = TString("base");
      sumNames[3] = TString("baseend");
      sumNames[4] = TString("accept");
      sumNames[5] = TString("total");
      sumNames[6] = TString("singlet");
      sumNames[7] = TString("dublet");
      sumNames[8] = TString("triplet");
      sumNames[9] = TString("ngood");

      for (int iv = 0; iv < SUMVARS; ++iv)
        tSummary->SetBranchAddress(sumNames[iv], &sumVars[iv]);

      for (int iset = 0; iset < NSETS; ++iset)
      {
        setTotal[iset] = 0;
        setRuns[iset] = 0;
      }

      // ratios

      if (tSummary)
        printf(" tSummary has  %lld  entries \n", tSummary->GetEntries());
      else
        printf(" no tSummary \n ");

      // all entries and fill the histograms
      for (Long64_t jent = 0; jent < tSummary->GetEntries(); jent++)
      {
        tSummary->GetEntry(jent);
        //printf(" tSummary entry %lld \n ",jent);
        for (int iv = 0; iv < SUMVARS; ++iv)
        {
          int jset = int(sumVars[ESET]);
          int nset = int(sumVars[NGOOD]);
          setTotal[jset] += int(sumVars[NGOOD]);
          ++setRuns[jset];
          //printf(" \t %i %s %f \n",iv, sumNames[iv].Data(), sumVars[iv]);
        }
      }

      for (int iset = 0; iset < NSETS; ++iset)
        printf(" set %i runs %i total events %i \n ", iset, setRuns[iset], setTotal[iset]);
    }

    TTree *tPre = NULL;
    inFile->GetObject("EvPre", tPre);
    enum
    {
      NPRE = 8
    };
    float preVars[NPRE];
    TString preNames[NPRE];
    preNames[0] = TString("run");
    preNames[1] = TString("set");
    preNames[2] = TString("flag");
    preNames[3] = TString("sum");
    preNames[4] = TString("singlet");
    preNames[5] = TString("triplet");
    preNames[6] = TString("late");
    preNames[7] = TString("latetime");

    for (int iv = 0; iv < NPRE; ++iv)
      tPre->SetBranchAddress(preNames[iv], &preVars[iv]);

    vector<double> prePass;
    prePass.resize(NPRE);
    vector<double> preNum;
    preNum.resize(NPRE);

    for (Long64_t jent = 0; jent < tPre->GetEntries(); jent++)
    {
      tPre->GetEntry(jent);
      int iset = preVars[1];
      int flag = int(preVars[2]);
      preNum[iset] = preNum[iset] + 1.;
      if (flag == 0)
        prePass[iset] = prePass[iset] + 1.;
    }

    vector<double> preFrac;
    preFrac.resize(NPRE);

    for (unsigned iset = 0; iset < preFrac.size(); ++iset)
    {
      preFrac[iset] = prePass[iset] / preNum[iset];
    }

    vector<double> vset;
    for (unsigned iset = 0; iset < preFrac.size(); ++iset)
      vset.push_back(double(iset));

    TGraph *gPass = new TGraph(NSETS, &vset[0], &preFrac[0]);

    TCanvas *canpre = new TCanvas("prePass", "prePass");
    canpre->SetGridy();
    gPass->SetName("prePass");
    gPass->SetTitle("pass fraction ");
    gPass->SetMarkerColor(kRed);
    gPass->SetMarkerStyle(22);
    gPass->SetMarkerSize(1.4);
    gPass->GetXaxis()->SetTitle("set");
    gPass->GetYaxis()->SetTitle("pass fraction");
    gPass->GetYaxis()->SetRangeUser(0., 1.);

    gPass->Draw("ap");

    TTree *tEvent = NULL;
    inFile->GetObject("Event", tEvent);
    float evVars[EVVARS];
    TString evNames[EVVARS];
    evNames[0] = TString("run");
    evNames[1] = TString("set");
    evNames[2] = TString("flag");
    evNames[3] = TString("sum");
    evNames[4] = TString("singlet");
    evNames[5] = TString("triplet");
    evNames[6] = TString("late");
    evNames[7] = TString("latetime");
    evNames[8] = TString("wfsinglet");
    evNames[9] = TString("wfmin");

    for (int iv = 0; iv < EVVARS; ++iv)
      tEvent->SetBranchAddress(evNames[iv], &evVars[iv]);

    printf(" total events %lld \n", tEvent->GetEntries());

    // all entries and fill the histograms
    double PPMError[NSETS];
    for (int iset = 0; iset < NSETS; ++iset)
    {
      ngood[iset] = 0;
      PPMError[iset] = 0;
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

    int thisRun = -1;
    for (Long64_t jent = 0; jent < tEvent->GetEntries(); jent++)
    {
      tEvent->GetEntry(jent);
      int iset = evVars[EVSET];
      int run = evVars[EVRUN] - 20000 + 30;
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

          //printf("run %i size %lu triplet ave %f  err %f   ",thisRun, tripleRun.size(), ave , err );
          tripleRun.clear();
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
      hLate[iset]->Fill(evVars[EVLATE]);
      hLateTime[iset]->Fill(evVars[EVLATETIME]);
      bool pass = -1.0 * evVars[EVSINGLET] > singleCut[iset];
      for (int iv = 0; iv < EVVARS; ++iv)
      {
        if (evVars[EVFLAG] == 0 && pass)
        {
          ++ngood[iset];
          single[iset].push_back(evVars[EVSINGLET]);
          triple[iset].push_back(evVars[EVTRIPLET]);
          tripleRun.push_back(-1.0 * evVars[EVTRIPLET]);
          singleRun.push_back(-1.0 * evVars[EVSINGLET]);
        }
      }
      hSingleSet[iset]->Fill(-1.0 * evVars[EVSINGLET]);
      hTripleSet[iset]->Fill(-1.0 * evVars[EVTRIPLET]);
      hRunTriplet->Fill(run, -1.0 * evVars[ETRIPLET]);
    }

    //for(unsigned j=0; j<vRun.size();++j) printf(" %i %i %f \n",j,int(vRun[j]),aveTripleRun[j]);

    printf(" nruns %lu %lu \n", vRun.size(), aveTripleRun.size());

    double smeanSet[NSETS];
    double smeanSetErr[NSETS];
    double tmeanSet[NSETS];
    double tmeanSetErr[NSETS];

    double setVector[NSETS] = {0, 1, 2, 3, 4, 5};
    double setVectorErr[NSETS] = {0, 0, 0, 0, 0, 0};

    for (int iset = 0; iset < NSETS; ++iset)
    {

      weightedMean(singleSet[iset], singleSetErr[iset], smeanSet[iset], smeanSetErr[iset]);
      weightedMean(tripleSet[iset], tripleSetErr[iset], tmeanSet[iset], tmeanSetErr[iset]);
    }

    TGraphErrors *gSingleSet = new TGraphErrors(NSETS, setVector, smeanSet, setVectorErr, smeanSetErr);
    TGraphErrors *gTripleSet = new TGraphErrors(NSETS, setVector, tmeanSet, setVectorErr, tmeanSetErr);

    TGraphErrors *gTripleRun = new TGraphErrors(vRun.size(), &vRun[0], &aveTripleRun[0], &vRunErr[0], &errTripleRun[0]);
    TCanvas *canRun = new TCanvas("triletRun", "tripletRun");
    gTripleRun->SetName("gTripleRun");
    gTripleRun->SetTitle("Triple by Run");
    gTripleRun->SetMarkerColor(kRed);
    gTripleRun->SetMarkerStyle(22);
    gTripleRun->SetMarkerSize(1.4);
    gTripleRun->GetXaxis()->SetTitle("run");
    gTripleRun->GetYaxis()->SetTitle("triplet yield");
    gTripleRun->Draw("ap");
    canRun->Print(".png");

    TCanvas *canSet = new TCanvas("bySet", "bySet");
    gSingleSet->SetName("gSingleSet");
    gSingleSet->SetTitle("single by set");
    gSingleSet->SetMarkerColor(kRed);
    gSingleSet->SetMarkerStyle(22);
    gSingleSet->SetMarkerSize(1.4);
    gSingleSet->GetXaxis()->SetTitle("set");
    gSingleSet->GetYaxis()->SetTitle("yield");
    gSingleSet->GetYaxis()->SetRangeUser(0, 100);

    gTripleSet->SetName("gTripleSet");
    gTripleSet->SetTitle("Triple by set");
    gTripleSet->SetMarkerColor(kBlue);
    gTripleSet->SetMarkerStyle(23);
    gTripleSet->SetMarkerSize(1.4);
    gTripleSet->GetXaxis()->SetTitle("set");
    gTripleSet->GetYaxis()->SetTitle("yield");
    gTripleSet->GetYaxis()->SetRangeUser(0, 100);
    gSingleSet->Draw("ap");
    gTripleSet->Draw("psame");
    canSet->Print(".png");

    TGraphErrors *gSingleRun = new TGraphErrors(vRun.size(), &vRun[0], &aveSingleRun[0], &vRunErr[0], &errSingleRun[0]);
    TCanvas *canSingletRun = new TCanvas("singleRun", "singleRun");
    gSingleRun->SetName("gSingleRun");
    gSingleRun->SetTitle("Single by Run");
    gSingleRun->SetMarkerColor(kRed);
    gSingleRun->SetMarkerStyle(22);
    gSingleRun->SetMarkerSize(1.4);
    gSingleRun->GetXaxis()->SetTitle("run");
    gSingleRun->GetYaxis()->SetTitle("singlet yield");
    gSingleRun->Draw("ap");
    canSingletRun->Print(".png");

    hSingleSet[5]->SetFillColor(kBlack);
    hSingleSet[5]->SetFillStyle(3001);

    TCanvas *canSingleSet = new TCanvas("canSingleSet", "canSingleSet");
    hSingleSet[0]->Draw("");
    for (int iset = 1; iset < NSETS; ++iset)
      hSingleSet[iset]->Draw("sames");
    canSingleSet->BuildLegend();
    canSingleSet->Print(".png");

    TCanvas *canTripleSet = new TCanvas("canTripleSet", "canTripleSet");
    hTripleSet[0]->Draw("");
    for (int iset = 1; iset < NSETS; ++iset)
      hTripleSet[iset]->Draw("sames");
    canTripleSet->BuildLegend();
    canTripleSet->Print(".png");

    double singlePeak1[NSETS];
    double singlePeak2[NSETS];

    for (int iset = 0; iset < NSETS; ++iset)
    {
      //singlePeak1[iset]=hSingleSet[iset]->GetBinLowEdge(hSingleSet[iset]->GetMaximumBin());
      double xlow = 0.0;
      double xhigh = 5.0;
      if (iset == 0)
      {
        xlow = -1.0;
        xhigh = 2.0;
      }
      singlePeak1[iset] = peakInRange(hSingleSet[iset], xlow, xhigh);
      singlePeak2[iset] = peakInRange(hSingleSet[iset], xhigh, 20);
      printf("set %i peak value  %f second %f \n", iset, singlePeak1[iset], singlePeak2[iset]);
    }

    printf(" means by set\n");
    for (int iset = 0; iset < NSETS; ++iset)
      printf(" set %i singlet %f +/- %f triplet %f +/- %f \n ", iset,
             smeanSet[iset], smeanSetErr[iset], tmeanSet[iset], tmeanSetErr[iset]);

    //normalize
    for (int i = 0; i < NSETS; ++i)
    {
      if (hInWave[i]) hWave[i] = rescale(hInWave[i],i);
    }
 
    for (int i = 0; i < NSETS; ++i)
    {
      if (hWave[i])
      {
        int istart = hWave[i]->FindBin(singletStart);
        int iend = hWave[i]->FindBin(singletEnd);
        if (i == 0)
          printf("single start %f bin %i end %f bin %i last %i \n", singletStart, istart, singletEnd, iend, hWave[i]->GetNbinsX());
        cout << " wave hist " << hWave[i]->GetName() << " nev " << ngood[i];
        cout << " singlet " << hWave[i]->Integral(istart, iend);
        cout << " triplet " << hWave[i]->Integral(iend, hWave[i]->GetNbinsX()) << endl;
      }
    }

    TCanvas *canLifeAll = new TCanvas("WaveALL", "WaveALL");
    //canLifeAll->SetLogy();
    gStyle->SetOptStat(0);
    for (int iset = 0; iset < NSETS; ++iset)
    {
      hWave[iset]->SetTitle(Form("set %i  ", iset));
      hWave[iset]->GetYaxis()->SetTitle("yield SPE/40 ns");
      hWave[iset]->GetYaxis()->SetRangeUser(-1.0E-3, 0.07);
      hWave[iset]->GetXaxis()->SetRangeUser(0.95, 1.1);
      hWave[iset]->SetMarkerStyle(setStyle[iset]);
      hWave[iset]->SetMarkerColor(setColor[iset]);
      hWave[iset]->SetLineColor(setColor[iset]);
      hWave[iset]->SetMarkerSize(0.7);
      if (iset == 0)
        hWave[iset]->Draw("");
      else
        hWave[iset]->Draw("same");
    }
    canLifeAll->BuildLegend();
    canLifeAll->Print(".png");

    return;

    double smean[NSETS];
    double srms[NSETS];
    double tmean[NSETS];
    double trms[NSETS];

    for (unsigned iset = 0; iset < single.size(); ++iset)
    {
      getStats(single[iset], smean[iset], srms[iset]);
      getStats(triple[iset], tmean[iset], trms[iset]);

      printf(" set %u %lu  %f %f %lu %f %f \n", iset, single[iset].size(), smean[iset], srms[iset], triple[iset].size(), tmean[iset], trms[iset]);
    }

    //ratios
    double smeanRatio[NSETS];
    double srmsRatio[NSETS];
    double tmeanRatio[NSETS];
    double trmsRatio[NSETS];
    for (unsigned iset = 0; iset < single.size(); ++iset)
    {
      ratioE(smean[iset], srms[iset], smean[0], srms[0], smeanRatio[iset], srmsRatio[iset]);
      ratioE(tmean[iset], trms[iset], tmean[0], trms[0], tmeanRatio[iset], trmsRatio[iset]);
    }

    TGraphErrors *gSinglet = new TGraphErrors(NSETS, PPM, smeanRatio, PPMError, srmsRatio);
    gSinglet->SetTitle("good singlet");

    TCanvas *cSinglet = new TCanvas("Singlet", "Singlet");
    cSinglet->SetGridx();
    cSinglet->SetGridy();
    gSinglet->SetName("SingletGraph");
    gSinglet->SetTitle("Singlet");
    gSinglet->SetMarkerColor(kRed);
    gSinglet->SetMarkerStyle(22);
    gSinglet->SetMarkerSize(1.4);
    gSinglet->GetXaxis()->SetTitle("dopant PPM");
    gSinglet->GetYaxis()->SetTitle("yield");
    gSinglet->GetYaxis()->SetRangeUser(0, 10);
    gSinglet->Draw("ap");
    //cSinglet->BuildLegend();
    cSinglet->Print(".png");

    TGraphErrors *gTriplet = new TGraphErrors(NSETS, PPM, tmeanRatio, PPMError, trmsRatio);
    gTriplet->SetTitle("good triplet");

    TCanvas *cTriplet = new TCanvas("Triplet", "Triplet");
    cTriplet->SetGridx();
    cTriplet->SetGridy();
    gTriplet->SetName("TripletGraph");
    gTriplet->SetTitle("Triplet");
    gTriplet->SetMarkerColor(kRed);
    gTriplet->SetMarkerStyle(22);
    gTriplet->SetMarkerSize(1.4);
    gTriplet->GetXaxis()->SetTitle("dopant PPM");
    gTriplet->GetYaxis()->SetTitle("yield");
    gTriplet->GetYaxis()->SetRangeUser(0, 10);
    gTriplet->Draw("ap");
    //cTriplet->BuildLegend();
    cTriplet->Print(".png");

    //
    fout->Write();
  }
