#include <iostream>
#include <fstream>
#include <numeric>
#include "compiled/BaconAnalysis.hh"


// time is in microseconds
using namespace TMath;
double microToNano = 1.0E3;
enum { NSETS = 6};
int ngood[NSETS];
int nall[NSETS];
TH1D *hInWave[NSETS];
TH1D* hWave[NSETS];
int setColor[NSETS];
int setStyle[NSETS];
double tripletEnd = 2.5;
TGraph *gModel[3];
TGraphErrors *gModelNorm[3];
TString modelGraphName[3];
TFile* fout;
enum
{
  NBITS =5
};
std::bitset<NBITS> rejectBits;

void getModelGraphs(TString modelFile)
{
  TFile *fmodel = new TFile(modelFile);
  cout << " \n \t getting model graphs from file " << modelFile << "\n" << endl;
  modelGraphName[0]=TString("RangeSingletModel");
  modelGraphName[1]=TString("RangeTripletModel");
  modelGraphName[2]=TString("TotalTripletModel");
  for(int i=0; i<3; ++i) {
    fmodel->GetObject(modelGraphName[i], gModel[i]);
    cout << "model graph "<< i << " " << gModel[i]->GetName() << endl;
    fout->Append(gModel[i]);
  }


}

void normModelGraphs(double* smean, double* smeanSetErr)
{
  // normalize to PPM zero in data
  double xval0,yval0;
  gModel[0]->GetPoint(0,xval0,yval0);
  double norm = smean[0]/yval0;

  for(int i=0; i<gModel[0]->GetN(); ++i) printf("\t point %i smean %f \n",i,smean[i]);
  for(int i=0; i<3; ++i) 
  {
    //gModelNorm[i] =  (TGraphErrors* gModel[i]->Clone(modelGraphName[i]+TString("-norm"));
    gModelNorm[i] =  new TGraphErrors( gModel[i]->GetN());
    gModelNorm[i]->SetName(Form("%s-norm",gModel[i]->GetName()));
    gModelNorm[i]->SetTitle(Form("%s",gModel[i]->GetName()));
    cout << "NNNNN norm model graph "<< i << " " << gModelNorm[i]->GetName() << " npoints " << gModel[i]->GetN()  << " norm " << norm << endl;

    for(int j=0; j< gModel[i]->GetN(); ++j) {
      double xval,yval;
      norm = smean[0]/yval0;
      if(i>0) norm = smean[j]/yval0;
      printf(" model %i point %i norm %f err %f \n",i,j,norm,smeanSetErr[j]);
      gModel[i]->GetPoint(j,xval,yval);
      gModelNorm[i]->SetPoint(j,xval,norm*yval);
      if(i>0) gModelNorm[i]->SetPointError(j,0.0,smeanSetErr[j]);
    }
    fout->Append(gModelNorm[i]);
  }
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

void eventAna()
{

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
  sumNames[8] = TString("ntot");
  sumNames[10] = TString("ngood");
  sumNames[11] = TString("over");
  sumNames[12] = TString("minbase");


  TString evNames[EVVARS];
  evNames[0] = TString("ev");
  evNames[1] = TString("run");
  evNames[2] = TString("set");
  evNames[3] = TString("flag");
  evNames[4] = TString("sum");
  evNames[5] = TString("singlet");
  evNames[6] = TString("triplet");
  evNames[7] = TString("late");
  evNames[8] = TString("latetime");
  evNames[9] = TString("wfsinglet");
  evNames[10] = TString("wfmin");


   
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

  TFile *inFile = new TFile("BaconAnalysisHighCut-1000.root", "READONLY");
  if (!inFile)
    return;

  fout = new TFile("eventAna.root", "RECREATE");

  TString modelFileName; modelFileName.Form("model-%.2f.root",tripletEnd);
  getModelGraphs(modelFileName);

  int setTotal[NSETS];
  int setGood[NSETS];
  int setRuns[NSETS];

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

  
  //
  TH2D* qRunTriplet= new TH2D("RunTriplet","run triplet",200,0,200,300,-100,200);
  //
  TH1D* hLate[NSETS];
  TH1D* hLateTime[NSETS];
  TH1D* hSingleSet[NSETS];
  TH1D* hTripleSet[NSETS];
  TH1D* hSingleSetN[NSETS];
  TH1D* hTripleSetN[NSETS];
  
  TString hname;
  TString htitle;
  for(int iset=0; iset<NSETS; ++iset) {
    hSingleSetN[iset] = NULL;
    hTripleSetN[iset] = NULL;
    hLate[iset] = NULL;
    hLateTime[iset] =  NULL;
    hSingleSet[iset] =  NULL;
    hTripleSet[iset] =  NULL;
    if (!hInWave[iset])
      continue;
    hname.Form("Late%i-%i-PPM", iset, int(PPM[iset]));
    htitle.Form("Late set %i, %i PPM",iset,int(PPM[iset]));
    hLate[iset] = new TH1D(hname,htitle,40,0.,.2);
    hLate[iset]->SetLineColor(pcolor[iset]);
    hname.Form("LateTime%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("LateTime set %i, %i PPM",iset,int(PPM[iset]));
    hLateTime[iset] = new TH1D(hname,htitle,200,0.,10000.);
    hLateTime[iset]->SetLineColor(pcolor[iset]);

    hname.Form("Singlet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Singlet yield set %i, %i PPM",iset,int(PPM[iset]));
    hSingleSet[iset] = new TH1D(hname,htitle,600,-10.,210.);
    hSingleSet[iset]->SetLineColor(pcolor[iset]);

    hname.Form("Triplet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Triplet yield set %i, %i PPM",iset,int(PPM[iset]));
    hTripleSet[iset] = new TH1D(hname,htitle,400,-100.,500.);
    hTripleSet[iset]->SetLineColor(pcolor[iset]);
  }


  TTree *tSummary = NULL;
  inFile->GetObject("Summary", tSummary);
  float_t sumVars[SUMVARS];
  if (tSummary)
    {
      printf(" sumVars \n");
      //TNtuple *ntSummary = new TNtuple("Summary", " Summary ", "run:set:base:baseend:accept:total:singlet:dublet:triplet:ngood");
    
      for (int iv = 0; iv < SUMVARS; ++iv)
        tSummary->SetBranchAddress(sumNames[iv], &sumVars[iv]);

      for (int iset = 0; iset < NSETS; ++iset)
      {
        setTotal[iset] = 0;
        setGood[iset] = 0;
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
        int jset = int(sumVars[ESET]);
        setTotal[jset] += int(sumVars[NTOT]);
        setGood[jset] += int(sumVars[NGOOD]);
        ++setRuns[jset];
        //printf(" \t %i %s %f \n",iv, sumNames[iv].Data(), sumVars[iv]);
      }

      for (int iset = 0; iset < NSETS; ++iset)
        printf(" set %i runs %i total events %i  good %i  \n ", iset, setRuns[iset], setTotal[iset],setGood[iset]);
    }

    fout->ls();
    return;

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

    printf(" preVars \n");

    for (int iv = 0; iv < NPRE; ++iv)
      if(tPre) tPre->SetBranchAddress(preNames[iv], &preVars[iv]);

    vector<double> prePass;
    prePass.resize(NPRE);
    vector<double> preNum;
    preNum.resize(NPRE);


    double passCount[NBITS][NSETS];
    for(int ib=0; ib<NBITS; ++ib) 
      for(int iset=0; iset<NSETS; ++iset) 
        passCount[ib][iset]=0;


    for (Long64_t jent = 0; jent < tPre->GetEntries(); jent++)
    {
      tPre->GetEntry(jent);
      int run = preVars[0];
      int iset = preVars[1];
      int flag = int(preVars[2]);
      rejectBits=flag;
      /* check
      printf(" jent %llu run %i flag %i ( ",jent,run,flag);
      for(int ib=0; ib<NBITS; ++ib) printf(" %i ",rejectBits.test(ib));
      printf(" ) \n");
      */
      for(int ib=0; ib<NBITS; ++ib) if(!rejectBits.test(ib)) passCount[ib][iset] = passCount[ib][iset]+1.;

      preNum[iset] = preNum[iset] + 1.;
      if (flag == 0)
        prePass[iset] = prePass[iset] + 1.;
    }


    vector<double> preFrac;
    preFrac.resize(NPRE);

    for (unsigned iset = 0; iset < preFrac.size(); ++iset)
    {
      if( preNum[iset]>0) preFrac[iset] = prePass[iset] / preNum[iset];
    }

    for(int ib=0; ib<NBITS; ++ib) {
      for (unsigned iset = 0; iset < NSETS; ++iset) {
        if( preNum[iset]>0) passCount[ib][iset] = passCount[ib][iset] / preNum[iset];
      }
    }


    vector<double> vset;
    for (unsigned iset = 0; iset < preFrac.size(); ++iset){
      vset.push_back(double(iset));
    }

    TGraph *gPass = new TGraph(NSETS, &vset[0], &preFrac[0]);

    TGraph *gPassBit[NBITS];
    for(int ib=0; ib<NBITS; ++ib)  gPassBit[ib] = new TGraph(NSETS,&vset[0], &passCount[ib][0]);


    TString gTitle;
    TMultiGraph* mgPass = new TMultiGraph();

    TCanvas *canpre = new TCanvas("prePass", "prePass");
    canpre->SetGridy();
    gPass->SetName("prePass");
    gPass->SetTitle("pass fraction ");
    gPass->SetMarkerColor(kBlack);
    gPass->SetMarkerStyle(22);
    gPass->SetMarkerSize(1.4);
    gPass->GetXaxis()->SetTitle("set");
    gPass->GetYaxis()->SetTitle("pass fraction");
    gPass->GetYaxis()->SetRangeUser(0., 1.);

    mgPass->Add(gPass);
    for(int ib=0; ib<NBITS; ++ib)  {
      gPassBit[ib]->SetName(Form("bitPass-%i",ib));
      gPassBit[ib]->SetTitle(Form("bit %i",ib));
      gPassBit[ib]->SetMarkerColor(pcolor[ib]);
      gPassBit[ib]->SetMarkerStyle(20+ib);
      gPassBit[ib]->GetYaxis()->SetRangeUser(0., 1.);
      mgPass->Add(gPassBit[ib]);
      fout->Add(gPassBit[ib]);
    }

    gTitle.Form("Cut Pass Fraction ; run set; pass fraction  ");
    mgPass->SetTitle(gTitle);
    mgPass->Draw("ap");
    gPass->GetYaxis()->SetRangeUser(0., 1.);
    gPad->BuildLegend();


    TTree *tEvent = NULL;
    inFile->GetObject("Event", tEvent);
    float evVars[EVVARS];
    
    printf(" evVars \n");

    for (int iv = 0; iv < EVVARS; ++iv)
      if(tEvent) tEvent->SetBranchAddress(evNames[iv], &evVars[iv]);

    printf(" total events %lld \n", tEvent->GetEntries());

    // all entries and fill the histograms
    double PPMError[NSETS];
    for (int iset = 0; iset < NSETS; ++iset)
    {
      ngood[iset] = 0;
      nall[iset] = 0;
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
      hLate[iset]->Fill(evVars[EVLATE]);
      hLateTime[iset]->Fill(evVars[EVLATETIME]);
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
      if(singleSet[iset].size()==0) continue;
      weightedMean(singleSet[iset], singleSetErr[iset], smeanSet[iset], smeanSetErr[iset]);
      weightedMean(tripleSet[iset], tripleSetErr[iset], tmeanSet[iset], tmeanSetErr[iset]);
    }

  
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


    for (int iset = 0; iset < NSETS; ++iset) {
      if (!hSingleSet[iset])
        continue;
      TString newName;
      newName.Form("%s-norm", hSingleSet[iset]->GetName());
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
      printf("set %i singlet peak value first %f second %f \n", iset, singlePeak1[iset], singlePeak2[iset]);
    }

    printf(" means by set\n");
    for (int iset = 0; iset < NSETS; ++iset)
      printf(" set %i singlet %f +/- %f triplet %f +/- %f \n ", iset,
             smeanSet[iset], smeanSetErr[iset], tmeanSet[iset], tmeanSetErr[iset]);

    //normalize
    for (int i = 0; i < NSETS; ++i)
    {
      if (hInWave[i]) {
        TString newName; newName.Form("%s-norm",hInWave[i]->GetName());
        hWave[i] = rescale(hInWave[i],ngood[i],newName);
      }
    }

    
    double tripletIntegral[NSETS];
    double tripletIntegralError[NSETS];

    for (int i = 0; i < NSETS; ++i)
    {
      if (hWave[i])
      {
        int istart = hWave[i]->FindBin(singletStart);
        int iend = hWave[i]->FindBin(singletEnd);
        int itrip =  hWave[i]->FindBin(tripletEnd);
        if (i == 0)
          printf("single start %f bin %i end %f bin %i last %i \n", singletStart, istart, singletEnd, iend, hWave[i]->GetNbinsX());
        cout << " wave hist " << hWave[i]->GetName() << " nev " << ngood[i];
        cout << " singlet " << hWave[i]->Integral(istart, iend);
        double tripInt =  hWave[i]->Integral(iend,itrip);
        tripletIntegral[i]=tripInt;
        tripletIntegralError[i]=sqrt(tripInt);
        cout << " triplet " << tripInt << endl;
      }
        cout << " " <<  endl;
    }

  
    TCanvas *canLifeAll = new TCanvas("NormalizedWaveALL", "Normalized Wave All");
    //canLifeAll->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    for (int iset = 0; iset < NSETS; ++iset)
    {
      if(! hWave[iset]) continue;
      if( hWave[iset]->GetEntries()==0) continue;
      hWave[iset]->SetTitle(Form("set %i  ", iset));
      hWave[iset]->GetYaxis()->SetTitle("yield SPE/40 ns");
      hWave[iset]->GetYaxis()->SetRangeUser(1.E-5, 0.65);
      hWave[iset]->GetXaxis()->SetRangeUser(0.95, 6.5);
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


     
    TGraphErrors *gSingleSet = new TGraphErrors(NSETS-1, PPM, smeanSet, setVectorErr, smeanSetErr);
    TGraphErrors *gTripleSet = new TGraphErrors(NSETS-1, PPM, tmeanSet, setVectorErr, tmeanSetErr);
    TGraphErrors *gTripleInt = new TGraphErrors(NSETS-1, PPM,tripletIntegral, PPMError,tripletIntegralError);

    TGraphErrors *gTripleRun = new TGraphErrors(vRun.size(), &vRun[0], &aveTripleRun[0], &vRunErr[0], &errTripleRun[0]);
    TCanvas *canRun = new TCanvas("tripletRun", "tripletRun");
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
    canSet->SetTickx();
    canSet->SetTicky();
    canSet->SetGridx();
    canSet->SetGridy();

    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gSingleSet->SetName("gSingleSet");

    gSingleSet->SetMarkerColor(kRed);
    gSingleSet->SetMarkerStyle(22);
    gSingleSet->SetMarkerSize(1.4);
    gSingleSet->SetTitle(Form("Singlet Integral to %0.2f",singletEnd));
    gSingleSet->GetXaxis()->SetTitle("PPM");
    gSingleSet->GetYaxis()->SetTitle("yield");
    gSingleSet->GetYaxis()->SetRangeUser(0, 100);

    gTripleInt->SetName("gTripleInt");
    gTripleInt->SetTitle(Form("Triplet Integral to %0.2f",tripletEnd));
    gTripleInt->SetMarkerColor(kBlue);
    gTripleInt->SetMarkerStyle(23);
    gTripleInt->SetMarkerSize(1.4);
    gTripleInt->GetXaxis()->SetTitle("PPM ");
    gTripleInt->GetYaxis()->SetTitle("yield");
    gTripleInt->GetYaxis()->SetRangeUser(0, 100);
    gSingleSet->Draw("ap");
    gTripleInt->Draw("psame");
    canSet->BuildLegend();
    canSet->Print(".png");


    // norm model graphs 
    normModelGraphs(smean,smeanSetErr);

    //ratios
    double smeanRatio[NSETS];
    double srmsRatio[NSETS];
    double tmeanRatio[NSETS];
    double trmsRatio[NSETS];
    for (unsigned iset = 0; iset < single.size(); ++iset)
    {
      ratioE(smean[iset], srms[iset], smean[0], srms[0], smeanRatio[iset], srmsRatio[iset]);
      ratioE(tmean[iset], trms[iset], tmean[0], trms[0], tmeanRatio[iset], trmsRatio[iset]);
      printf(" set %u singlet %f(%f) ratio  %f(%f)  triplet %f(%f) ratio  %f(%f) \n", 
           iset, smean[iset], srms[iset],smeanRatio[iset],srmsRatio[iset],
           tmean[iset], trms[iset],tmeanRatio[iset],trmsRatio[iset]);
    }

    TGraphErrors *gSinglet = new TGraphErrors(NSETS-1, PPM, smean, PPMError, srms);
    gSinglet->SetTitle("good singlet");

    TMultiGraph* mgSinglet = new TMultiGraph();


    TCanvas *cSinglet = new TCanvas("SingletSet", "SingletSet");
    cSinglet->SetGridx();
    cSinglet->SetGridy();
    gSinglet->SetName("SingletGraph");
    gSinglet->SetTitle("Singlet");
    gSinglet->SetTitle(Form("Singlet Integral %0.2f to %0.2f",singletStart,singletEnd));
    gSinglet->SetMarkerColor(kBlack);
    gSinglet->SetMarkerStyle(21);
    gSinglet->SetMarkerSize(1.4);
    gModelNorm[0]->SetMarkerStyle(23);
    gModelNorm[0]->SetMarkerColor(kRed);


    //gModelNorm[0]->SetMarkerColor(kBlack);
    //gModelNorm[0]->SetMarkerStyle(23);
    //gModelNorm[0]->SetMarkerSize(1.4);

    gSinglet->GetXaxis()->SetTitle("dopant PPM");
    gSinglet->GetYaxis()->SetTitle("yield");
    //gSinglet->GetYaxis()->SetRangeUser(0, 2);
    mgSinglet->Add(gSinglet);
    mgSinglet->Add(gModelNorm[0]);
    gTitle.Form("Singlet Light Yield  ; Xe PPM dopant ; Yield ");
    mgSinglet->SetTitle(gTitle);
    mgSinglet->Draw("ap");
    //cTriplet->BuildLegend();
    gPad->BuildLegend();
    cSinglet->Print(".png");

    TGraphErrors *gTriplet = new TGraphErrors(NSETS-1, PPM, tmean, PPMError, trms);
    gTriplet->SetTitle("good triplet");

    TMultiGraph* mgTriplet = new TMultiGraph();
   
    TCanvas *cTriplet = new TCanvas("TripletSet", "TripletSet");
    cTriplet->SetGridx();
    cTriplet->SetGridy();
    gTriplet->SetName(Form("TripletGraph integral %.2f to %.2f",singletEnd,tripletEnd));
    gTriplet->SetTitle("Triplet total integral");
    gTriplet->SetMarkerColor(kRed);
    gTriplet->SetMarkerStyle(22);
    gTriplet->SetMarkerSize(1.4);

    gTripleInt->SetTitle(Form("Triplet range integral %0.2f to %0.2f",singletEnd,tripletEnd));
    gTripleInt->SetMarkerColor(kBlack);
    gTripleInt->SetMarkerStyle(21);
    gTripleInt->SetMarkerSize(1.4);


    gTriplet->GetXaxis()->SetTitle("dopant PPM");
    gTriplet->GetYaxis()->SetTitle("yield");
    //gTriplet->GetYaxis()->SetRangeUser(0, 2);
    //mgTriplet->Add(gTriplet);
    mgTriplet->Add(gTripleInt);
    gModelNorm[2]->SetMarkerStyle(29);
    gModelNorm[2]->SetMarkerColor(kGreen);

    gModelNorm[1]->SetMarkerStyle(23);
    gModelNorm[1]->SetMarkerColor(kRed);
    mgTriplet->Add(gModelNorm[1]);
    mgTriplet->Add(gModelNorm[2]);
    gTitle.Form("Triplet Light Yield to %.2f ; Xe PPM dopant ; Yield ",tripletEnd);
    mgTriplet->SetTitle(gTitle);
    mgTriplet->Draw("ap");

    //cTriplet->BuildLegend();
    gPad->BuildLegend();
    cTriplet->Print(".png");

    //
    fout->Write();
  }
