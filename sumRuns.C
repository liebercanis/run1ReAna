#include <fstream>
#include <iostream>
#include <numeric>
#include "compiled/BaconAnalysis.hh"

// time is in microseconds
using namespace TMath;
double microToNano = 1.0E3;
enum { NSETS = 6};
enum { NPARS = 7};

int ngood[NSETS];
int nall[NSETS];
double fitpar[NSETS][7];
TH1D* hWave[NSETS];
TH1D* hSingletSet[NSETS];
TH1D* hTripletSet[NSETS];

TH1D* hTripletIntegralSet[NSETS];
TH1D* hModelIntegral[NSETS];
TH1D* hModelWave[NSETS];


TH1D* hSingletSetPass[NSETS];
TH1D* hTripletSetPass[NSETS];
TF1* LifeFit[NSETS];
ofstream textFile;




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





void getHistSum(TH1D* h, int istart, int iend, double& sum, double& err)
{
  sum=0;
  err=0;
  for(int ibin = istart; ibin < iend; ++ibin) {

    sum+= h->GetBinContent(ibin);
    err+= pow(h->GetBinError(ibin),2.);
  }
  err=sqrt(err);
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



void histNorm(int i)
{
  wdir->cd();
  TString hname;
  hname.Form("waveNorm-%i",runList[i]);
  TH1D *hNorm = (TH1D*) waveList[i]->Clone(hname);
  hNorm->SetTitle(Form("waveNorm-%i",runList[i]));
  //hNorm->Sumw2();
  double wevents =  double(waveList[i]->GetEntries())/double(waveList[i]->GetNbinsX());
  double wnorm  = aveEvents/wevents;

  hNorm->Reset();
  for (int ibin = 0; ibin < waveList[i]->GetNbinsX(); ++ibin)
  {
    double val = waveList[i]->GetBinContent(ibin) - allMinimum;
    hNorm->SetBinContent(ibin, val*wnorm);
    hNorm->SetBinContent(ibin, val*wnorm);
    hNorm->SetBinError(ibin, 0); // set later
  }

  waveNormed.push_back(hNorm);
  //fout->Append(hNorm);


  double baseline =  hNorm->Integral(0.,900.,"width")/900.;
  double singlet = hNorm->Integral(sStart,sEnd) - baseline*(sEnd-sStart);
  double triplet = hNorm->Integral(sEnd,tEnd) - baseline*(tEnd-sEnd);
  double tripletRange = hNorm->Integral(tStart,tTripletCut) - baseline*(tTripletCut-tStart);
  int set = getRunSet(fileNumber[i]);
  double cr = singlet/(singlet+triplet);
  printf(" HISTNORM %i set %i file %i %s  events %i singlet %.3E triplet %.3E  cr %.3f  baseline %.3E  \n",
      i,set, runList[i], hNorm->GetName(),  int(wevents), singlet, triplet,cr , baseline );
  setNumber.push_back(   double(set)*25+ runsForSet[set].size() );
  runsForSet[set].push_back(i);
  singletSum.push_back(singlet);
  tripletSum.push_back(triplet);
  hSingletSet[set]->Fill(singlet);
  hTripletSet[set]->Fill(tripletRange);
  bool pass = singlet>singletSetCut[set]&& triplet>tripletSetCut[set];

  if(pass) hSingletSetPass[set]->Fill(singlet);
  if(pass) hTripletSetPass[set]->Fill(tripletRange);


  ntRun->Fill(float(fileNumber[i]),float(set),float(pass),float(wevents),float(singlet),float(triplet),float(tripletRange));
  fout->cd();
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

unsigned long getRunWaves() 
{

  TDirectory* waveDir=NULL;
  wdir= fout->mkdir("waveDir");
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
  }
  return waveList.size();
}

void getModelGraphs(TString modelFile,TString tag)
{
  TFile *fmodel = new TFile(modelFile);
  cout << " \n \t getting model graphs from file " << modelFile << "\n" << endl;
  modelGraphName[0]=TString("RangeSingletModel")+tag;
  modelGraphName[1]=TString("RangeTripletModel")+tag;
  modelGraphName[2]=TString("TotalTripletModel")+tag;
  for(int i=0; i<3; ++i) {
    gModel[i] = NULL;
    fmodel->GetObject(modelGraphName[i], gModel[i]);
    if (!gModel[i]) {
      cout << " did not find " << modelGraphName[i] << endl;
      continue;
    }
    cout << "model graph " << i << " " << gModel[i]->GetName() << endl;
    fout->Append(gModel[i]);
  }

  for(int iset=0; iset<NSETS-1; ++iset) {
    fmodel->GetObject(Form("FitHist-%i",iset),hModelWave[iset]);
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
    if(!gModel[i]) continue;
    // gModelNorm[i] =  (TGraphErrors* gModel[i]->Clone(modelGraphName[i]+TString("-norm"));
    gModelNorm[i] = new TGraphErrors(gModel[i]->GetN());
    gModelNorm[i]->SetName(Form("%s-norm",gModel[i]->GetName()));
    gModelNorm[i]->SetTitle(Form("%s norm",gModel[i]->GetName()));
    cout << "NNNNN norm model graph " << i << " " << gModelNorm[i]->GetName() << " npoints " << gModel[i]->GetN() << " norm " << norm << " points " << gModel[i]->GetN() << " NNNN " << endl;

    for(int j=0; j<= gModel[i]->GetN(); ++j) {
      double xval,yval;
      norm = smean[0]/yval0;
      if(i>0) smean[j]/yval0;
      printf(" model %i point %i norm %E err %E \n",i,j,norm,smeanSetErr[j]);
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

void getStats(std::vector<double> v,double& mean, double& stdev)
{

  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  mean = sum / v.size();

  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(),
      std::bind2nd(std::minus<double>(), mean));
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  stdev = std::sqrt(sq_sum / v.size());
}

void writeHist(TH1D *h)
{
    textFile << "hist " << h->GetName() << "  title " << h->GetTitle() << " y axis title " << h->GetYaxis()->GetTitle() << endl;

    for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
        textFile << ibin << "  " << h->GetBinContent(ibin) << "  " << h->GetBinError(ibin) << endl;
}



void sumRuns()
{


  runsForSet.resize(NSETS);
  for(int iset=0; iset<NSETS; ++iset) runsForSet[iset].clear();

 
  singletSetCut[0]=1.8E4;
  singletSetCut[1]=1.8E4;
  singletSetCut[2]=1.8E4;
  singletSetCut[3]=1.8E4;
  singletSetCut[4]=1.8E4;
  singletSetCut[5]=1.8E4;

  
  tripletSetCut[0]=160.E3;
  tripletSetCut[1]=220.E3;
  tripletSetCut[2]=240.E3;
  tripletSetCut[3]=260.E3;
  tripletSetCut[4]=260.E3;
  tripletSetCut[5]=220.E3;



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



  double PPM[NSETS];
  for(int i=0; i<NSETS; ++i) PPM[i]= fitpar[i][2]; 


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

  
  fout = new TFile("sumRuns.root", "RECREATE");

  ntRun = new TNtuple("ntRun","run","run:set:pass:nev:singlet:triple:tRange");


  int setTotal[NSETS];
  int setGood[NSETS];
  int setRuns[NSETS];

  
  //
  TH2D* qRunTriplet= new TH2D("RunTriplet","run triplet",200,0,200,300,-100,200);
  //
  
  TString hname;
  TString htitle;
  fout->cd();
  for(int iset=0; iset<NSETS; ++iset) {
    hSingletSet[iset] =  NULL;
    hTripletSet[iset] =  NULL;

    hname.Form("Singlet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Singlet yield set %i, %i PPM",iset,int(PPM[iset]));
    hSingletSet[iset] = new TH1D(hname,htitle,100,5.E3,5.E4);
    hSingletSet[iset]->SetLineColor(setColor[iset]);

    hname.Form("Triplet%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Triplet yield set %i, %i PPM",iset,int(PPM[iset]));
    hTripletSet[iset] = new TH1D(hname,htitle,100,1.E5,5.E5);
    hTripletSet[iset]->SetLineColor(setColor[iset]);

    
    hname.Form("SingletPass%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Pass Singlet yield set %i, %i PPM",iset,int(PPM[iset]));
    hSingletSetPass[iset] = new TH1D(hname,htitle,500,5.E3,5.E4);
    hSingletSetPass[iset]->SetLineColor(setColor[iset]);

    hname.Form("TripletPass%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Pass Triplet yield set %i, %i PPM",iset,int(PPM[iset]));
    hTripletSetPass[iset] = new TH1D(hname,htitle,500,1.E4,5.E5);
    hTripletSetPass[iset]->SetLineColor(setColor[iset]);


    hname.Form("TripletIntegral%i-%i-PPM",iset,int(PPM[iset]));
    htitle.Form("Triplet Integral yield set %i, %i PPM",iset,int(PPM[iset]));
    hTripletIntegralSet[iset] = new TH1D(hname,htitle,nDigi,0,nDigi);
    hTripletIntegralSet[iset]->SetLineColor(setColor[iset]);
    hTripletIntegralSet[iset]->GetXaxis()->SetRangeUser(tStart,tTripletCut);



    hWave[iset] = new TH1D(Form("Wave%i", iset), Form("wave %i", iset), nDigi,0,nDigi);
    //hWave[iset]->Sumw2();
    hWave[iset]->SetLineColor(setColor[iset]);
    hWave[iset]->SetMarkerColor(setColor[iset]);
    hWave[iset]->SetMarkerStyle(7);

  }


  if (tSummary)
  {
    for (int iset = 0; iset < NSETS; ++iset)
    {
      setTotal[iset] = 0;
      setGood[iset] = 0;
      setRuns[iset] = 0;
    }
    printf(" tSummary has  %lld  entries \n", tSummary->GetEntries());

    // all entries and fill the histograms
    for (Long64_t jent = 0; jent < tSummary->GetEntries(); jent++)
    {
      tSummary->GetEntry(jent);
      int jset = int(sumVars[ESET]);
      setTotal[jset] += int(sumVars[NTOTAL]);
      setGood[jset] += int(sumVars[NGOOD]);
      ++setRuns[jset];
      fileNumber.push_back(int(sumVars[ERUN]));
      printf(" tSummary entry %lld  run %i set %i total %i good %i \n ",jent, int(sumVars[ERUN]), jset, int(sumVars[NTOTAL]), int(sumVars[NGOOD])  );

      //for(int iv=0; iv<SUMVARS; ++iv ) printf(" \t %i %s %f \n",iv, sumNames[iv].Data(), sumVars[iv]);
    }

    for (int iset = 0; iset < NSETS; ++iset)
      printf(" set %i runs %i total events %i  good events  %i  \n ", iset, setRuns[iset], setTotal[iset],setGood[iset]);
    }

    // pre events
       
    vector<double> prePass;
    prePass.resize(NSETS);
    vector<double> preNum;
    preNum.resize(NSETS);


    double passCount[NBITS][NSETS];
    for(int ib=0; ib<NBITS; ++ib) 
      for(int iset=0; iset<NSETS; ++iset) 
        passCount[ib][iset]=0;


    for (Long64_t jent = 0; jent < tEvPre->GetEntries(); jent++)
    {
      tEvPre->GetEntry(jent);
      int run = preVars[EVRUN];
      int iset = preVars[EVSET];
      int flag = int(preVars[EVFLAG]);
      rejectBits=flag;
      //printf(" jent %llu run %i flag %i ( ",jent,run,flag);
      //for(int ib=0; ib<NBITS; ++ib) printf(" %i ",rejectBits.test(ib));
      //printf(" ) \n");
      for(int ib=0; ib<NBITS; ++ib) if(!rejectBits.test(ib)) passCount[ib][iset] = passCount[ib][iset]+1.;

      preNum[iset] = preNum[iset] + 1.;
      if (flag == 0)
        prePass[iset] = prePass[iset] + 1.;
    }


    // events 
    int nInSet[NSETS]={0,0,0,0,0,0};

    int thisRun = -1;
    for (Long64_t jent = 0; jent < tEvent->GetEntries(); jent++)
    {
      tEvent->GetEntry(jent);
      int iset = evVars[EVSET];
      if(iset<0||iset>5) continue;
      int run = evVars[EVRUN];
      bool pass = evVars[EVSINGLET] > singleCut[iset];
      if(!pass) continue;
      if (evVars[EVFLAG] == 0 && pass)
      {
        ++ngood[iset];
        double rnumber = double(iset)*500. + double(nInSet[iset]);
        ++nInSet[iset];
        runNumber.push_back(rnumber);
        tripletRun.push_back(evVars[EVTRIPLET]);
        singletRun.push_back(evVars[EVSINGLET]);
        //hRunTriplet->Fill(run, evVars[ETRIPLET]);
      } 
    }

    vector<double> preFrac;
    preFrac.resize(NSETS);

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

    TGraph *gPass = new TGraph(NSETS-1, &vset[0], &preFrac[0]);

    TGraph *gPassBit[NBITS];
    for(int ib=0; ib<NBITS; ++ib)  gPassBit[ib] = new TGraph(NSETS-1,&vset[0], &passCount[ib][0]);


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
      gPassBit[ib]->SetMarkerColor(setColor[ib]);
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


    // all entries and fill the histograms
    double PPMError[NSETS];
    for (int iset = 0; iset < NSETS; ++iset)
    {
      ngood[iset] = 0;
      nall[iset] = 0;
      PPMError[iset] = 0;
    }

    
    
    
    printf(" >>>> number of runs in file %s = %lu <<<< \n", inFile->GetName(), getRunWaves());


    for(int iset=0; iset<NSETS; ++iset) setRuns[iset]=0;

    for(int i=0; i<waveList.size(); ++i) getHistStats(i);


    allMinimum = *std::min_element(wMinimum.begin(), wMinimum.end());


    aveEvents = std::accumulate(wEvents.begin(), wEvents.end(), 0)/double(wEvents.size());

    for(int i=0; i<waveList.size(); ++i)  histNorm(i);

    /*
    for(int i=0; i<runsForSet.size(); ++i) {
      printf(" set %i number runs %lu \n", i, runsForSet[i].size());
      for(unsigned srun=0; srun< runsForSet[i].size(); ++srun) printf("\t run %i \n",runsForSet[i][srun]);
    }
    */

    printf(" \t\t >>>>> waveNormed %lu ave events %f all min %f \n",waveNormed.size(), aveEvents,allMinimum);


    int theSet = -1;
    can=NULL;
    bool first=true;
    int theColor=kBlack;
    int setCount=0;
    for(unsigned i=0; i< waveList.size(); ++i) {
      int iset = getRunSet(runList[i]);
      if(iset!=theSet) {
        if(can) {
          printf(" finished set %i count %i name %s \n", theSet, setCount, cname.Data());
          can->BuildLegend();
          can->Print(".pdf");
        }
        theSet=iset;
        cname.Form("AllWavesForSet-%i",iset);
        can = new TCanvas(cname,cname);
        can->SetLogy();
        first=true;
        theColor=kBlack;
        setCount=0;
      }

      bool pass = singletSum[i]>singletSetCut[iset]&& tripletSum[i]>tripletSetCut[iset];
      printf(" %i file %i set %i count %i %s entries %.4E singlet %E triplet %E pass %i \n",
          i,fileNumber[i],iset, setCount, waveNormed[i]->GetName() , waveNormed[i]->GetEntries(),singletSum[i],tripletSum[i],int(pass) );

      if(!pass) continue;

      waveNormed[i]->SetLineColor(theColor+setCount++);
      if(first) {
        waveNormed[i]->Draw("");
        first=false;
      }
      else 
        waveNormed[i]->Draw("same");
    }

    if(can) {
      printf(" finished set %i count %i name %s \n", theSet, setCount, cname.Data());
      can->BuildLegend();
      can->Print(".pdf");
    }


    
    // calc and set bin errors
    double mean,sdev;
    int iwave=0;
    for(unsigned iwave=0 ; iwave < waveNormed.size(); ++iwave) {
      std::vector<double> vin;
      for(int ip=waveNormed[iwave]->GetNbinsX()-5000 ; ip<waveNormed[iwave]->GetNbinsX(); ++ip) vin.push_back(waveNormed[iwave]->GetBinContent(ip));
      getStats(vin,mean,sdev);     
      //printf(" i1=%i mean %f sdev %f \n",iwave,mean,sdev);
      for(unsigned ip=0; ip< waveNormed[iwave]->GetNbinsX();++ip)  waveNormed[iwave]->SetBinError(ip+1,sdev);
      vin.clear();
    }
    

    // sum waveforms over sets as weighted mean,error
    // loop over sets
    for (int iset = 0; iset<NSETS ; ++iset) {
      // loop over bins
      for (int jbin = 1; jbin < hWave[iset]->GetNbinsX()+1 ; ++jbin) {
        double sum=0;
        double wsum=0;
        // ave over runs in this set for this bin 
        for (int jrun = 0; jrun < runsForSet[iset].size() ; ++jrun) {
          int krun = runsForSet[iset][jrun];
          //if(jbin>1000&&jbin<1100) printf(" jbin %i run %i set %i \n",jbin,krun,iset);
          double wi=pow(waveNormed[krun]->GetBinError(jbin),-2.);
          sum+=  waveNormed[krun]->GetBinContent(jbin)*wi;
          wsum+= wi;
          //if(jbin>1000&&jbin<1100) printf("in set  %i bin %i  %f %f sum %f %f \n",iset,jbin,
            //waveNormed[jrun]->GetBinContent(jbin),waveNormed[jrun]->GetBinError(jbin),sum,wsum);
        }
        // 
        double mu = sum/wsum;
        double muError = 1.0/sqrt(wsum);
        hWave[iset]->SetBinContent(jbin,mu);
        hWave[iset]->SetBinError(jbin,muError);
        //if(jbin>1000&&jbin<1100) printf("set %i bin %i wsum %f %f mu %f %f \n",iset,jbin,sum,wsum, mu,muError);
      }
    }

    // write out waveforms
    textFile.open("run1SummedWaves.txt");
    for(int iset=0; iset<NSETS; ++iset)  writeHist(hWave[iset]);
    textFile.close();


    // plot summed by set
    for(int iset=0; iset< NSETS-1; ++iset) {
      cname.Form("SummedForSet-%i",iset);
      can = new TCanvas(cname,cname);
      can->SetLogy();
      hWave[iset]->Draw();
      can->Print(".pdf");
    }

    
    // plot summed all sets
    cname.Form("AllSetsSummed");
    can = new TCanvas(cname,cname);
    can->SetLogy();
    can->SetGridx(); can->SetGridy();

    gStyle->SetOptStat(0);
    for(int iset=0; iset< NSETS-1; ++iset) {
      hWave[iset]->SetTitle("");
      hWave[iset]->SetStats(0);
      hWave[iset]->GetYaxis()->SetRangeUser(4.0,1200.);
      if(iset==0) hWave[iset]->Draw();
      else hWave[iset]->Draw("same");
    }
    can->BuildLegend();
    can->Print(".pdf");


    
    TMultiGraph *mgBeforeCuts = new TMultiGraph();

    TGraph *gSinglet = new TGraphErrors(setNumber.size()-3, &setNumber[0], &singletSum[0]);
    TGraph *gTriplet = new TGraphErrors(setNumber.size()-3, &setNumber[0], &tripletSum[0]);
    TCanvas *cTripletBefore = new TCanvas("set-sums-base","set-sums-base");
    cTripletBefore->SetGridx(); cTripletBefore->SetGridy();
    gSinglet->SetName("gSinglet");
    gSinglet->SetTitle("singlet");
    gSinglet->SetMarkerColor(kRed);
    gSinglet->SetMarkerStyle(22);
    gSinglet->SetMarkerSize(.4);
    gSinglet->GetHistogram()->GetXaxis()->SetTitle("  set*20+sequential run");
    gSinglet->GetHistogram()->GetYaxis()->SetTitle(" singlet yield");

    gTriplet->SetName("gTriplet");
    gTriplet->SetTitle("triplet");
    gTriplet->SetMarkerColor(kBlack);
    gTriplet->SetMarkerStyle(21);
    gTriplet->SetMarkerSize(.4);
    gTriplet->GetHistogram()->GetXaxis()->SetTitle("  set*20+sequential run");
    gTriplet->GetHistogram()->GetYaxis()->SetTitle(" triplet yield");

    mgBeforeCuts->Add(gTriplet);
    mgBeforeCuts->Add(gSinglet);
    mgBeforeCuts->SetTitle("Light Yield; set*20+sequential run ; yield in range");
    mgBeforeCuts->Draw("ap");
    cTripletBefore->BuildLegend();
    cTripletBefore->Print(".pdf");


    fout->Write();
}
