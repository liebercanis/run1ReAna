#include <iostream>
#include <fstream>
#include "compiled/BaconAnalysis.hh"
// time is in microseconds
using namespace TMath;
static double singletFraction = 0.14;
static double tTriplet=2.100; 
static double tSinglet=5.E-3;
static double tXe = 20.0E-3;
static double kxe = 8.8E-2;
static double kplusZero  = 1.3E-1;
static double xstart=1.;
static double xstop=3.;
static double tailStart=7.0;
static double tailStop=10.0;
static double pmodtG = 3.48;
double k1Zero = kxe*131./40.;


static void lineExtrap(double a, double ae, double b, double be, double &z, double &ze) {
  z = -b/a;
  ze = pow(z*ae,2.)-z*ae*be + be*be;
  ze = sqrt(ze);
}

static void ratioE(double x, double xe, double y, double ye, double& r,double& re){
  r = x/y;
  re = r * sqrt( pow(xe/x,2.)+pow(ye/y,2.) );
}



static void weightAve(int n, double *x, double *xe, double &ave, double &err)
{
  ave = 0;
  err = 0;
  if (n < 1)
        return;
    for (int i = 0; i < n; ++i)
        printf(" weight Ave %i x %f xe %f \n", i, x[i], xe[i]);
    for (int i = 0; i < n; ++i)
    {
        double wi = 1. / pow(xe[i], 2.);
        err += wi;
        ave += wi * x[i];
    }
    ave = ave / err;
    err = 1. / sqrt(err);
    return;
}

static void aveVec(vector<double> x, vector<double> xe, double &a, double &ae)
{
    ae = 0;
    a = 0;
    if (x.size() < 1)
        return;
    //for(unsigned  i=0; i<x.size(); ++i) printf(" vec Ave %i x %f xe %f \n",i,x[i],xe[i]);
    for (unsigned i = 0; i < x.size(); ++i)
    {
        double xei = xe[i];
        if (xei == 0)
            xei = 1E-9;
        double wi = 1. / pow(xei, 2.);
        ae += wi;
        a += wi * x[i];
    }
    a = a / ae;
    ae = 1. / sqrt(ae);
}

// add quenching
// par binwidth norm PPM kplus type 
//     0        1     2   3     4 

static Double_t lightModel(Double_t *xx, Double_t *par)
{ 
  double x=xx[0]-xstart;
  double bw = par[0];
  double norm = par[1];
  double kx = k1Zero*par[2];
  double tau3 = par[3];
  double kp = par[4];
  int type = int(par[5]);
  double td = 1/kx;
  double tq = 1./(1./tau3+kp);
  double tr = 1./(kx + 1./tq);
  double alpha1 = bw*singletFraction*norm;
  double alpha3 = bw*(1.-singletFraction)*norm;
  double fs =  alpha1/tSinglet*Exp( -x/tSinglet);
  double f3 =  alpha3/tau3*Exp(-x/tr);
  double fx =  alpha3*pow(kx,2)*tq*( Exp( -x/td) - Exp( -x/tr) );
  double f=0;
  if (type == 0)
    f = fs;
  else if (type == 1)
    f = f3;
  else if (type == 2)
    f = fx;
  else if (type == 3)
    f = fs+f3+fx;
  return f;
}



class tbSFitAll
{
public:
    //double taut = 2100.;
    tbSFitAll(bool fix = true, double lifetime3 = tTriplet); //2.93); 3.48
    virtual ~tbSFitAll() { ; }
    TFile *inFile;
    enum
    {
        MAXSETS = 6
    };
    enum
    {
        MAXFITS = 6
    };
    enum
    {
        npars = 6 
    };

    Float_t sumVars[SUMVARS];
    int setTotal[MAXSETS];
    int setRuns[MAXSETS];


    TF1 *setFit(int iset, TH1D *hLife, double ppm);
    TF1 *tailFit(int iset, TH1D *hTail);
    void tbFit1(int iset, TH1D *hLife);
    TCanvas *fhistCan(TString title, TH1D *h, TF1 *fp);
    double rangeIntegral(TH1D *h, double xlow, double high);

    void writeHist(TH1D *h);
    ofstream textFile;

    TString fitTag;
    TF1 *fp[MAXSETS];
    TF1 *pmodelFit[MAXSETS][MAXSETS];
    TF1 *fmodelFit[MAXSETS][MAXFITS];
    TCanvas *canFmodel[MAXSETS];
    TCanvas *canFit[MAXSETS];

    TString runTag[MAXSETS];
    TString runRange[MAXSETS];
    TH1D *hMpvSet[MAXSETS];
    TH1D *hLifeSum[MAXSETS];
    TH1D *hLifeInt[MAXSETS];

    TH1D *hLifeBlah[MAXSETS];
    TH1D *hMpvSum[MAXSETS];
    int pcolor[MAXSETS];
    int pstyle[MAXSETS];
    int setColor[MAXSETS];
    int setStyle[MAXSETS];
    double runPPM[MAXSETS];
    double tauGSet[MAXSETS];
    double setChi[MAXSETS];
    double pmodtT;
    bool fixTriplet;

    std::vector<double> vset;
    std::vector<double> vsetErr;
    std::vector<double> va0;
    std::vector<double> va0Err;
    std::vector<double> vk1;
    std::vector<double> vk1Err;

    std::vector<double> vsetMu;
    std::vector<double> vsetErrMu;
    std::vector<double> va0Mu;
    std::vector<double> va0ErrMu;
    std::vector<double> vk1Mu;
    std::vector<double> vk1ErrMu;

    std::vector<double> vspeMu;
    std::vector<double> vspeErrMu;
    std::vector<double> vspe;
    std::vector<double> vspeErr;

    std::vector<double> vsetAll;
    std::vector<double> vsetAllErr;
    std::vector<double> vk1All;
    std::vector<double> vk1AllErr;

    std::vector<double> va0All;
    std::vector<double> va0AllErr;
    std::vector<double> va0AllMu;
    std::vector<double> va0AllErrMu;

    std::vector<double> vppmAll;
    std::vector<double> vppmAllErr;
    std::vector<double> vppmAllMu;
    std::vector<double> vppmAllErrMu;

    std::vector<double> vspeAll;
    std::vector<double> vspeAllErr;
    std::vector<double> vtautAll;
    std::vector<double> vtautAllErr;
    std::vector<double> vtaugAll;
    std::vector<double> vtaugAllErr;

    std::vector<double> vg0All;
    std::vector<double> vg0AllErr;

    std::vector<double> vgfracAll;
    std::vector<double> vgfracAllErr;

    std::vector<double> vXeFrac;
    std::vector<double> vXeFracErr;
    std::vector<double> vppm0;
    std::vector<double> vppm0Err;
  };

void tbSFitAll::writeHist(TH1D *h)
{

    textFile << " hist " << h->GetName() << "  title " << h->GetTitle() << " y axis title " << h->GetYaxis()->GetTitle() << endl;

    for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
        textFile << ibin << "  " << h->GetBinContent(ibin) << "  " << h->GetBinError(ibin) << endl;
}

double tbSFitAll::rangeIntegral(TH1D *h, double xlow, double xhigh)
{
    double sum = 0;
    int i1 = h->FindBin(xlow);
    int i2 = h->FindBin(xhigh);
    for (int ibin = i1; ibin < i2; ++ibin)
        sum += h->GetBinContent(ibin);
    return sum;
}

TCanvas *tbSFitAll::fhistCan(TString title, TH1D *h, TF1 *fp)
{
    TCanvas *can = new TCanvas(title, title);
    can->SetLogy();
    gStyle->SetOptFit();
    //if(fp) fp->SetLineColor(kTeal-6);
    //h->GetYaxis()->SetRangeUser(1E-1,1E3);
    h->Draw();
    if (fp)
        fp->Draw("sames");
    can->Print(".png");
    return can;
}

TF1 *tbSFitAll::tailFit(int iset, TH1D *hTail)
{
  TF1 *ftail=NULL;
  return ftail;
}

/*

TF1 *tbSFitAll::tailFit(int iset, TH1D *hTail)
{
    // fit exponential tail
    double tauG = 3.140;
    double binwidth = hTail->GetBinWidth(1);
    double tailIntegral = rangeIntegral(hTail, tailStart, tailStop);
    TF1 *ftail = new TF1(Form("fexpSet-%i", iset), fexp, singletEnd, xstop, 3);
    ftail->SetNpx(1000); // numb points for function
    ftail->SetParNames("lifetime", "integral", "binwidth");
    ftail->SetParameters(tauG, tailIntegral, binwidth);
    //ftail->FixParameter(0,3.2);
    ftail->FixParameter(2, binwidth);

    ftail->Print();
    for (int ii = 0; ii < 3; ++ii)
    {
        printf("\t  param %i %s %.3f +/- %.3f \n", ii, ftail->GetParName(ii), ftail->GetParameter(ii), ftail->GetParError(ii));
    }

    hTail->Fit(ftail, "RLE0+", "", tailStart, tailStop);
    return ftail;
}
*/
TF1 *tbSFitAll::setFit(int iset, TH1D *hLife, double ppm)
{
    Double_t binwidth = hLife->GetBinWidth(1);
    // ppm is by number 
    Double_t integral = rangeIntegral(hLife, singletEnd, xstop);
    double norm = integral; //-G0Par;
    int ifit = 3;
    double kplus = kplusZero;

    TF1* fp  = new TF1(Form("LfFit-%i-%0.f-PPM", iset, ppm),lightModel,singletEnd,tailStop,npars);
    fp->SetParName(0,"binw");
    fp->SetParName(1,"norm");
    fp->SetParName(2,"PPM");
    fp->SetParName(3,"tau3");
    fp->SetParName(4,"k+");
    fp->FixParameter(0, binwidth);
    fp->SetParameter(1, integral);
    fp->SetParameter(2, ppm);
    fp->SetParameter(3, tTriplet);
    fp->FixParameter(4, kplus);
    fp->FixParameter(5,ifit);
    fp->SetTitle(Form("LfFit-%i-%0.f-PPM", iset, ppm));
    fp->SetNpx(1000); // numb points for function
    /* do the fit here */
    TFitResultPtr fptr = hLife->Fit(fp, "LE0S+","",xstart,xstop);
    TMatrixDSym cov = fptr->GetCorrelationMatrix();
    printf(" correlation set %i cov(2,3) %f \n", iset, cov(2, 3));
    fptr->Print("V");
    fp->Print();
    fp->SetLineColor(pcolor[iset]);
    fp->SetMarkerStyle(pstyle[iset]);

    printf(" \n\n >>> setFit starting parameters \n");
    for (int ii = 0; ii < npars; ++ii)
    {
        printf("\t  param %i %s %.3f +/- %.3f \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
    }
    // fill arrays
    setChi[iset] = fp->GetChisquare();
    return fp;
}

// MAIN 
tbSFitAll::tbSFitAll(bool fix, double lifetime3)
{
    double markerSize = 0.5;
    
    for(int iset=0; iset<MAXSETS; ++iset) {
      hMpvSet[iset]=NULL;
      hLifeSum[iset]=NULL;
      hLifeInt[iset]=NULL;
      hLifeBlah[iset]=NULL;
      hMpvSum[iset]=NULL;
    }

    pmodtT = lifetime3;
    textFile.open("run1Data.txt");
    if (fix)
        fitTag.Form("fixedTaut=%.3f", lifetime3);
    else
        fitTag = TString("floatTaut");
    cout << " fit Tag = " << fitTag << endl;
    fixTriplet = fix;
    inFile = new TFile("BaconAnalysisHighCut-1000.root", "READONLY");
    if (!inFile)
        return;
    TFile *fout = new TFile("tbSFitOut.root", "RECREATE");

    runTag[0] = TString("00PPM-Ran");
    runTag[1] = TString("01PPM-Ran");
    runTag[2] = TString("02PPM-Ran");
    runTag[3] = TString("05PPM-Ran");
    runTag[4] = TString("10PPM-Ran");
    runTag[5] = TString("10PPM-Mu");

    runRange[0] = TString("03000-03020");
    runRange[1] = TString("20005-20020");
    runRange[2] = TString("20025-20040");
    runRange[3] = TString("20045-20060");
    runRange[4] = TString("20065-20080");
    //runRange[5]=TString("20200-20202");
    runRange[5] = TString("20220-20222");
    //runRange[5]=TString("20205-20215");

    runPPM[0] = 0;
    runPPM[1] = 1;
    runPPM[2] = 2;
    runPPM[3] = 5;
    runPPM[4] = 10;
    runPPM[5] = 10;

    tauGSet[0] = 3.140;
    tauGSet[1] = 3.004;
    tauGSet[2] = 3.004;
    tauGSet[3] = 3.004;
    tauGSet[4] = 3.004;
    tauGSet[5] = 3.140;
    for (int iset = 0; iset < MAXSETS; ++iset)
        tauGSet[iset] = pmodtG;

    for (int i = 0; i < MAXSETS; ++i)
    {
        pcolor[i] = i;
        pstyle[i] = 20 + i;
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

    pcolor[5] = kBlack;
    pcolor[4] = kBlue;
    pcolor[3] = kYellow - 2;
    pcolor[2] = kGreen;
    pcolor[1] = kMagenta + 2;
    pcolor[0] = kRed;

    TIter next(inFile->GetListOfKeys());
    TKey *key;
    TH1D *hlife = NULL;
    int iset3 = 0;
    while ((key = (TKey *)next()))
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1"))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();
        TString hname(h->GetName());
        cout << hname << endl;
        if (hname.Contains("Wave-"))
          hLifeBlah[iset3++] = h;
    }

    TTree* tSummary=NULL;

    inFile->GetObject("Summary",tSummary);
    //TNtuple *ntSummary = new TNtuple("Summary", " Summary ", "run:set:base:baseend:accept:total:singlet:dublet:triplet:ngood");
    TString sumNames[SUMVARS];
    sumNames[0]=  TString("run");
    sumNames[1]=  TString("set");
    sumNames[2]=  TString("base");
    sumNames[3]=  TString("baseend");
    sumNames[4]=  TString("accept");
    sumNames[5]=  TString("total");
    sumNames[6]=  TString("singlet");
    sumNames[7]=  TString("dublet");
    sumNames[8]=  TString("triplet");
    sumNames[9]=  TString("ngood");


    for(int iv=0; iv< SUMVARS; ++iv ) tSummary->SetBranchAddress(sumNames[iv],&sumVars[iv]);

    for(int iset=0; iset<MAXSETS; ++iset ) { 
      setTotal[iset]=0;
      setRuns[iset]=0;
    }


    printf(" got %i \n", iset3);
    if(tSummary) printf(" tSummary has  %lld  entries \n",tSummary->GetEntries());
    else printf(" no tSummary \n ");

    // all entries and fill the histograms
    for (Long64_t jent=0;jent<tSummary->GetEntries() ;jent ++) {
      tSummary->GetEntry(jent);
      //printf(" tSummary entry %lld \n ",jent);
      for(int iv=0; iv< SUMVARS; ++iv ) {
        int jset = int(sumVars[ESET]);
        int nset = int(sumVars[NGOOD]);
        setTotal[jset]+=  int(sumVars[NGOOD]);
        ++setRuns[jset];
        //printf(" \t %i %s %f \n",iv, sumNames[iv].Data(), sumVars[iv]); 
      }
    }

    for(int iset=0; iset<MAXSETS; ++iset ) printf(" set %i runs %i total events %i \n ", iset, setRuns[iset], setTotal[iset]);




    if(iset3!= MAXSETS) return;
    printf(" nbins is  %i \n",hLifeBlah[0]->GetNbinsX());


    for (int i = 0; i < MAXSETS; ++i)
    {
        hLifeSum[i] = (TH1D *)hLifeBlah[i]->Clone(Form("LifeToFitSet%i", i));
        hLifeInt[i] = (TH1D *)hLifeBlah[i]->Clone(Form("IntegralToFitSet%i", i));
        hLifeSum[i]->Reset();
        hLifeInt[i]->Reset();
        double yieldSum = 0;
        double yieldSumErr2 = 0;

        // shift trigger times
        int maxBin = hLifeBlah[i]->GetMaximumBin();
        int startBin = hLifeBlah[i]->FindBin(1.0);

        // normalize to total events
        double normFact = double(setTotal[0])/double(setTotal[i]);
        printf(" \t\t norm factor set %i = %f \n",i,normFact);
        for (int ibin = 0; ibin < hLifeBlah[i]->GetNbinsX(); ++ibin)
        {
          int jbin = ibin - maxBin+startBin;
          double val = hLifeBlah[i]->GetBinContent(ibin);
          hLifeSum[i]->SetBinContent(jbin,val*normFact);
          hLifeSum[i]->SetBinError(jbin, normFact*sqrt(val));
          if(val>1) {
            yieldSum += val;
            yieldSumErr2 += 1/val;
          }
          hLifeInt[i]->SetBinContent(jbin, yieldSum*normFact);
          hLifeInt[i]->SetBinError(jbin, normFact*sqrt(yieldSumErr2));
          //if(i==0&& val>1 ) printf(" bin %i val %f  sum %f err2 %f err %f \n",jbin,val,yieldSum,  yieldSumErr2 , sqrt(yieldSumErr2));
        }

        double spe = rangeIntegral(hLifeSum[i], singletEnd, tailStart);
        if (i == 5)
        {
            vspeMu.push_back(spe);
            vspeErrMu.push_back(sqrt(spe));
        }
        else
        {
            vspe.push_back(spe);
            vspeErr.push_back(sqrt(spe));
        }
        vspeAll.push_back(spe);
        vspeAllErr.push_back(sqrt(spe));

        printf(" set %i starting %.0f  normalized  %.0f \n", i, hLifeBlah[i]->Integral(), hLifeSum[i]->Integral());
    }


    // rbin
    printf(" rbinning 5 \n");
    for (int iset = 0; iset < MAXSETS; ++iset)
    {
        printf(" %i  %s ", iset,  hLifeSum[iset]->GetName());
        cout << " nbins " << hLifeSum[iset]->GetNbinsX() << " wid  " << hLifeSum[iset]->GetBinWidth(1) << " inte " << hLifeSum[iset]->Integral();
        hLifeSum[iset]->Rebin(5);
        cout << " nbins " << hLifeSum[iset]->GetNbinsX() << " wid  " << hLifeSum[iset]->GetBinWidth(1) << " inte " << hLifeSum[iset]->Integral() << endl;
    }

    TCanvas *canLifeAll = new TCanvas("LifeALL", "LifeALL");
    canLifeAll->SetLogy();
    gStyle->SetOptStat(0);
    hLifeSum[0]->Draw();
    double yup = 1E6;
    hLifeSum[0]->GetYaxis()->SetRangeUser(1E-1, yup);
    hLifeSum[0]->SetTitle(Form("sets %i ",0));
    for (int iset = 0; iset < MAXSETS; ++iset)
    {
        hLifeSum[iset]->SetTitle(Form("sets %i ",iset));
        hLifeSum[iset]->GetYaxis()->SetTitle("yield SPE/40 ns");
        hLifeSum[iset]->GetYaxis()->SetRangeUser(1E-1, yup);
        hLifeSum[iset]->Draw("sames");
        hLifeSum[iset]->SetMarkerStyle(setStyle[iset]);
        hLifeSum[iset]->SetMarkerColor(setColor[iset]);
        hLifeSum[iset]->SetLineColor(setColor[iset]);
        hLifeSum[iset]->SetMarkerSize(markerSize);
        writeHist(hLifeSum[iset]);
    }
    canLifeAll->BuildLegend();
    canLifeAll->Print(".png");

    textFile.close();

    TCanvas *canIntAll = new TCanvas("IntALL", "IntALL");
    gStyle->SetOptStat(0);
    yup = 5E6;
    for (int iset = 0; iset < MAXSETS; ++iset)
    {
      hLifeInt[iset]->SetTitle("SPE versus Time all sets ");
      hLifeInt[iset]->GetYaxis()->SetTitle("integral yield SPE/40 ns");
      hLifeInt[iset]->GetXaxis()->SetRangeUser(1, 4);
      hLifeInt[iset]->GetYaxis()->SetRangeUser(0,yup);
      hLifeInt[iset]->SetMarkerStyle(setStyle[iset]);
      hLifeInt[iset]->SetMarkerColor(setColor[iset]);
      hLifeInt[iset]->SetLineColor(setColor[iset]);
      hLifeInt[iset]->SetMarkerSize(markerSize);
      if(iset==0) hLifeInt[iset]->Draw("");
      else  hLifeInt[iset]->Draw("sames");
    }
    canIntAll->Print(".png");

    for (int iset = 0; iset < MAXSETS; ++iset)
      tbFit1(iset, hLifeSum[iset]);


    
    TGraphErrors *gTriplet0 = new TGraphErrors(vset.size()-1, &vset[0], &va0All[0], &vsetErr[0], &va0AllErr[0]);
    TGraphErrors *gTripletMu = new TGraphErrors(1,&vset[MAXSETS-1], &va0All[MAXSETS-1], &vsetErr[MAXSETS-1], &va0AllErr[MAXSETS-1]);
    TMultiGraph *mgTriplet = new TMultiGraph();

    TCanvas *cTriplet0 = new TCanvas("Triplet0","Triplet0");
    cTriplet0->SetGridx(); cTriplet0->SetGridy();
    gTriplet0->SetName("Triplet0-By-PPM");
    gTriplet0->SetTitle("Triplet0-By-PPM");
    gTriplet0->SetMarkerColor(kRed);
    gTriplet0->SetMarkerStyle(22);
    gTriplet0->SetMarkerSize(1.4);
    gTripletMu->SetMarkerColor(kBlue);
    gTripletMu->SetMarkerStyle(21);
    gTripletMu->SetMarkerSize(1.4);
    gTriplet0->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    gTriplet0->GetHistogram()->GetYaxis()->SetTitle(" fitted Triplet yield");
    mgTriplet->Add(gTriplet0,"p");
    mgTriplet->Add(gTripletMu,"p");
    mgTriplet->SetTitle("Light Yield; Xe PPM dopant ; fitted Triplet yield");
    mgTriplet->Draw("a");
    cTriplet0->Print(".png");

    
    TGraphErrors *gPPM = new TGraphErrors(vset.size()-1, &vset[0], &vppmAll[0], &vsetErr[0], &vppmAllErr[0]);
    TGraphErrors *gPPMMu = new TGraphErrors(1, &vset[MAXSETS-1], &vppmAll[MAXSETS-1], &vsetErr[MAXSETS-1], &vppmAllErr[MAXSETS-1]);
    TMultiGraph *mgppm = new TMultiGraph();



    TCanvas *cFitPPM = new TCanvas("FitPPM","FitPPM");
    cFitPPM->SetGridx(); cFitPPM->SetGridy();
    gPPM->SetName("FitPPM-By-PPM");
    gPPM->SetTitle("FitPPM-By-PPM");
    gPPM->SetMarkerColor(kRed);
    gPPM->SetMarkerStyle(22);
    gPPM->SetMarkerSize(1.4);
    gPPMMu->SetMarkerColor(kBlue);
    gPPMMu->SetMarkerStyle(21);
    gPPMMu->SetMarkerSize(1.4);
    gPPM->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    gPPM->GetHistogram()->GetYaxis()->SetTitle(" fitted PPM");
    mgppm->Add(gPPM,"p");
    mgppm->Add(gPPMMu,"p");
    mgppm->SetTitle("Fitted PPM; Xe PPM dopant ; PPM fitted");
    mgppm->Draw("a");
    cFitPPM->Print(".png");


     
    TGraphErrors *gtau3 = new TGraphErrors(vset.size()-1, &vset[0], &vtautAll[0], &vsetErr[0], &vtautAllErr[0]);
    TGraphErrors *gtau3Mu = new TGraphErrors(1, &vset[MAXSETS-1], &vtautAll[MAXSETS-1], &vsetErr[MAXSETS-1], &vtautAllErr[MAXSETS-1]);
    TMultiGraph *mgtau3 = new TMultiGraph();



    TCanvas *cFittau3 = new TCanvas("Fittau3","Fittau3");
    cFittau3->SetGridx(); cFittau3->SetGridy();
    gtau3->SetName("Fittau3-By-tau3");
    gtau3->SetTitle("Fittau3-By-tau3");
    gtau3->SetMarkerColor(kRed);
    gtau3->SetMarkerStyle(22);
    gtau3->SetMarkerSize(1.4);
    gtau3Mu->SetMarkerColor(kBlue);
    gtau3Mu->SetMarkerStyle(21);
    gtau3Mu->SetMarkerSize(1.4);
    gtau3->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    gtau3->GetHistogram()->GetYaxis()->SetTitle(" fitted triplet lifetime");
    mgtau3->Add(gtau3,"p");
    mgtau3->Add(gtau3Mu,"p");
    mgtau3->SetTitle("Fitted triplet lifetime ; Xe dopant PPM ; fitted triplet lifetime");
    mgtau3->Draw("a");
    cFittau3->Print(".png");


    for(int iset=0; iset<MAXSETS; ++iset) {
      printf("........ fitted parameters set %i\n", iset);
      for(int j=0; j<npars; ++j) printf(" par[%i][%i] = %f ;\n",iset,j,fp[iset]->GetParameter(j));
    }

    fout->Write();
}

void tbSFitAll::tbFit1(int iset, TH1D *hLife)
{

    double markerSize = 0.5;
    vset.push_back(runPPM[iset]);
    vsetErr.push_back(0);

    /***** fitting *****/
    Double_t binwidth = hLife->GetBinWidth(1);

    // pmodel
    double Single0 = 0;

    // fmodel
    double tauX = tXe;
    double tauA = pmodtT;
    double tauM = pmodtG;

    double integral = rangeIntegral(hLife, singletEnd, tailStart);

    // fit exponential tail
    double gnorm = 0;
    TH1D *hTail = (TH1D *)hLife->Clone(Form("TailSet%i", iset));
    double gfrac = 0;
    double G0Par = 0;
    double tauG = 3.140;
    double tauGErr = 0;
    TF1 *ftail = tailFit(iset, hTail);
    if (ftail)
    {
        //for(int ii=0; ii<3; ++ii) {
        //  printf("\t  param %i %s %.3f +/- %.3f \n",ii,ftail->GetParName(ii),ftail->GetParameter(ii),ftail->GetParError(ii));
        //}
        TString tailtitle;
        tailtitle.Form("TailFitToSet%i", iset);
        ftail->SetLineColor(kRed);
        ftail->SetLineWidth(2);
        fhistCan(tailtitle, hTail, ftail);
        //
        double sum = ftail->Integral(tailStart, tailStop) / hTail->GetBinWidth(1);
        tauG = ftail->GetParameter(0);
        tauGErr = ftail->GetParError(0);
        G0Par = ftail->Integral(singletEnd, xstop) / hTail->GetBinWidth(1);
        gfrac = G0Par / integral;
        printf(" >> ftail set %i tauG %.3f lifesum %.2E sum %.2E G0Par %.2E gfrac %.3f\n",
               iset, tauG, integral, sum, G0Par, gfrac);
    }
    vtaugAll.push_back(tauG);
    vtaugAllErr.push_back(tauGErr);
    vg0All.push_back(G0Par);
    vg0AllErr.push_back(sqrt(G0Par));

    fp[iset] = setFit(iset, hLife, runPPM[iset]);
   
    printf(" xxxx  %.f PPM \n", runPPM[iset]);

    va0All.push_back(fp[iset]->GetParameter(1));
    va0AllErr.push_back(fp[iset]->GetParError(1));
    vppmAll.push_back(fp[iset]->GetParameter(2));
    vppmAllErr.push_back(fp[iset]->GetParError(2));

    
    //fitted lifetime 
    vtautAll.push_back(fp[iset]->GetParameter(3));
    vtautAllErr.push_back(fp[iset]->GetParError(3));
    
    double g0frac;
    double g0fracErr;
    ratioE(fp[iset]->GetParameter(5), fp[iset]->GetParError(5), fp[iset]->GetParameter(1), fp[iset]->GetParError(1), g0frac, g0fracErr);
    vgfracAll.push_back(g0frac);
    vgfracAllErr.push_back(g0fracErr);

    double lightAll = fp[iset]->Integral(singletEnd, xstop) / binwidth;
    printf(" set %i lightAll %.3E \n", iset, lightAll);

    gStyle->SetOptFit(1111111);
    gStyle->cd();

    canFit[iset] = new TCanvas(Form("LifeFit-%.f-PPM-%s", runPPM[iset], runTag[iset].Data()), Form("LifeFit-%.f-PPM-%s", runPPM[iset], runTag[iset].Data()));
    canFit[iset]->SetLogy();
    hLife->SetTitle(Form("LifeFit set %i Xe %s ", iset + 1, runTag[iset].Data()));
    //hLife->GetYaxis()->SetRangeUser(1E-1,1E3);
    hLife->SetMarkerSize(0.5);
    hLife->SetLineColor(setColor[iset]);
    hLife->SetMarkerColor(setColor[iset]);
    hLife->Draw("p");
    fp[iset]->SetLineColor(kTeal - 6);
    fp[iset]->SetLineStyle(5);
    fp[iset]->SetLineWidth(4);
    fp[iset]->Draw("sames");
    if(ftail) ftail->Draw("sames");
    TPaveStats *st = (TPaveStats *)hLife->GetListOfFunctions()->FindObject("stats");
    gStyle->SetOptFit(); //for example
    canFit[iset]->Modified();
    canFit[iset]->Update();
    canFit[iset]->Print(".png");

    gStyle->SetOptFit(1111111);
    gStyle->cd();

    double ArSum = fp[iset]->Integral(singletEnd, xstop);
    double XeSum = fp[iset]->Integral(singletEnd, xstop);
    double totSum = ArSum + XeSum;
    double XeFrac = XeSum / totSum;
    double XeFracErr = sqrt(XeFrac * (1. - XeFrac) / totSum);
    //ratioE(XeSum,sqrt(XeSum), totSum,sqrt(totSum), XeFrac, XeFracErr);

    vXeFrac.push_back(XeFrac);
    vXeFracErr.push_back(XeFracErr);

 

}
