#include <iostream>
#include <fstream>
#include "compiled/BaconAnalysis.hh"
/*
** just segretto fit
** */
// time is in nanoseconds
using namespace TMath;

double tTripletS = 2100.0;

static double expGaus(double x, double tau)
{
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * Exp(arg1) * Erfc(arg2);
  return f;
}

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
}

// add quenching
// par binwidth norm PPM tau3 kplus type 
//     0        1     2   3     4     5
//
//
void getIntegrals(TH1D* hWave, double base, double& singInt,double& tripInt, double& tot, double& totTriplet ) 
{
  /*
  static double sStart = 900.;
  static double sEnd = 1020.;
  static double tStart = 1200;
  static double tTripletCut = 2300.;
  static double tEnd = 10000;
  */
  int istart = hWave->FindBin(sStart);
  int iend = hWave->FindBin(sEnd);
  int jstart = hWave->FindBin(tStart);
  int jend = hWave->FindBin(tTripletCut);
  int kend = hWave->FindBin(tEnd);
  double baseSum  = hWave->Integral(0,istart,"width") - hWave->GetBinWidth(1)*base*double(istart);
  singInt = hWave->Integral(istart, iend,"width") - hWave->GetBinWidth(1)*base*double(iend-istart);
  tripInt = hWave->Integral(jstart, jend,"width") - hWave->GetBinWidth(1)*base*double(jend-jstart);
  tot = hWave->Integral(istart, kend,"width") - hWave->GetBinWidth(1)*base*double(kend-istart);
  totTriplet = hWave->Integral(iend, hWave->GetNbinsX(),"width") - hWave->GetBinWidth(1)*base*double(hWave->GetNbinsX()-iend);

  printf(" getIntegrals %s base (%f) %f %f %f %f %f \n", hWave->GetName(),base, baseSum , singInt,tripInt,tot,totTriplet );

}


static Double_t lightModel(Double_t *xx, Double_t *par)
{
  double f = 0;
  double x = xx[0] - xTrigger;
  //if (x < 0)
    //return f;
  double bw = par[0];
  double norm = par[1];
  double kx = k1Zero*par[2];
  double tau3 = par[3];
  double kp = par[4]*kplusZero;
  double sfrac = par[5];
  int type = int(par[6]);
  double td = 1/kx;
  double tq = 1./(1./tau3+kp);
  double tr = 1./(kx + 1./tq);
  double alpha1 = bw*singletFraction*norm;
  double alpha3 = bw*(1.-singletFraction)*norm;
  double fs = alpha1 / tSinglet * expGaus(x, tSinglet);
  double f3 = alpha3 / tau3 * expGaus(x, tr);
  double fx = alpha3 * pow(kx, 2) * tq * (expGaus(x, td) - expGaus(x, tr));
  //double fs =  alpha1/tSinglet*Exp( -x/tSinglet);
  //double f3 =  alpha3/tau3*Exp(-x/tr);
  //double fx =  alpha3*pow(kx,2)*tq*( Exp( -x/td) - Exp( -x/tr) );
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



class tbTFitAll
{
public:
    //double taut = 2100.;
    tbTFitAll(bool fix = true, double lifetime3 = tTriplet); //2.93); 3.48
    virtual ~tbTFitAll() { ; }
    TFile *inFile;
    TFile *fout;
    enum
    {
        NSETS = 6
    };
    enum
    {
        MAXFITS = 6
    };
    enum
    {
        npars =7 
    };

    TString outTag;
    int setTotal[NSETS];
    int setGood[NSETS];
    double  setBase[NSETS];


    double fitPPM[NSETS];
    double fitPPMErr[NSETS];

    double mSingletIntegral[NSETS];
    double mTripletIntegral[NSETS];
    double mTotalIntegral[NSETS];
    double mTripletTotalIntegral[NSETS];

    double mSingletIntErr[NSETS];
    double mTripletIntErr[NSETS];
    double mTripletTotalIntErr[NSETS];
    double mTotalIntErr[NSETS];

    double baseline[NSETS];




    TF1 *setFitD(int iset, TH1D *hLife, double ppm);
    TF1 *tailFit(int iset, TH1D *hTail);
    void tbFit1(int iset, TH1D *hLife);
    TCanvas *fhistCan(TString title, TH1D *h, TF1 *fp);
    void rangeIntegral(TH1D *h, double xlow, double high, double& sum, double& err);
    void makeModelGraphs(TH1D** hList,TString tag);
    void makeDataGraphs(TH1D** hList,TString tag);


    void writeHist(TH1D *h);
    ofstream textFile;

    TString fitTag;
    TF1 *fp[NSETS];
    TF1 *fpD[NSETS];
    TF1 *pmodelFit[NSETS][NSETS];
    TF1 *fmodelFit[NSETS][MAXFITS];
    TCanvas *canFmodel[NSETS];
    TCanvas *canFit[NSETS];

    TString runTag[NSETS];
    TString runRange[NSETS];
    TH1D *hMpvSet[NSETS];
    TH1D *hLifeSet[NSETS];
    TH1D *hLifeInt[NSETS];
    TH1D* hFitHist[NSETS];
    TH1D* hFitDHist[NSETS];

    TH1D *hMpvSum[NSETS];
    TH1D *hBase[NSETS];
    int pstyle[NSETS];
    int setColor[NSETS];
    int setStyle[NSETS];
    double runPPM[NSETS];
    double tauGSet[NSETS];
    double setChi[NSETS];
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

    // graphs
    TGraphErrors *gSinglet;
    TGraphErrors *gTriplet;
    TGraphErrors *gTotTriplet;
    TGraphErrors *gTot;

    TGraphErrors *dSinglet;
    TGraphErrors *dTriplet;
    TGraphErrors *dTotTriplet;
    TGraphErrors *dTot;

};


void tbTFitAll::makeDataGraphs(TH1D** hList,TString tag)
{

  double  singInt, tripInt,  tot,  totTriplet;

  double dSingletIntegral[NSETS];
  double dTripletIntegral[NSETS];
  double dTotalIntegral[NSETS];
  double dTripletTotalIntegral[NSETS];

  double dSingletIntErr[NSETS];
  double dTripletIntErr[NSETS];
  double dTripletTotalIntErr[NSETS];
  double dTotalIntErr[NSETS];



  for(int iset=0; iset<NSETS; ++iset) {
    if(hList[iset]==NULL) continue;
    getIntegrals(hList[iset],baseline[iset],singInt,tripInt,tot,totTriplet );
    dSingletIntegral[iset] = singInt;
    dTripletIntegral[iset] = tripInt;
    dTotalIntegral[iset] = tot;
    dTripletTotalIntegral[iset] = totTriplet;

    dSingletIntErr[iset] = sqrt(singInt);
    dTripletIntErr[iset] = sqrt(tripInt);
    dTotalIntErr[iset] = sqrt(tot);
    dTripletTotalIntErr[iset] = sqrt(totTriplet);

  }

  fitPPM[0]=0.0;

   //ratios
    double tSingletIntegral[NSETS];
    double tSingletIntErr[NSETS];
    double tTripletIntegral[NSETS];
    double tTripletIntErr[NSETS];
    double tTotalIntegral[NSETS];
    double tTotalIntErr[NSETS];
    double tTripletTotalIntegral[NSETS];
    double tTripletTotalIntErr[NSETS];

     for (unsigned iset = 0; iset < NSETS; ++iset)
    {
      ratioE(dSingletIntegral[iset],dSingletIntErr[iset], dSingletIntegral[0],   dSingletIntErr[0],    
             tSingletIntegral[iset],tSingletIntErr[iset]);

      ratioE(dTripletIntegral[iset],  dTripletIntErr[iset], dSingletIntegral[iset], dSingletIntErr[iset],
             tTripletIntegral[iset],  tTripletIntErr[iset]);

      ratioE(dTotalIntegral[iset],  dTotalIntErr[iset], dSingletIntegral[iset], dSingletIntErr[iset],
             tTotalIntegral[iset],  tTotalIntErr[iset]);

      ratioE(dTripletTotalIntegral[iset], dTripletTotalIntErr[iset], dSingletIntegral[iset], dSingletIntErr[iset], 
             tTripletTotalIntegral[iset], tTripletTotalIntErr[iset]);
    }




  dSinglet = new TGraphErrors(NSETS-1, fitPPM, tSingletIntegral,fitPPMErr, tSingletIntErr);
  dSinglet->SetName(Form("RangeSingletModel%s",tag.Data()));
  dSinglet->SetTitle(Form("data singlet range integral %0.2f - %0.2f ", sStart, sEnd));
  dSinglet->SetMarkerColor(kRed);
  dSinglet->SetMarkerStyle(22);
  dSinglet->SetMarkerSize(1.4);

  dTriplet = new TGraphErrors(NSETS-1,fitPPM, tTripletIntegral,fitPPMErr, tTripletIntErr);
  dTriplet->SetName(Form("RangeTripletModel%s",tag.Data()));
  dTriplet->SetTitle(Form("data triplet range integral %0.2f - %0.2f ", tStart, tTripletCut));
  dTriplet->SetMarkerColor(kGreen);
  dTriplet->SetMarkerStyle(21);
  dTriplet->SetMarkerSize(1.4);

  dTotTriplet = new TGraphErrors(NSETS-1, fitPPM, tTripletTotalIntegral,fitPPMErr, tTripletTotalIntErr );
  dTotTriplet->SetName(Form("TotalTripletModel%s",tag.Data()));
  dTotTriplet->SetTitle(Form("data total triplet  %0.2f - %0.2f ", sEnd, tEnd));
  dTotTriplet->SetMarkerColor(kGreen + 2);
  dTotTriplet->SetMarkerStyle(21);
  dTotTriplet->SetMarkerSize(1.4);

  dTot = new TGraphErrors(NSETS-1, fitPPM, tTotalIntegral,fitPPMErr,tTotalIntErr );
  dTot->SetName(Form("TotalModel%s",tag.Data()));
  dTot->SetTitle(Form("data total integral %0.2f - %0.2f ", sStart, tEnd));
  dTot->SetMarkerColor(kBlue);
  dTot->SetMarkerStyle(21);
  dTot->SetMarkerSize(1.4);

 
}



void tbTFitAll::makeModelGraphs(TH1D** hList,TString tag)
{
  /*
  static double sStart = 900.;
  static double sEnd = 1020.;
  static double tStart = 1200;
  static double tTripletCut = 2300.;
  static double tEnd = 10000;
  */

  double  singInt, tripInt,  tot,  totTriplet;

  for(int iset=0; iset<NSETS; ++iset) {
    if(hList[iset]==NULL) continue;
    getIntegrals(hList[iset],baseline[iset],singInt,tripInt,tot,totTriplet );
    mSingletIntegral[iset] = singInt;
    mTripletIntegral[iset] = tripInt;
    mTotalIntegral[iset] = tot;
    mTripletTotalIntegral[iset] = totTriplet;

    mSingletIntErr[iset] = sqrt(singInt);
    mTripletIntErr[iset] = sqrt(tripInt);
    mTotalIntErr[iset] = sqrt(tot);
    mTripletTotalIntErr[iset] = sqrt(totTriplet);

  }

  fitPPM[0]=0.0;


    //ratios
    double tSingletIntegral[NSETS];
    double tSingletIntErr[NSETS];
    double tTripletIntegral[NSETS];
    double tTripletIntErr[NSETS];
    double tTotalIntegral[NSETS];
    double tTotalIntErr[NSETS];
    double tTripletTotalIntegral[NSETS];
    double tTripletTotalIntErr[NSETS];

    for (unsigned iset = 0; iset < NSETS; ++iset)
    {
      ratioE(mSingletIntegral[iset],mSingletIntErr[iset],     mSingletIntegral[0],   mSingletIntErr[0],    
             tSingletIntegral[iset],tSingletIntErr[iset]);

      ratioE(mTripletIntegral[iset],  mTripletIntErr[iset],  mSingletIntegral[iset], mSingletIntErr[iset],
             tTripletIntegral[iset],  tTripletIntErr[iset]);

      ratioE(mTotalIntegral[iset],  mTotalIntErr[iset],   mSingletIntegral[iset], mSingletIntErr[iset],
             tTotalIntegral[iset],  tTotalIntErr[iset]);

      ratioE(mTripletTotalIntegral[iset], mTripletTotalIntErr[iset], mSingletIntegral[iset], mSingletIntErr[iset], 
             tTripletTotalIntegral[iset], tTripletTotalIntErr[iset]);
    }


  gSinglet = new TGraphErrors(NSETS-1, fitPPM, tSingletIntegral,fitPPMErr, tSingletIntErr);
  gSinglet->SetName(Form("RangeSingletModel%s",tag.Data()));
  gSinglet->SetTitle(Form("model singlet range integral %0.2f - %0.2f ", sStart, sEnd));
  gSinglet->SetMarkerColor(kRed);
  gSinglet->SetMarkerStyle(30);
  gSinglet->SetMarkerSize(1.4);

  gTriplet = new TGraphErrors(NSETS-1,fitPPM, tTripletIntegral,fitPPMErr, tTripletIntErr);
  gTriplet->SetName(Form("RangeTripletModel%s",tag.Data()));
  gTriplet->SetTitle(Form("model triplet range integral %0.2f - %0.2f ", tStart,  tTripletCut));
  gTriplet->SetMarkerColor(kGreen);
  gTriplet->SetMarkerStyle(25);
  gTriplet->SetMarkerSize(1.4);

  gTotTriplet = new TGraphErrors(NSETS-1, fitPPM, tTripletTotalIntegral,fitPPMErr, tTripletTotalIntErr );
  gTotTriplet->SetName(Form("TotalTripletModel%s",tag.Data()));
  gTotTriplet->SetTitle(Form("model total triplet  %0.2f - %0.2f ", sEnd, tEnd));
  gTotTriplet->SetMarkerColor(kGreen + 2);
  gTotTriplet->SetMarkerStyle(27);
  gTotTriplet->SetMarkerSize(1.4);

  gTot = new TGraphErrors(NSETS-1, fitPPM, tTotalIntegral,fitPPMErr,tTotalIntErr );
  gTot->SetName(Form("TotalModel%s",tag.Data()));
  gTot->SetTitle(Form("model total integral %0.2f - %0.2f ", sStart, tEnd));
  gTot->SetMarkerColor(kBlue);
  gTot->SetMarkerStyle(26);
  gTot->SetMarkerSize(1.4);

}



void tbTFitAll::writeHist(TH1D *h)
{
    textFile << " hist " << h->GetName() << "  title " << h->GetTitle() << " y axis title " << h->GetYaxis()->GetTitle() << endl;

    for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
        textFile << ibin << "  " << h->GetBinContent(ibin) << "  " << h->GetBinError(ibin) << endl;
}

void tbTFitAll::rangeIntegral(TH1D *h, double xlow, double xhigh, double& sum, double& err)
{
    sum = 0;
    err = 0;
    int i1 = h->FindBin(xlow);
    int i2 = h->FindBin(xhigh);
    for (int ibin = i1; ibin < i2; ++ibin) {
        sum += h->GetBinContent(ibin);
        err += pow(h->GetBinContent(ibin),2.);
    }
    err = sqrt(err);
}

TCanvas *tbTFitAll::fhistCan(TString title, TH1D *h, TF1 *fp)
{
    TCanvas *can = new TCanvas(title, title);
    //can->SetLogy();
    gStyle->SetOptFit(0011);
    //if(fp) fp->SetLineColor(kTeal-6);
    //h->GetYaxis()->SetRangeUser(1E-1,1E3);
    h->GetXaxis()->SetRangeUser(.9,6);
    h->Draw();
    if (fp)
        fp->Draw("sames");
    can->Print(".pdf");
    return can;
}

TF1 *tbTFitAll::tailFit(int iset, TH1D *hTail)
{
  TF1 *ftail=NULL;
  return ftail;
}

/*

TF1 *tbTFitAll::tailFit(int iset, TH1D *hTail)
{
    // fit exponential tail
    double tauG = 3.140;
    double binwidth = hTail->GetBinWidth(1);
    double tailIntegral = rangeIntegral(hTail, tailStart, tailStop);
    TF1 *ftail = new TF1(Form("fexpSet-%i", iset), fexp, sEnd, tTripletCut, 3);
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

TF1 *tbTFitAll::setFitD(int iset, TH1D *hLife, double ppm)
{
    Double_t binwidth = hLife->GetBinWidth(1);
    // ppm is by number
    Double_t integral, integralErr;
    rangeIntegral(hLife, sEnd, tailStart, integral, integralErr);
    double norm = integral; //-G0Par;
    int ifit = 3;
    double kplus = 1;// kplusZero;

    //TF1* fp  = new TF1(Form("LifeFit-%i-%0.f-PPM", iset, ppm),lightModel,sEnd,tailStop,npars);

    TF1* fp  = new TF1(Form("LifeFit-%i-%0.f-PPM", iset, ppm),lightModel,xTrigger-100.,tailStart,npars);
    printf(" tbFit1 set %i fit range %f to %f k+ %E \n", iset,xTrigger-100.,tailStart,kplus);
    fp->SetParName(0,"binw");
    fp->SetParName(1,"norm");
    fp->SetParName(2,"PPM");
    fp->SetParName(3,"tau3");
    fp->SetParName(4,"k+");
    fp->SetParName(5,"sfrac");
    fp->FixParameter(0, binwidth);
    fp->SetParameter(1, integral);
    fp->SetParameter(2, ppm);
    fp->SetParLimits(2,0.,100.);
    //fp->FixParameter(2, ppm);
    fp->FixParameter(3, tTripletS);
    fp->FixParameter(4,kplus);
    fp->FixParameter(5, singletFraction);
    //fp->SetParLimits(5,0.,1.);
    fp->FixParameter(6,ifit);
    fp->SetTitle(Form("LfFit-%i-%0.f-PPM", iset, ppm));
    fp->SetNpx(1000); // numb points for function
    /* do the fit here */
    TFitResultPtr fptr = hLife->Fit(fp, "LE0S+", "", xTrigger, tailStart);
    TMatrixDSym cov = fptr->GetCorrelationMatrix();
    printf(" correlation set %i cov(2,3) %f \n", iset, cov(2, 3));
    fptr->Print("V");
    fp->Print();
    fp->SetLineColor(setColor[iset]);
    fp->SetMarkerStyle(setStyle[iset]);

    printf(" \n\n >>> setFitD fitted parameters set %i \n",iset);
    for (int ii = 0; ii < npars; ++ii)
    {
        printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
    }
    // fill arrays
    return fp;
}
void tbTFitAll::tbFit1(int iset, TH1D *hLife)
{
    double markerSize = 0.5;
    vset.push_back(runPPM[iset]);
    vsetErr.push_back(0);

    /***** fitting *****/
    Double_t binwidth = hLife->GetBinWidth(1);

    Double_t integral, integralErr;
    rangeIntegral(hLife, sEnd, tTripletCut, integral, integralErr);
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
        for(int ii=0; ii<3; ++ii) {
          printf("\t  param %i %s %.3f +/- %.3f \n",ii,ftail->GetParName(ii),ftail->GetParameter(ii),ftail->GetParError(ii));
        }
        TString tailtitle;
        tailtitle.Form("TailFitToSet%i", iset);
        ftail->SetLineColor(kRed);
        ftail->SetLineWidth(2);
        fhistCan(tailtitle, hTail, ftail);
        //
        double sum = ftail->Integral(tailStart, tailStop) / hTail->GetBinWidth(1);
        tauG = ftail->GetParameter(0);
        tauGErr = ftail->GetParError(0);
        G0Par = ftail->Integral(sEnd, tTripletCut) / hTail->GetBinWidth(1);
        gfrac = G0Par / integral;
        printf(" >> ftail set %i tauG %.3f lifesum %.2E sum %.2E G0Par %.2E gfrac %.3f\n",
               iset, tauG, integral, sum, G0Par, gfrac);
    }
    vtaugAll.push_back(tauG);
    vtaugAllErr.push_back(tauGErr);
    vg0All.push_back(G0Par);
    vg0AllErr.push_back(sqrt(G0Par));

    fpD[iset] = setFitD(iset, hLife, runPPM[iset]);

    // use default segreto model PPM in graphs
    setChi[iset] = fpD[iset]->GetChisquare();
    fitPPM[iset] = fpD[iset]->GetParameter(2);
    fitPPMErr[iset] = fpD[iset]->GetParError(2);


    //hFitHist[iset] =  (TH1D*) fp[iset]->GetHistogram();
    //hFitHist[iset]->SetName(Form("FitHist-%i",iset));
    //fout->Append(hFitHist[iset]);

    hFitHist[iset] = (TH1D*) fpD[iset]->GetHistogram();
    hFitHist[iset]->SetName(Form("FitHist-%i",iset));
    fout->Append(hFitDHist[iset]);

     

    // triplet integral
    va0All.push_back(fpD[iset]->GetParameter(1));
    va0AllErr.push_back(fpD[iset]->GetParError(1));
    vppmAll.push_back(fitPPM[iset]);
    vppmAllErr.push_back(fitPPMErr[iset]);

    
    //fitted lifetime 
    vtautAll.push_back(fpD[iset]->GetParameter(3));
    vtautAllErr.push_back(fpD[iset]->GetParError(3));

    vk1All.push_back(fpD[iset]->GetParameter(4));
    vk1AllErr.push_back(fpD[iset]->GetParError(4));

    double g0frac;
    double g0fracErr;
    ratioE(fpD[iset]->GetParameter(5), fpD[iset]->GetParError(5), fpD[iset]->GetParameter(1), fpD[iset]->GetParError(1), g0frac, g0fracErr);
    vgfracAll.push_back(g0frac);
    vgfracAllErr.push_back(g0fracErr);

    double lightAll = fpD[iset]->Integral(sEnd, tTripletCut) / binwidth;
    printf(" set %i lightAll %.3E \n", iset, lightAll);

    //gStyle->cd();
   
    TCanvas* canF = new TCanvas(Form("LifeFitLinear-%.f-PPM-%s", runPPM[iset], runTag[iset].Data()), Form("LifeFit-%.f-PPM-%s", runPPM[iset], runTag[iset].Data()));
    canFit[iset]=canF;
    gStyle->SetOptFit(11);
    gStyle->SetOptStat(1);

    canF->SetLogy(0);
    canF->SetGridx();
    canF->SetGridy();
    //hLife->SetTitle("");
    hLife->GetYaxis()->SetRangeUser(0.,2000.);
    hLife->GetXaxis()->SetRangeUser(900.,3000.);
    hLife->SetTitle(Form("LifeFit set %i Xe %.0f PPM ", iset, runPPM[iset]));
    hLife->GetXaxis()->SetTitle("time (nanoseconds)");
    //hLife->GetYaxis()->SetRangeUser(1E-1,1E3);
    hLife->SetMarkerSize(0.5);
    hLife->SetLineColor(setColor[iset]);
    hLife->SetMarkerColor(setColor[iset]);

    hLife->Draw("p");
    //fp[iset]->SetLineColor(kBlack);
    //fp[iset]->SetLineStyle(1);
    //fp[iset]->SetLineWidth(2);
    //fp[iset]->Draw("sames");
    fpD[iset]->SetLineColor(kBlack);
    fpD[iset]->SetLineStyle(1);
    fpD[iset]->SetLineWidth(2);
    fpD[iset]->Draw("sames");

    if(ftail) ftail->Draw("sames");
    TString ppmText;
    ppmText.Form("PPM (default)  %.2f +/- %.2f  ",fpD[iset]->GetParameter(2), fpD[iset]->GetParError(2));

    // add line
    canF->Update();
    TPaveText *stats = (TPaveText*) canFit[iset]->GetPrimitive("stats");
    if(stats) { 
      cout << " \t\t\t adding to stats " << ppmText << "to  " << hLife->GetName() << "  can " << canFit[iset]->GetName()  << endl;
      TList *listOfLines = stats->GetListOfLines();
      TText *tconst = stats->GetLineWith("Entries");
      listOfLines->Remove(tconst);
      tconst = stats->GetLineWith("Mean");
      listOfLines->Remove(tconst);
      tconst = stats->GetLineWith("Std Dev");
      listOfLines->Remove(tconst);
      stats->SetName("Mystats");
      //TLatex *myt = new TLatex(0,0,ppmText);
      //myt->SetTextFont(42);
      //myt->SetTextSize(0.06);
      //listOfLines->Add(myt);
      //stats->AddText(ppmText);
      //stats->SetTextColor(kBlue);
    } 

    hLife->SetStats(0);
    //canFit[iset]->Update();
    canF->Modified();
    canF->Print(".pdf");

    double ArSum = fpD[iset]->Integral(sEnd, tTripletCut);
    double XeSum = fpD[iset]->Integral(sEnd, tTripletCut);
    double totSum = ArSum + XeSum;
    double XeFrac = XeSum / totSum;
    double XeFracErr = sqrt(XeFrac * (1. - XeFrac) / totSum);
    //ratioE(XeSum,sqrt(XeSum), totSum,sqrt(totSum), XeFrac, XeFracErr);

    vXeFrac.push_back(XeFrac);
    vXeFracErr.push_back(XeFracErr);

}
// MAIN 
tbTFitAll::tbTFitAll(bool fix, double lifetime3)
{

    double markerSize = 0.5;
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

   

  for(int iset=0; iset<NSETS; ++iset) {
    hMpvSet[iset]=NULL;
    hLifeSet[iset]=NULL;
    hLifeInt[iset]=NULL;
    hFitHist[iset]=NULL;
    hFitDHist[iset]=NULL;
  }

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



    setColor[0] = kOrange + 4;
    setColor[1] = kBlue;
    setColor[2] = kRed;
    setColor[3] = kGreen - 2;
    setColor[4] = kCyan - 2;
    setColor[5] = kMagenta + 1;


    TString inFileName("sumRuns.root");
    inFile = new TFile(inFileName,"READONLY");
    if (!inFile->IsOpen()) {
      cout << " could not open " << inFileName << endl;  
      return;
    }
    printf(" opening file %s \n",inFile->GetName());
    fout = new TFile("tbTFitOut.root", "RECREATE");



    TIter next(inFile->GetListOfKeys());
    TKey *key;
    TH1D *hlife = NULL;
    int theSet = 0;
    while ((key = (TKey *)next()))
    {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1"))
        continue;
      TH1D *h = (TH1D *)key->ReadObj();
      TString hname(h->GetName());
      if (hname.Contains("Wave")){
        hLifeSet[theSet] = h;
        cout << " appending set " << theSet << "  " << hname << "  " << h->GetTitle() << endl;
        fout->Append(h);
        ++theSet;
      }
    }

    if(theSet<NSETS) {
      printf(" returning theSet = %i\n",theSet);
      return;
    }
    printf(" # sets %i nbins is  %i binwidth %f  \n",theSet, hLifeSet[0]->GetNbinsX(),  hLifeSet[0]->GetBinWidth(1));

    // baseline fits 
    for (int iset = 0; iset < NSETS-1 ; ++iset){
      //hLifeSet[iset]->Fit("pol1","","",0,900);
      //baseline[iset]=hLifeSet[iset]->GetFunction("pol1")->GetParameter(0);
      baseline[iset]= hLifeSet[iset]->Integral(0,900,"width")/900.;  
    }
    for (int iset = 0; iset < NSETS-1 ; ++iset){
      printf("set %i base %f \n",iset,baseline[iset]);
    }


    TCanvas *canLifeAll = new TCanvas("LifeALL", "LifeALL");
    canLifeAll->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    for (int iset = 0; iset < NSETS-1 ; ++iset)
    {
      hLifeSet[iset]->SetTitle("");
      hLifeSet[iset]->SetStats(1);
      hLifeSet[iset]->GetYaxis()->SetRangeUser(4.0,2.E3);
      hLifeSet[iset]->SetTitle(Form("sets %i ",iset));
      hLifeSet[iset]->GetYaxis()->SetTitle("yield SPE/40 ns");
      hLifeSet[iset]->GetXaxis()->SetTitle("time (nanoseconds)");
      hLifeSet[iset]->SetMarkerStyle(7);
      hLifeSet[iset]->SetMarkerSize(1);
      hLifeSet[iset]->SetMarkerColor(setColor[iset]);
      hLifeSet[iset]->SetLineColor(setColor[iset]);
      hLifeSet[iset]->SetMarkerSize(markerSize);
      if(iset==0) hLifeSet[iset]->Draw("");
      else hLifeSet[iset]->Draw("same");
    }
    canLifeAll->BuildLegend();
    canLifeAll->Print(".pdf");

    // singlet comparison
    TCanvas *canSingletAll = new TCanvas("SingletALL", "SingletALL");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();

    TH1D* hLifeSinglet[NSETS];
    canSingletAll->SetGridx(); canSingletAll->SetGridy();

    for (int iset = 0; iset < NSETS-1 ; ++iset) {
      hLifeSinglet[iset] = (TH1D*) hLifeSet[iset]->Clone(Form("LifeSinglet%i",iset));
      hLifeSinglet[iset]->SetMarkerSize(1.4);
      hLifeSinglet[iset]->SetLineWidth(2);

      hLifeSinglet[iset]->GetYaxis()->SetRangeUser(0.01,1.2E3);
      hLifeSinglet[iset]->GetXaxis()->SetRangeUser(980,1100);
      if(iset==0) hLifeSinglet[iset]->Draw("HISTL");
      else hLifeSinglet[iset]->Draw("HISTL same");
    }
    canSingletAll->BuildLegend();
    canSingletAll->Print(".pdf");




    for (int iset = 0; iset < NSETS; ++iset)
        tbFit1(iset, hLifeSet[iset]);


    printf(" XXXXXX makeModelGraphs   \n");   
    outTag=TString("segreto");
    makeModelGraphs(hFitHist,outTag);

    printf(" XXXXXX makeDataGraphs   \n");   
    makeDataGraphs(hLifeSet,"data");



    // integral plots
    int startBin = hLifeSet[0]->FindBin(900.);
    for (int i = 0; i < NSETS; ++i)
    {
      double yieldSum=0;
      double yieldSumErr2=0;
      hLifeInt[i] = (TH1D *)hLifeSet[i]->Clone(Form("IntegralToFitSet%i-%i", i, int(runPPM[i])));
      hLifeInt[i]->Reset("ICESM");
      hLifeInt[i]->SetTitle(Form("integral set %i PPM %i ", i, int(runPPM[i])));
      printf(" set %i %s bins %i \n",i,hLifeInt[i]->GetName(),hLifeInt[i]->GetNbinsX());
      for (int ibin = 0; ibin < hLifeInt[i]->GetNbinsX(); ++ibin)
      {
        hLifeInt[i]->SetBinContent(ibin,0);
        hLifeInt[i]->SetBinError(ibin, 0);
        if(ibin<startBin) continue;
        double val = hLifeSet[i]->GetBinContent(ibin);
        yieldSum += val;
        yieldSumErr2 += hLifeSet[i]->GetBinError(ibin);
        hLifeInt[i]->SetBinContent(ibin, yieldSum);
        hLifeInt[i]->SetBinError(ibin, sqrt(yieldSumErr2));
      }
    }


    TCanvas *canIntAll = new TCanvas("IntALL", "IntALL");
    gStyle->SetOptStat(0);
    canIntAll->SetGridx(); canIntAll->SetGridy();
    for (int iset = 0; iset < NSETS-1; ++iset)
    {
      hLifeInt[iset]->SetTitle("");
      hLifeInt[iset]->SetStats(1);
      hLifeInt[iset]->GetYaxis()->SetTitle("integral yield SPE/40 ns");
      hLifeInt[iset]->GetXaxis()->SetTitle("time (nanoseconds)");
      hLifeInt[iset]->GetXaxis()->SetRangeUser(startBin,6000);
      hLifeInt[iset]->GetYaxis()->SetRangeUser(0,350E3);
      hLifeInt[iset]->SetLineWidth(1);
      hLifeInt[iset]->SetMarkerStyle(7);
      hLifeInt[iset]->SetMarkerColor(setColor[iset]);
      hLifeInt[iset]->SetLineColor(setColor[iset]);
      hLifeInt[iset]->SetMarkerSize(markerSize);
      if(iset==0) hLifeInt[iset]->Draw("");
      else  hLifeInt[iset]->Draw("sames");
    }
    canIntAll->BuildLegend();
    canIntAll->Print(".pdf");


    TGraphErrors *gTriplet0 = new TGraphErrors(vset.size()-1, &vset[0], &va0All[0], &vsetErr[0], &va0AllErr[0]);
    TGraphErrors *gTripletMu = new TGraphErrors(1,&vset[NSETS-1], &va0All[NSETS-1], &vsetErr[NSETS-1], &va0AllErr[NSETS-1]);
    TMultiGraph *mgTriplet = new TMultiGraph();

    TCanvas *cTriplet0 = new TCanvas("Integral-Fit-Value","Integral-Fit-Value");
    cTriplet0->SetGridx(); cTriplet0->SetGridy();
    gTriplet0->SetName("TripletFit");
    gTriplet0->SetTitle("Triplet-Fit-By-PPM");
    gTripletMu->SetTitle("Triplet-Fit-By-PPM-Mu");
    gTriplet0->SetMarkerColor(kRed);
    gTriplet0->SetMarkerStyle(22);
    gTriplet0->SetMarkerSize(1.4);
    gTripletMu->SetName("TripletFitMu");
    gTripletMu->SetMarkerColor(kBlue);
    gTripletMu->SetMarkerStyle(21);
    gTripletMu->SetMarkerSize(1.4);
    gTriplet0->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    gTriplet0->GetHistogram()->GetYaxis()->SetTitle(" fitted Triplet yield");
    mgTriplet->Add(gTriplet0,"p");
    //mgTriplet->Add(gTripletMu,"p");
    mgTriplet->SetTitle("Total Light Yield; Xe PPM dopant ; fitted Triplet yield");
    mgTriplet->Draw("a");
    cTriplet0->Print(".pdf");

    
    TGraphErrors *gPPM = new TGraphErrors(vset.size()-1, &vset[0], &vppmAll[0], &vsetErr[0], &vppmAllErr[0]);

    TF1 *fitppm  = new TF1("fitppm", "[0]+[1]*x",1,10);
    fitppm->SetParName(0, "intercept");
    fitppm->SetParName(1, "slope");
    fitppm->SetLineColor(kBlack);
    TFitResultPtr fitppmResult=gPPM->Fit("fitppm","R");

    printf(" \n\n >>> fit parameters \n");
    for (int ii = 0; ii < 2; ++ii)
    {
      printf("\t  param %i %s %.3f +/- %.3f \n", ii, fitppm->GetParName(ii), fitppm->GetParameter(ii), fitppm->GetParError(ii));
    }



    TGraphErrors *gPPMMu = new TGraphErrors(1, &vset[NSETS-1], &vppmAll[NSETS-1], &vsetErr[NSETS-1], &vppmAllErr[NSETS-1]);
    TMultiGraph *mgppm = new TMultiGraph();



    TCanvas *cFitPPM = new TCanvas("FitPPM","FitPPM");
    cFitPPM->SetGridx(); cFitPPM->SetGridy();
    gPPM->SetName("FitPPM");
    gPPMMu->SetName("FitPPMMu");
    gPPM->SetTitle("FitPPM-By-PPM");
    gPPMMu->SetTitle("FitPPM-By-PPM-Mu");
    gPPM->SetMarkerColor(kRed);
    gPPM->SetMarkerStyle(22);
    gPPM->SetMarkerSize(1.4);
    gPPMMu->SetMarkerColor(kBlue);
    gPPMMu->SetMarkerStyle(25);
    gPPMMu->SetMarkerSize(1.4);
    gPPM->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    gPPM->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    mgppm->Add(gPPM,"p");
    //mgppm->Add(gPPMMu,"p");
    mgppm->SetTitle("Fitted PPM; Xe PPM dopant ; PPM fitted");
    mgppm->Draw("a");
    //mgppm->GetYaxis()->SetRangeUser(0,12);
    //mgppm->GetXaxis()->SetRangeUser(0,12);
    mgppm->SetTitle("Fitted PPM; Xe PPM dopant ; PPM fitted");
    mgppm->Draw("a");
    cFitPPM->Print(".pdf");


    
    /*/ kplus fit
    TGraphErrors *gk1 = new TGraphErrors(vset.size() - 1, &vset[0], &vk1All[0], &vsetErr[0], &vk1AllErr[0]);
    TGraphErrors *gk1Mu = new TGraphErrors(1, &vset[NSETS - 1], &vk1All[NSETS - 1], &vsetErr[NSETS - 1], &vk1AllErr[NSETS - 1]);
    TMultiGraph *mgk1 = new TMultiGraph();

    TCanvas *cFitk1 = new TCanvas("FitkPlus", "FitkPlus");
    cFitk1->SetGridx();
    cFitk1->SetGridy();
    gk1->SetName("FitkPlus");
    gk1Mu->SetName("FitkPlusMu");
    gk1->SetTitle("FitkPlus-By-kPlus");
    gk1Mu->SetTitle("FitkPlus-By-kPlus-Mu");
    gk1->SetMarkerColor(kRed);
    gk1->SetMarkerStyle(22);
    gk1->SetMarkerSize(1.4);
    gk1Mu->SetMarkerColor(kBlue);
    gk1Mu->SetMarkerStyle(21);
    gk1Mu->SetMarkerSize(1.4);
    gk1->GetHistogram()->GetXaxis()->SetTitle(" doped PPM");
    gk1->GetHistogram()->GetYaxis()->SetTitle(" fitted k+");
    mgk1->Add(gk1, "p");
    //mgk1->Add(gk1Mu, "p");
    mgk1->SetTitle("Fitted k+; Xe dopant PPM ; fitted triplet lifetime");
    mgk1->Draw("a");
    cFitk1->Print(".pdf");
    */

    for(int iset=0; iset<NSETS; ++iset) {
      printf("//........ fitted parameters set %i\n", iset);
      for(int j=0; j<npars; ++j) printf(" par[%i][%i] = %f ;\n",iset,j,fpD[iset]->GetParameter(j));
    }


    

    fout->Append(gTriplet0);
    fout->Append(gTripletMu);
    fout->Append(gPPM);
    fout->Append(gPPMMu);

   
   
    // comparison plot
    double microToNano=1.E3;
    TMultiGraph *mgAfterCuts = new TMultiGraph();
    TString canname;
    canname.Form("set-sums-pass-%s",outTag.Data());
    TCanvas *cTripletAfter = new TCanvas(canname,canname);
    cTripletAfter->SetGridx(); cTripletAfter->SetGridy();
    TString gTitle;
    gTitle.Form("Triplet Light Yield %.3f to %.3f ; Xe PPM dopant ; Yield ",tStart/microToNano,tTripletCut/microToNano);
    mgAfterCuts->SetTitle(gTitle);

    mgAfterCuts->Add(gSinglet);
    mgAfterCuts->Add(gTriplet);
    mgAfterCuts->Add(gTotTriplet);
    mgAfterCuts->Add(dSinglet);
    mgAfterCuts->Add(dTriplet);

    mgAfterCuts->Draw("ap");
    cTripletAfter->BuildLegend();
    cTripletAfter->Print(".pdf");
  
    printf(" %s %s %s %s %s  \n", 
        gSinglet->GetName(),
        gTriplet->GetName(),
        gTotTriplet->GetName(),
        dSinglet->GetName(),
        dTriplet->GetName());


    //fout->ls();
      
    fout->Append(gSinglet);
    fout->Append(gTriplet);
    fout->Append(gTotTriplet);
    fout->Append(gTot);
    fout->Append(dSinglet);
    fout->Append(dTriplet);


    fout->Write();
}


