// time is in microseconds
// model with absorption
#include "compiled/BaconAnalysis.hh"
//static double kplusZero = 1.3E-4;
//static double kxe = 8.8E-5;
using namespace TMath;

double aveEvents = 2607.841463;
double binwidth = 0.04;
double kabsorption=0.0;
enum
{
  NSETS = 5
};
double PPM[NSETS];
double TRIPLET[NSETS];
double NORM[NSETS];
double BASE[NSETS];

int setColor[NSETS];
int setStyle[NSETS];
enum
{
  NPARS = 7 
};

TFile *fout;



double getNorm(int iset)
{
  double ppm=PPM[iset];
  double t3=TRIPLET[iset];
  double alpha3 = (1. - singletFraction);
  double kx = k1Zero * ppm;
  double tq = 1. / (1. /t3 + kplusZero);
  double tr = 1. / (kx + 1. / tq);
  double norm = singletFraction + alpha3 * tr / t3 + alpha3 * kx / (kx + 1.0 / tq);
  printf(" set %i \t ... ppm %.2f %f %f norm %f \n",iset, ppm,alpha3*tr/t3,alpha3*kx/(kx+1.0/tq) ,norm);
  return norm;
}

static double expGaus(double x, double tau)
{
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * Exp(arg1) * Erfc(arg2);
  return f;
}

void printModel(int iset)
{
  double ppm=PPM[iset];
  double tau3=TRIPLET[iset];

  double bw = binwidth;
  double norm = 3.0E7;
  double kx = k1Zero *ppm;
  double kp = kplusZero;
  double kabs = kabsorption*kplusZero*ppm;

  int type = 3;
  double td = 1 / kx;
  double ka = 1/tau3+kx+kp;
  double alpha1 = bw * singletFraction * norm;
  double alpha3 = bw * (1. - singletFraction) * norm;
  printf("set %i ppm %f tau3 %f alpha1 %E alpha3 %E ratio %E \n", iset, ppm, tau3, alpha1, alpha3, alpha3 / alpha1);
  //printf();
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
  double kabs = par[6]*kplusZero*par[2];
  double k3 = 1./par[3];
  double kp = par[4]*kplusZero;
  double kT = k3+kx+kp+kabs;
  double k1 = 1/tSinglet;
  double kS = k1+kx+kp+kabs;
  double sfrac = par[5];
  int type = int(par[7]);
  double alpha1 = bw*singletFraction*norm;
  double alpha3 = bw*(1.-singletFraction)*norm;
  double tS = 1/kS;
  double fs = alpha1 / tSinglet * expGaus(x, tS);
  double tT = 1./kT;
  double f3 = alpha3 *k3* expGaus(x, tT);
  double tx = 1 / kx;
  double a1 = alpha1 *kx*(kx+kabs)/(kS -kx)  * ( expGaus(x, tx)- expGaus(x, tS) ) ;
  double a3 = alpha3 *kx*(kx+kabs)/(kT -kx) * (expGaus(x, tx) - expGaus(x, tT));

  double fx =  a1 + a3;
  fs = fs - a1;
  f3 = f3*k3/kT - a3;
  /*double td = 1/kx;
  double tq = 1./(k3+kp);
  double tr = 1./(kx + 1./tq);
  double fs = alpha1 / tSinglet * expGaus(x, tSinglet);
  double f3 = alpha3*k3* expGaus(x, tr);
  double fx = alpha3 * pow(kx, 2) * tq * (expGaus(x, td) - expGaus(x, tr));
  */
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

TF1 *setFit(int iset)
{
  double tTripletS = 1600;//2100.0;
  int ifit = 3;
  double ppm=PPM[iset];
  double tau3=TRIPLET[iset];
  double integral = NORM[iset];
  double kplus = 1;

  TF1 *fp = new TF1(Form("LifeFit-%i", iset), lightModel, sStart, tEnd, NPARS);
  printf(" tbFit1 set %i fit range %f to %f k+ %E \n", iset,xTrigger-100.,tailStart,kplus);
    fp->SetParName(0,"binw");
    fp->SetParName(1,"norm");
    fp->SetParName(2,"PPM");
    fp->SetParName(3,"tau3");
    fp->SetParName(4,"k+");
    fp->SetParName(5,"sfrac");
    fp->SetParName(6,"ka");
    fp->SetParName(7,"type");

    fp->FixParameter(0, binwidth);
    fp->SetParameter(1, integral);
    //fp->SetParameter(2, ppm);
    //fp->SetParLimits(2,0.1,100.);
    fp->FixParameter(2, ppm);
    fp->FixParameter(3, tTripletS);
    fp->FixParameter(4,kplus);
    //fp->SetParLimits(4,1.E-3,100.);
    fp->FixParameter(5, singletFraction);
    //fp->FixParameter(6,0);
    fp->SetParameter(6,0);
    fp->SetParLimits(6,1.E-4,100.);
    //fp->SetParLimits(5,0.,1.);
    fp->FixParameter(7,ifit);
    fp->SetTitle(Form("LfFit-%i-%0.f-PPM", iset, ppm));
    fp->SetNpx(1000); // numb points for function
    fp->Print();
    fp->SetLineColor(setColor[iset]);
    //fp->SetMarkerStyle(setStyle[iset]);

    printf(" \n\n >>> setFitD fitted parameters set %i \n",iset);
    for (int ii = 0; ii < NPARS; ++ii)
    {
        printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
    }
  return fp;
}

void cmodel(double kc=0.0)
{

  kabsorption = kc;
  printf(" kabsorption = %E \n",kc);

  PPM[0]=0;
  PPM[1]=1;
  PPM[2]=2;
  PPM[3]=5;
  PPM[4]=10;

  for(int iset=0; iset<NSETS; ++iset) {
    TRIPLET[iset]=1600.;
    NORM[iset]= 10000;
  }

  /*
  for(int iset=0; iset<NSETS; ++iset) {
    printf(" ... set %i \n", iset );
    for(int ipar=0; ipar<NPARS; ++ipar)
      printf("\t %i %f \n", ipar, fitpar[iset][ipar]);
  }
  */


  setColor[0] = kOrange + 4;
  setColor[1] = kBlue;
  setColor[2] = kRed;
  setColor[3] = kGreen - 2;
  setColor[4] = kCyan - 2;



  fout = new TFile(Form("gmodel-%0.f.root",TRIPLET[0]) , "RECREATE");
  // if(iset>=NSETS) return;
  TF1 *fp[NSETS];
  for (int iset = 0; iset < NSETS; ++iset)
  {
    fp[iset] = setFit(iset);
  }


  // histograms
  TH1D *hWave[NSETS];
  /*
  for (int iset = 0; iset < NSETS; ++iset)
  {
    hWave[iset] = new TH1D(Form("ModelWaveSet-%i", iset), Form("model wave set %i ppm %.0f tau3 %.0f ", iset,PPM[iset],TRIPLET[iset]), nDigi, 0, nDigi);
    hWave[iset]->SetLineColor(setColor[iset]);
  }
  */

  TString canTitle;
  canTitle.Form("lightModel");
  TCanvas *can = new TCanvas(canTitle, canTitle);
  fp[NSETS - 1]->GetHistogram()->SetXTitle("time nanoseconds");
  fp[NSETS - 1]->Draw();
  for (int iset = 0; iset < NSETS; ++iset)
    fp[iset]->Draw("same");

  can->BuildLegend();
  for (int iset = 0; iset < NSETS; ++iset) fout->Append(fp[iset]);
  gPad->SetLogy();

  for (int iset = 0; iset < NSETS; ++iset)
  {
    printModel(iset);
  }

  /*
  for (int iset = 0; iset < NSETS; ++iset){
    for (int ibin = 0; ibin < hWave[iset]->GetNbinsX(); ++ibin) {
      double x = hWave[iset]->GetBinLowEdge(ibin) + 0.5 * hWave[iset]->GetBinWidth(ibin);
      if (x < sStart)  hWave[iset]->SetBinContent(ibin, BASE[iset]);
      else hWave[iset]->SetBinContent(ibin, fp[iset]->Eval(x)+BASE[iset]);
    }
    fout->Append(hWave[iset]);
  }
  */

  for (int iset = 0; iset < NSETS; ++iset) {
    hWave[iset] = (TH1D*) fp[iset]->GetHistogram();
    hWave[iset]->SetName(Form("ModelWaveSet-%i",iset));
    hWave[iset]->SetLineColor(setColor[iset]);
    fout->Append(hWave[iset]);
  }

  double singletIntegral[NSETS];
  double tripletIntegral[NSETS];
  double totalIntegral[NSETS];
  double tripletTotalIntegral[NSETS];
  double totalNorm[NSETS];

  double binWidth = hWave[0]->GetBinWidth(0);

  for (int i = 0; i < NSETS; ++i)
  {
    int istart = hWave[i]->FindBin(sStart);
    int iend = hWave[i]->FindBin(sEnd);
    int jstart = hWave[i]->FindBin(tStart);
    int jend = hWave[i]->FindBin(tTripletCut);
    int kend = hWave[i]->FindBin(tEnd);
    if (i == 0)
      printf(" istart %i %f iend %i %f jstart %i %f jend %i %f kend %i %f \n",
             istart, hWave[i]->GetBinLowEdge(istart),
             iend, hWave[i]->GetBinLowEdge(iend),
             jstart, hWave[i]->GetBinLowEdge(jstart),
             jend, hWave[i]->GetBinLowEdge(jend),
             kend, hWave[i]->GetBinLowEdge(kend));

    double singInt = hWave[i]->Integral(istart, iend);
    double tripInt = hWave[i]->Integral(jstart, jend);
    double tot = hWave[i]->Integral(istart, kend);
    double totTriplet = hWave[i]->Integral(iend, hWave[i]->GetNbinsX());
    singletIntegral[i] = singInt;
    tripletIntegral[i] = tripInt;
    totalIntegral[i] = tot;
    tripletTotalIntegral[i] = totTriplet;
    totalNorm[i] = getNorm(i);
  }


  for(int i=0; i < NSETS; ++i) printf("%i singlet %f triplet %f total %f triplet %f total norm %f \n",
      i, singletIntegral[i],  tripletIntegral[i], totalIntegral[i],tripletTotalIntegral[i], totalNorm[i]   );




  TMultiGraph *mg = new TMultiGraph();

  TGraph *gSinglet = new TGraph(NSETS, PPM, singletIntegral);
  gSinglet->SetName("RangeSingletModel");
  gSinglet->SetTitle(Form("model singlet range integral %0.2f - %0.2f ", sStart, sEnd));
  gSinglet->SetMarkerColor(kRed);
  gSinglet->SetMarkerStyle(30);
  gSinglet->SetMarkerSize(1.4);

  TGraph *gTriplet = new TGraph(NSETS, PPM, tripletIntegral);
  gTriplet->SetName("RangeTripletModel");
  gTriplet->SetTitle(Form("model triplet range integral %0.2f - %0.2f ", sEnd, tEnd));
  gTriplet->SetMarkerColor(kGreen);
  gTriplet->SetMarkerStyle(25);
  gTriplet->SetMarkerSize(1.4);

  TGraph *gTotTriplet = new TGraph(NSETS, PPM, tripletTotalIntegral);
  gTotTriplet->SetName("TotalTripletModel");
  gTotTriplet->SetTitle(Form("model total triplet  %0.2f - %0.2f ", sEnd, hWave[0]->GetBinLowEdge(nDigi)));
  gTotTriplet->SetMarkerColor(kGreen + 2);
  gTotTriplet->SetMarkerStyle(27);
  gTotTriplet->SetMarkerSize(1.4);

  TGraph *gTot = new TGraph(NSETS, PPM, totalIntegral);
  gTot->SetName("RangeTotalModel");
  gTot->SetTitle(Form("model total range integral %0.2f - %0.2f ", sStart, tEnd));
  gTot->SetMarkerColor(kBlue);
  gTot->SetMarkerStyle(26);
  gTot->SetMarkerSize(1.4);

  TGraph *gTotNorm = new TGraph(NSETS, PPM, totalNorm);
  gTotNorm->SetName("TotalModel");
  gTotNorm->SetTitle("model total norm");
  gTotNorm->SetMarkerColor(kBlack);
  gTotNorm->SetMarkerStyle(24);
  gTotNorm->SetMarkerSize(1.4);

  mg->Add(gSinglet);
  mg->Add(gTotTriplet);
  mg->Add(gTriplet);
  mg->Add(gTot);
  mg->Add(gTotNorm);

  TCanvas *cLight = new TCanvas(Form("IntegralsBySet-%.2f", tEnd), "IntegralsBySet");
  cLight->SetGridx();
  cLight->SetGridy();
  TString gTitle;
  gTitle.Form("Model Light Yield  ; Xe PPM dopant ; Yield ");
  mg->SetTitle(gTitle);
  mg->Draw("ap");
  //gl->Draw("Csame");
  cLight->BuildLegend();
  cLight->Print(".png");

  printf(" results for tEnd = %.3f \n", tEnd);
  for (int i = 0; i < NSETS; ++i)
    printf("%i %.3E %.3E %.3E norm %.3f triplet ratio %.3f\n",
           i, singletIntegral[i], tripletIntegral[i], totalIntegral[i], totalNorm[i], tripletIntegral[i] / tripletIntegral[0]);

  canTitle.Form("lightModel");
  can = new TCanvas(canTitle, canTitle);
  fp[NSETS - 1]->Draw();
  for (int iset = 0; iset < NSETS; ++iset)
    fp[iset]->Draw("same");
  for (int iset = 0; iset < NSETS; ++iset)
    fout->Append(fp[iset]);
  gPad->SetLogy();


  /*
  TString cname;
  for(int iset=0; iset< NSETS-1; ++iset) {
    cname.Form("WaveGraphSet-%i",iset);
    can = new TCanvas(cname,cname);
    hWave[iset]->Draw("");
    fp[iset]->Draw("same");
  }
  */


  canTitle;
  canTitle.Form("WaveModel");
  TCanvas *canw = new TCanvas(canTitle, canTitle);
  // hWave[NSETS-1]->Draw();
  hWave[4]->GetXaxis()->SetRangeUser(980,6000);
  //hWave[4]->GetYaxis()->SetRangeUser(1E-2,1E2);
  hWave[4]->Draw("");
  for (int iset = 0; iset < NSETS; ++iset) {
    hWave[iset]->GetXaxis()->SetRangeUser(980,6000);
    //hWave[iset]->GetYaxis()->SetRangeUser(1E-2,1E2);
    hWave[iset]->Draw("same");
    //fp[iset]->Draw("same");
  }
  gPad->BuildLegend();
  gPad->SetLogy();

  // for(int iset=0; iset<NSETS; ++iset) {
  //   printModel(iset);
  // }


  fout->Append(gSinglet);
  fout->Append(gTriplet);
  fout->Append(gTotTriplet);
  fout->Append(gTot);
  fout->Append(gTotNorm);


  fout->Write();
}


