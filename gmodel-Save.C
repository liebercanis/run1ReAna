// time is in microseconds
#include "compiled/BaconAnalysis.hh"
using namespace TMath;

double aveEvents = 2607.841463;
double binwidth = 0.04;
enum
{
  NSETS = 5
};
double PPM[NSETS];
double TRIPLET[NSETS];
double NORM[NSETS];
double BASE[NSETS];

int setColor[NSETS];
enum
{
  NPARS = 7 
};

double fitpar[NSETS][NPARS];
TFile *fout;


void getParam(double par[NSETS][NPARS]){

//........ fitted parameters set 0
 par[0][0] = 1.000000 ;
 par[0][1] = 119980.397087 ;
 par[0][2] = -0.451276 ;
 //par[0][2] = 0. ;
 par[0][3] = 2538.172504 ;
 par[0][4] = 0.000130 ;
 par[0][5] = 0.140000 ;
 par[0][6] = 3.000000 ;
//........ fitted parameters set 1
 par[1][0] = 1.000000 ;
 par[1][1] = 225748.883783 ;
 par[1][2] = 3.804387 ;
 par[1][3] = 4072.803926 ;
 par[1][4] = 0.000130 ;
 par[1][5] = 0.140000 ;
 par[1][6] = 3.000000 ;
//........ fitted parameters set 2
 par[2][0] = 1.000000 ;
 par[2][1] = 242531.322443 ;
 par[2][2] = 4.841292 ;
 par[2][3] = 4489.702524 ;
 par[2][4] = 0.000130 ;
 par[2][5] = 0.140000 ;
 par[2][6] = 3.000000 ;
///........ fitted parameters set 3
 par[3][0] = 1.000000 ;
 par[3][1] = 250798.185081 ;
 par[3][2] = 7.142214 ;
 par[3][3] = 4448.561633 ;
 par[3][4] = 0.000130 ;
 par[3][5] = 0.140000 ;
 par[3][6] = 3.000000 ;
//........ fitted parameters set 4
 par[4][0] = 1.000000 ;
 par[4][1] = 282960.706190 ;
 par[4][2] = 10.923587 ;
 par[4][3] = 4748.391581 ;
 par[4][4] = 0.000130 ;
 par[4][5] = 0.140000 ;
 par[4][6] = 3.000000 ;
//........ fitted parameters set 5
 //par[5][0] = 1.000000 ;
 //par[5][1] = 288354.739170 ;
 //par[5][2] = 10.573530 ;
 //par[5][3] = 3798.823502 ;
 //par[5][4] = 0.000130 ;
 //par[5][5] = 0.140000 ;
 //par[5][6] = 3.000000 ;
 /* baseline from tbTFitAll.cc
    p0                        =      15.0874   +/-   0.00772584  
    p0                        =       12.393   +/-   0.0554649   
    p0                        =      12.7556   +/-   0.0181344   
    p0                        =      13.0048   +/-   0.0107551   
    p0                        =       13.259   +/-   0.0140966  
    */
  BASE[0]=15.0874;
  BASE[1]=12.393;
  BASE[2]=12.7556;
  BASE[3]=13.0048;
  BASE[4]=13.259;
}

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
  int type = 3;
  double td = 1 / kx;
  double tq = 1. / (1. / tau3 + kp);
  double tr = 1. / (kx + 1. / tq);
  double alpha1 = bw * singletFraction * norm;
  double alpha3 = bw * (1. - singletFraction) * norm;
  printf("set %i ppm %f tau3 %f alpha1 %E alpha3 %E ratio %E \n", iset, ppm, tau3, alpha1, alpha3, alpha3 / alpha1);
}

static Double_t lightModel(Double_t *xx, Double_t *para)
{
  double f = 0;
  double x = xx[0] - xTrigger;
  double bw = para[0];
  double norm = para[1];
  double kx = k1Zero * para[2];
  double tau3 = para[3];
  double kp = para[4];
  int type = int(para[5]);
  double td = 1 / kx;
  double tq = 1. / (1. / tau3 + kp);
  double tr = 1. / (kx + 1. / tq);
  double alpha1 = bw * singletFraction * norm;
  double alpha3 = bw * (1. - singletFraction) * norm;
  double fs = alpha1 / tSinglet * expGaus(x, tSinglet);
  double f3 = alpha3 / tau3 * expGaus(x, tr);
  double fx = alpha3 * pow(kx, 2) * tq * (expGaus(x, td) - expGaus(x, tr));
  if (type == 0)
    f = fs;
  else if (type == 1)
    f = f3;
  else if (type == 2)
    f = fx;
  else if (type == 3)
    f = fs + f3 + fx;
  return f;
}

TF1 *setFit(int iset)
{

  int type = 3;
  double ppm=PPM[iset];
  double tau3=TRIPLET[iset];
  double norm = NORM[iset];

  TF1 *fp = new TF1(Form("LifeFit-%i", iset), lightModel, sStart, tEnd, NPARS);
  fp->SetParName(0, "binw");
  fp->SetParName(1, "norm");
  fp->SetParName(2, "PPM");
  fp->SetParName(3, "tau3");
  fp->SetParName(4, "k+");
  fp->SetParName(5, "type");
  fp->FixParameter(0, binwidth);
  fp->SetParameter(1, norm);
  fp->SetParameter(2, ppm);
  fp->SetParameter(3, tau3);
  fp->FixParameter(4, kplusZero);
  fp->FixParameter(5, type);
  fp->SetTitle(Form("LifFit-%i-%0.2f-PPM", iset, PPM[iset]));
  fp->SetNpx(1000); // numb points for function
  fp->SetLineColor(setColor[iset]);

  printf(" \n\n >>> setFit starting parameters \n");
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.3f +/- %.3f \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }
  //fout->Append(fp);
  return fp;
}

void gmodel()
{

  getParam(fitpar);

  for(int iset=0; iset<NSETS; ++iset) {
    PPM[iset]=fitpar[iset][2];
    TRIPLET[iset]=fitpar[iset][3];
    //TRIPLET[iset]=2100.;
    NORM[iset]= fitpar[iset][1];
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
  for (int iset = 0; iset < NSETS; ++iset)
  {
    hWave[iset] = new TH1D(Form("ModelWaveSet-%i", iset), Form("model wave set %i ppm %.0f tau3 %.0f ", iset,PPM[iset],TRIPLET[iset]), nDigi, 0, nDigi);
    hWave[iset]->SetLineColor(setColor[iset]);
  }

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

  for (int iset = 0; iset < NSETS; ++iset){
    for (int ibin = 0; ibin < hWave[iset]->GetNbinsX(); ++ibin) {
      double x = hWave[iset]->GetBinLowEdge(ibin) + 0.5 * hWave[iset]->GetBinWidth(ibin);
      if (x < sStart)  hWave[iset]->SetBinContent(ibin, BASE[iset]);
      else hWave[iset]->SetBinContent(ibin, fp[iset]->Eval(x)+BASE[iset]);
    }
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
  for (int iset = 0; iset < NSETS; ++iset) {
    //hWave[iset]->GetXaxis()->SetRangeUser(980,10000);
    if(iset==0) hWave[iset]->Draw("");
    else hWave[iset]->Draw("same");
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


