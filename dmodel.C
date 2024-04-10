// time is in microseconds
// model with absorption
#include "compiled/BaconAnalysis.hh"

using namespace TMath;

double ab = 0.5;
double aveEvents = 2607.841463;
double binwidth = 0.04;
enum
{
  NSETS = 5
};
double PPM[NSETS];
double TRIPLET= 1600;
double NORM[NSETS];
double BASE[NSETS];

int setColor[NSETS];
int setStyle[NSETS];
enum
{
  NPARS = 8 
};

TFile *fout;

double getNorm(int iset) {
  return 1.;
}


static double expGaus(double x, double tau)
{
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * Exp(arg1) * Erfc(arg2);
  return f;
}

static Double_t lightModel(Double_t *xx, Double_t *par)
{
  for(int i=0; i<NPARS; ++i ) printf(" par %i %f \n",i,par[i]);
  double x = xx[0] - xTrigger;
  double bw =  par[0];
  double norm = par[1];
  double ppm=  par[2];
  double tTrip = par[3];
  double kp = par[4]*kplusZero;
  double alpha1 = par[5]*bw*norm;
  double alpha3 = (1.-par[5])*bw*norm;
  double ab= par[6];
  int type = int(par[7]);
  //
  double kx = k1Zero*ppm;
  double lS =  1./tSinglet;
  double lT = 1/ tTrip;
  double l1 = 1./tSinglet + kp + kx;
  double l3 = 1./tTrip  + kp + kx;
  double lX = 1./tXe;
  double dX = lX - kx;
  double dS = lX - lS;
  double dT = lX - lT;

  printf("ppm %.1f kx %.3E kp %.3E lS %.3E lT %.3E lX %.3E\n",ppm,kx,kp,lS,lT,lX);

  double c1 = kx + ab*lS;
  double c3 = kx + ab*lT;
  double tKx = 1 / kx;
  double t1 = 1./l1;
  double t3 = 1./l3;

  double f=0;
  double fs = (1.-ab)*alpha1 / tSinglet * expGaus(x, t1);
  double ft = (1.-ab)*alpha3/tTrip * expGaus(x, t3);
  double fm = alpha1*c1/(l1-kx)*(  expGaus(x, tKx)  - expGaus(x, t1) ) + alpha3*c3/(l3-kx)*(  expGaus(x, tKx)  - expGaus(x, t3) );
  double x1 = c1*kx*alpha1/(l1-kx) * (  (expGaus(x, tKx)    -  expGaus(x, tXe))/(lX - kx)  - (expGaus(x, t1)     -  expGaus(x, tXe))/(lX - l1) );
  double x3 = c3*kx*alpha3/(l3-kx) * (  (expGaus(x, tKx) -  expGaus(x, tXe))/(lX - kx)  -  (expGaus(x, t3)  -  expGaus(x, tXe))/(lX - l3));

  printf(" %E %E %E %E %E  \n"  ,fs,ft,fm,x1,x3);

  double fx = x1 + x3;

  if (type == 0)
    f = fs;
  else if (type == 1)
    f = ft;
  else if (type == 2)
    f = fm;
  else if (type == 3)
    f = fx;
  else if (type == 4)
    f = fs+ft+fx;

  return f;


}

TF1 *setFit(int iset, int ifit)
{
  double ppm=PPM[iset];
  double integral = NORM[iset];
  double kplus = 1;

  double kx = k1Zero*ppm;
  double kp = ppm*kplusZero;
  double lS =  1./tSinglet;
  double lT = 1/  tTriplet;
  double l1 = 1./tSinglet + kp + kx;
  double l3 = 1./ tTriplet;  + kp + kx;
  double lX = 1./tXe;
  double dX = lX - kx;
  double d1 = lX - l1;
  double d3 = lX - l3;

  printf(" set %i ppm %.1f kx %.3E kp %.3E lS %.3E lT %.3E lX %.3E\n",iset,ppm,kx,kp,lS,lT,lX);
  printf(" dX %E d1 %E d3 %E  \n",dX,d1,d3);


  TF1 *fp = new TF1(Form("LifeFit-type-%i-set-%i",ifit,iset), lightModel, xTrigger-20., 3000, NPARS);
  printf(" set %i fit range %f to %f k+ %E \n", iset,xTrigger-20.,tailStart,kplus);
  fp->SetParName(0,"binw");
  fp->SetParName(1,"norm");
  fp->SetParName(2,"PPM");
  fp->SetParName(3,"tau3");
  fp->SetParName(4,"kp");
  fp->SetParName(5,"sfrac");
  fp->SetParName(6,"ab");
  fp->SetParName(7,"type");

  fp->FixParameter(0, binwidth);
  fp->SetParameter(1, integral);
  fp->FixParameter(2, ppm);
  fp->FixParameter(3, tTriplet);
  fp->FixParameter(4,kplus);
  fp->FixParameter(5, singletFraction);
  fp->SetParameter(6,ab);
  fp->SetParLimits(6,1.E-4,100.);
  fp->FixParameter(7,ifit);
  fp->SetTitle(Form("LfFit-type-%i-set-%i-%0.1f-PPM-ab-%.2f", ifit, iset, ppm,ab));
  fp->SetNpx(1000); // numb points for function
  fp->Print();
  fp->SetLineColor(setColor[iset]);
  //fp->SetMarkerStyle(setStyle[iset]);

  printf(" \n\n >>> setFitD fitted parameters set %i \n",iset);
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }

  cout << "defined function " << fp->GetName() << "  " <<  fp->GetTitle() << endl;
  return fp;
}

void dmodel(double absorption = 0.5)
{
  ab = absorption;

  printf(" absorption = %E \n",ab);

  PPM[0]=0.1;
  PPM[1]=1;
  PPM[2]=2;
  PPM[3]=5;
  PPM[4]=10;

  for(int iset=0; iset<NSETS; ++iset) {
    NORM[iset]= aveEvents;
  }


  setColor[0] = kOrange + 4;
  setColor[1] = kBlue;
  setColor[2] = kRed;
  setColor[3] = kGreen - 2;
  setColor[4] = kCyan - 2;



  fout = new TFile(Form("gmodel-ab-%0.3f.root",ab), "RECREATE");
  // if(iset>=NSETS) return;
  TF1 *fp0[NSETS];
  TF1 *fp1[NSETS];
  TF1 *fp2[NSETS];
  TF1 *fp3[NSETS];
  TF1 *fp4[NSETS];


  for (int iset = 0; iset < NSETS; ++iset)
  {
    printf(" set %i \n",iset);
    fp0[iset] = setFit(iset,0);
    fout->Append(fp0[iset]);
    fp1[iset] = setFit(iset,1);
    fout->Append(fp1[iset]);
    fp2[iset] = setFit(iset,2);
    fout->Append(fp2[iset]);
    fp3[iset] = setFit(iset,3);
    fout->Append(fp3[iset]);
    fp4[iset] = setFit(iset,4);
    fout->Append(fp4[iset]);
  }
 

  fout->Write();

  TString canTitle;
  canTitle.Form("lightModel-ab-%.2f",ab);
  TCanvas* can = new TCanvas(canTitle, canTitle);

  fp4[NSETS - 1]->GetYaxis()->SetRangeUser(1.E-3,1E1);
  fp4[NSETS - 1]->Draw();
  for (int iset = 0; iset < NSETS; ++iset){
    fp4[iset]->GetYaxis()->SetRangeUser(1.E-3,1E1);
    fp4[iset]->Draw("same");
  }
  gPad->SetLogy();
  can->BuildLegend(1,1,1,1,canTitle);
  can->Modified(); can->Update();
  can->Print(".pdf");

  
  canTitle.Form("Mixed-ab-%.2f",ab);
  can = new TCanvas(canTitle, canTitle);

  fp2[NSETS - 1]->GetYaxis()->SetRangeUser(1.E-3,1E2);
  fp2[NSETS - 1]->Draw();
  for (int iset = 0; iset < NSETS; ++iset){
    fp2[iset]->GetYaxis()->SetRangeUser(1.E-3,1E2);
    fp2[iset]->Draw("same");
  }
  gPad->SetLogy();
  can->BuildLegend(1,1,1,1,canTitle);
  can->Modified(); can->Update();
  can->Print(".pdf");
  return;




  // histograms
  TH1D *hWave[NSETS];

  canTitle.Form("lightModel");
  can = new TCanvas(canTitle, canTitle);
  fp4[NSETS - 1]->GetHistogram()->SetXTitle("time nanoseconds");
  fp4[NSETS - 1]->Draw();
  for (int iset = 0; iset < NSETS; ++iset)
    fp4[iset]->Draw("same");

  can->BuildLegend();
  gPad->SetLogy();

  for (int iset = 0; iset < NSETS; ++iset) {
    hWave[iset] = (TH1D*) fp4[iset]->GetHistogram();
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

  

  canTitle;
  canTitle.Form("WaveModel-ab-%.2f",ab);
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


