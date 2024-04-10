// time is in microseconds
// model with absorption
#include "modelFit.hh"
using namespace TMath;

modelFit* model[NTYPES][NSETS];
TF1* fit[NTYPES][NSETS];
TH1D *hWave[NSETS];
TH1D *hWaveInt[NSETS];


TFile *fout;

double getNorm(int iset) {
  return 1.;
}


static double absFunc(double *xx, double *par) {
  double c = xx[0];
  double a0 = 0.62;
  double c0 = par[0];
  if(c>0.1)  return a0;
  double a = 1.-1.*Exp(-c/c0);
  return a;
};

void makeWaves() {
  for (int iset = 0; iset < NSETS; ++iset) {
    hWave[iset] = (TH1D*) fit[4][iset]->GetHistogram();
    hWave[iset]->SetName(Form("ModelWaveSet-%i",iset));
    hWave[iset]->SetLineColor(model[4][iset]->setColor[iset]);
    fout->Append(hWave[iset]);
  }

  // integral plots
  int startBin = hWave[0]->FindBin(900.);
  for (int i = 0; i < NSETS; ++i)
  {
    double runPPM = fit[4][i]->GetParameter(2);
    double runKp = fit[4][i]->GetParameter(4);
    double runKprime = fit[4][i]->GetParameter(7);
    double yieldSum=0;
    double yieldSumErr=0;
    TString htitle; htitle.Form("IntegralToFitSet%i-PPM%.2f-kP%.2f-kPrime-%.2f", i,runPPM,runKp,runKprime);
    hWaveInt[i] = (TH1D *)hWave[i]->Clone(htitle);
    hWaveInt[i]->Reset("ICESM");
    hWaveInt[i]->SetTitle(Form("integral set %i PPM %i ", i, int(runPPM)));
    printf(" set %i %s bins %i \n",i,hWaveInt[i]->GetName(),hWaveInt[i]->GetNbinsX());
    for (int ibin = 0; ibin < hWaveInt[i]->GetNbinsX(); ++ibin)
    {
      hWaveInt[i]->SetBinContent(ibin,0);
      hWaveInt[i]->SetBinError(ibin, 0);
      if(ibin<startBin) continue;
      double val = hWave[i]->GetBinContent(ibin);
      yieldSum += val;
      yieldSumErr = hWave[i]->GetBinError(ibin);
      hWaveInt[i]->SetBinContent(ibin, yieldSum);
      hWaveInt[i]->SetBinError(ibin, yieldSumErr);
    }
  }


}

void plotSet(int iset) {


  TString canTitle;
  TCanvas *can; 
  canTitle.Form("ModelSet-%if",iset);
  can = new TCanvas(canTitle, canTitle);
  can->SetLogy(1);

  for (int i = 0; i < NTYPES; ++i){
    fit[i][iset]->GetYaxis()->SetRangeUser(1.E-3,1E3);
    fit[i][iset]->SetLineColor(model[i][iset]->setColor[i]);
    if(i==0) fit[i][iset]->Draw("");
    else fit[i][iset]->Draw("same");
  }
  gPad->SetLogy();
  can->BuildLegend(.3,.2,.3,.2,canTitle);
  can->Modified(); can->Update();
  can->Print(".pdf");

}



void runModel(double absorption = 0.62, double kprime =1.0,  double kp=1.0)
{

  /*
    param 0 binw 4.000E-02 +/- 0.000E+00 
	  param 1 norm 1.000E+04 +/- 0.000E+00 
	  param 2 PPM 1.000E+03 +/- 0.000E+00 
	  param 3 tau3 1.600E+03 +/- 0.000E+00 
	  param 4 kp 1.000E+00 +/- 0.000E+00 
	  param 5 sfrac 1.400E-01 +/- 0.000E+00 
	  param 6 ab 0.000E+00 +/- 0.000E+00 
	  param 7 kprime 0.000E+00 +/- 0.000E+00 
	  param 8 type 4.000E+00 +/- 0.000E+00 
  ** 
  */

  TF1 *fabs = new TF1("fabs",absFunc,0,1,1);   // create TF1 class.
  fabs->SetParName(0,"C0");
  fabs->SetParameter(0,.103);
  fabs->SetNpx(1000000); 
  TCanvas* canabs = new TCanvas("absFitFunc","absFitFunc");
  fabs->Draw("");
  canabs->Print(".pdf");
   
  TString fileTitle;
  fileTitle.Form("integral-light-model-kp-%.2f-kPrime-%.2f-Abs-%.2f.root",kp,kprime,absorption);
  cout <<  "\t ....  " << fileTitle << endl;
  fout = new TFile(fileTitle, "RECREATE");
  // if(iset>=NSETS) return
  int fitColor[NTYPES] = {kYellow, kBlue, kRed, kGreen, kBlack};
  TString typeName[NTYPES];
  typeName[0]=TString("singlet");
  typeName[1]=TString("triplet");
  typeName[2]=TString("mixed");
  typeName[3]=TString("xenon");
  typeName[4]=TString("total");
  for (int ifit = 0; ifit < NTYPES; ++ifit){
    for (int iset = 0; iset < NSETS; ++iset){
      printf(" fit %i set %i \n",ifit,iset);
      model[ifit][iset]= new modelFit(ifit,iset);
      fit[ifit][iset] = model[ifit][iset]->fp;
      double sppm =  model[ifit][iset]->XPPM[iset];
      //double sppm =  model[ifit][iset]->PPM[iset];
      double afactor = absorption;
      if(sppm<1&&absorption>0.001)  afactor = fabs->Eval( sppm );
      fit[ifit][iset]->SetParameter(2, sppm);
      fit[ifit][iset]->SetParameter(4,kp);
      //absorption = fabs->Eval(sppm);
      fit[ifit][iset]->SetParameter(6,afactor);
      fit[ifit][iset]->SetParameter(7,kprime);
      //fit[ifit][iset]->SetParameter(9,10.E3);

      fit[ifit][iset]->SetTitle(Form("%s-PPM-%0.1f",typeName[ifit].Data() , sppm ));
      //model[ifit][iset]->show();
      printf(" fit %i set %i ppm %.2f abs %.4f \n",ifit,iset,sppm, afactor);
      fout->Append(fit[ifit][iset]);
    }
  }
  model[4][NSETS-1]->show();

  double xval[NTYPES][NSETS];
  double yval[NTYPES][NSETS];


  for (int iset = 0; iset < NSETS; ++iset) plotSet(iset);

  double sStart = 900.;
  for (int ifit = 0; ifit < NTYPES; ++ifit){
    for (int iset = 0; iset < NSETS; ++iset) {
      double binw = fit[ifit][iset]->GetParameter(0);
      //xval[ifit][iset] =fit[ifit][iset]->GetParameter(2);
      xval[ifit][iset] = fit[ifit][iset]->GetParameter(2);
      yval[ifit][iset]= fit[ifit][iset]->Integral(sStart,xMax)/binw;

      printf(" total light integral set %i fit %i norm %f int %f \n",iset,ifit,xval[ifit][iset],yval[ifit][iset]);
    }
  }

  TGraph *gInt[NTYPES];

  for (int ifit = 0; ifit < NTYPES; ++ifit) 
    gInt[ifit] = new TGraph(NSETS,&xval[ifit][0],&yval[ifit][0]);

  TString canTitle;
  double runKp = fit[4][NSETS-1]->GetParameter(4);
  double ab = fit[4][NSETS-1]->GetParameter(6);
  double runKprime = fit[4][NSETS-1]->GetParameter(7);



  canTitle.Form("integral-light-model-kp-%.2f-kPrime-%.2f-abs-%.2f",runKp,runKprime,ab);
  TCanvas *cani = new TCanvas(canTitle, canTitle);

  TMultiGraph *mgInt = new TMultiGraph();
  for (int ifit = 0; ifit < NTYPES; ++ifit){ 
    gInt[ifit]->SetTitle(canTitle);
    gInt[ifit]->SetMarkerStyle(20+ifit);
    if(ifit==4)  gInt[ifit]->SetMarkerStyle(29);
    gInt[ifit]->SetMarkerSize(1.5);
    gInt[ifit]->SetMarkerColor(fitColor[ifit]);
    gInt[ifit]->SetTitle(typeName[ifit]);
    gInt[ifit]->GetYaxis()->SetTitle(" total light yield ");
    gInt[ifit]->GetXaxis()->SetTitle(" dopant ");
    mgInt->Add(gInt[ifit]);
  }
  mgInt->SetTitle("yield; Xe dopant PPM ; yield ");
  mgInt->Draw("ap");
  gPad->SetLogy(0);
  gPad->SetLogx(1);
  gPad->SetGridx();
  gPad->SetGridy();
  cani->BuildLegend();
  cani->Print(".pdf");
  //cani->Update();

  printf(" set %i int @ %.2f PPM %.0f \n",0,xval[4][0], yval[4][0] );
  printf(" set %i int @ %.2f PPM %.0f \n",4,xval[4][4], yval[4][4] );


  TCanvas *can1; 
  canTitle.Form("Mixed-ab-%.2f",absorption);
  can1 = new TCanvas(canTitle, canTitle);

  for (int iset = 0; iset < NSETS; ++iset){
    fit[2][iset]->GetYaxis()->SetRangeUser(1.E-2,1E3);
    if(iset==0) fit[2][iset]->Draw("");
    else fit[2][iset]->Draw("same");
  }
  gPad->SetLogy();
  //By default the legend is automatically placed with width = x1= x2 = 0.3 and height = y1= y2 = 0.21. 
  //x1,y1,x2,y2	The TLegend coordinates 
  can1->BuildLegend(.3,.2,.3,.2,canTitle);
  can1->Modified(); can1->Update();
  can1->Print(".pdf");

  
 
  makeWaves();

  canTitle.Form("IntegralLight-ab-%.2f-kP-%.2f-kPrime-%.2f", ab,runKp,runKprime);
  TCanvas *can3; 
  can3 = new TCanvas(canTitle, canTitle);
  gStyle->SetOptTitle(0); //this will disable the title for all coming histograms or

  for (int iset = 0; iset < NSETS; ++iset){
    hWaveInt[iset]->GetYaxis()->SetRangeUser(1.E-2,1.1E4);
    if(iset==0) hWaveInt[iset]->Draw("");
    else hWaveInt[iset]->Draw("same");
  }
  gPad->SetLogy();
  can3->BuildLegend(.3,.2,.3,.2,canTitle);
  can3->Modified(); can3->Update();
  can3->Print(".pdf");

  canTitle.Form("lightModel-ab-%.2f",absorption);
  TCanvas *can2; 
  can2 = new TCanvas(canTitle, canTitle);
  for (int iset = 0; iset < NSETS; ++iset){
    fit[4][iset]->GetYaxis()->SetRangeUser(1.E-2,1E2);
    fit[4][iset]->SetLineColor(model[4][iset]->setColor[iset]);
    //fit[3][iset]->SetLineColor(model[3][iset]->setColor[iset]);

    if(iset==0) fit[4][iset]->Draw("");
    else fit[4][iset]->Draw("same");
    //fit[3][iset]->Draw("same");

  }
  gPad->SetLogy();
  can2->BuildLegend(1,1,1,1,canTitle);
  can2->Modified(); can2->Update();
  can2->Print(".pdf");



  return;
  fout->Write();
}

