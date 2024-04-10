static double tres = 5.4; //ns

static double expGaus(double* xx, double* par)
{
  double tau = par[0];
  double x = xx[0];
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2)/tau;
  return f;
}

void drawExpGaus()
{
  double tau = 1600.;
  TF1* fp = new TF1(Form("ExpGaus-%.2f-%.2f",tau,tres),expGaus,0.,10*tau,1);
  fp->SetParName(0,"tau");
  fp->SetParameter(0,tau);
  fp->Print();

  TString canTitle;
  TCanvas *can1; 
  canTitle.Form("ExpGaus%.2f",tau);
  can1 = new TCanvas(canTitle, canTitle);
  fp->Draw("");
  gPad->SetLogy();
  //By default the legend is automatically placed with width = x1= x2 = 0.3 and height = y1= y2 = 0.21. 
  //x1,y1,x2,y2	The TLegend coordinates 
  can1->BuildLegend(.3,.2,.3,.2,canTitle);
  can1->Modified(); can1->Update();
  can1->Print(".pdf");

  printf(" integral %f %f \n",fp->Integral(0,10.*tau),fp->Integral(-10.*tres,10.*tau) );

}
