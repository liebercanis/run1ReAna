static double tres = 7.86; //ns

static double expG(double* xx, double* par)
{
  double tau = par[0];
  double x = xx[0];
  double f =  TMath::Exp(-x/tau)/tau;
  return f;
}

void drawExp()
{
  double tau = 1600.;
  TF1* fp = new TF1(Form("Exp-%.2f",tau),expG,0.,10*tau,1);
  fp->SetParName(0,"tau");
  fp->SetParameter(0,tau);
  fp->Print();

  TString canTitle;
  TCanvas *can1; 
  canTitle.Form("Exp%.0f",tau);
  can1 = new TCanvas(canTitle, canTitle);
  fp->Draw("");
  gPad->SetLogy();
  //By default the legend is automatically placed with width = x1= x2 = 0.3 and height = y1= y2 = 0.21. 
  //x1,y1,x2,y2	The TLegend coordinates 
  can1->BuildLegend(.3,.2,.3,.2,canTitle);
  can1->Modified(); can1->Update();
  can1->Print(".pdf");

  printf(" integral %f \n",fp->Integral(0,10.*tau) );

}
