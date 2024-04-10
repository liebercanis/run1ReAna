// time is in microseconds
#include "compiled/BaconAnalysis.hh"
using namespace TMath;

double aveEvents = 2607.841463;
double binwidth = 0.04;
enum
{
  NSETS = 5
};
double PPM[NSETS] = {0, 1, 2, 5, 10};
int setColor[NSETS];
enum
{
  npars = 6
};
TFile *fout;

enum {MAXSETS=5};
double PPM[MAXSETS]={0,1,2,5,10};
int pcolor[MAXSETS] = {kRed,kMagenta+2,kGreen, kYellow-2,kBlue};
TFile *fout;
double par[NSETS][npars+1];



void printModel(int iset)
{ 
  double bw = par[iset][0];
  double norm = par[iset][1];
  double kx = k1Zero*par[iset][2];
  double tau3 = par[iset][3];
  double kp = par[iset][4];
  int type = int(par[iset][5]);
  double td = 1/kx;
  double tq = 1./(1./tau3+kp);
  double tr = 1./(kx + 1./tq);
  double alpha1 = bw*singletFraction*norm;
  double alpha3 = bw*(1.-singletFraction)*norm;

  printf("set %i tau3 %f alpha1 %E alpha3 %E ratio %E \n",iset,tau3,alpha1,alpha3,alpha3/alpha1);
}



static Double_t lightModel(Double_t *xx, Double_t *para)
{ 
  double x=xx[0]-xstart;
  double bw = para[0];
  double norm = para[1];
  double kx = k1Zero*para[2];
  double tau3 = para[3];
  double kp = para[4];

  int type = int(para[5]);
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


TF1* setFit(int iset, double norm=1)
{

  int type = 3;
  double ppm=PPM[iset]; 
  double kplus = kplusZero;

  TF1* fp  = new TF1(Form("LfFit-%i-%0.f-PPM", iset, ppm),lightModel, sStart, tEnd, npars);
  fp->SetParName(0,"binw");
  fp->SetParName(1,"norm");
  fp->SetParName(2,"PPM");
  fp->SetParName(3,"tau3");
  fp->SetParName(4,"k+");
  fp->SetParName(5,"type");
  fp->FixParameter(0, par[iset][0]);
  fp->SetParameter(1, par[iset][1]);
  fp->SetParameter(2, par[iset][2]);
  fp->SetParameter(3, par[iset][3]);
  fp->FixParameter(4, par[iset][4]);
  fp->FixParameter(5, par[iset][5]);
  fp->SetTitle(Form("LfFit-%i-%0.f-PPM", iset, ppm));
  fp->SetNpx(1000); // numb points for function
  fp->SetLineColor(pcolor[iset]);

  printf(" \n\n >>> setFit starting parameters \n");
  for (int ii = 0; ii < npars; ++ii)
  {
    printf("\t  param %i %s %.3f +/- %.3f \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }
  fout->Append(fp);
  return fp;
}



void model()
{
  //........ fitted parameters set 0
  par[0][0] = 1.000000 ;
  par[0][1] = 106755.581074 ;
  par[0][2] = -0.511235 ;
  par[0][3] = 2100.000000 ;
  par[0][4] = 0.000130 ;
  par[0][5] = 0.140000 ;
  par[0][6] = 3.000000 ;
  //........ fitted parameters set 1
  par[1][0] = 1.000000 ;
  par[1][1] = 219597.924134 ;
  par[1][2] = 3.342483 ;
  par[1][3] = 2100.000000 ;
  par[1][4] = 0.000130 ;
  par[1][5] = 0.140000 ;
  par[1][6] = 3.000000 ;
  //........ fitted parameters set 2
  par[2][0] = 1.000000 ;
  par[2][1] = 244210.401667 ;
  par[2][2] = 4.055385 ;
  par[2][3] = 2100.000000 ;
  par[2][4] = 0.000130 ;
  par[2][5] = 0.140000 ;
  par[2][6] = 3.000000 ;
  //........ fitted parameters set 3
  par[3][0] = 1.000000 ;
  par[3][1] = 253277.558390 ;
  par[3][2] = 6.250083 ;
  par[3][3] = 2100.000000 ;
  par[3][4] = 0.000130 ;
  par[3][5] = 0.140000 ;
  par[3][6] = 3.000000 ;
  //........ fitted parameters set 4
  par[4][0] = 1.000000 ;
  par[4][1] = 283584.755090 ;
  par[4][2] = 10.061776 ;
  par[4][3] = 2100.000000 ;
  par[4][4] = 0.000130 ;
  par[4][5] = 0.140000 ;
  par[4][6] = 3.000000 ;


  fout = new TFile("model.root","RECREATE");
  //if(iset>=MAXSETS) return;
  TF1* fp[MAXSETS];
  for(int iset=0; iset<MAXSETS; ++iset) {
    fp[iset] = setFit(iset);
  }

  TString canTitle;
  canTitle.Form("lightModel");
  TCanvas *can = new TCanvas(canTitle,canTitle);
  fp[MAXSETS-1]->Draw();
  for(int iset=0; iset<MAXSETS; ++iset) fp[iset]->Draw("same");
  for(int iset=0; iset<MAXSETS; ++iset) fout->Append(fp[iset]);
  gPad->SetLogy();

  fout->Write();

  for(int iset=0; iset<MAXSETS; ++iset) {
    printModel(iset);
  }

}
