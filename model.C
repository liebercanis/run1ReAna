// time is in microseconds
#include "compiled/BaconAnalysis.hh"
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
double binwidth = 0.04;
enum {MAXSETS=5};
double PPM[MAXSETS]={0,1,2,5,10};
int pcolor[MAXSETS] = {kRed,kMagenta+2,kGreen, kYellow-2,kBlue};
enum {npars=6};
TFile *fout;
double par[6][6];




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

  TF1* fp  = new TF1(Form("LfFit-%i-%0.f-PPM", iset, ppm),lightModel,xstart,xstop,npars);
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
  par[0][0] = 0.005000 ;
  par[0][1] = 1010853.295568 ;
  par[0][2] = 0.705848 ;
  par[0][3] = 1.911888 ;
  par[0][4] = 0.130000 ;
  par[0][5] = 3.000000 ;
  par[1][0] = 0.005000 ;
  par[1][1] = 2556929.017753 ;
  par[1][2] = 3.910191 ;
  par[1][3] = 3.364905 ;
  par[1][4] = 0.130000 ;
  par[1][5] = 3.000000 ;
  par[2][0] = 0.005000 ;
  par[2][1] = 3075164.849148 ;
  par[2][2] = 4.902713 ;
  par[2][3] = 3.876612 ;
  par[2][4] = 0.130000 ;
  par[2][5] = 3.000000 ;
  par[3][0] = 0.005000 ;
  par[3][1] = 3229307.580083 ;
  par[3][2] = 7.810860 ;
  par[3][3] = 4.264110 ;
  par[3][4] = 0.130000 ;
  par[3][5] = 3.000000 ;
  par[4][0] = 0.005000 ;
  par[4][1] = 3673658.847956 ;
  par[4][2] = 12.620146 ;
  par[4][3] = 8.837865 ;
  par[4][4] = 0.130000 ;
  par[4][5] = 3.000000 ;
  par[5][0] = 0.005000 ;
  par[5][1] = 4764881.683593 ;
  par[5][2] = 12.600111 ;
  par[5][3] = 8.188048 ;
  par[5][4] = 0.130000 ;
  par[5][5] = 3.000000 ;

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
