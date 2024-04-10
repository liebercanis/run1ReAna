// time is in microseconds
using namespace TMath;
static double singletFraction = 0.14;
static double tTriplet=2.100; 
static double tSinglet=5.E-3;
static double tXe = 20.0E-3;
static double kxe = 8.8E-2;
static double kplusZero  = 1.3E-1;
static double xstart=1.;
static double xstop=10.;
static double tailStart=7.0;
static double tailStop=10.0;
static double pmodtG = 3.48;
double k1Zero = kxe*131./40.;
enum {MAXSETS=5};
double PPM[MAXSETS]={0,1,2,5,10};
int pcolor[MAXSETS] = {kRed,kMagenta+2,kGreen, kYellow-2,kBlue};
double tripletFit[MAXSETS];
double ppmFit[MAXSETS];
enum {npars=6};
TFile *fout;
double par[MAXSETS][npars];
TH1D* hWave[MAXSETS];
int nDigi = 10000;
double singletStart = 0.9;
double singletEnd = 1.02;
double tripletEnd = 2.5;

enum {NDATA = 6};
TString dataGraphName[NDATA];
TGraph *gData[NDATA];


static Double_t lightModel(Double_t *xx, Double_t *para)
{ 
  double x=xx[0]-xstart;
  if(x<0.) return 0.;
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


double getNorm(double ppm){
  double alpha3 = (1.-singletFraction);
  double kx = k1Zero*ppm;
  double tq = 1./(1./tTriplet+kplusZero);
  double tr = 1./(kx + 1./tq);
  double norm = singletFraction+alpha3*tr/tTriplet+alpha3*kx/(kx+1.0/tq);
  //printf(" \t ... ppm %.2f %f %f norm %f \n",ppm,alpha3*tr/tTriplet,alpha3*kx/(kx+1.0/tq) ,norm);
  return norm;
}


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



TF1* setFit(int iset, double norm=1)
{

  int type = 3;
  double ppm=PPM[iset]; 
  double kplus = kplusZero;

  TF1* fp  = new TF1(Form("LModel-%i-%0.f-PPM", iset, ppm),lightModel,0.0,xstop,npars);
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

  printf(" \n\n >>> setFit starting parameters set %i  \n",iset);
  for (int ii = 0; ii < npars; ++ii)
  {
    printf("\t  param %i %s %.3f +/- %.3f \n", ii, fp->GetParName(ii),fp->GetParameter(ii), fp->GetParError(ii));
  }
  fout->Append(fp);
  return fp;
}

//TripletFit	TripletFitMu Fittau3 Fittau3Mu FitPPM FitPPMMu

void getDataGraphs(TString dataFile)
{
  TFile *fdata = new TFile(dataFile);
  dataGraphName[0]=TString("TripletFit");
  dataGraphName[1]=TString("TripletFitMu");

  dataGraphName[2]=TString("Fittau3");
  dataGraphName[3]=TString("Fittau3Mu");

  dataGraphName[4]=TString("FitPPM");
  dataGraphName[5]=TString("FitPPMMu");


  for(int i=0; i<NDATA; ++i) {
    fdata->GetObject(dataGraphName[i], gData[i]);
    cout << "model graph "<< i << " " << gData[i]->GetName() << endl;
    fout->Append(gData[i]);
  }

  // get triplet 
  for(int i=0; i<MAXSETS; ++i) {
    double xval,yval;
    gData[2]->GetPoint(i,xval,yval);
    tripletFit[i]=yval;
    printf(" set %i PPM %f tau3 %f \n",i,xval,yval);
  }

  
  // get triplet 
  for(int i=0; i<MAXSETS; ++i) {
    double xval,yval;
    gData[4]->GetPoint(i,xval,yval);
    ppmFit[i]=yval;
    printf(" set %i PPM %f fit PPM %f \n",i,xval,yval);
  }

}



void model()
{
  TString modelFileName; modelFileName.Form("model-%.2f.root",tripletEnd);
  fout = new TFile(modelFileName,"RECREATE");
  printf(" the model \n");

  getDataGraphs(TString("tbSFitOut.root"));

   
  // this is theory 
  for(int iset=0; iset<MAXSETS; ++iset) {
    //PPM[iset]= PPM[iset] + 1.0;
    hWave[iset] = new TH1D(Form("WaveSet-%i",iset) ,Form("model wave set %i ",iset), nDigi,0,10.0); // time in microseconds
    hWave[iset]->SetLineColor(pcolor[iset]);
    par[iset][0]=hWave[iset]->GetBinWidth(0);  
    par[iset][1]=1.0;
    par[iset][2]=PPM[iset];
    //par[iset][3]=tripletFit[iset];
    par[iset][3]=tTriplet;
    par[iset][4]=kplusZero;
    par[iset][5]=3;
  }


  int max=1000;
  std::vector<double> ltot;
  std::vector<double> xppm;

  double maxppm=100.0;
  for(int i =0; i <=max; ++i) {
    double pi = double(i)/double(max)*maxppm;
    double ni = getNorm(pi);
    xppm.push_back(pi);
    ltot.push_back(ni);
    //printf("... %i %.0f %f ppm %f  norm %f %f \n",i,pi, double(i)/double(max),pi,ni,maxppm);
  }

  TGraph *gl = new TGraph(xppm.size(), &xppm[0],&ltot[0]);
  TCanvas *canl = new TCanvas("ltot","ltot");
  canl->SetGridx();
  canl->SetGridy();

  gl->SetTitle("total light");
  //gl->SetMarkerSize(0.007);
  gl->SetMarkerStyle(22);
  gl->SetLineWidth(3);
  gl->SetLineColor(kRed);
  gl->GetHistogram()->GetXaxis()->SetTitle("PPM");
  gl->GetHistogram()->GetYaxis()->SetTitle("triplet light");
  gl->Draw("ac");


  //if(iset>=MAXSETS) return;
  TF1* fp[MAXSETS];
  for(int iset=0; iset<MAXSETS; ++iset) {
    fp[iset] = setFit(iset);
    for(int ibin=0; ibin< hWave[iset]->GetNbinsX(); ++ibin ) {
      double x =  hWave[iset]->GetBinLowEdge(ibin)+0.5*hWave[iset]->GetBinWidth(ibin);
      if(x<xstart) continue;
      hWave[iset]->SetBinContent(ibin,fp[iset]->Eval(x));
     // if(iset==0) printf(" %i x %f val %f \n",ibin,x,fp[iset]->Eval(x));
    }
  }

  double singletIntegral[MAXSETS];
  double tripletIntegral[MAXSETS];
  double totalIntegral[MAXSETS];
  double tripletTotalIntegral[MAXSETS];
  double totalNorm[MAXSETS];

  double binWidth = hWave[0]->GetBinWidth(0);

  for (int i = 0; i < MAXSETS; ++i)
  {
    int istart = hWave[i]->FindBin(singletStart);
    int iend = hWave[i]->FindBin(singletEnd);
    int itrip =  hWave[i]->FindBin(tripletEnd);
    if(i==0) printf(" istart %i %f iend %i %f itrip %i %f \n",
          istart,hWave[i]->GetBinLowEdge(istart),
          iend,hWave[i]->GetBinLowEdge(iend),
          itrip,hWave[i]->GetBinLowEdge(itrip)
        );
    double singInt =  hWave[i]->Integral(istart, iend);
    double tripInt =  hWave[i]->Integral(iend,itrip);
    double tot =  hWave[i]->Integral(istart,itrip);
    double totTriplet =  hWave[i]->Integral(iend,hWave[i]->GetNbinsX());
    singletIntegral[i]=singInt;
    tripletIntegral[i]=tripInt;
    totalIntegral[i] = tot;
    tripletTotalIntegral[i] = totTriplet;
    totalNorm[i]=getNorm(PPM[i]);
  }

  /*
  double snorm = singletIntegral[0];
  double tnorm = tripletIntegral[0];
  for (int i = 0; i < MAXSETS; ++i) {
    printf("set %i ppm %f  %f trip %f ... ",i, PPM[i], singletIntegral[i], tripletIntegral[i]);
    singletIntegral[i]=singletIntegral[i]/snorm;
    tripletIntegral[i]=tripletIntegral[i]/tnorm;
    printf("set %i ppm %f  %f trip %f \n ",i, PPM[i], singletIntegral[i], tripletIntegral[i]);
  }
  */

    TMultiGraph* mg = new TMultiGraph();

    TGraph *gSinglet = new TGraph(MAXSETS, PPM, singletIntegral);
    gSinglet->SetName("RangeSingletModel");
    gSinglet->SetTitle(Form("model singlet range integral %0.2f - %0.2f ",singletStart,singletEnd));
    gSinglet->SetMarkerColor(kRed);
    gSinglet->SetMarkerStyle(30);
    gSinglet->SetMarkerSize(1.4);

    TGraph *gTriplet = new TGraph(MAXSETS, PPM, tripletIntegral);
    gTriplet->SetName("RangeTripletModel");
    gTriplet->SetTitle(Form("model triplet range integral %0.2f - %0.2f ",singletEnd,tripletEnd));
    gTriplet->SetMarkerColor(kGreen);
    gTriplet->SetMarkerStyle(25);
    gTriplet->SetMarkerSize(1.4);


    TGraph *gTotTriplet = new TGraph(MAXSETS, PPM, tripletTotalIntegral);
    gTotTriplet->SetName("TotalTripletModel");
    gTotTriplet->SetTitle(Form("model total triplet  %0.2f - %0.2f ",singletEnd,  hWave[0]->GetBinLowEdge(nDigi)));
    gTotTriplet->SetMarkerColor(kGreen+2);
    gTotTriplet->SetMarkerStyle(27);
    gTotTriplet->SetMarkerSize(1.4);


    TGraph *gTot = new TGraph(MAXSETS, PPM, totalIntegral);
    gTot->SetName("RangeTotalModel");
    gTot->SetTitle(Form("model total range integral %0.2f - %0.2f ",singletStart,tripletEnd));
    gTot->SetMarkerColor(kBlue);
    gTot->SetMarkerStyle(26);
    gTot->SetMarkerSize(1.4);

    
    TGraph *gTotNorm = new TGraph(MAXSETS, PPM, totalNorm);
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


    TCanvas *cLight  = new TCanvas(Form("IntegralsBySet-%.2f",tripletEnd), "IntegralsBySet");
    cLight->SetGridx(); cLight->SetGridy();
    TString gTitle; gTitle.Form("Model Light Yield  ; Xe PPM dopant ; Yield ");
    mg->SetTitle(gTitle);
    mg->Draw("ap");
    gl->Draw("Csame");
    cLight->BuildLegend();
    cLight->Print(".png");

    printf(" results for tripletEnd = %.3f \n", tripletEnd);
    for(int i=0; i<MAXSETS; ++i) printf("%i %.3f %.3f %.3f norm %.3f triplet ratio %.3f\n",
        i, singletIntegral[i],tripletIntegral[i],totalIntegral[i],totalNorm[i], tripletIntegral[i]/tripletIntegral[0] );




  TString canTitle;
  canTitle.Form("lightModel");
  TCanvas *can = new TCanvas(canTitle,canTitle);
  fp[MAXSETS-1]->Draw();
  for(int iset=0; iset<MAXSETS; ++iset) fp[iset]->Draw("same");
  for(int iset=0; iset<MAXSETS; ++iset) fout->Append(fp[iset]);
  gPad->SetLogy();

  
  canTitle;
  canTitle.Form("WaveModel");
  TCanvas *canw = new TCanvas(canTitle,canTitle);
  //hWave[MAXSETS-1]->Draw();
  for(int iset=0; iset<MAXSETS; ++iset) hWave[iset]->Draw("same");
  gPad->BuildLegend();
  gPad->SetLogy();


  //for(int iset=0; iset<MAXSETS; ++iset) {
  //  printModel(iset);
  //}

  fout->Append(gSinglet);
  fout->Append(gTriplet);
  fout->Append(gTotTriplet);
  fout->Append(gTot);
  fout->Append(gTotNorm);



  fout->Write();
}
