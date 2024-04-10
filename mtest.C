using namespace TMath;
static double singletFraction = 0.14;
static double tTriplet=2.100; 
static double tSinglet=5.E-3;
static double tXe = 20.0E-3;
static double kxe = 8.8E-2;
static double kplusZero  = 1.3E-1;
double k1Zero = kxe*131./40.;
double binwidth = 5.0E-3;

double getNorm(double ppm){
  double alpha3 = 1.0; //binwidth*(1.-singletFraction);
  double kx = k1Zero*ppm;
  double tq = 1./(1./tTriplet+kplusZero);
  double tr = 1./(kx + 1./tq);
  double norm =  alpha3*tr/tTriplet+alpha3*kx/(kx+1.0/tq);
  //printf(" \t ... ppm %.2f %f %f norm %f \n",ppm,alpha3*tr/tTriplet,alpha3*kx/(kx+1.0/tq) ,norm);
  return norm;
}

void mtest(){
  int max=100;
  std::vector<double> ltot;
  std::vector<double> xppm;

  double maxppm=20.0;
  for(int i =0; i <=max; ++i) {
    double pi = double(i)/double(max)*maxppm;
    double ni = getNorm(pi);
    xppm.push_back(pi);
    ltot.push_back(ni);
    printf("... %i %.0f %f ppm %f  norm %f %f \n",i,pi, double(i)/double(max),pi,ni,maxppm);
  }

   
  TGraph *gl = new TGraph(xppm.size(), &xppm[0],&ltot[0]);
  new TCanvas("ltot","ltot");
  gl->SetTitle("triplet light");
  //gl->SetMarkerSize(0.007);
  gl->SetMarkerStyle(22);
  gl->GetHistogram()->GetXaxis()->SetTitle("PPM");
  gl->GetHistogram()->GetYaxis()->SetTitle("triplet light");
  gl->Draw("apc");


}

