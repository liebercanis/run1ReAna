using namespace TMath;

void noise(int nnoise = 1000)
{
  
  TFile *fout = new TFile("noiseTest.root","recreate");
  double s=10;
  int nsig=1000000;
  int ntot = nsig+nnoise;

  TRandom3 *ran = new TRandom3();

  TH1D* hWave = new TH1D("wave",Form("sum-nsig-%i-nev-%i",nsig,ntot),1000,0,10);
  TH1D* hSig = (TH1D*) hWave->Clone("sig");
  TH1D* hNoise = (TH1D*) hWave->Clone("noise");
  hNoise->SetTitle(Form("noise-%i",nnoise));
  hSig->SetTitle(Form("signal-%i",nsig));
  hNoise->SetLineColor(kRed);
  //hNoise->SetFillColor(kRed);
  //hNoise->SetFillStyle(3315);
  hWave->Reset();
  hSig->Reset();
  hNoise->Reset();

  double tau = 2.1/1000.;
  for(int j=0; j<nsig; ++j) {
    //double val = s*(ran->Exp(tau));
    double val = s*(-tau*Log(1.-ran->Rndm()) );
    hWave->Fill(val);
    hSig->Fill(val);
  }

  int nbins = hWave->GetNbinsX();

  for(int iev = 0; iev<nnoise; ++iev) {
      for(int ibin = 0; ibin<nbins; ++ibin) {
        double val = ran->Gaus();
        //val= ran->Gaus();
        hWave->SetBinContent(ibin,  hWave->GetBinContent(ibin) + val);
        hNoise->SetBinContent(ibin,  hNoise->GetBinContent(ibin) + val);
      }
  }

  if(ntot>0){ 
    hWave->Scale(1./double(ntot));
    hNoise->Scale(1./double(ntot));
    hSig->Scale(1./double(ntot));

  }
  //hWave->SetFillColor(kBlue);
  hNoise->SetMarkerColor(kRed);
  hNoise->SetMarkerStyle(20);
  hNoise->SetMarkerSize(.5);
  hWave->SetMarkerColor(kBlack);
  hWave->SetMarkerStyle(21);
  hWave->SetMarkerSize(.5);
  hSig->SetLineColor(kGreen);
  hSig->SetMarkerColor(kGreen);
  hSig->SetMarkerStyle(22);
  hSig->SetMarkerSize(.5);



  //hWave->SetFillStyle(3351);

  //hWave->GetYaxis()->SetRangeUser(-1.5*s,1.5*s);

  TCanvas* can = new TCanvas(Form("noiseTest-nsig-%i-nev-%i",nsig,ntot) ,Form("noiseTest-%i",ntot));
  hWave->Draw("histp");
  hNoise->Draw("same");
  hSig->Draw("same");
  hWave->Draw("samep");
  gStyle->SetOptLogy();
  can->BuildLegend();

  can->Print(".png");

  fout->Write();

}

