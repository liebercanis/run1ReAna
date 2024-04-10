enum
{
  NSETS = 5
};

TH1D* hWave0[NSETS];
TH1D* hWave1[NSETS];



void compModel() 
{

  TFile *fout = new TFile("compModel.root","RECREATE");

  TFile *f0 = new TFile("gmodel-2100.root","readonly");
  TFile *f1 = new TFile("gmodel-2538.root","readonly");

  for(int i=0; i<NSETS; ++i) {
    TString hname; hname.Form("WaveSet-%i",i);
    f0->GetObject(hname,hWave0[i]);
    f1->GetObject(hname,hWave1[i]);
    hWave0[i]->SetName(Form("WaveSet-%i-file0",i));
    hWave1[i]->SetName(Form("WaveSet-%i-file1",i));
    fout->Append(hWave0[i]);
    fout->Append(hWave1[i]);
  }


  fout->ls();
  fout->Write();


}
