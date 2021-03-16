#include "anaPulses.hh"
//

anaPulses::anaPulses(Int_t irun){
  irun = 1;
  //TString outFileName ; outFileName.Form("recirc%i/PulseAnalysisRecirc_%i.root",irun,irun);
  //TString outFileName ; outFileName.Form("PulseAnalysisLAr_1_%i.root",irun);


  TString dir = "/home/nmcfadden/RooT/PMT/rootData/";

  //TString fileName; fileName.Form("rootData/baconRun_100kohms_%i.root",irun);
  //LED run 
  //TString fileName = "baconRun_runNov2018_1_0.root";
  //TString fileName = "cosmic_1_11.27.2018_10000_events.root";
  //TString fileName = "run_3.root";
  //TString fileName; fileName.Form("rootData/recirc%i.root",irun);
  //TString fileName; fileName.Form("rootData/baconRun_%i.root",irun);
  //TString fileName; fileName.Form("rootData/LAr_%i.root",irun);
  for(int z = 3; z <= 31; z++){
  if(z == 10) continue;
  TString fileName = TString("run_") +to_string(z) +TString(".root");

  TString outFileName = TString("Out_")+fileName; 
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());
  
  TNtuple *ntuplePulse = new TNtuple("ntuplePulse","ntuplePulse","irun:ientry:pmtNum:nhits:charge:startTime:peakWidth:T0:vMax:vMaxTime");
  TNtuple *ntupleEvent = new TNtuple("ntupleEvent","ntupleEvent","irun:ientry:pmtNum:Sdev:baseline:integral:deltaT:deltaV");
  
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(dir+fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  else 
    printf("  found file %s \n",fileName.Data() );

  // get pmtTree from file 
  fin->GetObject("pmtTree",pmtTree);
  Long64_t nentries = pmtTree->GetEntries();
  cout << " number of entries is " << nentries << endl;

  // set up memory for reading
  pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);


  //switch to output file
  outFile->cd();

  
  for(Long64_t ientry = 0; ientry < nentries; ientry++){
    pmtTree->GetEntry(ientry);
    signal.resize(1);derivative.resize(signal.size());integral.resize(derivative.size());

    //if( ientry == 100) break;
   
    /*
    Int_t skipBin = 5;
    if(skipBin < ientry ) break;
    //if( skipBin != ientry) continue;
    */

    NEvents = pmtEvent->time.size();
    Double_t deltaT = (pmtEvent->time[NEvents -1] - pmtEvent->time[0])/(NEvents-1);

    if(ientry == 0){
      hSum[0] = new TH1D("Sum_pmt0","sum_pmt0",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //hSum[1] = new TH1D("Sum_pmt1","sum_pmt1",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hSumBaseline[0] = new TH1D("SumBaseline_pmt0","sum_pmt0",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //hSumBaseline[1] = new TH1D("SumBaseline_pmt1","sum_pmt1",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    }
    signal[0] = pmtEvent->volt1;
    //signal[1] = pmtEvent->volt2;
    //TF1 * fPoly[signal.size()];
    TF1 * fBaseSag[signal.size()];
    TH1D * hBaseFit[signal.size()];
    TH1D * hPeakFinding[signal.size()];
    Int_t nDer = 10, nInt =10;
    Double_t sum = 0;
    for(int j = 0; j < signal.size();j++){
      hSignal[j] = new TH1D(TString("signal_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hIntegral[j] = new TH1D(TString("Integral_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //fPoly[j] = new TF1(TString("Poly_baselineFit")+to_string(j),"pol2(0)",pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hBaseFit[j] = new TH1D(TString("baseFit_")+to_string(j)+TString("_")+to_string(ientry),TString("baseFit")+to_string(j),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hPeakFinding[j] = new TH1D(TString("PeakFinding_")+to_string(j)+TString("_")+to_string(ientry),TString("PeakFinding")+to_string(j),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      /*
      for(int i = 0; i < signal[j].size(); i++){
        hBaseFit[j]->SetBinContent(i,signal[j][i]);
      }
      */
      
      //Int_t window = 500;
      //std::vector<Double_t> baselineFit = BaselineSubtraction(signal[j],window,j,ientry);

      //hBaseFit[j]->Fit(fPoly[j],"Q");
      std::vector<Double_t> integral;
      for(int i = 0; i < signal[j].size(); i++){
        //signal[j][i] = signal[j][i] - fPoly[j]->Eval(pmtEvent->time[i]);
        //signal[j][i] = signal[j][i] - baselineFit[i];
        hSignal[j]->SetBinContent(i+1,signal[j][i]);
        hSum[j]->SetBinContent(i+1,hSum[j]->GetBinContent(i+1)+signal[j][i]);
        sum +=  signal[j][i];
        //hIntegral[j]->SetBinContent(i+1,sum*deltaT);
        integral.push_back(sum);
      }
      signal[j] = Derivative(integral,nInt);
      hDerivative[j] = new TH1D(TString("derivative_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    }
    for(int j = 0; j < signal.size(); j++){
      //fBaseSag[j] = new TF1(TString("Baseline_Sag")+to_string(j),"[0]*TMath::Exp(-x/[1])*TMath::Sin(2*(x-[2])/([1]))",pmtEvent->time[hSignal[j]->GetMaximumBin()],4e-6); 
      //fBaseSag[j]->SetParameter(2,0);
      fBaseSag[j] = new TF1(TString("Baseline_Sag")+to_string(j),"[0]*TMath::Exp(-x/[1])",pmtEvent->time[hSignal[j]->GetMinimumBin()],4e-6); 
      fBaseSag[j]->SetParameter(0,-hSignal[j]->GetMinimum()/12.5);
      fBaseSag[j]->SetParameter(1,3e-6);
      //rezero sum after baseline
      sum = 0;
      Double_t sigMin = hSignal[j]->GetMinimum();
      Int_t sigMinBin = hSignal[j]->GetMinimumBin();//FindBin(sigMin);
      for(int i = 0; i < signal[j].size(); i++){
        signal[j][i] /= (double) nInt;
        hBaseFit[j]->SetBinContent(i+1,signal[j][i]);
        if(pmtEvent->time[sigMinBin] < pmtEvent->time[i]){
          signal[j][i] -= fBaseSag[j]->Eval(pmtEvent->time[i]);
        }
        sum += signal[j][i];
        hIntegral[j]->SetBinContent(i+1,sum*deltaT);
        hSignal[j]->SetBinContent(i+1,signal[j][i]);
        hPeakFinding[j]->SetBinContent(i+1,signal[j][i]);
      }
      derivative[j] = Derivative(signal[j],nDer);
      for(int i = 0; i < derivative[j].size(); i++){
        hDerivative[j]->SetBinContent(i+1,derivative[j][i]);
      }
    }
    //Peak Finding
    std::vector< std::vector<Int_t> > peakTime;peakTime.resize(signal.size());
    for(int i = 0; i < signal.size(); i++){
      std::vector<Int_t> pTime = PeakFinding(signal[i],derivative[i],i,ientry);
      if(pTime.size() == 0) continue;
      //cout<<pTime.size()<<" "<<ientry<<endl;
      //Cull the small pulses
      for(int j = 0; j< pTime.size() -1 ;j +=2){
        Int_t startBin = pTime[j],stopBin = pTime[j+1];
        Double_t peakWidth = std::fabs(pmtEvent->time[stopBin] - pmtEvent->time[startBin] );
        //cout<<"\t peakWidth "<<peakWidth<<endl;
        if(peakWidth > minPeakWidth && peakWidth < maxPeakWidth) {
          peakTime[i].push_back(pTime[j]);
          peakTime[i].push_back(pTime[j+1]);
        }
      }
    }
    //cout<<"baseline "<<endl;
    //Calculate baseline and fill ntuple
    for(int j = 0; j < signal.size(); j++){
      ///*
      std::pair<Double_t, Double_t> baseLine_Sdev;
      if(peakTime[j].empty())
        baseLine_Sdev = BaselineSubtraction(signal[j],j,ientry);
      else
        baseLine_Sdev = BaselineSubtraction(signal[j],peakTime[j],j,ientry);
      //Double_t baseline = baseLine_Sdev.first;
      Double_t baseline = CalculateMean(signal[j],0,500);
      //Double_t Sdev = baseLine_Sdev.second;
      Double_t Sdev =  CalculateSdev(signal[j],0,500,baseline);
      ntupleEvent->Fill(irun,ientry,j,Sdev,baseline,-sum,deltaT,Sdev*deltaT);
      //*/
      //hBaseline[j] = new TH1D(TString("Baseline_Channel_")+to_string(j)+TString("_event_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //hSignal[j] = new TH1D(TString("signal_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      TString sTitle = TString("");
      /*
      //Double_t baselineConstant = BaselineSubtraction(signal[j],j,ientry).first;
      //Int_t window = 60;
      //std::vector<Double_t> baselineFit = BaselineSubtraction(signal[j],window,j,ientry);
      for(int i = 0; i < signal[j].size(); i++){
        Int_t timeBin = hSignal[j]->FindBin(pmtEvent->time[i]);
        //hSignal[j]->SetBinContent(timeBin,-signal[j][i]);
        //signal[j][i] = signal[j][i]- baselineFit[i];
        //signal[j][i] = signal[j][i]- baselineConstant;
        //signal[j][i] = signal[j][i]- baseline;
        hSignal[j]->SetBinContent(timeBin,-signal[j][i]);
      }
      
      for(int i = 0;i < signal[j].size(); i++){
        Int_t timeBin = hSignal[j]->FindBin(pmtEvent->time[i]);
        hSignal[j]->SetBinContent(timeBin,-signal[j][i]);
        //invert signal to positive pulses
        signal[j][i] = signal[j][i]- baselineConstant;
        //signal[j][i] = signal[j][i]- baselineFit[i];
        hSum[j]        ->SetBinContent(timeBin,-signal[j][i]+hSum[j]->GetBinContent(timeBin));
        //hSumBaseline[j]->SetBinContent(timeBin,-baselineFit[i]+hSumBaseline[j]->GetBinContent(timeBin));
      }
      */
      if(peakTime[j].empty()) continue;
      //cout<<"ntupleFill"<<endl;
      Double_t T0=0;
      //cout<<"Event "<<ientry<<endl;
      for(int i = 0; i < peakTime[j].size()-1; i += 2){
        Double_t vMax = 0,charge = 0,vMaxTime = 0,startTime,stopTime,peakWidth;
        Int_t startBin = peakTime[j][i],stopBin = peakTime[j][i+1];
        //Int_t startBin = hSignal[j]->FindBin(400e-9);
        //Int_t stopBin  = hSignal[j]->FindBin(500e-9);
        startTime = pmtEvent->time[startBin];
        if(i == 0) T0 = startTime;
        stopTime = pmtEvent->time[stopBin];
        peakWidth = stopTime - startTime;
        if(peakWidth < 0) cout<<"Negative pulse width "<<peakWidth<<" run "<<z<<" entry " <<ientry<<endl;
        for(int k = startBin; k < stopBin; k++){
          charge += signal[j][k];
          if(std::fabs(signal[j][k]) > std::fabs(vMax) ){
            vMax = signal[j][k];
            vMaxTime = pmtEvent->time[k];
          }
          hPeakFinding[j]->SetBinContent(k,0);
        }
        //cout<<"\t Start/stop "<<startTime<<"/"<<stopTime<<", width "<<peakWidth<<", vMax "<<vMax<<", charge "<<-1e9*deltaT*charge<<endl;
        ntuplePulse->Fill(irun,ientry,j,peakTime[j].size()/2,-charge*deltaT,startTime,peakWidth,T0,-vMax,vMaxTime,Sdev,baseline);
        sTitle += TString("Charge_")+to_string(-1e9*deltaT*charge)+TString("_start/stop_")+to_string(1.e6*startTime)+TString("/")+to_string(1.e6*stopTime)+TString("_vMax_")+to_string(vMax)+TString("_");
      }
      hSignal[j]->SetTitle(sTitle);
    }

    if(ientry%1000 == 0) cout<<"completed "<<ientry<<" events"<<endl;

    //###################//
    //End of Event Clean-Up//
    //###################//

    if(ientry > NHistograms){
      for(int i = 0; i < signal.size(); i++){
        delete hSignal[i];
        delete hDerivative[i];
        //delete fPoly[i];
        delete fBaseSag[i];
        delete hBaseFit[i];
        delete hIntegral[i];
        delete hPeakFinding[i];
//        delete hBaseline[i];
        //delete hFilter[i];
      }
    }
    signal.clear();
    derivative.clear();
  }
  //ntuplePulse->Write();
  outFile->Write();
  }
  cout<<"you did it! your code ran with crashing"<<endl;
}

std::vector<Double_t> anaPulses::Derivative(std::vector<Double_t> sig,Int_t N){
  std::vector<Double_t> diff;
  if(N%2 != 0) N++;
  
  for(int i = 0; i < sig.size(); i++){
    Double_t front= 0, back = 0;
    Double_t nFront = 0, nBack = 0;

    for(int j = i;j<i+N;j++){
      if(j >= sig.size()) continue;
      front += sig[j];
      //cout<<front<<" "<<sig[j]<<" "<<j<<endl;
      nFront++;
    }
    for(int j = i;j>i-N;j--){
      if(j < 0) continue;
      back  += sig[j];
      nBack++;
    }
    if(nFront == 0 || nBack == 0){
      diff.push_back(0);
      cout<<"nFront = "<<nFront<<", nBack = "<<nBack<<endl;
    }
    //Double subtraction leaving tiny remaining difference
    else if( log(std::fabs(front/nFront-back/nBack))  < -20){
      diff.push_back(0);
    }
    else
      diff.push_back(front/nFront-back/nBack);
  }
  
  return diff;
}
std::vector<Double_t> anaPulses::Integral(std::vector<Double_t> sig,Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> integ;
  Double_t sum = 0;
  for(int i = 0; i < sig.size(); i++){
    sum += sig[i];
    integ.push_back(sum);
  }
  return integ;
}
std::vector<Int_t> anaPulses::PeakFinding(std::vector<Double_t> sig, std::vector<Double_t> diff,Int_t pmtNum,Int_t ientry){
  std::vector<Int_t> peakTime;
  /*
  std::vector<Double_t> sort = diff;
  std::sort(sort.begin(),sort.end());
  Double_t RMS = 5*std::fabs(sort[sort.size()*.68]);
  hDerivative[pmtNum]->SetTitle(TString("Threshold_")+to_string(RMS));
  */
  //Int_t rmsBin = hDerivative[pmtNum]->FindBin(150e-9);
  Int_t rmsBin = 500;//sig.size();
  //2.25 for Scint Runs
  //3.5 for LED runs
  Double_t loRMS = 3.0*CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  //2.25 for Scint Runs
  //2.5 for LED runs
  Double_t hiRMS = 2.5*CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  hDerivative[pmtNum]->SetTitle(TString("loThreshold_")+to_string(loRMS)+TString("hiThreshold_")+to_string(hiRMS)+TString("_mean_")+to_string(CalculateMean(diff)) );
  Bool_t low1 = false,low2 = false,high1 = false, high2 = false, ll = false, hh = false;

  Int_t LowStartBin = -1,LowStopBin = -1,HighStartBin = -1,HighStopBin = -1;
  for(int i = 1; i < diff.size(); i++){
    if(diff[i] <= -loRMS && diff[i-1] > -loRMS){
      LowStartBin = i;
      low1 =  true;
      //cout<<pmtEvent->time[LowStartBin]<<"  "<<diff[i]<<" "<<diff[i-1]<<endl;
    }
    if(diff[i] >= -loRMS && diff[i-1] < -loRMS ){
      LowStopBin = i;
      low2 = true;
    }

    if(diff[i] >= hiRMS && diff[i-1] < hiRMS ){
      HighStartBin = i;
      high1 = true;
    }
    if(diff[i] <= hiRMS && diff[i-1] > hiRMS ){
      HighStopBin = i;
      high2 = true;
    }
    
    if(low1 && low2) {
      //cout<<"low "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<" "<<diff[i]<<" "<<-RMS<<endl;
      //cout<<"\thi "<<pmtEvent->time[HighStartBin]<<" "<<pmtEvent->time[HighStopBin]<<endl;
      peakTime.push_back(LowStartBin);
      peakTime.push_back(LowStopBin);
      low1 = false; low2 = false;
      if(hh)
        ll = false;
      else
        ll = true;
      //LowStartBin = -1;LowStopBin = -1;
    }
    
    if(high1 && high2){
      //cout<<"hi "<<pmtEvent->time[HighStartBin]<<" "<<pmtEvent->time[HighStopBin]<<" "<<diff[i]<<" "<<RMS<<endl;
      //cout<<"\tlow "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<endl;
      peakTime.push_back(HighStartBin);
      peakTime.push_back(HighStopBin);
      high1 = false;high2 = false;
      if(ll)
        hh = true;
      else 
        hh = false;
      //HighStartBin = -1;HighStopBin = -1;
    }
    
    if( ll && hh) {
      //cout<<"llhh "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[HighStopBin]<<endl;
      peakTime.erase(peakTime.end() - 4,peakTime.end() );
      peakTime.push_back(LowStartBin);
      peakTime.push_back(HighStopBin);
      low1 = false; low2 = false;
      high1 = false;high2 = false;
      ll = false;hh = false;
      //LowStartBin = -1;LowStopBin = -1;HighStartBin = -1;HighStopBin = -1;
    }
   
  }

  return peakTime;
}

//calc baseline and exclude region where peaks are found
std::pair<Double_t,Double_t> anaPulses::BaselineSubtraction(std::vector<Double_t> sig, std::vector<Int_t> weight, Int_t pmtNum, Int_t ientry){
  Int_t startBin = 0, stopBin = 0;
  
  for(int i = 0; i < weight.size() - 1; i += 2){
    startBin = weight[i];
    stopBin = weight[i+1];
    for(int k = 0; k< sig.size(); k++){
      if(k >= startBin && k <= stopBin)
        sig.erase(sig.begin()+k);
    } 
  }
  //std::sort(sig.begin(),sig.end());

  return std::make_pair(CalculateMean(sig),CalculateSdev(sig,0,sig.size(),CalculateMean(sig)) );
}
//calc moving bassline
//P. Funk, G. Funk ect...
std::vector<Double_t> anaPulses::BaselineSubtraction(std::vector<Double_t> sig, Int_t weight, Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> baseline;
  std::vector<Double_t> gaus;
  //weight should be even
  if(weight%2 == 1) weight++;

  Double_t sigma = weight;
  Double_t norm = 0;
  for(int i = -weight; i <= weight; i++){
    Double_t val = exp(-i*i/(2*sigma*sigma));
//    gaus.push_back((1/sqrt(sigma*sigma*3.1415*2))*exp(-i*i/(2*sigma*sigma)));
    gaus.push_back(val);
    norm += val;
  }
  std::vector<Double_t> median;
  //for(int i = sig.size() - 1; i >= 0; i--){
  for(int i = 0; i < sig.size(); i++){
    std::vector<Double_t> localMedian;
    for(int j = i-weight; j <= i+weight; j++){
      if(j < 0) continue;
      if(j >= sig.size() ) continue;
      localMedian.push_back(sig[j]);  
    }
    //sorts vector from lowest to highest
    std::sort(localMedian.begin(),localMedian.end());
    //middle value is the median of the window per Local Median method
    //tuned value to 0.8 to better match long/large peak baselines
    //N. McFadden 
    median.push_back(localMedian[(int)localMedian.size()*.5]);
    //median.push_back(localMedian[localMedian.size()*.2]);
  }
  for(int i = 0; i < median.size(); i++){
  //for(int i = sig.size() - 1; i >= 0; i--){
    Double_t val = 0;
      for(int j = i-weight; j <= i+weight; j++){
        if(j < 0) continue;
        if(j >= median.size() ) continue;
        val += median[j]*gaus[j+weight-i]/norm;
        //cout<<"\t"<<median[j]<<" "<<gaus[j+weight-i]/norm<<" "<<i<<" "<<j<<" "<<j+weight-i<<endl;;
      }
      baseline.push_back(val);
  }
  //hBaseline[pmtNum] = new TH1D(TString("Baseline_Channel_")+to_string(pmtNum)+TString("_event_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
  for(int i = 0; i < baseline.size();i++){
    //Int_t timeBin = hBaseline[pmtNum]->FindBin(pmtEvent->time[i]);
    //hBaseline[pmtNum]->SetBinContent(timeBin,-baseline[i]);
  }
  return baseline;
}
//calculate baseline of whole waveform
std::pair<Double_t,Double_t> anaPulses::BaselineSubtraction(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry){
  std::sort(sig.begin(),sig.end());
  return std::make_pair(sig[sig.size()/2],sig[sig.size()*.68]);
}

// seems like it is some type of smoothing filter...
std::vector<Double_t> anaPulses::TrigFilter(std::vector<Double_t> sig, Int_t N, Int_t pmtNum, Int_t ientry){
  
  if(N%2 != 0) N++;
  std::vector<Double_t> vec(sig.size(),0); 
  for(int i = 0; i < sig.size(); i++){
    Double_t topSum = 0,botSum = 0;
    for(int j = i-N/2; j <=i+N/2; j++){
      if(j < 0) continue;
      if(j >= sig.size() ) continue;
      topSum += sig[j]*(1+std::cos((j-i)*2*M_PI/N));
      botSum += (1+std::cos((j-i)*2*M_PI/N));
    }
    if(botSum == 0) continue;
    vec[i] = topSum/botSum;
  }

  return vec;
}

Double_t anaPulses::CalculateMean(std::vector<Double_t> vec){
  Double_t mean = 0;
  for(int i = 0; i< vec.size(); i++) mean += vec[i];
  return mean/vec.size();
}

Double_t anaPulses::CalculateMean(std::vector<Int_t> vec){
  Double_t mean = 0;
  for(int i = 0; i< vec.size(); i++) mean += vec[i];
  return mean/vec.size();
}

Double_t anaPulses::CalculateMean(std::vector<Double_t> vec,Int_t startBin, Int_t stopBin){
  Double_t mean = 0;
  Double_t N = (stopBin-startBin);
  for(int i = startBin; i< stopBin; i++) mean += vec[i];
  return mean/N;
}

Double_t anaPulses::CalculateSdev(std::vector<Double_t> vec,Int_t startBin, Int_t stopBin, Double_t mean){
  Double_t sdev = 0;
  for(int i = startBin; i < stopBin; i++){
    sdev += (mean-vec[i])*(mean-vec[i]);
  }
  return std::sqrt(sdev/(stopBin-startBin));
}

std::vector<Double_t> anaPulses::BubbleSort(std::vector<Double_t> A){
  int i, j, N = A.size();

  for (i = 0; i < N; i++){
    for (j = i + 1; j < N; j++){
      if (A[j] < A[i]){
        Double_t buff;
        buff = A[i];
        A[i] = A[j];
        A[j] = buff;
      }
    }
  }
  return A;
}

std::vector<Double_t> anaPulses::TrapFilter(std::vector<Double_t> sig, Int_t ramp, Int_t flat,Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> filter;
  for(int i = 0; i < sig.size(); i++){
    Double_t val1 = 0,val2 = 0;
    for(int j = i; j < ramp+i; j++){
      if(j >= sig.size()) continue;
      val1 += sig[j];
    }
    val1 /= ramp;
    for(int j = i +ramp +flat; j < i+ 2*ramp + flat; j++){
      if(j >= sig.size()) continue;
      val2 += sig[j];
    }
    val2 /= ramp;
    filter.push_back(val2-val1); 
  }
  /*
  hFilter[pmtNum] = new TH1D(TString("Trap_Filter_Channel_")+to_string(pmtNum)+TString("_event_")+to_string(ientry)+TString("_flat_")+to_string(flat)+TString("_ramp_")+to_string(ramp),"",filter.size(),0,filter.size());
  for(int i = 0; i < filter.size(); i++){
    hFilter[pmtNum]->SetBinContent(i,filter[i]);
  }
  */

  return filter;
}
