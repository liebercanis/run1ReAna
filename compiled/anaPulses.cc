#include "anaPulses.hh"
//
std::vector<std::complex<double>> anaPulses::FFT(std::vector<double> vin)
{
  std::vector<std::complex<double>> VectorComplex;
  int nsamples = (int)vin.size();
  for (int is = 0; is < nsamples; ++is)
    fFFT->SetPoint(is, vin[is]);
  fFFT->Transform(); //

  std::vector<Double_t> realVec, imVec;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im); //.real or .imag accessors
    VectorComplex.push_back(c);
  }
  return VectorComplex;
}

std::vector<double> anaPulses::inverseFFT(std::vector<std::complex<double>> vin)
{
  std::vector<Double_t> vout;
  int nsamples = (int)vin.size();
  for (int is = 0; is < nsamples; ++is)
  {
    fInverseFFT->SetPoint(is, vin[is].real(), vin[is].imag());
  }
  fInverseFFT->Transform();
  double norm = 0;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl = fInverseFFT->GetPointReal(i);
    norm += rl;
    vout.push_back(rl);
  }

  //for(unsigned j=0; j<vout.size(); ++j) vout[j] *= sum/norm;
  for (unsigned j = 0; j < vout.size(); ++j)
    vout[j] /= double(nsamples);

  return vout;
}

bool anaPulses::initDecon()
{
  bool rc = false;
  TFile *fin = NULL;
  fin = new TFile("decon.root", "readonly");
  vResponseFFT.clear();
  if (!fin)
    return rc;
  if (!fin->IsOpen())
    return rc;
  // find summed hist
  TList *list = fin->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  TGraph *gr = NULL;
  while ((key = (TKey *)next()))
  {
    obj = key->ReadObj();
    if (obj->InheritsFrom("TGraph"))
    {
      TString hname = obj->GetName();
      if (hname.Contains("complexShapeFFT"))
      {
        gr = (TGraph *)obj->Clone();
        break;
      }
    }
  }

  if (!gr)
    return rc;

  double real, imag;
  for (int i = 0; i < gr->GetN(); ++i)
  {
    gr->GetPoint(i, real, imag);
    std::complex<double> c(real, imag); //.real or .imag accessors
    vResponseFFT.push_back(c);
  }

  rc = true;
  printf(" initDecon success points %ld  \n", vResponseFFT.size());

  return rc;
}

std::vector<Double_t> anaPulses::decon(std::vector<Double_t> inputWave)
{
  std::vector<Double_t> outWave;
  outWave.resize(inputWave.size());

  // make deconvolution FFT for Wave
  std::vector<std::complex<double>> vWaveFFT;
  vWaveFFT = FFT(inputWave);

  std::vector<std::complex<double>> vDWaveFFT;
  vDWaveFFT.resize(vWaveFFT.size());
  for (unsigned iw = 0; iw < vWaveFFT.size(); ++iw)
  {
    std::complex<double> d = vWaveFFT[iw];
    std::complex<double> g = vResponseFFT[iw];
    std::complex<double> gstar = std::conj(g);
    double c = std::abs(d);
    double c2 = std::abs(g);
    if (c2 > 0)
      vDWaveFFT[iw] = gstar * d / (g * gstar + complexSmall);
    else
      vDWaveFFT[iw] = d;
  }
  outWave = inverseFFT(vDWaveFFT);

  return outWave;
}

anaPulses::anaPulses(TString tag, Int_t maxEvents)
{
  printf(" starting anaPulses tag %s \n", tag.Data());
<<<<<<< HEAD
  if (!initDecon())
  {
    printf(" !!!! initDecon failed, abort!  \n");
    return;
  }

  // spe fit func
  fPulseBaseline = new TF1("fPulseBaseline", "[0]+x*[1]", 0., 90.);
  fPulseBaseline->SetParName(0, "const");
  fPulseBaseline->SetParName(1, "slope");
  fPulseBaseline->SetParameters(0., 0.);

=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
  theTag = tag;

  for (int num = 0; num < NPMT; ++num)
  {
    pulseFirstShapeNorm[num] = 0;
    pulseShapeNorm[num] = 0;
  }

  // open ouput file and make some histograms

  TString fileName;
  fileName.Form("/data1/bacon/rootData/DS2/%s.root", tag.Data());
  printf(" looking for file %s\n", fileName.Data());
  TFile *fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
  {
    printf(" couldnt open file %s\n", fileName.Data());
    return;
  }
  else
    printf("  found file %s \n", fileName.Data());
  // get pmtTree from file
  TTree *pmtTree = NULL;
  fin->GetObject("pmtTree", pmtTree);
  Long64_t nentries = pmtTree->GetEntries();

  // set up memory for reading
  pmtEvent = new TPmtEvent();
  Int_t iset = pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
  if (iset != 0)
  {
    printf(" couldnt find pmtEvent in file  %s\n", fileName.Data());
  }

  //switch to output file
  TString outfileName;
  outfileName.Form("run1RootData/anaPulsesDecon_%s_%i.root", tag.Data(), maxEvents);
  printf(" opening output file %s \n", outfileName.Data());
  outfile = new TFile(outfileName, "recreate");
  pDir = outfile->mkdir("pulseHist");
<<<<<<< HEAD
  evDir = outfile->mkdir("eventHist");
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
  outfile->cd();
  tbTree = new TTree("TBacon", "TBacon data");
  bEvent = new TBaconEvent();
  tbTree->Branch("bevent", &bEvent);

  /** from pmtSim **/
  ntuplePulse = new TNtuple("ntuplePulse", "ntuplePulse", "ientry:pmtNum:nhits:charge:startTime:peakWidth:T0:V0:C0:vMax:vMaxTime");
  //TNtuple *ntupleEvent = new TNtuple("ntupleEvent","ntupleEvent","irun:ientry:pmtNum:Sdev:baseline:integral:deltaT:deltaV:vMax:tMax");
  ntupleEvent = new TNtuple("ntupleEvent", "ntupleEvent", "irun:ientry:pmtNum:Sdev:baseline:integral:deltaT:tMax:vMax:cMax:totalCharge:windowSum:nHits");
  hWindowSum[0] = new TH1D("WindowSum_0", "WindowSum_0", 500, -1e-11, 5e-9);
  hWindowSum[1] = new TH1D("WindowSum_1", "WindowSum_1", 500, -1e-11, 5e-9);
  hSumWaveform[0] = new TH1D("SumWaveform0", "SumWaveform0", 10000, 0, 10e-6);
  hSumWaveform[1] = new TH1D("SumWaveform1", "SumWaveform1", 10000, 0, 10e-6);
  processedTree = new TTree("processedTree", "processedTree");
  processedEvent = new TPmtEvent();
  processedTree->Branch("processedEvent", &processedEvent);
  pmtRun = new TPmtRun();
  processedTree->Branch("pmtRun", &pmtRun);

  hPulseCharge = new TH1D("PulseCharge", "PulseOne", 1000, 0, 100);
  pDir->cd();
  hPulseOne = new TH1D("PulseOne", "PulseOne", 400, 0, 400);
<<<<<<< HEAD
  //hPulseOne->Sumw2();
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
  outfile->cd();

  // set up memory for reading
  cout << "Event branch set" << endl;
  /** ***/

  // loop over entries
  if (maxEvents > 0)
    nentries = maxEvents;
  printf(" STARTING RUN %s with  events  %lld of %lld  \n", tag.Data(), nentries, pmtTree->GetEntries());
  eventsSaved = 0;
  maxEventsSaved = 100;
  int gotPMT = 0;
  for (Long64_t ientry = 0; ientry < nentries; ientry++)
  {
    bEvent->clear();
    pmtTree->GetEntry(ientry);
    nSamples = pmtEvent->time.size();

    if (ientry == 0)
    {
      printf(" .... events %lld samples %i PMT0 %zu PMT1 %zu \n", pmtTree->GetEntries(), nSamples, pmtEvent->volt1.size(), pmtEvent->volt2.size());
      // initialize fft
      fFFT = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
      fInverseFFT = TVirtualFFT::FFT(1, &nSamples, "C2R M K");
    }
    if (ientry % 100 == 0)
      cout << "\t \t .... " << ientry << " events" << endl;

    gotPMT = 0;
    if (pmtEvent->volt1.size() > 0)
      ++gotPMT;
    if (gotPMT < 1)
      break;
    // define pmt signal histograms
    Double_t pmtXLow = pmtEvent->time[0] * microSec;
    Double_t pmtXHigh = pmtEvent->time[nSamples - 1] * microSec;
    // define histos
    if (ientry == 0)
    {
      TString hname;
      for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
      {
        hname.Form("WaveSum%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hWaveSum[ipmt] = new TH1D(hname, hname, nSamples, 0, nSamples);
        hname.Form("PulseSum%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hPulseSum[ipmt] = (TH1D *)hPulseOne->Clone(hname);
<<<<<<< HEAD
        hname.Form("PulseSumUn%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hPulseSumUn[ipmt] = (TH1D *)hPulseOne->Clone(hname);
        hname.Form("PulseFirstSum%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hPulseFirstSum[ipmt] = new TH1D(hname, hname, 1500, 0, 1500);
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
        //hPulseSum[ipmt] = new TH1D(hname, hname, 1500, 0, 1500);
        hname.Form("Sum_pmt%i", ipmt);
        hSum[ipmt] = new TH1D(hname, hname, nSamples, pmtEvent->time[ipmt], pmtEvent->time[nSamples - 1]);
        hname.Form("SumBaseline_pmt%i", ipmt);
        hSumBaseline[ipmt] = new TH1D(hname, hname, nSamples, pmtEvent->time[ipmt], pmtEvent->time[nSamples - 1]);
        outfile->cd();
      }
      evDir->cd();
      for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
      {
        hname.Form("WaveRawOne%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hWaveRawOne[ipmt] = new TH1D(hname, hname, nSamples, pmtEvent->time[ipmt], pmtEvent->time[nSamples - 1]);
        evDir->cd();
        hname.Form("EventOne%i-Ev%lld-%s", ipmt, nentries, tag.Data());
        hWaveOne[ipmt] = new TH1D(hname, hname, nSamples, pmtEvent->time[ipmt], pmtEvent->time[nSamples - 1]);
      }
      outfile->cd();
      outfile->ls();
    }

    // event histograms
    for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
    {
      hWaveOne[ipmt]->Reset();
      hWaveRawOne[ipmt]->Reset();
    }

    /*for (unsigned isample = 0; isample < pmtEvent->volt1.size(); isample++)
    {
      Double_t volt0 = pmtEvent->volt1[isample];
      Double_t digi0 = double(volt0);
      hWaveRawOne[0]->SetBinContent(isample + 1, digi0);
    }
    */
    //hWaveRawOne[0]->SetEntries(double(pmtEvent->volt1.size()));

    /* do analysis here */
    anaEntry(ientry);
    fillBaconEvent(ientry);
<<<<<<< HEAD
    // save/clone  event histograms
    if (eventsSaved < maxEventsSaved)
    {
      ++eventsSaved;
      evDir->cd();
      for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
      {
        TString hEvName;
        hEvName.Form("Wave%i-event%lld-%s", ipmt, ientry, tag.Data());
        hWave[ipmt] = (TH1D *)hWaveOne[ipmt]->Clone(hEvName);
        hEvName.Form("WaveRaw%i-event%lld-%s", ipmt, ientry, tag.Data());
        hWaveRaw[ipmt] = (TH1D *)hWaveRawOne[ipmt]->Clone(hEvName);
      }
      //printf(" event %lld  saved %i  %s entries in list %i in histo  %f\n ", ientry, eventsSaved, hWave[0]->GetName(), evDir->GetList()->GetEntries(),hWaveRaw[0]->GetEntries());
    }
    else
    {
      for (int j = 0; j < signal.size(); j++)
      {
        delete hSignal[j];
        delete hDerivative[j];
        delete hPeakFinding[j];
      }
    }
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
  }
  outfile->cd();

<<<<<<< HEAD
  printf(" END RUN %s  %lld pulses %lld \n", tag.Data(), nentries, ntuplePulse->GetEntries());
  printf(" pulse norm first %i late %i  \n", pulseFirstShapeNorm[0], pulseShapeNorm[0]);
  printf(" events %i saved \n ", eventsSaved);
  delete hPulseOne;
  for (int ipmt = 0; ipmt < gotPMT; ++ipmt)
  {
    delete hWaveOne[ipmt];
    delete hWaveRawOne[ipmt];
  } // normalize
=======
  //outfile->ls();
  delete hPulseOne;
  // normalize
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
  TH1D *hres = hPulseSum[0];
  if (pulseShapeNorm[0])
  {
    for (int ibin = 0; ibin < hres->GetNbinsX(); ++ibin)
      hres->SetBinContent(ibin, hres->GetBinContent(ibin) / double(pulseShapeNorm[0]));
  }
<<<<<<< HEAD
  hres = hPulseSumUn[0];
  if (pulseShapeNorm[0])
  {
    for (int ibin = 0; ibin < hres->GetNbinsX(); ++ibin)
      hres->SetBinContent(ibin, hres->GetBinContent(ibin) / double(pulseShapeNorm[0]));
  }
  hres = hPulseFirstSum[0];
  if (pulseFirstShapeNorm[0])
  {
    for (int ibin = 0; ibin < hres->GetNbinsX(); ++ibin)
      hres->SetBinContent(ibin, hres->GetBinContent(ibin) / double(pulseFirstShapeNorm[0]));
  }
  outfile->ls();
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
  outfile->Write();
}

void anaPulses::fillBaconEvent(Long64_t ientry)
{
  Int_t jpmt = 0;
  bEvent->clear();
  bEvent->event = ientry;
  bEvent->run = atoi(theTag(theTag.Last('_') + 1, theTag.Length() - theTag.Last('_')).Data());
  bEvent->npmt = 1;
  bEvent->totQ = pmtRun->totalCharge[jpmt];
  bEvent->T0 = pmtRun->T0[jpmt];
  bEvent->totalCharge = pmtRun->totalCharge[jpmt];
  bEvent->tMax = pmtRun->tMax[jpmt];
  bEvent->vMax = pmtRun->vMax[jpmt];
  bEvent->cMax = pmtRun->cMax[jpmt];
  bEvent->baseline = pmtRun->baseline[jpmt];
  bEvent->sDev = pmtRun->sDev[jpmt];
  Int_t nPulses = pmtRun->charge[jpmt].size();
  bEvent->npulse = nPulses;
  bEvent->spe = 4.60e-11;
  bEvent->muVmax = pmtRun->vMax[jpmt];
  bEvent->F40 = 0;

  // list of pulse start times in ns
  std::vector<int> hitStartTime;
  std::vector<int> hitEndTime;
  for (int ip = 0; ip < nPulses; ip++)
  {
    hitStartTime.push_back(int(pmtRun->startTimes[jpmt][ip] * 1E9));
    hitEndTime.push_back(int(pmtRun->startTimes[jpmt][ip] * 1E9) + int(pmtRun->peakWidths[jpmt][ip] * 1E9));
  }

  // fill puleses
  for (int ip = 0; ip < nPulses; ip++)
  {
    Double_t charge = pmtRun->charge[jpmt][ip];
    Double_t startTime = pmtRun->startTimes[jpmt][ip];
    Double_t peakWidth = pmtRun->peakWidths[jpmt][ip];
    Double_t peakHeight = pmtRun->peakHeights[jpmt][ip];
    Double_t peakTime = pmtRun->peakTimes[jpmt][ip];

    TPulse thePulse;
    thePulse.istart = int(startTime * 1E9);
    thePulse.iend = int(startTime * 1E9) + int(peakWidth * 1E9);
    thePulse.time = startTime;
    thePulse.tpeak = peakTime;
    thePulse.pwidth = peakWidth;
    thePulse.peak = peakHeight;
    thePulse.q = charge;
<<<<<<< HEAD
    int pLength = 100; // add a bit to each end of pulse
=======
    int pLength = 50; // add a bit to each end of pulse
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
    //thePulse.qerr=phitQErr;
    if (thePulse.iend - thePulse.istart < 1)
      continue;
    // fill pulse digi
    for (unsigned idigi = TMath::Max(0, thePulse.istart - pLength); idigi < TMath::Min(int(signal[jpmt].size()), thePulse.iend + pLength); ++idigi)
      thePulse.digi.push_back(signal[jpmt][idigi]);
    //printf(" xxxx iev %ldd ip %i %E %i %i digi %u \n", ientry, ip, startTime, thePulse.istart, thePulse.iend, thePulse.digi.size());
    //printf(" xxx ipeak %i size %u  \n", ipeak, thePulse.digi.size());
    int phistCount = pDir->GetList()->GetEntries() - 1;
    // find peak
    int ipeak = 0;
    double dmax = 0;
    for (unsigned idigi = 0; idigi < thePulse.digi.size(); ++idigi)
    {
      double di = -1.0 * (thePulse.digi[idigi]);
      if (di > dmax)
      {
        dmax = di;
        ipeak = idigi;
      }
    }
    int pOffSet = hPulseOne->GetNbinsX() / 2;
    bool isolated = false;
    int isoTime = 100;
    int timeCut = 6000;
<<<<<<< HEAD
    double qcut = 0.0564 * 3.0; //0.2; from fit to noise peak
=======
    double qcut = 0.2;
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
    if (ip == 0 && hitStartTime[ip + 1] - hitEndTime[ip] > isoTime)
      isolated = true;
    if (ip < hitStartTime.size() - 1 && hitStartTime[ip] - hitEndTime[ip - 1] > isoTime && hitStartTime[ip + 1] - hitEndTime[ip] > isoTime)
      isolated = true;
    if (ip == hitStartTime.size() - 1 && hitStartTime[ip] - hitEndTime[ip - 1])
      isolated = true;

    //printf(" pulse %i qp %f hist %i start %i q %f \n", ip, dmax, phistCount, thePulse.istart, thePulse.q * -1.0E9);
    if (phistCount < 500 && thePulse.istart > timeCut && thePulse.q * 1E9 > qcut && isolated)
    {
      TString hpname;
      hpname.Form("PulseOne-%lld-%i-%ld-%i", ientry, ip, thePulse.digi.size(), int(thePulse.q * 1E12));
      TString htitle;
      htitle.Form("event %lld hit# %i length %ld q %i", ientry, ip, thePulse.digi.size(), int(thePulse.q * 1E12));
      pDir->cd();
      TH1D *hPOne = (TH1D *)hPulseOne->Clone(hpname);
      hPOne->Reset();
      hPOne->SetTitle(htitle);
      outfile->cd();
<<<<<<< HEAD
      double fitMin = pOffSet - pLength - thePulse.pwidth * 0.5E9;
      double fitMax = pOffSet - thePulse.pwidth * 0.5E9;
      //printf(" %i width %f fit from %f to %f \n", ip, thePulse.pwidth * 1E9, fitMin, fitMax);
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
      for (unsigned idigi = 0; idigi < thePulse.digi.size(); ++idigi)
      {
        int pbin = idigi - ipeak + pOffSet;
        //printf("\t pulse %i peak %i qp %f digi %i pbin %i q %f \n", ip, ipeak, dmax, idigi, pbin, thePulse.digi[idigi] * -1.0);
        hPOne->SetBinContent(pbin, thePulse.digi[idigi]);
      }
<<<<<<< HEAD
      fPulseBaseline->SetParameters(0., 0.);
      hPOne->Fit(fPulseBaseline, "QOLM", "", fitMin, fitMax); //Q for quiet O do not draw L like W minuit
      if (fPulseBaseline->IsValid())
      {
        thePulse.baseline = fPulseBaseline->GetParameter(0);
        thePulse.baselineErr = fPulseBaseline->GetParError(0);
        thePulse.slope = fPulseBaseline->GetParameter(1);
        thePulse.slopeErr = fPulseBaseline->GetParError(1);
        thePulse.baseChi = fPulseBaseline->GetChisquare() / double(fPulseBaseline->GetNDF());
      }
      if (!fPulseBaseline->IsValid() || thePulse.baseChi > 10)
        printf(" \t !!! pulse fit fails   %lld ip const %lu min %f max %f base %f chi %f \n", ientry, ip, thePulse.digi.size(), fitMin, fitMax, thePulse.baseline, thePulse.baseChi);
      /* subtract fit
      if (fPulseBaseline->IsValid())
      {
        for (unsigned idigi = 0; idigi < thePulse.digi.size(); ++idigi)
        {
          int pbin = idigi - ipeak + pOffSet;
          //printf("\t pulse %i peak %i qp %f digi %i pbin %i q %f \n", ip, ipeak, dmax, idigi, pbin, thePulse.digi[idigi] * -1.0);
          hPOne->SetBinContent(pbin, hPOne->GetBinContent(pbin) - thePulse.baseline);
        }
      }*/
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
    }

    //printf(" ,,,,, %s ipeak %i size %u  (%i,%i) q %E peak %E sum %E \n", hpname.Data(), ip, thePulse.digi.size(), thePulse.istart, thePulse.iend, charge, thePulse.digi[0], hPOne->Integral(0, hPOne->GetNbinsX()));
    /** fill summed pulse **/
    hPulseCharge->Fill(thePulse.q * 1E9);
    if (thePulse.istart > timeCut && thePulse.q * 1E9 > qcut && isolated)
    {
      ++pulseShapeNorm[jpmt];

      for (unsigned idigi = 0; idigi < thePulse.digi.size(); ++idigi)
      {
        int pbin = idigi - ipeak + pOffSet;
<<<<<<< HEAD
        hPulseSum[jpmt]->SetBinContent(pbin, hPulseSum[jpmt]->GetBinContent(pbin) + thePulse.digi[idigi] - thePulse.baseline);
        hPulseSumUn[jpmt]->SetBinContent(pbin, hPulseSumUn[jpmt]->GetBinContent(pbin) + thePulse.digi[idigi]);
      }
    } // end do pulse
    if (thePulse.istart > 1000)
    {
      ++pulseFirstShapeNorm[jpmt];

      for (unsigned idigi = 0; idigi < thePulse.digi.size(); ++idigi)
      {
        int pbin = idigi - ipeak + pOffSet;
        double val = thePulse.digi[idigi];
        hPulseFirstSum[jpmt]->SetBinContent(pbin, hPulseFirstSum[jpmt]->GetBinContent(pbin) + val);
      }
    } // end do pulse

=======
        double val = thePulse.digi[idigi];
        hPulseSum[jpmt]->SetBinContent(pbin, hPulseSum[jpmt]->GetBinContent(pbin) + val);
      }
    } // end do pulse
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
    // add the pulse
    bEvent->hits.push_back(thePulse);
  } // end pulse loop
  tbTree->Fill();
}
void anaPulses::anaEntry(Long64_t ientry)
{
  processedEvent->event = ientry;
  signal.clear();
  signal.resize(1);
  derivative.clear();
  derivative.resize(signal.size());
  integral.resize(derivative.size());
  signal[0] = pmtEvent->volt1;

  deltaT = pmtEvent->time[1] - pmtEvent->time[0];
  if (ientry == 0)
    cout << ientry << " time bin size is " << deltaT << endl;

  //Setting up PMT Run class
  //Class contains all data needs for light analysis
  //
  pmtRun = new TPmtRun();
  pmtRun->nPMTs = signal.size();
  pmtRun->entry = ientry;
  pmtRun->deltaT = deltaT;
  pmtRun->resize();
  //if(ientry > 99) break;
  peakFindingDebug = false;
  /*
    //TCanvas * c0 = new TCanvas("Signal","Signal");
    //c0->cd();
    Int_t skipBin = 21;//4212;
    if( ientry == skipBin+1){
      if(eventInfo != NULL)
        eventInfo->Write();
      break;
    }
    peakFindingDebug = true;
    NHistograms = skipBin;
    //if(skipBin < ientry ) break;
    if( skipBin != ientry) continue;
    hSum[0] = new TH1D(TString("Sum_")+to_string(0),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
    hSumFresh[0] = new TH1D(TString("SumVirgin_")+to_string(0),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
    */

  Int_t nInt = 1;
  Double_t sum = 0;
  if (peakFindingDebug)
    printf("creating event histograms event %lld \n", ientry);
  TString hName;
  for (int j = 0; j < signal.size(); j++)
  {
    evDir->cd();
    hName.Form("Signal_%i_EV_%lld", j, ientry);
    hSignal[j] = new TH1D(hName, hName, nSamples, pmtEvent->time[0], pmtEvent->time[nSamples - 1]);
    hName.Form("Signalb_%i_EV_%lld", j, ientry);
    hSignalb[j] = new TH1D(hName, hName, nSamples, pmtEvent->time[0], pmtEvent->time[nSamples - 1]);
    outfile->cd();
    if (ientry == 0)
    {
      hSum[j] = new TH1D(TString("SumModified_") + to_string(j), "", nSamples, pmtEvent->time[0], pmtEvent->time[nSamples - 1]);
      hSumFresh[j] = new TH1D(TString("SumVirgin_") + to_string(j), "", nSamples, pmtEvent->time[0], pmtEvent->time[nSamples - 1]);
    }
    hPeakFinding[j] = new TH1D(TString("PeakFinding_") + to_string(j) + TString("_Ev_") + to_string(ientry), TString("PeakFinding") + to_string(ientry), nSamples, pmtEvent->time[0], pmtEvent->time[nSamples - 1]);
    //low pass filter
    //signal[j] = LowPassFrequency(signal[j], 20e6 /*MHz*/,deltaT);
    //Trap filter
    //signal[j] = TrapFilter(integral,40,5,j,ientry);
    //DownSampler does not work
    //signal[j] = DownSampler(signal[j],5);
    //Notch Filter
    //Double_t notchFreq = 5e-9;
    //Double_t BW = 0.05*deltaT;
    //signal[j] = NotchFilter(signal[j],notchFreq,BW,deltaT);

    derivative[j] = Derivative(signal[j], nDer);
    hDerivative[j] = new TH1D(TString("derivative_") + to_string(j) + TString("_Ev_") + to_string(ientry), "", nSamples, pmtEvent->time[0], pmtEvent->time[nSamples - 1]);
  }
  sum = 0;

  //cout << "Filling event histograms " << hSignal[0]->GetName() << endl;
  for (int j = 0; j < derivative.size(); j++)
  {
    Double_t baseline = CalculateMean(signal[j], 0, signal[j].size() * 0.05); //500);
    for (int i = 0; i < derivative[j].size(); i++)
    {
      //signal[j][i] /= (double) nInt;
      sum += signal[j][i];
      signal[j][i] -= baseline;
      hSignal[j]->SetBinContent(i + 1, signal[j][i]);
      hWaveSum[j]->SetBinContent(i + 1, signal[j][i] + hWaveSum[j]->GetBinContent(i + 1));
      hPeakFinding[j]->SetBinContent(i + 1, 0);
      hDerivative[j]->SetBinContent(i + 1, derivative[j][i]);
      //hSumWaveform[j]->SetBinContent(i+1,hSumWaveform[j]->GetBinContent(i+1)- signal[j][i]);
    }
  }

  if (peakFindingDebug)
    hSignal[0]->Draw();
  //cout<<"setting derivative"<<endl;
  //Peak Finding
  std::vector<std::vector<Int_t>> peakTime;
  peakTime.resize(signal.size());
  std::vector<std::vector<Double_t>> peakTimeDoubleStart;
  peakTimeDoubleStart.resize(signal.size());
  std::vector<std::vector<Double_t>> peakVoltsStart;
  peakVoltsStart.resize(signal.size());
  std::vector<std::vector<Double_t>> peakTimeDoubleStop;
  peakTimeDoubleStop.resize(signal.size());
  std::vector<std::vector<Double_t>> peakVoltsStop;
  peakVoltsStop.resize(signal.size());
  Int_t maxPeakWidth = 0;
  if (peakFindingDebug)
    cout << "starting peak finding" << endl;
  for (int i = 0; i < signal.size(); i++)
  {
    //cout<<"peakFinding"<<" derSize "<<derivative[i].size()<<" "<<ientry<<endl;
    std::vector<Int_t> pTime = PeakFinding(signal[i], derivative[i], i, ientry);
    if (pTime.size() == 0)
      continue;
    //cout<<"Culling the small and weak"<<endl;
    //Cull the small pulses
    for (int j = 0; j < pTime.size() - 1; j += 2)
    {
      Int_t startBin = pTime[j], stopBin = pTime[j + 1];
      Double_t peakWidth = pmtEvent->time[stopBin] - pmtEvent->time[startBin];
      if ((stopBin - startBin) > maxPeakWidth)
        maxPeakWidth = (stopBin - startBin);

      if (peakFindingDebug)
        cout << "start " << startBin << " " << pmtEvent->time[startBin] << ", stop " << stopBin << " " << pmtEvent->time[stopBin] << endl;
      //cout<<"\t peakWidth "<<peakWidth<<endl;
      if (peakWidth > minPeakWidth && peakWidth < maxPeakWidth)
      {
        peakTime[i].push_back(pTime[j]);
        peakVoltsStart[i].push_back(signal[i][pTime[j]]);
        peakTimeDoubleStart[i].push_back(pmtEvent->time[startBin]);
        //cout<<signal[i][pTime[j]]<<" "<< pTime[j]<<endl;
        peakTime[i].push_back(pTime[j + 1]);
        peakVoltsStop[i].push_back(signal[i][pTime[j + 1]]);
        peakTimeDoubleStop[i].push_back(pmtEvent->time[stopBin]);
        //cout<<"\t"<<signal[i][pTime[j+1]]<<" "<< pTime[j+1]<<endl;
      }
    }
  }
  TGraph *tFindingStart = new TGraph(peakTime[0].size(), &(peakTimeDoubleStart[0][0]), &(peakVoltsStart[0][0]));
  if (peakFindingDebug)
  {
    tFindingStart->SetMarkerColor(3);
    tFindingStart->Draw("same*");
  }
  TGraph *tFindingStop = new TGraph(peakTime[0].size(), &(peakTimeDoubleStop[0][0]), &(peakVoltsStop[0][0]));
  if (peakFindingDebug)
  {
    tFindingStop->SetMarkerStyle(2);
    tFindingStop->SetMarkerColor(2);
    tFindingStop->Draw("same*");
  }
  if (peakFindingDebug)
    cout << "finished peakFinding" << endl;
  //2N + 1 moving window for WMA
  Int_t movingWindow;
  if (signal[0].size() == 10000)
    movingWindow = 2000;
  else if (signal[0].size() == 500)
    movingWindow = 100;
  //if(signal[0].size() == 500)
  //movingWindow = 10;

  std::vector<std::vector<Double_t>> BaSeLiNe;
  BaSeLiNe.resize(signal.size());
  Double_t waveformFraction = 0.60;
  //if(dsNum == 1) for(int j = 0; j < signal.size();j++) BaSeLiNe[j] = BaselineWMA(signal[j],peakTime[j],ientry,movingWindow,j);
  for (int j = 0; j < signal.size(); j++)
    BaSeLiNe[j] = BaselineWMA(signal[j], peakTime[j], ientry, movingWindow, j, waveformFraction);
  //BaSeLiNe[0].resize(0);
  //BaSeLiNe[0].resize(signal[0].size());
  if (peakFindingDebug)
    cout << "finished Baseline calculations" << endl;

  ///*
  Double_t windowSum = 0;
  for (int j = 0; j < signal.size(); j++)
  {
    Int_t startWindowBin = hSignal[j]->FindBin(50e-9);
    Int_t stopWindowBin = hSignal[j]->FindBin(100e-9);
    Double_t baseline = CalculateMean(signal[j], 0, signal[j].size() * 0.05); //500);
    for (int k = 0; k < signal[j].size(); k++)
    {
      //if(dsNum == 2) BaSeLiNe[j].push_back(baseline);
      if (BaSeLiNe[j].size() > 0)
      {
        //hSignal[j]->SetBinContent(k+1,signal[j][k]-BaSeLiNe[j][k]);
        hSignalb[j]->SetBinContent(k + 1, signal[j][k]);
        hSumWaveform[j]->SetBinContent(k + 1, hSumWaveform[j]->GetBinContent(k + 1) - signal[j][k] + BaSeLiNe[j][k]);
        //hSignal[j]->SetBinContent(k+1,signal[j][k]);
        if (j == 0)
          processedEvent->volt1.push_back(signal[j][k] - BaSeLiNe[j][k]);
        if (j == 1)
          processedEvent->volt2.push_back(signal[j][k] - BaSeLiNe[j][k]);
      }
      else
      {
        hSignalb[j]->SetBinContent(k + 1, signal[j][k]);
        if (j == 0)
          processedEvent->volt1.push_back(signal[j][k] - baseline);
        if (j == 1)
          processedEvent->volt2.push_back(signal[j][k] - baseline);
      }
      processedEvent->time.push_back(pmtEvent->time[k]);
      if (k >= startWindowBin && k < stopWindowBin)
      {
        if (BaSeLiNe[j].size() > 0)
          windowSum += (signal[j][k] - BaSeLiNe[j][k]) * deltaT;
        else
          windowSum += (signal[j][k] - baseline) * deltaT;
      }
    }
    hWindowSum[j]->Fill(-windowSum);
  }
  if (peakFindingDebug)
    cout << "finished windowSum" << endl;
  //*/

  /*
    for(int i = 0; i < signal[0].size();i++){
      if(peakTime[0].empty()) continue;
      hSum[0]->SetBinContent(i+1,hSum[0]->GetBinContent(i+1)-signal[0][i]+BaSeLiNe[i]);
      hSumFresh[0]->SetBinContent(i+1,hSumFresh[0]->GetBinContent(i+1)-signal[0][i]);
    }
    */
  if (peakFindingDebug)
    cout << "Filling Ntuple" << endl;
  //filling ntuple
  for (int j = 0; j < signal.size(); j++)
  {
    Double_t vMaxEvent = 0, tMaxEvent = 0, cMaxEvent = 0;
    TString sTitle = TString("");
    Double_t baseline = CalculateMean(signal[j], 0, signal[j].size() * 0.05); //500);
    Double_t Sdev = CalculateSdev(signal[j], 0, signal[j].size() * 0.05, baseline);

    if (peakTime[j].empty())
      continue;
    Double_t T0 = 0, V0 = 0, C0 = 0;
    Double_t nFound = 0, nMiss = 0, nNoise = 0, totalCharge = 0;
    if (peakFindingDebug)
      cout << "entry " << ientry << endl;
    for (int i = 0; i < peakTime[j].size() - 1; i += 2)
    {
      Double_t vMax = 0, charge = 0, vMaxTime = 0, cMax = 0, startTime, stopTime, peakWidth;
      Int_t startBin = peakTime[j][i], stopBin = peakTime[j][i + 1];
      //Int_t startBin = hSignal[j]->FindBin(400e-9);
      //Int_t stopBin  = hSignal[j]->FindBin(500e-9);
      startTime = pmtEvent->time[startBin];
      stopTime = pmtEvent->time[stopBin];
      peakWidth = stopTime - startTime;
      if (peakFindingDebug)
        cout << "Looping over pulse start/stop " << startTime << "/" << stopTime << endl;
      for (int k = startBin; k < stopBin; k++)
      {
        //if(signal[j][k]-BaSeLiNe[k] > 0) continue;
        charge += signal[j][k] - BaSeLiNe[j][k];
        hSum[0]->SetBinContent(k + 1, hSum[0]->GetBinContent(k + 1) - signal[j][k] + BaSeLiNe[j][i]);
        hSumFresh[0]->SetBinContent(k + 1, hSumFresh[0]->GetBinContent(k + 1) - signal[j][k]);
        if (std::fabs(signal[j][k] - BaSeLiNe[j][k]) > std::fabs(vMax))
        {
          vMax = signal[j][k] - BaSeLiNe[j][k];
          vMaxTime = pmtEvent->time[k];
          cMax = charge;
        }
        //hPeakFinding[j]->SetBinContent(k+1,signal[j][k]-BaSeLiNe[j][k]);
        hPeakFinding[j]->SetBinContent(k + 1, signal[j][k]);
      }
      if (peakFindingDebug)
        cout << "Finished charge calculation" << endl;
      if (i == 0)
      {
        T0 = startTime;
        V0 = vMax;
        C0 = -charge * deltaT;
      }
      ntuplePulse->Fill(ientry, j, peakTime[j].size() / 2, -charge * deltaT, startTime, peakWidth, T0, V0, C0, -vMax, vMaxTime);

      pmtRun->charge[j].push_back(-charge * deltaT);
      pmtRun->startTimes[j].push_back(startTime);
      pmtRun->peakWidths[j].push_back(peakWidth);
      pmtRun->peakHeights[j].push_back(-vMax);
      pmtRun->peakTimes[j].push_back(vMaxTime);

      //sTitle += TString("Charge_")+to_string(-1e9*deltaT*charge)+TString("_start/stop_")+to_string(1.e6*startTime)+TString("/")+to_string(1.e6*stopTime)+TString("_vMax_")+to_string(-vMax)+TString("_");
      if (std::fabs(vMax) > std::fabs(vMaxEvent))
      {
        vMaxEvent = vMax;
        tMaxEvent = vMaxTime;
        cMaxEvent = cMax * deltaT;
      }

      bool doPulse = true;
      int pmtNum = j;
      std::vector<double> digi = pmtEvent->volt1;
      int ndigi = digi.size();
<<<<<<< HEAD
      /*** do deconvolution here ****/
      std::vector<Double_t> deconWave = decon(signal[j]);
      for (int i = 0; i < deconWave.size(); i++)
      {
        hWaveRawOne[j]->SetBinContent(i + 1, signal[j][i]);
        hWaveOne[j]->SetBinContent(i + 1, deconWave[i]);
      }
=======
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
    }
    //end of pmtNum loop
    sTitle += TString("tMax_") + to_string(tMaxEvent * 1e6) + TString("_vMax_") + to_string(vMaxEvent) + TString("_totalCharge_") + to_string(totalCharge);

    if (peakFindingDebug)
      cout << "Filling data ntupleEvent" << endl;
    ntupleEvent->Fill(0, ientry, j, Sdev, baseline, -sum, deltaT, tMaxEvent, -vMaxEvent, -cMaxEvent, -totalCharge * deltaT, 0, peakTime[j].size() / 2.); //-vMaxEvent,tMaxEvent);
    if (peakFindingDebug)
      cout << "Filling pmtRun class" << endl;

    pmtRun->T0[j] = T0;
    pmtRun->totalCharge[j] = -totalCharge * deltaT;
    pmtRun->tMax[j] = tMaxEvent;
    pmtRun->vMax[j] = -vMaxEvent;
    pmtRun->cMax[j] = -cMaxEvent;
    pmtRun->baseline[j] = baseline;
    pmtRun->sDev[j] = Sdev;
    //for(int i = 0; i < simNoise[j].size() && z == 0; i++) sTitle +=TString("_") + to_string(1e6*simNoise[j][i]);
    //for(int i = 0; i < simNoise[j].size() && isSim; i++) sTitle +=TString("_") + to_string(1e6*simNoise[j][i]);
    hSignal[j]->SetTitle(sTitle);
  }
  if (peakFindingDebug)
    cout << "finished filling ntuple" << endl;
  hPeakFinding[0]->SetLineColor(2);
  if (signal.size() == 2)
    hPeakFinding[1]->SetLineColor(2);
  if (peakFindingDebug)
    hPeakFinding[1]->Draw("same");
  //###################//
  //End of Event Clean-Up//
  //###################//
  if (peakFindingDebug)
  {
    cout << "End of Event Clean-Up" << endl;
  }

  processedTree->Fill();
  //pmtRun->Write();
  processedEvent->clear();
<<<<<<< HEAD
=======
  if (ientry > NHistograms)
  {
    for (int i = 0; i < signal.size(); i++)
    {
      delete hSignal[i];
      delete hDerivative[i];
      delete hPeakFinding[i];
    }
  }
>>>>>>> f0dfb5122b6d53a2719b5e05c077401073376aeb
}

std::vector<Double_t> anaPulses::DownSampler(std::vector<Double_t> sig, Int_t N)
{

  if (N % 2 == 0)
    N + 1;
  std::vector<Double_t> downSample;
  cout << "starting down sampler " << endl;
  for (int i = 0; i < sig.size();)
  {
    Int_t arr[] = {N, i, (Int_t)sig.size() - 1 - i};
    Int_t window = *std::min_element(arr, arr + 3);
    Double_t val = 0;
    for (int j = 1; j < window; j++)
    {
      val += (sig[i + j] + sig[i - j]) / window;
    }
    downSample.push_back(window);
    if (window == 0)
      i++;
    else
      i += window;
  }
  cout << downSample.size() << endl;
  return downSample;
}

std::vector<Double_t> anaPulses::Derivative(std::vector<Double_t> sig, Int_t N)
{
  std::vector<Double_t> diff;
  //diff.push_back(0);
  Double_t leftSum = 0, rightSum = 0;
  //if(N<2) N = 2;
  for (int i = 0; i < sig.size(); i++)
  {
    Double_t der = 0;
    Int_t arr[] = {N, i, (Int_t)sig.size() - 1 - i};
    Int_t window = *std::min_element(arr, arr + 3);
    if (i <= N && 2 * i <= sig.size() - 1)
    {
      leftSum = leftSum + sig[i - 1];
      rightSum = rightSum - sig[i] + sig[2 * i - 1] + sig[2 * i];
    }
    else if (i > N && i + N <= sig.size() - 1)
    {
      leftSum = leftSum + sig[i - 1] - sig[i - 1 - N];
      rightSum = rightSum - sig[i] + sig[i + N];
    }
    else if (i + N > sig.size() - 1 && 2 * i > sig.size() - 1)
    {
      leftSum = leftSum + sig[i - 1] - sig[2 * i - sig.size() - 1] - sig[2 * i - sig.size()];
      rightSum = rightSum - sig[i];
    }
    else
    {
      cout << "you done messed the if else statements" << endl;
    }
    diff.push_back((rightSum - leftSum) / (1 + window));
    /*
    for(int j = 1; j < window ; j++){
      der += (sig[i+j] - sig[i-j]);
    }
    diff.push_back(der);
    */
  }
  return diff;
}

Double_t anaPulses::RMSCalculator(std::vector<Double_t> vec, Int_t ientry)
{
  Double_t rms = 0, sum = 0;
  TH1D *h1 = new TH1D(TString("rmsHistogram_") + to_string(ientry), TString("rmsHistogram_") + to_string(ientry), 500, -BubbleSort(vec)[vec.size() - 1], BubbleSort(vec)[vec.size() - 1]);
  for (int i = 0; i < vec.size(); i++)
  {
    h1->Fill(vec[i]);
  }
  Int_t centerBin = h1->GetNbinsX() / 2;
  h1->SetBinContent(centerBin, sqrt(h1->GetBinContent(centerBin) * (h1->GetBinContent(centerBin + 1) + h1->GetBinContent(centerBin - 1)) / 2.));
  for (int i = 0; i < h1->GetNbinsX(); i++)
  {
    h1->SetBinContent(i + 1, TMath::Exp(h1->GetBinContent(i + 1) / vec.size()) - 1);
  }
  TF1 *fGaus = new TF1("fit", "gaus(0)");
  h1->Fit(fGaus, "Q");
  rms = fGaus->GetParameter(2);
  if (ientry > NHistograms)
    delete h1;
  return rms;
}
//Weighted Moving Average
std::vector<Double_t> anaPulses::BaselineWMA(std::vector<Double_t> sig, std::vector<Int_t> peaks, Int_t ientry, Int_t N, Int_t pmtNum, Double_t waveformFraction)
{
  std::vector<Double_t> Baseline;
  if (peaks.size() == 0)
  {
    return Baseline;
  }
  /*
  //subsection of the signal for large events
  std::vector<Double_t> tempSig;
  if(waveformFraction != 0){
    //copy(sig.at(int(sig.size()*waveformFraction)),sig.end(),back_inserter(tempSig));
    //tempSig.assign(sig.at(int(sig.size()*waveformFraction)),sig.end());
    int spot = sig.size()*waveformFraction;
    tempSig.assign(sig.begin() + spot,sig.end());
    sig.clear();
    sig.resize(tempSig.size());
    sig = tempSig;
  }
  */
  Int_t ipeak = 1;
  Int_t startPeak = peaks[0];
  Int_t stopPeak = -1; //peaks[ipeak+1];
  std::vector<Double_t> weight;
  Int_t windowCounter = 0;
  //set peak weights
  for (int i = 0; i < sig.size(); i++)
  {
    if (i <= startPeak && i >= stopPeak)
    {
      weight.push_back(startPeak - stopPeak - 1);
      if (windowCounter > N)
        N = windowCounter;
      windowCounter = 0;
    }
    else
    {
      weight.push_back(1.e-6);
      windowCounter++;
    }
    if (i == peaks[ipeak])
    {
      ipeak += 2;
      if (ipeak - 1 == peaks.size())
      {
        startPeak = sig.size();
      }
      else
        startPeak = peaks[ipeak - 1];
      stopPeak = peaks[ipeak - 2];
    }
    //cout<<"Bin "<<i<<", startPeak "<<startPeak<<", stopPeak "<<stopPeak<<", weight "<<weight[i]<<", window "<<window[i]<<", ipeak "<<ipeak<<endl;
  }
  //if(N < 10)
  //N*=1.5;
  //  N += 10;
  //cout<<N<<endl;
  TH1D *f1 = new TH1D(TString("WMA_") + to_string(pmtNum) + TString("_") + to_string(ientry), TString("WMA") + to_string(ientry) + TString("_N_") + to_string(N), sig.size(), 0, pmtEvent->time[nSamples - 1]);
  /*
  //calculate WMA
  for(int i = 0; i < sig.size(); i++){
    Int_t hiWindow = std::min(i+N,(Int_t)sig.size() -1);
    Int_t loWindow = std::max(0,i-N);
    Double_t top = 0, bot = 0;
    f2->SetBinContent(i+1,weight[i]);
    //cout<<"loWindow "<<loWindow<<", hiWindow "<<hiWindow<<endl;
    for(int j = loWindow;j<hiWindow;j++){
      top +=sig[j]*weight[j]*(1+cos((j-i)*pi/N));
      bot +=weight[j]*(1+cos((j-i)*pi/N));
      //cout<<"\t top "<<top<<", bot "<<bot<<endl;
    }
    Baseline.push_back(top/bot);
    f1->SetBinContent(i+1,top/bot);
  }
  */
  ///*
  //obnoxious way of making things faster
  //i.e. run times does not depend on how big the window, N, is.
  //It only depends on total number of points
  Double_t K_sw = 0;
  Double_t C_sw = 0;
  Double_t S_sw = 0;
  Double_t K_w = 0;
  Double_t C_w = 0;
  Double_t S_w = 0;
  Double_t kc = cos(pi / N);
  Double_t ks = sin(pi / N);
  Int_t resetBin = 100, lowBinMax = 15;
  Double_t lowBinWeight = N;
  if (sig.size() == 500)
    resetBin = 20;
  for (int i = 0; i < sig.size(); i++)
  {
    Int_t hiWindow = std::min(i + N, (Int_t)sig.size() - 1);
    Int_t loWindow = std::max(0, i - N);
    if (weight[i] <= lowBinMax && weight[i] >= 1 && sig[i] > 0)
      weight[i] *= lowBinWeight;
    if (i % resetBin == 0)
    {
      K_sw = 0;
      C_sw = 0;
      S_sw = 0;
      K_w = 0;
      C_w = 0;
      S_w = 0;
      for (int j = loWindow; j < hiWindow; j++)
      {
        if (weight[j] <= lowBinMax && weight[j] >= 1 && sig[j] > 0)
          weight[j] *= lowBinWeight;
        K_sw += sig[j] * weight[j];
        C_sw += sig[j] * weight[j] * cos((j - i) * pi / N);
        S_sw += sig[j] * weight[j] * sin((j - i) * pi / N);
        K_w += weight[j];
        C_w += weight[j] * cos((j - i) * pi / N);
        S_w += weight[j] * sin((j - i) * pi / N);
      }
    }
    else
    {
      if (i - N <= 0 && i + N <= sig.size() - 1)
      {
        K_sw = K_sw + sig[i + N] * weight[i + N];
        Double_t tempC_sw = C_sw;
        C_sw = kc * C_sw + ks * S_sw - sig[i + N] * weight[i + N];
        S_sw = kc * S_sw - ks * tempC_sw;

        K_w = K_w + weight[i + N];
        Double_t tempC_w = C_w;
        C_w = kc * C_w + ks * S_w - weight[i + N];
        S_w = kc * S_w - ks * tempC_w;
      }
      else if (i - N > 0 && i + N <= sig.size() - 1)
      {
        K_sw = K_sw - sig[i - 1 - N] * weight[i - 1 - N] + sig[i + N] * weight[i + N];
        Double_t tempC_sw = C_sw;
        C_sw = kc * C_sw + ks * S_sw + kc * sig[i - 1 - N] * weight[i - 1 - N] - sig[i + N] * weight[i + N];
        S_sw = kc * S_sw - ks * tempC_sw - ks * sig[i - 1 - N] * weight[i - 1 - N];

        K_w = K_w - weight[i - 1 - N] + weight[i + N];
        Double_t tempC_w = C_w;
        C_w = kc * C_w + ks * S_w + kc * weight[i - 1 - N] - weight[i + N];
        S_w = kc * S_w - ks * tempC_w - ks * weight[i - 1 - N];
      }
      else if (i - N > 0 && i + N > sig.size() - 1)
      {
        K_sw = K_sw - sig[i - 1 - N] * weight[i - 1 - N];
        Double_t tempC_sw = C_sw;
        C_sw = kc * C_sw + ks * S_sw + kc * sig[i - 1 - N] * weight[i - 1 - N];
        S_sw = kc * S_sw - ks * tempC_sw - ks * sig[i - 1 - N] * weight[i - 1 - N];

        K_w = K_w - weight[i - 1 - N];
        Double_t tempC_w = C_w;
        C_w = kc * C_w + ks * S_w + kc * weight[i - 1 - N];
        S_w = kc * S_w - ks * tempC_w - ks * weight[i - 1 - N];
      }
      else if (i - N <= 0 && i + N > sig.size() - 1)
      {
        K_sw = K_sw;
        Double_t tempC_sw = C_sw;
        C_sw = kc * C_sw + ks * S_sw;
        S_sw = kc * S_sw - ks * tempC_sw;

        K_w = K_w - weight[i - 1 - N];
        Double_t tempC_w = C_w;
        C_w = kc * C_w + ks * S_w;
        S_w = kc * S_w - ks * tempC_w;
      }
      else
      {
        cout << "you messed up the logic on the moving baseline" << endl;
      }
    }
    //cout<<"ientry"<<ientry<<" "<<i<<"...K_sw "<<K_sw<<", C_sw "<<C_sw<<", S_sw "<<S_sw<<", K_w "<<K_w<<", C_w "<<C_w<<"..."<< (K_sw+C_sw)/(K_w+C_w)<<endl;
    f1->SetBinContent(i + 1, (K_sw + C_sw) / (K_w + C_w));
    Baseline.push_back((K_sw + C_sw) / (K_w + C_w));
  }
  //*/
  f1->SetLineColor(6);
  if (peakFindingDebug && pmtNum == 1)
    f1->Draw("same");
  //f1->Write();
  if (ientry > NHistograms)
  {
    delete f1;
  }

  return Baseline;
}

std::vector<Double_t> anaPulses::Integral(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry)
{
  std::vector<Double_t> integ;
  Double_t sum = 0;
  for (int i = 0; i < sig.size(); i++)
  {
    sum += sig[i];
    integ.push_back(sum);
  }
  return integ;
}
std::vector<Int_t> anaPulses::PeakFinding(std::vector<Double_t> sig, std::vector<Double_t> diff, Int_t pmtNum, Int_t ientry)
{
  std::vector<Int_t> peakTime;
  std::vector<Int_t> peakTimeTrimmed;
  //-1 ll, 0, llhh, 1 hh
  std::vector<Int_t> peakFlag;

  Int_t rmsBin = sig.size() * 0.05; //500;
  //Double_t rmsFit = RMSCalculator(diff,ientry);
  Double_t rmsMean = CalculateMean(diff);
  Double_t rmsSdevlo = CalculateSdev(diff, 0, rmsBin, rmsMean);
  Double_t rmsSdevhi = CalculateSdev(diff, diff.size() - rmsBin, diff.size() - 1, rmsMean);
  Double_t rmsSdevAll = CalculateSdev(diff, 0, diff.size() - 1, rmsMean);
  Double_t arr[] = {rmsSdevlo, rmsSdevhi, rmsSdevAll};
  Double_t rmsSdev = *std::min_element(arr, arr + 3);
  //Double_t rmsSdev = CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  Double_t rms = rmsSdev;     // min(rmsFit,rmsSdev);
  Double_t loRMS = 3.5 * rms; //CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  //if midRMS is larger than hiRMS than the peakfind breaks
  Double_t midRMS = 2.5 * rms;
  bool midRMSFlag = false;
  Double_t hiRMS = 3.5 * rms; //CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  //cout<<"RMS = "<<rms<<endl;
  hDerivative[pmtNum]->SetTitle(TString("loThreshold_") + to_string(loRMS) + TString("_hiThreshold_") + to_string(hiRMS) + TString("_midThresh_") + to_string(midRMS));
  hDerivative[pmtNum]->SetLineColor(8);
  Bool_t low1 = false, low2 = false, mid1 = false, mid2 = false, high1 = false, high2 = false, ll = false, hh = false, mm = false;

  Int_t LowStartBin = -1, LowStopBin = -1, MidStartBin = -1, MidStopBin = -1, HighStartBin = -1, HighStopBin = -1;
  Int_t zeroCrossing = 0;
  for (int i = 1; i < diff.size(); i++)
  {

    //sometimes noise will trigger a lolo then a long time later a hihi will
    //be found resulting in a long noise peak that forms a lohi
    //if baseline crosses zero more than onces between lo hi, don't count it
    //resets if a hi && lo were found
    if (((diff[i] < 0 && diff[i - 1] > 0) || (diff[i] > 0 && diff[i - 1] < 0)))
    {
      zeroCrossing++;
    }

    if (diff[i] >= midRMS && diff[i - 1] < midRMS)
    {
      MidStartBin = i;
      mid1 = true;
    }
    if (diff[i] <= midRMS && diff[i - 1] > midRMS)
    {
      MidStopBin = i;
      mid2 = true;
    }
    if (diff[i] >= hiRMS && diff[i - 1] < hiRMS)
    {
      HighStartBin = i;
      high1 = true;
    }
    if (diff[i] <= hiRMS && diff[i - 1] > hiRMS)
    {
      HighStopBin = i;
      high2 = true;
    }

    if (high1 && high2)
    {
      if (peakFindingDebug)
      {
        cout << "hi " << pmtEvent->time[HighStartBin] << " " << pmtEvent->time[HighStopBin] << " " << diff[i] << " " << rms << endl;
        cout << "\tlow " << pmtEvent->time[LowStartBin] << " " << pmtEvent->time[LowStopBin] << endl;
        cout << "\tmid " << pmtEvent->time[MidStartBin] << " " << pmtEvent->time[MidStopBin] << endl;
      }
      peakTime.push_back(HighStartBin);
      peakTime.push_back(HighStopBin);
      peakFlag.push_back(1);
      high1 = false;
      high2 = false;
      if (ll)
        hh = true;
      else
        hh = false;
    }
    else if (mid1 && mid2 && ll && midRMSFlag)
    {
      if (peakFindingDebug)
      {
        cout << "mid " << pmtEvent->time[MidStartBin] << " " << pmtEvent->time[MidStopBin] << " " << diff[i] << " " << rms << endl;
        cout << "\tlow " << pmtEvent->time[LowStartBin] << " " << pmtEvent->time[LowStopBin] << endl;
        cout << "\thi " << pmtEvent->time[HighStartBin] << " " << pmtEvent->time[HighStopBin] << endl;
      }
      peakTime.push_back(MidStartBin);
      peakTime.push_back(MidStopBin);
      peakFlag.push_back(2);
      mid1 = false;
      mid2 = false;
      if (ll)
        mm = true;
      else
        mm = false;
    }

    if (ll && hh)
    {
      if (zeroCrossing > 1)
      {
        zeroCrossing = 0;
        ll = false;
        hh = false;
        low1 = false;
        low2 = false;
        continue;
      }
      if (peakFindingDebug)
        cout << "llhh " << pmtEvent->time[LowStartBin] << " " << pmtEvent->time[HighStopBin] << " zeroCrossing " << zeroCrossing << endl;
      peakFlag.erase(peakFlag.end() - 2, peakFlag.end());
      peakFlag.push_back(0);
      peakTime.erase(peakTime.end() - 4, peakTime.end());
      peakTime.push_back(LowStartBin);
      peakTime.push_back(HighStopBin);
      low1 = false;
      low2 = false;
      high1 = false;
      high2 = false;
      ll = false;
      hh = false;
    }
    else if (ll && mm)
    {
      if (zeroCrossing > 1)
      {
        zeroCrossing = 0;
        ll = false;
        mm = false;
        low1 = false;
        low2 = false;
        mid1 = false;
        mid2 = false;
        continue;
      }
      if (peakFindingDebug)
        cout << "llmm " << pmtEvent->time[LowStartBin] << " " << pmtEvent->time[MidStopBin] << " zeroCrossing " << zeroCrossing << endl;
      peakFlag.erase(peakFlag.end() - 2, peakFlag.end());
      peakFlag.push_back(2);
      peakTime.erase(peakTime.end() - 4, peakTime.end());
      peakTime.push_back(LowStartBin);
      peakTime.push_back(MidStopBin);
      low1 = false;
      low2 = false;
      mid1 = false;
      mid2 = false;
      ll = false;
      mm = false;
    }

    //Check low crossings last in order to avoid the next pulse from over lapping with the last pulse,
    //Only a problem when binSize >= 2e-9
    if (diff[i] <= -loRMS && diff[i - 1] > -loRMS)
    {
      LowStartBin = i;
      low1 = true;
    }
    if (diff[i] >= -loRMS && diff[i - 1] < -loRMS)
    {
      LowStopBin = i;
      low2 = true;
    }

    if (low1 && low2)
    {
      if (peakFindingDebug)
      {
        cout << "low " << pmtEvent->time[LowStartBin] << " " << pmtEvent->time[LowStopBin] << " " << diff[i] << " " << -rms << endl;
        cout << "\thi " << pmtEvent->time[HighStartBin] << " " << pmtEvent->time[HighStopBin] << endl;
        cout << "\tmid " << pmtEvent->time[MidStartBin] << " " << pmtEvent->time[MidStopBin] << endl;
      }

      peakTime.push_back(LowStartBin);
      peakTime.push_back(LowStopBin);
      peakFlag.push_back(-1);
      low1 = false;
      low2 = false;
      //Make sure mid peaks are not found before low peaks
      mid1 = false;
      mid2 = false;
      if (hh || mm)
        ll = false;
      else
      {
        ll = true;
        //reset zero crossing if a new lolo peak is found
        zeroCrossing = 0;
      }
    }
  }
  //find local minimum around found peaks
  //also merge peaks
  for (int i = 0; i < peakTime.size(); i += 2)
  {
    Int_t lowBin = peakTime[i];
    Int_t hiBin = peakTime[i + 1];
    Double_t startTime = pmtEvent->time[lowBin];
    Double_t stopTime = pmtEvent->time[hiBin];
    /*
    if( peakFlag[i/2] == 1 && (hiBin-lowBin)*deltaT < 1.5*minPeakWidth){
      if(peakFindingDebug) cout<<"Removing UP-TYPE small pulse start "<<pmtEvent->time[lowBin]<<", stop "<<pmtEvent->time[hiBin]<<endl;
      continue;
    }
    */
    if ((hiBin - lowBin) * deltaT < minPeakWidth)
    {
      if (peakFindingDebug)
        cout << "Removing small pulse start " << pmtEvent->time[lowBin] << ", stop " << pmtEvent->time[hiBin] << endl;
      continue;
    }
    //*/
    Double_t peakThresh = 1 * rms;
    Int_t shiftBin = sig.size() * 1e-3;
    if (shiftBin == 0)
      shiftBin = 1;
    shiftBin = 1;
    ///*
    for (int k = lowBin; k > 1; k--)
    {
      if (((diff[k] < 0 && diff[k - 1] > 0) || (diff[k] > 0 && diff[k - 1] < 0)))
      {
        zeroCrossing++;
      }
      //#-1 && #0 crossing
      if ((peakFlag[i / 2] == 0 || peakFlag[i / 2] == -1) && (diff[k] < -peakThresh && diff[k - 1] > -peakThresh))
      {
        if (peakFindingDebug)
          cout << "extending low bin from " << pmtEvent->time[lowBin] << ", to " << pmtEvent->time[k] << ", peakFlag " << peakFlag[i / 2] << endl;
        if (peakFindingDebug)
          cout << "\tdiff[k] " << diff[k] << ", diff[k-1] " << diff[k - 1] << " peakThresh " << peakThresh << endl;
        if (sig[k] > 0)
          lowBin = k + shiftBin; //10;
        else
          lowBin = k;
        break;
      }
      //#1 && #2 crossing
      else if ((peakFlag[i / 2] >= 1) && ((diff[k] < -peakThresh && diff[k - 1] > -peakThresh)))
      {
        if (peakFindingDebug)
          cout << "extending low bin from " << pmtEvent->time[lowBin] << ", to " << pmtEvent->time[k] << ", peakFlag " << peakFlag[i / 2] << endl;
        if (peakFindingDebug)
          cout << "\tdiff[k] " << diff[k] << ", diff[k-1] " << diff[k - 1] << " peakThresh " << peakThresh << endl;
        if (sig[k] > 0)
          lowBin = k + shiftBin; //10;
        else
          lowBin = k;
        break;
      }
    }
    //keep maximum for baseline fitting
    for (int k = hiBin; k < diff.size() - 1; k++)
    {
      //#4 crossing
      if ((diff[k + 1] < peakThresh && diff[k] > peakThresh))
      {
        if (peakFindingDebug)
          cout << "extending hi bin from " << pmtEvent->time[hiBin] << ", to " << pmtEvent->time[k] << ", peakFlag " << peakFlag[i / 2] << endl;
        if (peakFindingDebug)
          cout << "\tdiff[k+1] " << diff[k + 1] << ", diff[k] " << diff[k] << " peakThresh " << peakThresh << endl;
        if (sig[k] > 0)
          hiBin = k - shiftBin; //10;
        else
          hiBin = k;
        break;
      }
    }
    //*/
    //cull identical pulses
    if (peakTimeTrimmed.size() > 0 && lowBin == peakTimeTrimmed[peakTimeTrimmed.size() - 2] && hiBin == peakTimeTrimmed[peakTimeTrimmed.size() - 1])
    {
      if (peakFindingDebug)
        cout << "Culling identical pulses..(" << pmtEvent->time[lowBin] << "," << pmtEvent->time[hiBin] << ")..and.." << pmtEvent->time[peakTime[i - 2]] << "," << pmtEvent->time[peakTime[i - 1]] << endl;
      continue;
    }
    //cull pulses inside other pulses
    if (peakTimeTrimmed.size() > 0 && lowBin < peakTimeTrimmed[peakTimeTrimmed.size() - 1])
    {
      if (peakFindingDebug)
        cout << "Pulse (" << pmtEvent->time[lowBin] << "," << pmtEvent->time[hiBin] << ") starts inside this pulse (" << pmtEvent->time[peakTimeTrimmed[peakTimeTrimmed.size() - 2]] << "," << pmtEvent->time[peakTimeTrimmed[peakTimeTrimmed.size() - 1]] << ")" << endl;
      //if current found pulse starts at the end of previous pulse,
      //move the start time to where the previous pulse ended
      if (sig[peakTimeTrimmed[peakTimeTrimmed.size() - 1]] > 0)
        lowBin = peakTimeTrimmed[peakTimeTrimmed.size() - 1] + shiftBin;
      else
        lowBin = peakTimeTrimmed[peakTimeTrimmed.size() - 1];
      //continue;
    }

    if (peakFindingDebug)
      cout << "#### Pushing back lowBinTime " << pmtEvent->time[lowBin] << ", highBin " << pmtEvent->time[hiBin] << " ####" << endl;
    //if( peakFlag[i/2] == 2) cout<<"ientry ... "<<ientry<<" Mid peakfinding...lowBinTime "<<pmtEvent->time[lowBin]<<", midBin "<<pmtEvent->time[hiBin]<<" ####"<<endl;

    peakTimeTrimmed.push_back(lowBin);
    peakTimeTrimmed.push_back(hiBin);
  }

  return peakTimeTrimmed;
  //return peakTime;
}

//calc baseline and exclude region where peaks are found
std::pair<Double_t, Double_t> anaPulses::BaselineSubtraction(std::vector<Double_t> sig, std::vector<Int_t> weight, Int_t pmtNum, Int_t ientry)
{
  Int_t startBin = 0, stopBin = 0;

  for (int i = 0; i < weight.size() - 1; i += 2)
  {
    startBin = weight[i];
    stopBin = weight[i + 1];
    for (int k = 0; k < sig.size(); k++)
    {
      if (k >= startBin && k <= stopBin)
        sig.erase(sig.begin() + k);
    }
  }
  //std::sort(sig.begin(),sig.end());

  return std::make_pair(CalculateMean(sig), CalculateSdev(sig, 0, sig.size(), CalculateMean(sig)));
}
//calc moving bassline
//P. Funk, G. Funk ect...
std::vector<Double_t> anaPulses::BaselineSubtraction(std::vector<Double_t> sig, Int_t weight, Int_t pmtNum, Int_t ientry)
{
  std::vector<Double_t> baseline;
  std::vector<Double_t> gaus;
  //weight should be even
  if (weight % 2 == 1)
    weight++;

  Double_t sigma = weight;
  Double_t norm = 0;
  for (int i = -weight; i <= weight; i++)
  {
    Double_t val = exp(-i * i / (2 * sigma * sigma));
    //    gaus.push_back((1/sqrt(sigma*sigma*3.1415*2))*exp(-i*i/(2*sigma*sigma)));
    gaus.push_back(val);
    norm += val;
  }
  std::vector<Double_t> median;
  //for(int i = sig.size() - 1; i >= 0; i--){
  for (int i = 0; i < sig.size(); i++)
  {
    std::vector<Double_t> localMedian;
    for (int j = i - weight; j <= i + weight; j++)
    {
      if (j < 0)
        continue;
      if (j >= sig.size())
        continue;
      localMedian.push_back(sig[j]);
    }
    //sorts vector from lowest to highest
    std::sort(localMedian.begin(), localMedian.end());
    //middle value is the median of the window per Local Median method
    //tuned value to 0.8 to better match long/large peak baselines
    //N. McFadden
    median.push_back(localMedian[(int)localMedian.size() * .8]);
    //median.push_back(localMedian[localMedian.size()*.2]);
  }
  for (int i = 0; i < median.size(); i++)
  {
    //for(int i = sig.size() - 1; i >= 0; i--){
    Double_t val = 0;
    for (int j = i - weight; j <= i + weight; j++)
    {
      if (j < 0)
        continue;
      if (j >= median.size())
        continue;
      val += median[j] * gaus[j + weight - i] / norm;
      //cout<<"\t"<<median[j]<<" "<<gaus[j+weight-i]/norm<<" "<<i<<" "<<j<<" "<<j+weight-i<<endl;;
    }
    baseline.push_back(val);
  }
  //hBaseline[pmtNum] = new TH1D(TString("Baseline_Channel_")+to_string(pmtNum)+TString("_event_")+to_string(ientry),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
  for (int i = 0; i < baseline.size(); i++)
  {
    //Int_t timeBin = hBaseline[pmtNum]->FindBin(pmtEvent->time[i]);
    //hBaseline[pmtNum]->SetBinContent(timeBin,-baseline[i]);
  }
  return baseline;
}
//calculate baseline of whole waveform
std::pair<Double_t, Double_t> anaPulses::BaselineSubtraction(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry)
{
  std::sort(sig.begin(), sig.end());
  return std::make_pair(sig[sig.size() / 2], sig[sig.size() * .68]);
}

// seems like it is some type of smoothing filter...
std::vector<Double_t> anaPulses::TrigFilter(std::vector<Double_t> sig, Int_t N, Int_t pmtNum, Int_t ientry)
{

  if (N % 2 != 0)
    N++;
  std::vector<Double_t> vec(sig.size(), 0);
  for (int i = 0; i < sig.size(); i++)
  {
    Double_t topSum = 0, botSum = 0;
    Int_t window = min(i, N);
    for (int j = i - window / 2; j <= i + window / 2; j++)
    {
      if (j < 0)
        continue;
      if (j >= sig.size())
        continue;
      topSum += sig[j] * (1 + std::cos((j - i) * 2 * pi / N));
      botSum += (1 + std::cos((j - i) * 2 * pi / N));
    }
    if (botSum == 0)
      continue;
    vec[i] = topSum / botSum;
  }

  return vec;
}

Double_t anaPulses::CalculateMean(std::vector<Double_t> vec)
{
  Double_t mean = 0;
  for (int i = 0; i < vec.size(); i++)
    mean += vec[i];
  return mean / vec.size();
}

Double_t anaPulses::CalculateMean(std::vector<Int_t> vec)
{
  Double_t mean = 0;
  for (int i = 0; i < vec.size(); i++)
    mean += vec[i];
  return mean / vec.size();
}
Double_t anaPulses::CalculateMean(std::vector<Double_t> vec, Int_t startBin, Int_t stopBin)
{
  Double_t mean = 0;
  Double_t N = (stopBin - startBin);
  for (int i = startBin; i < stopBin; i++)
    mean += vec[i];
  return mean / N;
}

Double_t anaPulses::CalculateSdev(std::vector<Double_t> vec, Int_t startBin, Int_t stopBin, Double_t mean)
{
  Double_t sdev = 0;
  for (int i = startBin; i < stopBin; i++)
  {
    sdev += (mean - vec[i]) * (mean - vec[i]);
  }
  return std::sqrt(sdev / (stopBin - startBin));
}

std::vector<Double_t> anaPulses::BubbleSort(std::vector<Double_t> A)
{
  int i, j, N = A.size();

  for (i = 0; i < N; i++)
  {
    for (j = i + 1; j < N; j++)
    {
      if (A[j] < A[i])
      {
        Double_t buff;
        buff = A[i];
        A[i] = A[j];
        A[j] = buff;
      }
    }
  }
  return A;
}

std::vector<Double_t> anaPulses::TrapFilter(std::vector<Double_t> sig, Int_t ramp, Int_t flat, Int_t pmtNum, Int_t ientry)
{
  std::vector<Double_t> filter;
  if (flat % 2 != 0)
    flat++;
  for (int i = 0; i < sig.size(); i++)
  {
    Double_t val1 = 0, val2 = 0;
    Int_t minArr[] = {ramp, i, (Int_t)sig.size() - 1 - i};
    Int_t min = *std::min_element(minArr, minArr + 3);
    if (i < ramp + flat || i > sig.size() - ramp - flat - 1)
    {
      for (int j = 1; j < min; j++)
      {
        val1 += sig[i - j];
        val2 += sig[i + j];
      }
    }

    else
    {
      for (int j = 0; j < min; j++)
      {
        val1 += sig[i - j - flat];
        val2 += sig[i + j + flat];
      }
    }
    filter.push_back((val2 - val1) / (min + 1));
  }
  /*
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
  */
  /*
  hFilter[pmtNum] = new TH1D(TString("Trap_Filter_Channel_")+to_string(pmtNum)+TString("_event_")+to_string(ientry)+TString("_flat_")+to_string(flat)+TString("_ramp_")+to_string(ramp),"",filter.size(),0,filter.size());
  for(int i = 0; i < filter.size(); i++){
    hFilter[pmtNum]->SetBinContent(i,filter[i]);
  }
  */

  return filter;
}

std::vector<Double_t> anaPulses::LowPassFrequency(std::vector<Double_t> input, Double_t cutoff, Double_t sampleRate)
{

  std::vector<Double_t> output;
  output.resize(input.size());
  output[0] = input[0];
  Double_t RC = 1.0 / (cutoff * 2 * 3.1415);
  Double_t dt = sampleRate;
  Double_t alpha = dt / (RC + dt);
  for (int i = 1; i < input.size(); i++)
  {
    output[i] = output[i - 1] + (alpha * (input[i] - output[i - 1]));
  }
  return output;
}

std::vector<Double_t> anaPulses::NotchFilter(std::vector<Double_t> x, Double_t notch, Double_t bandWidth, Double_t sampleRate)
{
  int n = x.size();
  std::vector<Double_t> y;
  y.resize(n);

  //Both of these are expressed as a fraction of the sampling frequency, and therefore must be between 0 and 0.5
  double BW = bandWidth / sampleRate; // bandwith is expressed in a fraction of the sample rate
  double f = sampleRate / notch;      // not that this is inverse because we are dealing with smallest sample times
  if (f > 0.5 || BW > 0.5)
  {
    cout << "Notch filter parameters are choosen incorrectly... \n BW = " << BW << "...f = " << f << "see http://dspguide.com/ch19/3.htm" << endl;
    return x;
  }
  double R = 1 - 3 * BW;
  double K = (1 - 2 * R * std::cos(2 * 3.1415 * f) + R * R) / (2 - 2 * std::cos(2 * 3.1415 * f));

  double a0 = K;
  double a1 = -2 * K * std::cos(2 * 3.1415 * f);
  double a2 = K;
  double b1 = 2 * R * std::cos(2 * 3.1415 * f);
  double b2 = -R * R;

  static float x_2 = 0.0f; // delayed x, y samples
  static float x_1 = 0.0f;
  static float y_2 = 0.0f;
  static float y_1 = 0.0f;

  for (int i = 0; i < n; ++i)
  {
    y[i] = a0 * x[i] + a1 * x_1 + a2 * x_2 // IIR difference equation
           + b1 * y_1 + b2 * y_2;
    x_2 = x_1; // shift delayed x, y samples
    x_1 = x[i];
    y_2 = y_1;
    y_1 = y[i];
  }

  return y;
}
