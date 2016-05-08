#include "SeaveyDataHandler.h"
#include "TH2D.h"
#include "FFTtools.h"
#include "ProgressBar.h"
#include "TGraph2D.h"


// Silly globals...
const Int_t numAnts = 4;
const Int_t antNums[numAnts] = {30, 29, 26, 25};

// const Int_t azRange = 45;
// const Int_t deltaAz = 15;
// const Int_t elRange = 45;
// const Int_t deltaEl = 15;

// const Int_t nAnglesAz = 2*(azRange/deltaAz)+1;
// const Int_t nAnglesEl = 2*(elRange/deltaEl)+1;  

int main(){

  SeaveyDataHandler sData(8192*2);

  // Channel 1 is aligned, channel 4 is cross-pol
  Int_t channel = 1;

  const Double_t timeBefore = 5e-9; //2e-9;
  const Double_t timeAfter = 30e-9;

  // const Double_t timeBefore = 2e-9;
  // const Double_t timeAfter = 5e-9;

  TFile* fOut = new TFile("makeAntennaInformationPlots.root", "recreate");

  ProgressBar p(1);

  // List of plots I want...
  // Gain as a function of frequency (takes account of geometry of setup + Friis correction).
  // Group delay as a function of frequency on boresight.
  // Gain as a function of off axis delay.
  // Off axis delay -> Can't get this parameter because we're not rotating about the phase center.


  // list of steps for antenna information happyness...
  // to get the gain we need to first remove the effect of the pulse going through all the cables.
  // get all waveforms noise subtracted...


  // Reproducing the python analyis steps here...
  // get wave, noise subtract, window
  // remove copol cable responses (CRs.py)
  //      do some padding here... so that deltaF is the same for both waveforms
  //      fft the wave
  //      divide out the cable response
  //      divide out the pulse frequencies  
  
  for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){

    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    
    for(int antInd=0; antInd < numAnts; antInd++){

      Int_t antNum = antNums[antInd];
      
      TGraph* gr0 = sData.getBoresightGraphFromTFile(antNum, channel, pol);

      // const Double_t dt = gr0->GetX()[1] - gr0->GetX()[0];
      // const Double_t dtInterp = 0.1*dt;

      // gr0->SetTitle("Before");
      // gr0->SetName("grBefore");  
      // gr0->Write();
      sData.doNoiseSubtraction(gr0, antNum, channel, pol);

      // gr0->SetTitle("Noise subbed");
      // gr0->SetName("grNoiseSubbed");
      // gr0->Write();
      sData.windowPulse(gr0, timeBefore, timeAfter);
      // gr0->SetTitle("Final");
      // gr0->SetName("grFinal");

      TString gr0Name = TString::Format("gr%d_%d_%d_%d", antNum, 0, 0, polInd);
      gr0->SetName(gr0Name);
      gr0->Write();


      TGraph* gr0ps = FFTtools::makePowerSpectrum(gr0);
      gr0ps->SetName(gr0Name + "ps");
      gr0ps->Write();
      delete gr0ps;
      
      std::cout << gr0->GetN() << std::endl;
      
      FFTWComplex* seaveyToSeavey = sData.removeCopolResponse(gr0);
      std::cout << gr0->GetN() << std::endl;

      std::vector<Double_t> ps;
      std::vector<Double_t> freqs;

      for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){	
	ps.push_back(seaveyToSeavey[freqInd].getAbsSq());
	freqs.push_back(freqInd*sData.deltaF);
      }
      delete [] seaveyToSeavey;
      TGraph* gr0DecoPs = new TGraph(ps.size(), &freqs[0], &ps[0]);
      gr0DecoPs->SetName(gr0Name + "DecoPs");
      gr0DecoPs->Write();
      delete gr0DecoPs;
    }
  }


  std::vector<Double_t> pulsePower;
  std::vector<Double_t> freqs;
  for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){	
    pulsePower.push_back(sData.fftwComplexPulseFreqs[freqInd].getAbsSq());
    freqs.push_back(freqInd*sData.deltaF);
  }

  TGraph* grPulseFreqs = new TGraph(sData.numFreqs, &freqs[0], &pulsePower[0]);
  grPulseFreqs->SetName("grPulseFreqsPs");
  grPulseFreqs->Write();


  std::vector<Double_t> copolPower;  
  for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){	
    copolPower.push_back(sData.fftwComplexCopolCableResponse[freqInd].getAbsSq());
    freqs.push_back(freqInd*sData.deltaF);
  }

  TGraph* grCopolFreqs = new TGraph(sData.numFreqs, &freqs[0], &copolPower[0]);
  grCopolFreqs->SetName("grCopolFreqsPs");
  grCopolFreqs->Write();

  
  fOut->Write();
  fOut->Close();
}



// void makeMoreSightGainVsFrequencyPlots(const SeaveyDataHandler* sDataPtr){

  
  

// }
