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

// speed of light
const Double_t c = 299792458;
const Double_t separationMeters = 8.89; //m
const Double_t faceToPhaseCenter = 0.2; // 20cm


std::vector<Double_t> friisCorrection(const std::vector<Double_t>& ps, const std::vector<Double_t>& freqs);
inline void dBConversion(std::vector<double>& ps){
  for(UInt_t i=0; i < ps.size(); i++){
    // ps[i] = 10*TMath::Log10(ps.at(i)/10);
    ps[i] = 10*TMath::Log10(ps.at(i)/10);    
  }

}

inline void convertXaxisFromHzToMHz(TGraph* gr){
  for(int i=0; i < gr->GetN(); i++){
    gr->GetX()[i] *= 1e-6;
  }
}

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

  const int numSeaveyPoints = 5;
  Double_t seaveyNumsVPol[numSeaveyPoints] = {6.3, 7.7, 9.5, 9.0, 12.5};
  Double_t seaveyNumsHPol[numSeaveyPoints] = {6.0, 8.1, 10.1, 8.0, 12.8};
  Double_t seaveyFreqsMHz[numSeaveyPoints] = {200, 450, 700, 950, 1200};  

  TGraph* grSeaveyNumsHPol = new TGraph(numSeaveyPoints, seaveyFreqsMHz, seaveyNumsHPol);
  TGraph* grSeaveyNumsVPol = new TGraph(numSeaveyPoints, seaveyFreqsMHz, seaveyNumsVPol);

  grSeaveyNumsVPol->SetName("grSeaveyNumsVPol");
  grSeaveyNumsVPol->Write();
  
  grSeaveyNumsHPol->SetName("grSeaveyNumsHPol");
  grSeaveyNumsHPol->Write();
  
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

  const Double_t minFreq = 70e6;
  const Double_t maxFreq = 1425e6;

  std::vector<double> pulserPowerNoCables;
  std::vector<double> freqs;  
  for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
    double powerTimesTime = sData.fftwComplexPulseFreqs[freqInd].getAbsSq();
    double power = powerTimesTime/sData.dtCables;
    pulserPowerNoCables.push_back(power);
    freqs.push_back(freqInd*sData.deltaF);
  }

  
  for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){

    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    
    for(int antInd=0; antInd < numAnts; antInd++){
      if(antInd > 0) {
	continue;
      }
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
      convertXaxisFromHzToMHz(gr0ps);      
      gr0ps->Write();
      delete gr0ps;
      
      std::cout << gr0->GetN() << std::endl;
      double dt = gr0->GetX()[1] - gr0->GetX()[0];
      
      FFTWComplex* seaveyToSeavey = sData.removeCopolResponse(gr0);
      std::cout << gr0->GetN() << std::endl;

      std::vector<Double_t> ps;
      std::vector<Double_t> phaseResponse;
      std::vector<Double_t> groupDelay;
      std::vector<Double_t> theseFreqs;

      for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
	Double_t f = freqInd*sData.deltaF;
	if(f >= minFreq && f < maxFreq){
	  Double_t power = seaveyToSeavey[freqInd].getAbsSq();
	  power /= dt;

	  ps.push_back(power);	  
	  theseFreqs.push_back(f);

	  phaseResponse.push_back(seaveyToSeavey[freqInd].getPhase());
	}
      }
      
      std::vector<Double_t> gain = friisCorrection(ps, theseFreqs);
      
      dBConversion(gain);
      dBConversion(ps);
      
      TGraph* gr0Gain = new TGraph(gain.size(), &theseFreqs[0], &gain[0]);
      gr0Gain->SetName(gr0Name + "Gain");
      convertXaxisFromHzToMHz(gr0Gain);
      gr0Gain->Write();
      delete gr0Gain;

      TGraph* gr0Phase = new TGraph(phaseResponse.size(), &theseFreqs[0], &phaseResponse[0]);
      gr0Phase->SetName(gr0Name + "Phase");

      TGraph* gr0Group = (TGraph*) gr0Phase->Clone(gr0Name + "Group");
      
      for(int i=0; i < gr0Group->GetN()-1; i++){
	Double_t y0 = gr0Group->GetY()[i];
	Double_t y1 = gr0Group->GetY()[i+1];	
	while(y1 < y0){
	  y1 += TMath::TwoPi();
	}
	Double_t dy = y1 - y0;
	Double_t dx = gr0Group->GetX()[i+1] - gr0Group->GetX()[i];
	gr0Group->GetY()[i] = dy/dx;
	// gr0Group->GetY()[i] -= (separationMeters + 2*faceToPhaseCenter)/c;
      }
      gr0Group->RemovePoint(gr0Group->GetN()-1);

      convertXaxisFromHzToMHz(gr0Phase);
      convertXaxisFromHzToMHz(gr0Group);
      gr0Phase->Write();      
      gr0Group->Write();
      delete gr0Group;
      delete gr0Phase;
      

      TGraph* gr0DecoPs = new TGraph(ps.size(), &theseFreqs[0], &ps[0]);
      gr0DecoPs->SetName(gr0Name + "DecoPs");
      convertXaxisFromHzToMHz(gr0DecoPs);
      gr0DecoPs->Write();
      delete gr0DecoPs;
      
      delete [] seaveyToSeavey;      
    }
  }

  TGraph* grPulseFreqs = new TGraph(sData.numFreqs, &freqs[0], &pulserPowerNoCables[0]);
  grPulseFreqs->SetName("grPulseFreqsPs");
  grPulseFreqs->Write();

  std::vector<Double_t> copolPower;
  for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
    copolPower.push_back(sData.fftwComplexCopolCableResponse[freqInd].getAbs());
  }

  TGraph* grCopolFreqs = new TGraph(sData.numFreqs, &freqs[0], &copolPower[0]);
  grCopolFreqs->SetName("grCopolFreqsPs");
  grCopolFreqs->Write();
  
  fOut->Write();
  fOut->Close();
}



std::vector<Double_t> friisCorrection(const std::vector<Double_t>& ps, const std::vector<Double_t>& freqs){

  const double partFactor = 4*TMath::Pi()*(separationMeters + 2*faceToPhaseCenter)/c;

  std::vector<Double_t> gain(ps.size(), 0);
  
  // separationFactors = [(f*4*math.pi*distMeters/c)**2 for f in freqsMHz]
 
  for(UInt_t i=0; i < ps.size(); i++){

    Double_t f = freqs.at(i);
    Double_t friisFactor = (f*f*partFactor*partFactor);
    
    // ps[i] *= friisFactor
    gain.at(i) = TMath::Sqrt(ps[i]*friisFactor);
  }
  return gain;
}

// void makeMoreSightGainVsFrequencyPlots(const SeaveyDataHandler* sDataPtr){



  
  

// }
