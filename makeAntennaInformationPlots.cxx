#include "SeaveyDataHandler.h"
#include "TH2D.h"
#include "FFTtools.h"
#include "ProgressBar.h"
#include "TGraph2D.h"


// Silly globals...
// const Int_t numAnts = 4;
// const Int_t antNums[numAnts] = {30, 29, 26, 25};

// const Int_t azRange = 45;
// const Int_t deltaAz = 15;
// const Int_t elRange = 45;
// const Int_t deltaEl = 15;

// const Int_t nAnglesAz = 2*(azRange/deltaAz)+1;
// const Int_t nAnglesEl = 2*(elRange/deltaEl)+1;

// speed of light
const Double_t c = 299792458;
const Double_t faceToFaceSeparationMeters = 8.89; //m
const Double_t faceToPhaseCenter = 0.2; // 20cm
const double phaseCenterToPhaseCenterSeparationMeters = faceToFaceSeparationMeters + 2*faceToPhaseCenter;
// const Double_t phaseCenterToPhaseCenterTime = phaseCenterToPhaseCenterSeparationMeters/c;


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
  sData.kPrintWarnings = 1;

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

  // const Double_t minFreq = 70e6;
  const Double_t minFreq = 160e6;  
  // const Double_t maxFreq = 1425e6;
  const Double_t maxFreq = 1570e6;
  // const Double_t maxFreq = 1570e6;    
  const double maxFreqs[AnitaPol::kNotAPol] = {1570e6, 1430e6};
  
  std::vector<double> pulserPowerNoCables;
  std::vector<double> freqs;  
  for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
    double powerTimesTime = sData.fftwComplexPulseFreqs[freqInd].getAbsSq();
    double power = powerTimesTime/sData.dtCables;
    pulserPowerNoCables.push_back(power);
    freqs.push_back(freqInd*sData.deltaF);
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

  // for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
  for(int polInd = AnitaPol::kVertical; polInd >= AnitaPol::kHorizontal; polInd--){

    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    
    // for(int antInd=0; antInd < numAnts; antInd++){
    // for(int ant=0; ant < 52; ant++){
    // for(int antInd=0; antInd < 5; antInd++){
    for(int antInd=0; antInd < 51; antInd++){
      // for(int antInd=0; antInd < numAnts; antInd++){
      //   if(antInd > 0) {
      // 	continue;
      //   }
      // Int_t antNum = antNums[antInd];
      int antNum = antInd + 1;

      TGraph* gr0 = sData.getBoresightGraphFromTFile(antNum, channel, pol);
      std::cerr << gr0 << "\t" << antNum << std::endl;
      if(gr0==NULL){
	continue;
      }
      sData.doNoiseSubtraction(gr0, antNum, channel, pol);
      sData.windowPulse(gr0, timeBefore, timeAfter);

      TString gr0Name = TString::Format("gr%d_%d_%d_%d", antNum, 0, 0, polInd);
      gr0->SetName(gr0Name);
      gr0->Write();

      // TGraph* gr0ps = FFTtools::makePowerSpectrum(gr0);
      // gr0ps->SetName(gr0Name + "ps");
      // gr0ps->Write();

      // FFTWComplex* rawFFT = sData.doNormalizedFFT(gr0);
      // TGraph* grRawPhase = new TGraph();
      // for(int samp=0; samp < gr0ps->GetN(); samp++){
      // 	grRawPhase->SetPoint(grRawPhase->GetN(), gr0ps->GetX()[samp], rawFFT[samp].getPhase());
      // }
      // grRawPhase->SetName(gr0Name + "RawPhase");
      // grRawPhase->Write();

      // TGraph* grRawGroup = (TGraph*) grRawPhase->Clone(gr0Name + "RawGroup");
      // double dxRaw = TMath::TwoPi()*(grRawGroup->GetX()[1]-grRawGroup->GetX()[0]);
      // for(int i=0; i < grRawGroup->GetN()-1; i++){
      // 	double y1 = grRawPhase->GetY()[i+1];
      // 	double y0 = grRawPhase->GetY()[i];
      // 	if(y1 > y0){
      // 	  y1 -= TMath::TwoPi();
      // 	}
      // 	double dy = y1 - y0;
      // 	grRawGroup->GetY()[i] = -dy/dxRaw;
      // 	grRawGroup->GetY()[i] *= 1e9;
      // }
      // grRawGroup->Write();

      // std::cerr << "gr0Dt = " << gr0->GetX()[1] - gr0->GetX()[0] << std::endl;
      // std::cerr << "gr0N = " << gr0->GetN() << std::endl;
      // std::cerr << "dfRaw = " << gr0ps->GetX()[1] - gr0ps->GetX()[0] << std::endl;

      // convertXaxisFromHzToMHz(gr0ps);
      
      // delete gr0ps;
      // delete grRawPhase;
      // delete grRawGroup;
      // delete [] rawFFT;
      
      // std::cout << gr0->GetN() << std::endl;
      double dt = gr0->GetX()[1] - gr0->GetX()[0];
      
      FFTWComplex* seaveyToSeavey = sData.removeCopolResponse(gr0);
      std::cout << gr0->GetN() << std::endl;

      std::vector<Double_t> ps;
      std::vector<Double_t> phaseResponse;
      std::vector<Double_t> theseFreqs;

      FFTWComplex* impulseRespFreq = new FFTWComplex[sData.numFreqs];
      
      for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
	Double_t f = freqInd*sData.deltaF;
	if(f >= minFreq && f < maxFreq){
	  Double_t power = seaveyToSeavey[freqInd].getAbsSq();
	  power /= dt;

	  ps.push_back(power);
	  theseFreqs.push_back(f);

	  Double_t phase = seaveyToSeavey[freqInd].getPhase();
	  if(phase < -TMath::Pi()){
	    phase += TMath::TwoPi();
	  }
	  else if(phase >= TMath::Pi()){
	    phase -= TMath::TwoPi();
	  }

	  double mag = seaveyToSeavey[freqInd].getAbs();
	  // impulseRespFreq[freqInd].setMagPhase(mag, 0);
	  impulseRespFreq[freqInd].setMagPhase(mag, phase/2);
	  phaseResponse.push_back(phase);
	}
	else{
	  impulseRespFreq[freqInd].re = 0;
	  impulseRespFreq[freqInd].im = 0;
	}
      }
      
      // std::vector<Double_t> gain = ps;
      // double minGainInBand = 1e9;
      // double maxGainInBand = -1e9;
      // for(UInt_t i=0; i < gain.size(); i++){
      // 	gain.at(i)/=2; // half
      // }
      // dBConversion(gain);
      // for(UInt_t i=0; i < gain.size(); i++){
      // 	if(theseFreqs.at(i) > 200e6 && theseFreqs.at(i) < 1200e6){
      // 	  if(gain.at(i) < minGainInBand){
      // 	    minGainInBand = gain.at(i);
      // 	  }
      // 	  if(gain.at(i) > maxGainInBand){
      // 	    maxGainInBand = gain.at(i);
      // 	  }
      // 	}
      // }
      // for(UInt_t i=0; i < gain.size(); i++){
      // 	gain.at(i) -= maxGainInBand;
      // 	std::cerr << theseFreqs[i]/1e6 << "\t" << gain.at(i) << std::endl;
      // }
      // minGainInBand -= maxGainInBand;

      // Double_t lowEdge = -1000;
      // Double_t highEdge = 1e15;
      // for(UInt_t i=0; i < gain.size(); i++){
      // 	if(gain.at(i) < minGainInBand - 3 && theseFreqs.at(i) < 200e6 && theseFreqs.at(i) >= lowEdge){
      // 	  lowEdge = theseFreqs.at(i);
      // 	}
      // }
      // for(UInt_t i=gain.size()-1; i != 0; i--){
      // 	std::cerr << theseFreqs.at(i)/1e6 << "\t" << gain.at(i) << "\t" << minGainInBand << std::endl;
      //   if(gain.at(i) < minGainInBand - 3 && theseFreqs.at(i) > 1200e6){
      // 	  highEdge = theseFreqs.at(i);
      // 	}
      // 	  if(theseFreqs.at(i) < 1200e6){
      // 	    break;
      // 	  }
      // }
      
      // // std::cerr << doneHighEdge << "\t" << minGainInBand - maxGainInBand << std::endl;
      // std::cerr << "the band = " << lowEdge/1e6 << " - " << highEdge/1e6 << "MHz" << std::endl;
      
      // TGraph* gr0Gain = new TGraph(gain.size(), &theseFreqs[0], &gain[0]);
      // gr0Gain->SetName(gr0Name + "Gain");
      // convertXaxisFromHzToMHz(gr0Gain);
      // gr0Gain->Write();
      // delete gr0Gain;
      
      std::vector<Double_t> gain_dBi = friisCorrection(ps, theseFreqs);
      
      dBConversion(gain_dBi);
      dBConversion(ps);
      
      // TGraph* grSeaveyToSeavey = sData.doNormalizedInvFFT(gr0->GetN(), seaveyToSeavey, sData.deltaF);
      // grSeaveyToSeavey->SetName(gr0Name+ "SeaveyToSeavey");
      // grSeaveyToSeavey->Write();
      // delete grSeaveyToSeavey;
      
      TGraph* gr0GaindBi = new TGraph(gain_dBi.size(), &theseFreqs[0], &gain_dBi[0]);
      gr0GaindBi->SetName(gr0Name + "GaindBi");
      convertXaxisFromHzToMHz(gr0GaindBi);
      gr0GaindBi->Write();
      delete gr0GaindBi;

      TGraph* gr0Phase = new TGraph(phaseResponse.size(), &theseFreqs[0], &phaseResponse[0]);
      gr0Phase->SetName(gr0Name + "Phase");

      // TGraph* gr0PhaseUnwrapped = (TGraph*) gr0Phase->Clone(gr0Name + "PhaseUnwrapped");

      // if(gr0PhaseUnwrapped->GetN() > 1){
      // 	for(int i=0; i < gr0PhaseUnwrapped->GetN()-1; i++){
      // 	  Double_t y0 = gr0PhaseUnwrapped->GetY()[i];
      // 	  Double_t y1 = gr0PhaseUnwrapped->GetY()[i+1];

      // 	  if(y1 - y0 > TMath::Pi()){
      // 	    y1 -= TMath::TwoPi();
      // 	  }
      // 	  else if(y1 - y0 < -TMath::Pi()){
      // 	    y1 += TMath::TwoPi();	  
      // 	  }	  
      // 	}
      // }

      // deltaF is about 1.2MHz
      // this means that deltaOmega = 7.7 MRads per second
      // if the phase didn't change at all, the group delay would be zero
      // if it changes by 1 radian then dPhi/dOmega = 1.30380e-07
      // if it changes by pi then the group delay is dPhi/dOmega = 8.19200e-07
      // you can make the group delay larger by having abs(dPhi) > TwoPi.
      // but you can't make it smaller unless the phase variation is TINY...
      // does the interpolation do something strage here?
      
      TGraph* gr0Group = (TGraph*) gr0Phase->Clone(gr0Name + "Group");

      Double_t dx = TMath::TwoPi()*(theseFreqs[1]-theseFreqs[0]);
      for(int i=0; i < gr0Group->GetN()-1; i++){
	Double_t y0 = gr0Phase->GetY()[i];
	Double_t y1 = gr0Phase->GetY()[i+1];
	
	if(y1 - y0 > TMath::Pi()){
	  y1 -= TMath::TwoPi();
	}
	else if(y1 - y0 <= -TMath::Pi()){
	  y1 += TMath::TwoPi();	  
	}
	Double_t dy = -(y1 - y0);
	gr0Group->GetY()[i] = dy/dx;
      }
      gr0Group->RemovePoint(gr0Group->GetN()-1);

      // now divide by two to account for the fact that you've gone through two antennas...      
      for(int i=0; i < gr0Group->GetN()-1; i++){
	gr0Group->GetY()[i]/=2;
      }
      
      // and finally, I think this is how this should work.
      // I go through my assumed group delay to generate a new set of phases...
      // is this just the same as dividing the magnitude by two? Not sure at this point...
      // grGroup was generated doing [i+1] - [i], so if you have [i], you get [i+1].
      // assume an initial phase of pi...


      double f0 = gr0Group->GetX()[0];
      int freqInd0 = TMath::Nint(f0/sData.deltaF);
      impulseRespFreq[freqInd0].setMagPhase(impulseRespFreq[freqInd0].getAbs(), TMath::Pi()/2);
					    
      for(int i=0; i < gr0Group->GetN(); i++){
	Int_t freqInd = freqInd0 + i;	
	double thisPhase = impulseRespFreq[freqInd].getPhase();

	double minusDyByDx = gr0Group->GetY()[i];
	double thisDx = sData.deltaF*TMath::TwoPi();
	double deltaPhase = -1*minusDyByDx*thisDx;
	double nextPhase = thisPhase + deltaPhase;

	// impulseRespFreq[freqInd+1].setMagPhase(impulseRespFreq[freqInd+1].getAbs(), nextPhase);
	impulseRespFreq[freqInd+1].setMagPhase(impulseRespFreq[freqInd+1].getAbs(), nextPhase);	
      }
      

      TGraph* grAssumedImpulseResponse = sData.doNormalizedInvFFT(sData.nPointsCables,
								  impulseRespFreq,
								  sData.deltaF);
      grAssumedImpulseResponse->SetName(gr0Name + "ImpulseResponse");
      grAssumedImpulseResponse->Write();
      delete grAssumedImpulseResponse;
      delete [] impulseRespFreq;
	    
      
      // convertXaxisFromHzToMHz(gr0Phase);
      convertXaxisFromHzToMHz(gr0Group);
      // gr0Phase->Write();      

      // convert from s to ns
      for(int i=0; i < gr0Group->GetN(); i++){
	gr0Group->GetY()[i] *= 1e9; 
      }
      gr0Group->Write();      
      delete gr0Group;
      // delete gr0Phase;



      // TGraph* gr0DecoPs = new TGraph(ps.size(), &theseFreqs[0], &ps[0]);
      // gr0DecoPs->SetName(gr0Name + "DecoPs");
      // convertXaxisFromHzToMHz(gr0DecoPs);
      // gr0DecoPs->Write();
      // delete gr0DecoPs;
      
      delete [] seaveyToSeavey;
      delete gr0;
    }

    
  }

  
  fOut->Write();
  fOut->Close();
}



std::vector<Double_t> friisCorrection(const std::vector<Double_t>& ps, const std::vector<Double_t>& freqs){

  const double partFactor = 4*TMath::Pi()*(phaseCenterToPhaseCenterSeparationMeters)/c;

  std::vector<Double_t> gain_dBi(ps.size(), 0);
  
  // separationFactors = [(f*4*math.pi*distMeters/c)**2 for f in freqsMHz]
 
  for(UInt_t i=0; i < ps.size(); i++){

    Double_t f = freqs.at(i);
    Double_t friisFactor = (f*f*partFactor*partFactor);
    
    // ps[i] *= friisFactor
    gain_dBi.at(i) = TMath::Sqrt(ps[i]*friisFactor);
  }
  return gain_dBi;
}




  
  

// }
