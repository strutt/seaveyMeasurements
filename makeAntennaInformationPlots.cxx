#include "SeaveyDataHandler.h"
#include "TH2D.h"
#include "FFTtools.h"
#include "ProgressBar.h"
#include "TGraph2D.h"
#include "TGraphPolar.h"


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
const Double_t faceToFaceSeparationMeters = 8.89; //m
const Double_t faceToPhaseCenter = 0.2; // 20cm
const double phaseCenterToPhaseCenterSeparationMeters = faceToFaceSeparationMeters + 2*faceToPhaseCenter;
// const Double_t phaseCenterToPhaseCenterTime = phaseCenterToPhaseCenterSeparationMeters/c;

const double partFactor = 4*TMath::Pi()*(phaseCenterToPhaseCenterSeparationMeters)/c;


std::vector<Double_t> friisCorrection(const std::vector<Double_t>& ps,
				      const std::vector<Double_t>& freqs,
				      AnitaPol::AnitaPol_t pol);

inline void dBConversion(std::vector<double>& ps){
  for(UInt_t i=0; i < ps.size(); i++){
    // ps[i] = 10*TMath::Log10(ps.at(i)/10);
    // ps[i] = 10*TMath::Log10(ps.at(i)/10);
    ps[i] = 10*TMath::Log10(ps.at(i));        
  }
}

inline void undBConversion(std::vector<double>& ps){
  for(UInt_t i=0; i < ps.size(); i++){
    // ps[i] = 10*TMath::Log10(ps.at(i)/10);
    // ps[i] = 10*TMath::Log10(ps.at(i)/10);

    // double temp1 = ps.at(i)/10;
    // double temp2 = TMath::Log10(temp1);
    // double temp3 = 10*temp2;
    // ps[i] = temp3;
    double temp1 = ps[i]/10;
    double temp2 = pow(10, temp1);
    double temp3 = temp2*10;
    ps[i] = temp3;    
  }
}

inline void convertXaxisFromHzToMHz(TGraph* gr){
  for(int i=0; i < gr->GetN(); i++){
    gr->GetX()[i] *= 1e-6;
  }
}

std::vector<double> meanVPolGain;
std::vector<double> meanVPolGroupDelay;

int main(){

  SeaveyDataHandler sData(8192*2);
  sData.kPrintWarnings = 1;
  
  // Channel 1 is aligned, channel 4 is cross-pol
  Int_t channel = 1;

  // const Double_t timeBefore = 5e-9; //2e-9;
  // const Double_t timeAfter = 30e-9;
  const Double_t timeBefore = 4e-9;
  const Double_t timeAfter = 31e-9;

  // const Double_t timeBefore = 2e-9;
  // const Double_t timeAfter = 5e-9;

  TFile* fOut = new TFile("makeAntennaInformationPlots.root", "recreate");


  TGraph* grPulsePs = new TGraph();
  grPulsePs->Set(sData.numFreqs); // alloc arrays in TGraph
  std::vector<double> pulsePs(sData.numFreqs, 0);
  for(int f=0; f < sData.numFreqs; f++){
    pulsePs.at(f) = sData.fftwComplexPulseFreqs[f].getAbsSq();
  }
  dBConversion(pulsePs);
  for(int f=0; f < sData.numFreqs; f++){  
    grPulsePs->SetPoint(f, f*sData.deltaF, pulsePs.at(f));
  }
  grPulsePs->SetName("grPulsePs");
  grPulsePs->SetTitle("Picosecond pulser frequency content; Frequency (Hz); Power (dB)");  
  grPulsePs->Write();
  delete grPulsePs;

  
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
  // relative off axis delay?

  TGraph* grMeanGaindBi[AnitaPol::kNotAPol] = {NULL};
  TGraph* grMeanGroupDelay[AnitaPol::kNotAPol] = {NULL};
  for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){  
    grMeanGaindBi[polInd] = new TGraph();
    TString grName = polInd == AnitaPol::kHorizontal ? "grMeanGaindBiHPol" : "grMeanGaindBiVPol";
    TString grTitle = polInd == AnitaPol::kHorizontal ? "Mean HPol Gain" : "Mean VPol Gain";
    grTitle += "; Frequency (MHz); Gain (dBi)";
    grMeanGaindBi[polInd]->SetName(grName);
    grMeanGaindBi[polInd]->SetTitle(grTitle);


    grMeanGroupDelay[polInd] = new TGraph();
    grName = polInd == AnitaPol::kHorizontal ? "grMeanGroupDelayHPol" : "grMeanGroupDelayVPol";
    grTitle = polInd == AnitaPol::kHorizontal ? "Mean HPol Group Delay" : "Mean VPol Group Delay";
    grTitle += "; Frequency (MHz); Group delay (ns)";
    grMeanGroupDelay[polInd]->SetName(grName);
    grMeanGroupDelay[polInd]->SetTitle(grTitle);
  }


  // const Double_t minFreq = 160e6;
  const Double_t minFreq = 150e6;  
  // const double maxFreqs[AnitaPol::kNotAPol] = {1570e6, 1430e6};
  const double maxFreqs[AnitaPol::kNotAPol] = {1570e6, 1570e6};  
  
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

  const int numAntInds = 51;
  
  for(int polInd = AnitaPol::kVertical; polInd >= AnitaPol::kHorizontal; polInd--){

    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    
    for(int antInd=0; antInd < numAntInds; antInd++){
      int antNum = antInd + 1;

      TGraph* gr0 = sData.getBoresightGraphFromTFile(antNum, channel, pol);
      if(gr0==NULL){
	std::cerr << antNum << " was skipped!!!!" << std::endl;
	continue;
      }
      sData.doNoiseSubtraction(gr0, antNum, channel, pol);
      sData.windowPulse(gr0, timeBefore, timeAfter);

      TString gr0Name = TString::Format("gr%d_%d_%d_%d", antNum, 0, 0, polInd);
      gr0->SetName(gr0Name);
      gr0->Write();

      // double dt = gr0->GetX()[1] - gr0->GetX()[0];
      
      // FFTWComplex* seaveyToSeavey = sData.removeCopolResponse(gr0);
      FFTWComplex* seaveyToSeavey = sData.doNormalizedFFT(gr0);

      // std::cerr << seaveyToSeavey[1].re << "\t" << seaveyToSeavey[1].im << "\t"
      // 		<< seaveyToSeavey[1].getAbsSq() << std::endl;
      
      // std::cout << gr0->GetN() << std::endl;

      std::vector<Double_t> ps;
      std::vector<Double_t> phaseResponse;
      std::vector<Double_t> theseFreqs;
      std::vector<Double_t> gain_dBi;
      
      FFTWComplex* impulseRespFreq = new FFTWComplex[sData.numFreqs];

      for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
	Double_t f = freqInd*sData.deltaF;
	if(f >= minFreq && f < maxFreqs[polInd]){
	  Double_t power = seaveyToSeavey[freqInd].getAbsSq();
	  // power /= dt;

	  ps.push_back(power);
	  theseFreqs.push_back(f);

	  // double PrOverPt = power/(sData.fftwComplexPsPulserCopolFast[freqInd].getAbsSq()/sData.dtCables);
	  double PrOverPt = power/(sData.fftwComplexPsPulserCopolFast[freqInd].getAbsSq());

	  double GrGt = PrOverPt*(partFactor*partFactor*f*f);
	  // double gainSquared = PrOverPt;

	  if(polInd==((int)AnitaPol::kVertical)){
	    gain_dBi.push_back(TMath::Sqrt(GrGt));
	  }
	  else{
	    const int i = gain_dBi.size();
	    gain_dBi.push_back(GrGt/meanVPolGain.at(i));
	  }
	  // gain_dBi.push_back(TMath::Sqrt(gainSquared));

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

      // std::vector<Double_t> gain_dBi = friisCorrection(ps, theseFreqs, AnitaPol::AnitaPol_t (polInd));

      // std::vector<double> gain_dBi = getGain(&sData, seaveyToSeavey, minFreq, maxFreqs[polInd]);
      
      if(polInd==AnitaPol::kVertical){
	if(meanVPolGain.size()==0){
	  meanVPolGain = gain_dBi;
	}
	else{
	  for(UInt_t i=0; i < gain_dBi.size(); i++){
	    meanVPolGain.at(i) += gain_dBi.at(i);
	  }
	}
      }
      
      dBConversion(gain_dBi);
      dBConversion(ps);

      TGraph* gr0GaindBi = new TGraph(gain_dBi.size(), &theseFreqs[0], &gain_dBi[0]);
      gr0GaindBi->SetName(gr0Name + "GaindBi");
      convertXaxisFromHzToMHz(gr0GaindBi);
      gr0GaindBi->Write();

      if(grMeanGaindBi[polInd]->GetN()==0){
	for(int i=0; i < gr0GaindBi->GetN(); i++){
	  grMeanGaindBi[polInd]->SetPoint(i, gr0GaindBi->GetX()[i], gr0GaindBi->GetY()[i]);
	}
      }
      else{
	for(int i=0; i < gr0GaindBi->GetN(); i++){	
	  grMeanGaindBi[polInd]->GetY()[i] += gr0GaindBi->GetY()[i];
	}
      }
      delete gr0GaindBi;

      TGraph* gr0Gain = new TGraph(ps.size(), &theseFreqs[0], &ps[0]);
      gr0Gain->SetName(gr0Name + "Gain");
      convertXaxisFromHzToMHz(gr0Gain);
      gr0Gain->Write();
      delete gr0Gain;

      TGraph* gr0Phase = new TGraph(phaseResponse.size(), &theseFreqs[0], &phaseResponse[0]);
      gr0Phase->SetName(gr0Name + "Phase");

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
      if(polInd==AnitaPol::kVertical){
	for(int i=0; i < gr0Group->GetN(); i++){
	  gr0Group->GetY()[i]/=2;
	}	
	if(meanVPolGroupDelay.size()==0){
	  for(int i=0; i < gr0Group->GetN(); i++){
	    meanVPolGroupDelay.push_back(gr0Group->GetY()[i]);
	  }		  	
	}
	else{
	  for(Int_t i=0; i < gr0Group->GetN(); i++){
	    meanVPolGroupDelay.at(i) += gr0Group->GetY()[i];
	  }
	}
      }
      else{
	for(int i=0; i < gr0Group->GetN(); i++){
	  gr0Group->GetY()[i] -= meanVPolGroupDelay.at(i);
	}
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

	impulseRespFreq[freqInd+1].setMagPhase(impulseRespFreq[freqInd+1].getAbs(), nextPhase);	
      }
      
      TGraph* grAssumedImpulseResponse = sData.doNormalizedInvFFT(sData.nPointsCables,
								  impulseRespFreq,
								  sData.deltaF);
      grAssumedImpulseResponse->SetName(gr0Name + "ImpulseResponse");
      grAssumedImpulseResponse->Write();
      delete grAssumedImpulseResponse;
      delete [] impulseRespFreq;
	    
      convertXaxisFromHzToMHz(gr0Group);

      // convert from s to ns
      for(int i=0; i < gr0Group->GetN(); i++){
	gr0Group->GetY()[i] *= 1e9;
      }
      gr0Group->Write();
      
      // add to average
      if(grMeanGroupDelay[polInd]->GetN()==0){
	for(int i=0; i < gr0Group->GetN(); i++){
	  grMeanGroupDelay[polInd]->SetPoint(i, gr0Group->GetX()[i], gr0Group->GetY()[i]);
	}
      }
      else{
	for(int i=0; i < gr0Group->GetN(); i++){
	  grMeanGroupDelay[polInd]->GetY()[i] += gr0Group->GetY()[i];
	}
      }

      delete gr0Group;
      
      delete [] seaveyToSeavey;
      delete gr0;
    }

    for(int i=0; i < grMeanGaindBi[polInd]->GetN(); i++){
      grMeanGaindBi[polInd]->GetY()[i]/=numAntInds;
    }
    grMeanGaindBi[polInd]->Write();

    for(UInt_t i=0; i < meanVPolGain.size(); i++){
      meanVPolGain.at(i)/=numAntInds;
    }
    for(UInt_t i=0; i < meanVPolGroupDelay.size(); i++){
      meanVPolGroupDelay.at(i)/=numAntInds;
    }
    
    for(int i=0; i < grMeanGroupDelay[polInd]->GetN(); i++){
      grMeanGroupDelay[polInd]->GetY()[i]/=numAntInds;
    }
    grMeanGroupDelay[polInd]->Write();
  }

  for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    delete grMeanGaindBi[polInd];
    delete grMeanGroupDelay[polInd];
  }
  
  for(int polInd = AnitaPol::kVertical; polInd >= AnitaPol::kHorizontal; polInd--){

    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    
    for(int antInd=0; antInd < numAnts; antInd++){
      
      int antNumber = antNums[antInd];

      TGraph2D* grMeanGain2D = new TGraph2D();
      TGraphPolar* grAz = new TGraphPolar();
      TGraphPolar* grEl = new TGraphPolar();      
      
      double theBoresightVal = 0;
      double minVal = DBL_MAX;
      for(int az=-45; az <=45; az+=15){
	for(int el=-45; el <=45; el+=15){
	  if(antNumber==26 && el==30){
	    continue;
	  }

	  
	  sData.kPrintWarnings = 0;
	  TGraph* gr0 = sData.getOffAxisGraphFromTFile(antNumber, channel, pol, az, el);
	  sData.kPrintWarnings = 1;
	  if(gr0==NULL){
	    // std::cerr << antNum << " was skipped!!!!" << std::endl;
	    continue;
	  }

	  
	  sData.doNoiseSubtraction(gr0, antNumber, channel, pol);
	  sData.windowPulse(gr0, timeBefore, timeAfter);

	  if(el < 0){
	    for(int samp=0; samp < gr0->GetN(); samp++){
	      gr0->GetY()[samp] *= -1;
	    }
	  }

	  TString gr0Name = TString::Format("grOff%d_%d_%d_%d", antNumber, az, el, polInd);

	  gr0->SetName(gr0Name);
	  gr0->Write();

	  double dt = gr0->GetX()[1] - gr0->GetX()[0];

	  FFTWComplex* seaveyToSeavey = sData.removeCopolResponse(gr0);
	  std::vector<Double_t> ps;
	  std::vector<Double_t> phaseResponse;
	  std::vector<Double_t> theseFreqs;

	  FFTWComplex* impulseRespFreq = new FFTWComplex[sData.numFreqs];
      
	  for(int freqInd=0; freqInd < sData.numFreqs; freqInd++){
	    Double_t f = freqInd*sData.deltaF;
	    if(f >= minFreq && f < maxFreqs[polInd]){
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
      
      
	  std::vector<Double_t> gain_dBi = friisCorrection(ps, theseFreqs, AnitaPol::AnitaPol_t (polInd));
	  if(polInd==AnitaPol::kVertical){
	    if(meanVPolGain.size()==0){
	      meanVPolGain = gain_dBi;
	    }
	    else{
	      for(UInt_t i=0; i < gain_dBi.size(); i++){
	  	meanVPolGain.at(i) += gain_dBi.at(i);
	      }
	    }
	  }


      
	  dBConversion(gain_dBi);
	  dBConversion(ps);

	  TGraph* gr0GaindBi = new TGraph(gain_dBi.size(), &theseFreqs[0], &gain_dBi[0]);
	  gr0GaindBi->SetName(gr0Name + "GaindBi");
	  convertXaxisFromHzToMHz(gr0GaindBi);
	  gr0GaindBi->Write();
	  delete gr0GaindBi;


	  TGraph* gr0Gain = new TGraph(ps.size(), &theseFreqs[0], &ps[0]);
	  gr0Gain->SetName(gr0Name + "Gain");
	  convertXaxisFromHzToMHz(gr0Gain);
	  gr0Gain->Write();
	  delete gr0Gain;


	  double mean=0;
	  int numSamps = 0;
	  for(int samp=0; samp < gr0Gain->GetN(); samp++){
	    // std::cout << gr0Gain->GetX()[samp] << ", ";
	    if(gr0Gain->GetX()[samp] >= 200 && gr0Gain->GetX()[samp] < 1200){
	      mean += gr0Gain->GetY()[samp];
	      numSamps++;
	    }
	  }
	  // std::cout << std::endl << mean << "\t" << numSamps << std::endl;

	  mean /= numSamps;

	  // std::cout << antNumber << "\t" << az << "\t" << el << "\t" << mean << "\t" << std::endl;

	  if(az==0 && el==0){
	    theBoresightVal = mean;
	  }
	  if(el==0){
	    grAz->SetPoint(grAz->GetN(), az*TMath::DegToRad(), mean);
	  }
	  if(az==0){
	    grEl->SetPoint(grEl->GetN(), el*TMath::DegToRad(), mean);
	  }

	  if(mean < minVal){
	    minVal = mean;
	    // minVal = 0;	    
	  }
	  
	  
	  grMeanGain2D->SetPoint(grMeanGain2D->GetN(), az, el, mean);
	  
	  delete gr0;
	}
      }
      
      TString gr2DName = TString::Format("grMeanGain2D_%d_%d", antNumber, polInd);
      grMeanGain2D->SetName(gr2DName);

      for(int samp=0; samp < grMeanGain2D->GetN(); samp++){
	// grMeanGain2D->GetZ()[samp] -= theBoresightVal;
	grMeanGain2D->GetZ()[samp] -= minVal;
      }
      grMeanGain2D->Write();
      delete grMeanGain2D;

      TString grAzName = TString::Format("grAz_%d_%d", antNumber, polInd);
      grAz->SetName(grAzName);

      for(int samp=0; samp < grAz->GetN(); samp++){
	// grAz->GetY()[samp] -= theBoresightVal;
	grAz->GetY()[samp] -= minVal;	
      }
      grAz->Write();
      delete grAz;

      TString grElName = TString::Format("grEl_%d_%d", antNumber, polInd);
      grEl->SetName(grElName);

      for(int samp=0; samp < grEl->GetN(); samp++){
	// grEl->GetY()[samp] -= theBoresightVal;
	grEl->GetY()[samp] -= minVal;
      }
      grEl->Write();
      delete grEl;
    }
  }
  
  fOut->Write();
  fOut->Close();
}



std::vector<Double_t> friisCorrection(const std::vector<Double_t>& ps,
				      const std::vector<Double_t>& freqs,
				      AnitaPol::AnitaPol_t pol){

  std::vector<Double_t> gain_dBi(ps.size(), 0);
  
  // separationFactors = [(f*4*math.pi*distMeters/c)**2 for f in freqsMHz]
 
  for(UInt_t i=0; i < ps.size(); i++){

    Double_t f = freqs.at(i);
    Double_t friisFactor = (f*f*partFactor*partFactor);
    
    // pt/pr = GrGt(f*f*partFactor*partFactor)
    // for vpol Gr==Gt so do sqrt((pt/pr)/(f*f*partFactor*partFactor))
    // for hpol Gr==Gt so do (pt/pr)/(f*f*partFactor*partFactor)/Gt (which is the vpol Gr==Gt)
    if(pol==AnitaPol::kVertical){
      gain_dBi.at(i) = TMath::Sqrt(ps[i]*friisFactor);
    }
    else{
      gain_dBi.at(i) = ps[i]*friisFactor/meanVPolGain.at(i);
    }    
  }
  return gain_dBi;
}
