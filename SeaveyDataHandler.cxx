#include "SeaveyDataHandler.h"



SeaveyDataHandler::SeaveyDataHandler(){
  zeroPointers();
  makeCableResponses(0);
}


SeaveyDataHandler::SeaveyDataHandler(int padLength){
  zeroPointers();
  makeCableResponses(padLength);
}




void SeaveyDataHandler::zeroPointers(){
  
  kFile = TFile::Open("seaveyDataAsTGraphs.root");
  kPrintWarnings = false;
  if(kFile==NULL){
    std::cerr << "Warning! Unable to open TFile with TGraphs!" << std::endl;
  }

  grPsPulserDirectFast = NULL;
  grPsPulserCopolFast5ft = NULL;
  grPsPulserXpolFast5ft = NULL;
  grPsPulserCopolFast = NULL;
  grPsPulserXpolFast = NULL;

  fftwComplexPsPulserDirectFast = NULL;
  fftwComplexPsPulserCopolFast5ft = NULL;
  fftwComplexPsPulserXpolFast5ft = NULL;
  fftwComplexPsPulserCopolFast = NULL;
  fftwComplexPsPulserXpolFast = NULL;

  fftwComplexCopolCableResponse = NULL;
  fftwComplexXpolCableResponse = NULL;
  fftwComplexPulseFreqs = NULL;
  

}
SeaveyDataHandler::~SeaveyDataHandler(){


  if(grPsPulserDirectFast != NULL){
    delete grPsPulserDirectFast;
    grPsPulserDirectFast = NULL;    
  }
  if(grPsPulserCopolFast5ft != NULL){
    delete grPsPulserCopolFast5ft;
    grPsPulserCopolFast5ft = NULL;    
  }
  if(grPsPulserXpolFast5ft != NULL){
    delete grPsPulserXpolFast5ft;
    grPsPulserXpolFast5ft = NULL;    
  }
  if(grPsPulserCopolFast != NULL){
    delete grPsPulserCopolFast;
    grPsPulserCopolFast = NULL;    
  }
  if(grPsPulserXpolFast != NULL){
    delete grPsPulserXpolFast;
    grPsPulserXpolFast = NULL;    
  }  
  
  if(fftwComplexPsPulserDirectFast!=NULL){
    delete [] fftwComplexPsPulserDirectFast;
    fftwComplexPsPulserDirectFast = NULL;
  }
  if(fftwComplexPsPulserCopolFast5ft!=NULL){
    delete [] fftwComplexPsPulserCopolFast5ft;
    fftwComplexPsPulserCopolFast5ft = NULL;
  }
  if(fftwComplexPsPulserXpolFast5ft!=NULL){
    delete [] fftwComplexPsPulserXpolFast5ft;
    fftwComplexPsPulserXpolFast5ft = NULL;
  }
  if(fftwComplexPsPulserCopolFast!=NULL){
    delete [] fftwComplexPsPulserCopolFast;
    fftwComplexPsPulserCopolFast = NULL;
  }
  if(fftwComplexPsPulserXpolFast!=NULL){
    delete [] fftwComplexPsPulserXpolFast;
    fftwComplexPsPulserXpolFast = NULL;
  }

  if(fftwComplexCopolCableResponse != NULL){
    delete [] fftwComplexCopolCableResponse;
    fftwComplexCopolCableResponse = NULL;
  }
  if(fftwComplexXpolCableResponse != NULL){
    delete [] fftwComplexXpolCableResponse;
    fftwComplexXpolCableResponse = NULL;
  }

  if(fftwComplexPulseFreqs != NULL){
    delete [] fftwComplexPulseFreqs;
    fftwComplexPulseFreqs = NULL;
  }
  
  
  if(kFile!=NULL){
    kFile->Close();
  }

  
}


void SeaveyDataHandler::removeAttenuationTimeDomain(TGraph* gr, Double_t attendB){

    // # e.g. 20dB in power is 10dB in voltage
    // atten_v = atten_dB/2.
    // mult = pow(10, atten_v/10.)
    // vals = [v*mult for v in vals]
    // return vals

  // 20dB in power is 10dB in voltage...
  Double_t attenVolts = attendB/2;
  Double_t factor = pow(10, attenVolts/2);
  
  for(int samp=0; samp < gr->GetN(); samp++){
    gr->GetY()[samp] *= factor;
  }  
}





/** 
 * @brief Removes the factor of n picked up by a forwards and then inverse FFT on the forwards step.
 * 
 * @param n is the number of samples in the y array
 * @param y pointer to an array of length n
 * 
 * @return pointer to the FFTWComplex array, normalized such that the magntitude of each element is scaled by 1/n.
 */
FFTWComplex* SeaveyDataHandler::doNormalizedFFT(Int_t n, Double_t* y){

  FFTWComplex* theFFT = FFTtools::doFFT(n, y);
  Int_t nf = (n/2)+1;
  for(int freqInd=0; freqInd < nf; freqInd++){
    double phase = theFFT[freqInd].getPhase();
    double mag = TMath::Sqrt(theFFT[freqInd].getAbs());
    double newMag = mag/n;
    theFFT[freqInd].setMagPhase(newMag, phase);
  }
  return theFFT;
}



void SeaveyDataHandler::makeCableResponses(int padLength){

  grPsPulserDirectFast = (TGraph*) kFile->Get("gr_ps_pulser_direct_fast_Ch1");
  grPsPulserCopolFast5ft = (TGraph*) kFile->Get("gr_ps_pulser_copol_fast_5ft_Ch1");
  grPsPulserXpolFast5ft = (TGraph*) kFile->Get("gr_ps_pulser_xpol_fast_5ft_Ch1");
  grPsPulserCopolFast = (TGraph*) kFile->Get("gr_ps_pulser_copol_fast_Ch1");
  grPsPulserXpolFast = (TGraph*) kFile->Get("gr_ps_pulser_xpol_fast_Ch1");

  // put into an array for looping
  const int numCalibGraphs = 5;
  TGraph* theGrs[numCalibGraphs] = {grPsPulserDirectFast,
				    grPsPulserCopolFast5ft,
				    grPsPulserXpolFast5ft,
				    grPsPulserCopolFast,
				    grPsPulserXpolFast};

  dtCables = theGrs[0]->GetX()[1] - theGrs[0]->GetX()[0];

  if(padLength > 0){
    for(int i=0; i < numCalibGraphs; i++){
      if(padLength > theGrs[i]->GetN()){
	// std::cerr << "padding..." << std::endl;
	while(theGrs[i]->GetN() != padLength){
	  theGrs[i]->SetPoint(theGrs[i]->GetN(), theGrs[i]->GetX()[0] + dtCables*theGrs[i]->GetN(), 0);
	}
      }
      else if(padLength < theGrs[i]->GetN()){
	std::cerr << "Removing points... you should probably check this is ok!" << std::endl;
	while(theGrs[i]->GetN() != padLength){
	  theGrs[i]->RemovePoint(theGrs[i]->GetN()-1);
	}
      }
    }
  }
  
  // When we took these data we added 20dB of attenuation to the ps pulser.
  // Here we remove that attenuation  

  nPointsCables = theGrs[0]->GetN();
  deltaF = 1./(dtCables*nPointsCables);
  numFreqs = (nPointsCables/2) + 1;
  
  for(int i=0; i < numCalibGraphs; i++){
    removeAttenuationTimeDomain(theGrs[i], 20);

    Double_t thisDt = theGrs[i]->GetX()[1] - theGrs[i]->GetX()[0];
    Int_t thisN = theGrs[i]->GetN();

    if(thisDt!=dtCables || nPointsCables != thisN){
      std::cerr << theGrs[i]->GetName() << "\t" << thisN << "\t" << thisDt << std::endl;
    }
  }
  
  // Convert to frequency domain.
  fftwComplexPsPulserDirectFast = doNormalizedFFT(nPointsCables, grPsPulserDirectFast->GetY());
  fftwComplexPsPulserCopolFast5ft = doNormalizedFFT(nPointsCables, grPsPulserCopolFast5ft->GetY());
  fftwComplexPsPulserXpolFast5ft = doNormalizedFFT(nPointsCables, grPsPulserXpolFast5ft->GetY());
  fftwComplexPsPulserCopolFast = doNormalizedFFT(nPointsCables, grPsPulserCopolFast->GetY());
  fftwComplexPsPulserXpolFast = doNormalizedFFT(nPointsCables, grPsPulserXpolFast->GetY());


  // now make the Copol cable response
  fftwComplexCopolCableResponse = new FFTWComplex[numFreqs];
  for(int freqInd = 0; freqInd < numFreqs; freqInd++){
    fftwComplexCopolCableResponse[freqInd] = fftwComplexPsPulserCopolFast5ft[freqInd]/fftwComplexPsPulserDirectFast[freqInd];
  }

  // now make the XPol cable response
  fftwComplexXpolCableResponse = new FFTWComplex[numFreqs];
  for(int freqInd = 0; freqInd < numFreqs; freqInd++){
    fftwComplexXpolCableResponse[freqInd] = fftwComplexPsPulserXpolFast5ft[freqInd]/fftwComplexPsPulserDirectFast[freqInd];
  }

  // now get the raw pulser frequency content
  fftwComplexPulseFreqs = new FFTWComplex[numFreqs];
  for(int freqInd = 0; freqInd < numFreqs; freqInd++){
    fftwComplexPulseFreqs[freqInd] = fftwComplexPsPulserCopolFast[freqInd]/fftwComplexCopolCableResponse[freqInd];

    if(fftwComplexCopolCableResponse[freqInd].getAbsSq() == 0){
      std::cerr << freqInd << "\t" << fftwComplexCopolCableResponse[freqInd].re << "\t" << fftwComplexCopolCableResponse[freqInd].im << std::endl;
    }
  }
}






FFTWComplex* SeaveyDataHandler::removeCopolResponse(TGraph* gr){
    // def removeResponse(self, wave, dtNs, polKey):
    //     # In order to have the same df, we need to have 
    //     # N_1*dt_1 == N_2*dt_2. If that's not the case, we 
    //     # need to take action.

    //     numZerosToPad = len(self.waves[polKey])*self.dts[polKey]/dtNs - len(wave)
    //     assert numZerosToPad == int(numZerosToPad) # If it's not an integer, things will be tricky.
    //     if numZerosToPad >=0:
    //         for zeroInd in xrange(int(numZerosToPad)):
    //             wave.append(0)
    //     else:
    //         print 'Warning! Deleting things! You should probably check this is OK!'
    //         for zeroInd in xrange(int(abs(numZerosToPad))):
    //             wave.pop()

    //     # So now the frequencies should be the same
    //     assert 1e3/(len(wave)*dtNs) == 1e3/(len(self.waves[polKey])*self.dts[polKey])

    //     fft_wave = doNormalizedFFT(wave, dtNs)
    //     withoutCables = deconvolveFreqToFreq(fft_wave, self.responses[polKey])
    //     justSeaveyToSeavey = deconvolveFreqToFreq(withoutCables, self.pulseFreqs)
    //     return justSeaveyToSeavey

  Int_t n = gr->GetN();
  Double_t dt = gr->GetX()[1] - gr->GetX()[0];

  // need same df, which means I need n1*dt1 = n2*dt2...
  Double_t nPadDouble = (nPointsCables*dtCables)/dt;

  Int_t nPad = TMath::Nint(nPadDouble);

  std::cerr << n << "\t" << dt << "\t" << nPointsCables << "\t" << dtCables << "\t" << nPad << std::endl;

  std::cerr << nPad*dt << "\t" << nPointsCables*dtCables << std::endl;
  
  if(TMath::Abs(nPad - nPadDouble) > 0.00001){
    std::cerr << "Warning! I can't figure out how to pad the input to match frequency domains... in " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
  }

  if(nPad > n){
    // std::cerr << "padding..." << std::endl;
    while(gr->GetN() != nPad){
      gr->SetPoint(gr->GetN(), gr->GetX()[0] + dt*gr->GetN(), 0);
    }
  }
  else if(nPad < n){
    std::cerr << "Removing points... you should probably check this is ok!" << std::endl;
    while(gr->GetN() != nPad){
      gr->RemovePoint(gr->GetN()-1);
    }
  }

  FFTWComplex* fftGr = doNormalizedFFT(gr->GetN(), gr->GetY());
  // Int_t nFFT = (gr->GetN()/2) + 1;

  FFTWComplex* deconvolved = new FFTWComplex[numFreqs];

  for(int freqInd=0; freqInd < numFreqs; freqInd++){
    // std::cerr << freqInd << " (" << fftGr[freqInd].re << ", " << fftGr[freqInd].im << ")";    
    deconvolved[freqInd] = fftGr[freqInd]/fftwComplexCopolCableResponse[freqInd];
    // std::cerr << freqInd << " (" << fftGr[freqInd].re << ", " << fftGr[freqInd].im << ")";
    // std::cerr << std::endl;    
  }

  // withoutCables = deconvolveFreqToFreq(fft_wave, self.responses[polKey])
  // justSeaveyToSeavey = deconvolveFreqToFreq(withoutCables, self.pulseFreqs)

  for(int freqInd=0; freqInd < numFreqs; freqInd++){
    // std::cerr << freqInd << " (" << fftGr[freqInd].re << ", " << fftGr[freqInd].im << ")";    
    deconvolved[freqInd] /= fftwComplexPulseFreqs[freqInd];
    // std::cerr << freqInd << " (" << fftGr[freqInd].re << ", " << fftGr[freqInd].im << ")";
    // std::cerr << std::endl;
  }

  delete [] fftGr;
  
  return deconvolved;  
}














TGraph* SeaveyDataHandler::getGraphFromCsvFile(const TString& fileName){

  std::ifstream inFile(fileName.Data());
  const Int_t numLinesInHeader = 6;
  const Int_t maxLineSize = 1024;
  
  for(int line=0; line < numLinesInHeader; line++){
    char lineText[maxLineSize];
    inFile.getline(lineText,maxLineSize-1);
  }

  std::vector<Double_t> tVals;
  std::vector<Double_t> vVals;
  
  while(!inFile.eof()){
    char lineText[1024];
    inFile.getline(lineText,1023);
    TString lt = TString::Format("%s", lineText);

    // Fuck you strings and fuck you ROOT
    TObjArray* strVals = lt.Tokenize(",");
    for(int strInd = 0; strInd < strVals->GetEntries(); strInd++){
      TString theVal = ((TObjString*) strVals->At(strInd))->GetString();
      if(strInd==0){
	tVals.push_back(atof(theVal.Data()));
      }
      else{
	vVals.push_back(atof(theVal.Data()));
      }
    }
  }

  TGraph* gr = new TGraph(tVals.size(), &tVals[0], &vVals[0]);
  return gr;
}


TGraph* SeaveyDataHandler::getBoresightGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol){

  TString polStr = pol == AnitaPol::kHorizontal ? "hpol" : "vpol";
  
  TString grName = TString::Format("gr_rxp%d_%s_Ch%d",
				   antNumber,
				   polStr.Data(),
				   channel);

  TGraph* gr = (TGraph*) kFile->Get(grName);

  if(gr){
    TString title = TString::Format("RXP %d %s ch. %d", antNumber, polStr.Data(), channel);
    gr->SetTitle(title);
  }
  else{
    if(kPrintWarnings){    
      std::cerr << "Warning! Couldn't find TGraph named " << grName.Data() << " in " << kFile->GetName() << std::endl;
    }
  }

  return gr;
}


TGraph* SeaveyDataHandler::getAttenuatedGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol){

  TString polStr = pol == AnitaPol::kHorizontal ? "hpol" : "vpol";
  
  TString grName = TString::Format("gr_rxp%d_%s_60dBatten_Ch%d",
				   antNumber,
				   polStr.Data(),
				   channel);

  TGraph* gr = (TGraph*) kFile->Get(grName);

  if(gr){
    TString title = TString::Format("60db Attenuated RXP %d %s ch. %d", antNumber, polStr.Data(), channel);
    gr->SetTitle(title);  
  }
  else{
    if(kPrintWarnings){    
      std::cerr << "Warning! Couldn't find TGraph named " << grName.Data() << " in " << kFile->GetName() << std::endl;
    }
  }
  
  return gr;
}


TGraph* SeaveyDataHandler::getOffAxisGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol,
						  Int_t az, Int_t el){

  if(az==0 && el==0){
    return getBoresightGraphFromTFile(antNumber, channel, pol);
  }

  TString polStr = pol == AnitaPol::kHorizontal ? "hpol" : "vpol";

  // There should be a better way to do this... but I want it quickly
  TString elStr;
  if(el >=0 && el < 10){
    elStr = TString::Format("+0%d", el);
  }
  else if(el < 0 && el > -10){
    elStr = TString::Format("%+0d", el);
  }
  else{
    elStr = TString::Format("%+d", el);
  }

  TString azStr;
  if(az >=0 && az < 10){
    azStr = TString::Format("+0%d", az);
  }
  else if(az < 0 && az > -10){
    azStr = TString::Format("%d", az);    
  }
  else{
    azStr = TString::Format("%+d", az);
  }
  
  TString grName = TString::Format("gr_rxp%d_%s_%saz%sel_Ch%d",
				   antNumber,
				   polStr.Data(),
				   azStr.Data(),
				   elStr.Data(),
				   channel);

  TGraph* gr = (TGraph*) kFile->Get(grName);

  if(gr){
    TString title = TString::Format("RXP %d %s ch. %d Off-axis at Az %d El %d Degrees", antNumber, polStr.Data(), channel, az, el);
    gr->SetTitle(title);  
  }
  else{
    if(kPrintWarnings){    
      std::cerr << "Warning! Couldn't find TGraph named " << grName.Data() << " in " << kFile->GetName() << std::endl;
    }
  }
  return gr;
}


TGraph* SeaveyDataHandler::getOffAxisOffGraphFromTFile(Int_t antNumber, Int_t channel,
						     AnitaPol::AnitaPol_t pol,
						     Int_t az, Int_t el){

  TString polStr = pol == AnitaPol::kHorizontal ? "hpol" : "vpol";

  // There should be a better way to do this... but I want it quickly
  TString elStr;
  if(el >=0 && el < 10){
    elStr = TString::Format("+0%d", el);
  }
  else if(el < 0 && el > -10){
    elStr = TString::Format("%+0d", el);
  }
  else{
    elStr = TString::Format("%+d", el);
  }

  TString azStr;
  if(az >=0 && az < 10){
    azStr = TString::Format("+0%d", az);
  }
  else if(az < 0 && az > -10){
    azStr = TString::Format("%d", az);    
  }
  else{
    azStr = TString::Format("%+d", az);
  }
  
  TString grName = TString::Format("gr_rxp%d_%s_%saz%sel_off_Ch%d",
				   antNumber,
				   polStr.Data(),
				   azStr.Data(),
				   elStr.Data(),
				   channel);

  TGraph* gr = (TGraph*) kFile->Get(grName);

  if(gr){
    TString title = TString::Format("RXP %d %s ch. %d Off-axis Off at Az %d El %d Degrees", antNumber, polStr.Data(), channel, az, el);
    gr->SetTitle(title);  
  }
  else{
    if(kPrintWarnings){
      std::cerr << "Warning! Couldn't find TGraph named " << grName.Data() << " in " << kFile->GetName() << std::endl;
    }
  }
  return gr;
}


void SeaveyDataHandler::doNoiseSubtraction(TGraph* gr, Int_t antNumber, Int_t channel,
					   AnitaPol::AnitaPol_t  pol){
    
  // Let's assume everything is well behaved for now...
  TGraph* grNoise = getAttenuatedGraphFromTFile(antNumber, channel, pol);

  for(int samp=0; samp < gr->GetN(); samp++){
    gr->GetY()[samp] -= grNoise->GetY()[samp];
  }

  delete grNoise;

}


void SeaveyDataHandler::windowPulse(TGraph* gr, Double_t timeBeforePeak, Double_t timeAfterPeak){
    
  // Let's assume everything is well behaved for now...
  Int_t peakSamp = FFTtools::getPeakBin(gr);
  Double_t peakTime = gr->GetX()[peakSamp];
  
  for(int samp=0; samp < gr->GetN(); samp++){

    if(peakTime - gr->GetX()[samp] > timeBeforePeak || gr->GetX()[samp] - peakTime  > timeAfterPeak){
      gr->GetY()[samp] = 0;
    }
  }
}
