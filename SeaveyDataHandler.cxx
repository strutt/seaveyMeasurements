#include "SeaveyDataHandler.h"

SeaveyDataHandler::SeaveyDataHandler(){
  kFile = TFile::Open("seaveyDataAsTGraphs.root");
  kPrintWarnings = false;
  if(kFile==NULL){
    std::cerr << "Warning! Unable to open TFile with TGraphs!" << std::endl;
  }
}


SeaveyDataHandler::~SeaveyDataHandler(){
  if(kFile!=NULL){
    kFile->Close();
  }
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
