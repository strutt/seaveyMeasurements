/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Reads CSV files from the scope and gives you TGraphs
***********************************************************************************************************/

#ifndef SEAVEYCSVPARSER_H
#define SEAVEYCSVPARSER_H

// ROOT things
#include "TGraph.h"
#include "TObjString.h"
#include "TString.h"
#include "TObjArray.h"
#include "TFile.h"

// standard c++ things
#include <iostream>
#include <fstream>

#include "AnitaConventions.h"
#include "FFTtools.h"

class SeaveyDataHandler{
  
public:
  SeaveyDataHandler();
  ~SeaveyDataHandler();
  
  static TGraph* getGraphFromCsvFile(const TString& fileName);
  
  TGraph* getBoresightGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol);
  TGraph* getAttenuatedGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol);
  TGraph* getOffAxisGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol, Int_t az, Int_t el);
  TGraph* getOffAxisOffGraphFromTFile(Int_t antNumber, Int_t channel, AnitaPol::AnitaPol_t pol, Int_t az, Int_t el);

  void doNoiseSubtraction(TGraph* gr, Int_t antNumber, Int_t channel,
			  AnitaPol::AnitaPol_t  pol);
  void windowPulse(TGraph* gr, Double_t timeBeforePeak, Double_t timeAfterPeak);
  
  Bool_t kPrintWarnings;
  Bool_t kSubtractNoiseHistogram;
  
private:
  TFile* kFile;

};

#endif

