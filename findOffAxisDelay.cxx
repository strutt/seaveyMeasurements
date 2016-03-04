#include "SeaveyDataHandler.h"
#include "TH2D.h"
#include "FFTtools.h"

int main(){

  SeaveyDataHandler sData;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  Int_t antNum = 30;


  // Channel 1 is aligned, channel 4 is cross-pol
  Int_t channel = 1;

  const Double_t timeBefore = 5e-9;
  // const Double_t timeAfter = 50e-9;
  // const Double_t timeAfter = 30e-9;
  const Double_t timeAfter = 20e-9;

  TFile* fOut = new TFile("findOffAxisDelayPlots.root", "recreate");
  
  TGraph* gr0 = sData.getBoresightGraphFromTFile(antNum, 1, pol);

  const Double_t dt = gr0->GetX()[1] - gr0->GetX()[0];
  const Double_t dtInterp = 0.1*dt;

  gr0->SetTitle("Before");
  gr0->SetName("grBefore");  
  gr0->Write();
  sData.doNoiseSubtraction(gr0, antNum, channel, pol);

  gr0->SetTitle("Noise subbed");
  gr0->SetName("grNoiseSubbed");  
  gr0->Write();  
  sData.windowPulse(gr0, timeBefore, timeAfter);
  gr0->SetTitle("Final");
  gr0->SetName("grFinal");  
  gr0->Write();  

  
  const Int_t azRange = 45;
  const Int_t deltaAz = 15;
  const Int_t elRange = 45;
  const Int_t deltaEl = 15;

  TH2D* hGroupDelay = new TH2D("hGroupDelay", "Group delay; Azimuth (Degrees); Elevation (Degrees); Group delay (ns)", 100, -50, 50, 100, -50, 50);

  for(int az=-azRange; az <= azRange; az += deltaAz){
    for(int el=-elRange; el <= elRange; el += deltaEl){    
      TGraph* gr = sData.getOffAxisGraphFromTFile(antNum, channel, pol, az, el);

      if(gr){
	sData.doNoiseSubtraction(gr, antNum, channel, pol);
	sData.windowPulse(gr, timeBefore, timeAfter);
	
	// TGraph* grCor = FFTtools::getCorrelationGraph(gr, gr0);
	TGraph* grCor = FFTtools::getInterpolatedCorrelationGraph(gr, gr0, dtInterp);

	Int_t peakSamp = FFTtools::getPeakBin(grCor);
	const Double_t delay = grCor->GetX()[peakSamp];

	// std::cout << az << "\t" << el << "\t" << delay << std::endl;

	hGroupDelay->Fill(az, el, 1e9*delay);

	TString name = TString::Format("gr%d_%d_%d", antNum, az, el);
	gr->SetName(name);
	gr->Write();

	TString corName = TString::Format("grCor%d_%d_%d", antNum, az, el);
	grCor->SetName(corName);
	grCor->Write();
	
	delete gr;
	delete grCor;
      }
    }
  }
    
  fOut->Write();
  fOut->Close();
}
