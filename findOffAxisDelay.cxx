#include "SeaveyDataHandler.h"
#include "TH2D.h"
#include "FFTtools.h"
#include "ProgressBar.h"
#include "TGraph2D.h"

int main(){

  SeaveyDataHandler sData;

  // Channel 1 is aligned, channel 4 is cross-pol
  Int_t channel = 1;

  const Double_t timeBefore = 5e-9; //2e-9;
  const Double_t timeAfter = 30e-9;

  // const Double_t timeBefore = 2e-9;
  // const Double_t timeAfter = 5e-9;

  TFile* fOut = new TFile("findOffAxisDelayPlots.root", "recreate");

  const Int_t numAnts = 4;
  const Int_t antNums[numAnts] = {30, 29, 26, 25};

  const Int_t azRange = 45;
  const Int_t deltaAz = 15;
  const Int_t elRange = 45;
  const Int_t deltaEl = 15;

  const Int_t nAnglesAz = 2*(azRange/deltaAz)+1;
  const Int_t nAnglesEl = 2*(elRange/deltaEl)+1;  
  ProgressBar p(nAnglesAz*nAnglesEl*numAnts);
  
  for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){

    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    
    for(int antInd=0; antInd < numAnts; antInd++){

      Int_t antNum = antNums[antInd];
      
      TGraph* gr0 = sData.getBoresightGraphFromTFile(antNum, channel, pol);

      const Double_t dt = gr0->GetX()[1] - gr0->GetX()[0];
      const Double_t dtInterp = 0.1*dt;

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

  

      // TH2D* hGroupDelay = new TH2D("hGroupDelay", "Group delay; Azimuth (Degrees); Elevation (Degrees); Group delay (ns)", 100, -50, 50, 100, -50, 50);

      // const int numBins = 100;
      // const int histAngleRange = 50;
      
      TString gr2dName = TString::Format("grGroupDelay_%d_%d", antNum, polInd);
      TString gr2dTitle = TString::Format("Measured Group delay RXP %d ", antNum);
      gr2dTitle += pol == AnitaPol::kHorizontal ? "HPOL" : "VPOL";
      gr2dTitle += "; Azimuth (Degrees); Elevation (Degrees); Group delay (ns)";

      // TH2D* hGroupDelay = new TH2D(gr2dName, gr2dTitle,
      // 				   numBins, -histAngleRange, histAngleRange,
      // 				   numBins, -histAngleRange, histAngleRange);
      TGraph2D* grGroupDelay = new TGraph2D();
      grGroupDelay->SetName(gr2dName);
      grGroupDelay->SetTitle(gr2dTitle);
      
      for(int az=-azRange; az <= azRange; az += deltaAz){
	for(int el=-elRange; el <= elRange; el += deltaEl){    
	  TGraph* gr = sData.getOffAxisGraphFromTFile(antNum, channel, pol, az, el);

	  if(gr){
	    sData.doNoiseSubtraction(gr, antNum, channel, pol);


	    // // Hack for 26 HPOL +30 el, where the noise from the pulse box is bigger than the peak!!!
	    // const Double_t approxPulseTime = 0.15e-6; // Actually around 0.16e-6 so this is conservative
	    // for(int samp=0; samp < gr->GetN(); samp++){
	    //   if(gr->GetX()[samp] < approxPulseTime){
	    // 	gr->GetY()[samp] = 0;
	    //   }
	    // }
	    sData.windowPulse(gr, timeBefore, timeAfter);

	    // invert negative elevations to account for making measurement with flipped antenna
	    if(el < 0){
	      for(int samp=0; samp < gr0->GetN(); samp++){
		gr->GetY()[samp]*=-1;
	      }
	    }


	    
	    // TGraph* grCor = FFTtools::getCorrelationGraph(gr, gr0);
	    TGraph* grCor = FFTtools::getInterpolatedCorrelationGraph(gr, gr0, dtInterp);

	    Int_t peakSamp = FFTtools::getPeakBin(grCor);
	    const Double_t delay = grCor->GetX()[peakSamp];

	    // std::cout << az << "\t" << el << "\t" << delay << std::endl;

	    grGroupDelay->SetPoint(grGroupDelay->GetN(), az, el, 1e9*delay);

	    TString name = TString::Format("gr%d_%d_%d_%d", antNum, az, el, polInd);
	    gr->SetName(name);
	    gr->Write();

	    for(int samp=0; samp < gr->GetN(); samp++){
	      gr->GetX()[samp] -= delay;
	    }
	    name += "_shifted";
	    gr->SetName(name);
	    gr->Write();	    

	    TString corName = TString::Format("grCor%d_%d_%d_%d", antNum, az, el, polInd);
	    grCor->SetName(corName);
	    grCor->Write();
	
	    delete gr;
	    delete grCor;
	  }
	  p++;
	}
      }
      grGroupDelay->Write();
      delete grGroupDelay;
    }
  }  
  fOut->Write();
  fOut->Close();
}
