#include "SeaveyDataHandler.h"
#include "TFile.h"

int main(){

  TFile* fOut = new TFile("seaveyDataAsTGraphs.root", "recreate");

  TString dataDir = "seaveyDataPalestine2014/S21s/";

  const Int_t prefixLengthToRemove = 13;
  const Int_t suffixLengthToRemove = 4;

  std::ifstream listOfFiles("listOfCsvFiles.txt");

  Int_t numFiles = 0;
  while(!listOfFiles.eof()){
    char fileName[1024];
    listOfFiles >> fileName;
    std::cout << numFiles << "\t" << fileName << std::endl;

    TString fileName2 = TString::Format("%s", fileName);
    TGraph* gr = SeaveyDataHandler::getGraphFromCsvFile(dataDir + fileName2);

    Int_t newStrLength = fileName2.Length();
    newStrLength = newStrLength - suffixLengthToRemove - prefixLengthToRemove;
    TString grName = TString(fileName2(prefixLengthToRemove, newStrLength));
    grName = TString("gr") + grName;
    std::cerr << grName.Data() << std::endl;
    gr->SetName(grName);
    gr->Write();

    numFiles++;    
  }

 
  fOut->Write();
  fOut->Close();
}
