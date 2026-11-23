#ifndef __EcalCellPar__
#define  __EcalCellPar__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

struct Params {
  Double_t tdc0 = 1, tdc1 = 0;
  Double_t adc0 = 0, adc1 = 1, adc2 = 0;
  Double_t twc0 = 0, twc1 = 0, twc2 = 0;
  Double_t twc3 = 0, twc4 = 0;
};


istream& operator >> (istream & is, Params & par){
  if (is) is >> par.tdc0 >> par.tdc1 >> par.adc0 >> par.adc1 >> par.adc2 >> par.twc0 >> par.twc1 >> par.twc2;
  return is;
}

ostream& operator << (ostream & is, Params & par){
  if (is) {
    is << setw(15) << par.tdc0 << "\t"
    << setw(15) << par.tdc1 << "\t"
    << setw(15) << par.adc0 << "\t"
    << setw(15) << par.adc1 << "\t"
    << setw(15) << par.adc2 << "\t"
    << setw(15) << par.twc0 << "\t"
    << setw(15) << par.twc1 << "\t"
    << setw(15) << par.twc2;
  }
  return is;
  }


  struct EcalCellPar {
  EcalCellPar(){ }
  EcalCellPar(const Params & _oldParams, const Params & _newParams) : oldParams(_oldParams), newParams(_newParams){ }
  Params oldParams, newParams;
  };


void readParams(const TString &inputFile, EcalCellPar ecalCellPars[6][163], Option_t * option = "old")
{
  TString opt = option;
  opt.ToLower();
  ifstream inFile(inputFile);
  if (inFile.fail()) {
    cout << "-E- main : File " << inputFile << " does not exist!" << endl;
    return;
  }
  string line;
  do {
    getline(inFile, line);
  } while (line != "[EmcCalPar]");
  Params params;
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      inFile >> sector >> cell >> params;
      cout<<"Sector "<<sector<<" Cell "<<cell<<" "<<params<<endl;
      if (opt.Contains("old")) ecalCellPars[sector][cell].oldParams = params;
      ecalCellPars[sector][cell].newParams = params;
      if (inFile.eof()) break;
    }
  }
}

void writeParams(const TString &outputFile, EcalCellPar ecalCellPars[6][163])
{
  ofstream outFile(outputFile);
  if (outFile.fail()) {
    cout << "-E- main : File " << outputFile << " cannot be written!" << endl;
    return;
  }
  outFile << "##############################################################################\n"
          << "##############################################################################\n"
          << "# Calibration parameters for the Emc\n"
          << "# Format:\n"
          << "# sector(0...5)   position(0...162)   TdcSlope  TdcOffset AdcOffset1\n"
          << "# AdcOffset2 AdcSlope   Twc0  Twc1  Twc2 Twc3 Twc4\n"
          << "##############################################################################\n"
          << "[EmcCalPar]\n";
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      outFile << sector << "\t" << cell << "\t" << ecalCellPars[sector][cell].newParams << "\n";
    }
  }
}

#endif /* ifndef __EcalCellPar__ */
