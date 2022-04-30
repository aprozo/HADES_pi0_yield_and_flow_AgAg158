

template <typename Cell>
using Array = vector<vector<Cell> >;

struct Params {
  Double_t tdc0 = 1, tdc1 = 0;
  Double_t adc0 = 0, adc1 = 1, adc2 = 0;
  Double_t twc0 = 0, twc1 = 0, twc2 = 0;
};

class EcalCellPar {
public:
  EcalCellPar(const Int_t &newSector, const Int_t &newCell):
    sector(newSector),
    cell(newCell)
  {}
  Params oldParams, newParams;
private:
  Int_t sector, cell;
}; // end of a class

istream& operator >> (istream &is, Params &par)
{
  if (is) {
    is >> par.tdc0 >> par.tdc1 >> par.adc0 >> par.adc1 >> par.adc2 >> par.twc0 >> par.twc1 >> par.twc2;
  }
  return is;
}

Double_t twc3 = 0, twc4 = 0;
ostream& operator << (ostream &is, Params &par)
{
  if (is) {
    is << setw(10) << par.tdc0 << "\t"
       << setw(10)
       << par.tdc1 << "\t"
       << setw(10)
       << par.adc0 << "\t"
       << setw(10)
       << par.adc1 << "\t"
       << setw(10)
       << par.adc2 << "\t"
       << setw(10)
       << par.twc0 << "\t"
       << setw(10)
       << par.twc1 << "\t"
       << setw(10)
       << par.twc2;
  }
  return is;
}

void readParams(const TString &inputFile, Array<EcalCellPar> &ecalCells)
{
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
      ecalCells[sector][cell].oldParams = params;
      ecalCells[sector][cell].newParams = params;
      // cout<<sector<<"\t"<<cell<<"\t"<<ecalCells[sector][cell].oldParams<<endl;
      if (inFile.eof())
        break;
    }
  }
}

void readNewParams(const TString &inputFile, Array<EcalCellPar> &ecalCells)
{
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
      ecalCells[sector][cell].newParams = params;
      //cout<<sector<<"\t"<<cell<<"\t"<<ecalCells[sector][cell].oldParams<<endl;
      if (inFile.eof())
        break;
    }
  }
}

void writeParams(const TString &outputFile, Array<EcalCellPar> &ecalCells)
{
  ofstream outFile(outputFile);

  if (outFile.fail()) {
    cout << "-E- main : File " << outputFile << " cannot be written!" << endl;
    return;
  }
  outFile << "##############################################################################\n"
          << "##############continue;################################################################\n"
          << "# Calibration parameters for the Emc\n"
          << "# Format:\n"
          << "# sector(0...5)   position(0...162)   TdcSlope  TdcOffset AdcOffset1\n"
          << "# AdcOffset2 AdcSlope   Twc0  Twc1  Twc2 Twc3 Twc4\n"
          << "##############################################################################\n"
          << "[EmcCalPar]\n";
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      outFile << sector << "\t" << cell << "\t" << ecalCells[sector][cell].newParams << "\n";
    }
  }
}

Array<EcalCellPar> initialize()
{
  Array<EcalCellPar> ecalCells;

  ecalCells.reserve(6);
  for (Int_t sector = 0; sector < 6; sector++) {
    vector<EcalCellPar> oneSector;
    oneSector.reserve(163);
    for (Int_t cell = 0; cell < 163; cell++) {
      EcalCellPar a(sector, cell);
      oneSector.push_back(a);
    }
    ecalCells.push_back(oneSector);
  }
  return ecalCells;
  // gStyle->SetOptStat("n");
}
