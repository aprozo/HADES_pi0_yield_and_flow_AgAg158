#include "TFile.h"
#include "TH1.h"

struct cut {
  Float_t left, right;
  Bool_t withinBorders(const Float_t & beta)
  {
    if ((left) <= beta && (beta <= right)) return kTRUE;
    else return kFALSE;
  }
};

typedef struct {
  Float_t sector, cell, LeftBeamTime, RightBeamTime, LeftBeta, RightBeta;
} BetaBorders;


class BetaCut {
public:
  BetaCut(TString inFile = "cellBordersBeta.root")
  {
    TChain * t = new TChain("t");

    t->Add(inFile);
    BetaBorders in;
    t->SetBranchAddress("BetaBorders", &in.sector);
    Long_t nList = t->GetEntries();
    for (Int_t iList = 0; iList < nList; ++iList) {
      t->GetEntry(iList);
      Int_t sector  = (Int_t) in.sector;
      Int_t cell    = (Int_t) in.cell;
      Int_t beamBin = (Int_t) 200 * (in.LeftBeamTime - 67) / 21;
      bBeta[sector][cell][beamBin].left  = in.LeftBeta - 0.01;
      bBeta[sector][cell][beamBin].right = in.RightBeta + 0.01;
    }
  }

  Bool_t isGoodBeta(const PID & pid)
  {
    Int_t beamBin = (Int_t) 200 * (pid.beamTime - 67) / 21;
    Int_t sector  = (Int_t) pid.sector;
    Int_t cell    = (Int_t) pid.cell;

    if (bBeta[sector][cell][beamBin].withinBorders(pid.beta)) return kTRUE;
    else return kFALSE;
  };

private:
  cut bBeta[6][163][200];
};
//
//
// TFile * betaCut = new TFile("betaEnergyCut.root");
// TH1D * mean=(TH1D *)betaCut->Get("hEnergyBetaNoCuts_1");
// TH1D * sigma=(TH1D *)betaCut->Get("hEnergyBetaNoCuts_2");
//
// Bool_t isGoodBetaEn(Float_t beta,Float_t energy){
//         if (energy < 300) energy = 300;
//         if (energy > 1500) energy = 1490;
//         Float_t meanVal = mean->GetBinContent(mean->FindBin(energy));
//         Float_t twoSigmaVal = sigma->GetBinContent(sigma->FindBin(energy));
//         return (abs(beta  - meanVal) <twoSigmaVal);
// }
//
// Float_t getMeanBeta(Float_t energy){
//         return mean->GetBinContent(mean->FindBin(energy));
// }
//
// Float_t get2SigmaBeta(Float_t energy){
//         return sigma->GetBinContent(sigma->FindBin(energy));
// }
