#include "TTree.h"
#include "TChain.h"
#include "TObject.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TString.h"
#include "TSystem.h"
#include "TF1.h"
#include <iostream>

// #include "exp_setup_cells.h" // some useful functions for ECAL cell operations


using namespace std;


void fitHistsBetaBeamEnergy()
{
  const Int_t nBeamBins = 3;
  const Double_t xBeamBins[nBeamBins + 1] = { 0., 74.4, 82.5, 100 };
  Bool_t isSim = kFALSE;


  TFile * in;

  if (isSim) in = new TFile("photonsBetaSim.root", "read");
  else in = new TFile("photonsBeta3beamBins.root", "read");
  TFile * out;
  if (isSim) out = new TFile("betaBeamMeanSigmaSim.root", "recreate");
  else out = new TFile("betaBeamMeanSigma.root", "recreate");


  out->cd();


  TH2F * betaEnergyMean[nBeamBins];
  TH2F * betaEnergySigma[nBeamBins];
  TH1F * averageMeanBeam [nBeamBins];
  TH1F * averageSigmaBeam [nBeamBins];


  for (Int_t iBeamBin = 0; iBeamBin < nBeamBins; iBeamBin++) {
    betaEnergyMean[iBeamBin] = new TH2F(Form("beam%ibetaEnergyMean",
        iBeamBin), ";Energy, MeV;sector*200+cell ", 5, 50, 800, 1200, 0, 1200);
    betaEnergySigma[iBeamBin] = new TH2F(Form("beam%ibetaEnergySigma",
        iBeamBin), ";Energy, MeV;sector*200+cell ", 5, 50, 800, 1200, 0, 1200);
    averageMeanBeam [iBeamBin]  = new TH1F(Form("beam%iaverageMeanBeam", iBeamBin), ";cell", 1200, 0, 1200);
    averageSigmaBeam [iBeamBin] = new TH1F(Form("beam%iaverageSigmaBeam", iBeamBin), ";cell", 1200, 0, 1200);

    for (Int_t sector = 0; sector < 6; sector++) {
      for (Int_t cell = 0; cell < 163; cell++) {
        std::cout << "beam " << iBeamBin << "sector " << sector << " cell " << cell << '\n';


        TH2F * betaEnergy = (TH2F *) in->Get(Form("sector%icell%ibeam%ibetaEnergy", sector, cell, iBeamBin));
        betaEnergy->RebinY(2);

        if (betaEnergy->GetEntries() == 0) continue;

        TH1D * sumBetaEnergy = betaEnergy->ProjectionY("sum", 2, betaEnergy->GetNbinsX());

        TF1 * ff  = new TF1("ff", "gaus + pol3(3)", 0.9, 1.3);
        TF1 * ff1 = new TF1("ff1", "gaus", 0.9, 1.3);
        ff1->SetLineColor(kBlue);
        ff->SetParameters(sumBetaEnergy->GetMaximum(), 1., 0.05, 0., 0., 0.);
        if (isSim) { ff->SetRange(0.95, 1.05); ff1->SetRange(0.95, 1.05); }
        sumBetaEnergy->Fit("ff", "QRM");
        ff1->FixParameter(0, ff->GetParameter(0));
        ff1->FixParameter(1, ff->GetParameter(1));
        ff1->FixParameter(2, ff->GetParameter(2));
        sumBetaEnergy->Fit("ff1", "+Q");
        Float_t averageMean  = ff->GetParameter(1);
        Float_t averageSigma = ff->GetParameter(2);


        Float_t p0 = ff->GetParameter(3);
        Float_t p1 = ff->GetParameter(4);
        Float_t p2 = ff->GetParameter(5);

        if (averageMean < 0.9 || averageMean > 1.2 || averageSigma > 0.2 || averageSigma < 0.0) {
          averageMean  = 0.;
          averageSigma = 0.;
        }

        for (Int_t bin = betaEnergy->GetNbinsX(); bin >= 1; bin--) {
          TH1D * projY = betaEnergy->ProjectionY((TString) betaEnergy->GetName() + Form("energy%i", bin), bin, bin);
          if (projY->GetEntries() == 0) continue;
          ff->SetParameters(projY->GetMaximum(), averageMean, averageSigma, 0., 0., 0.);
          projY->Fit("ff", "QRM");
          ff1->FixParameter(0, ff->GetParameter(0));
          ff1->FixParameter(1, ff->GetParameter(1));
          ff1->FixParameter(2, ff->GetParameter(2));
          projY->Fit("ff1", "+Q");
          Float_t mean  = ff->GetParameter(1);
          Float_t sigma = ff->GetParameter(2);
          projY->Write();


          if (mean < 0.90 || mean > 1.2 || sigma > 0.2) {
            mean  = averageSigma;
            sigma = averageSigma;
          }

          betaEnergyMean[iBeamBin]->SetBinContent(bin, sector * 200 + cell + 1, mean);
          betaEnergySigma[iBeamBin]->SetBinContent(bin, sector * 200 + cell + 1, sigma);
        }// loop over  energy bin
        averageMeanBeam [iBeamBin]->SetBinContent(sector * 200 + cell + 1, averageMean);
        averageSigmaBeam [iBeamBin]->SetBinContent(sector * 200 + cell + 1, averageSigma);
        delete ff;
        delete ff1;
      }// loop over  cells
    }// loop over  sectors
    betaEnergyMean[iBeamBin]->Write();
    betaEnergySigma[iBeamBin]->Write();
    averageMeanBeam[iBeamBin]->Write();
    averageSigmaBeam[iBeamBin]->Write();
  }// loop over  beam bin
} // fitHists
