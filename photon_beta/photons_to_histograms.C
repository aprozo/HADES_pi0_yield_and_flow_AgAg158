#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCutG.h"
#include "htool.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObjString.h"
#include <iostream>
#include <map>
#include <stdio.h>
#include <iomanip>
#include "exp_setup_cells.h"

using namespace std;


void printProgress(Long64_t ct, Long64_t total, Long64_t step, TString comment)
{
  if (total <= 0) {
    cout << "Error::printProgress(), total is <=0 !" << endl;
    return;
  }
  if (ct < 0) {
    cout << "Error::printProgress(), ct is < 0 !" << endl;
    return;
  }
  if (ct == 0) cout << comment.Data() << flush;
  Long64_t frac   = (Long64_t) ((( ct / (Float_t) total) ) * 100);
  Long64_t frac_1 = (Long64_t) ((((ct - 1) / (Float_t) total) ) * 100);

  if (frac % step == 0 && (frac != frac_1) ) {
    if (frac_1 / step > 0) cout << "\b\b\b\b\b" << flush;
    cout << setw(3) << frac << " %" << flush;
    if (frac + step >= 100) cout << endl;
  }
}

Bool_t addFile(TString infile, TChain * fChain, TH1D * Counter_cent)
{
  TObjArray * elements = fChain->GetListOfFiles();
  Int_t nfiles         = elements->GetEntries();
  TFile * temp;

  if (gSystem->AccessPathName(infile.Data()) == 0) {
    std::cout << std::setw(4) << nfiles << " add file : " << infile.Data() << std::endl;
    fChain->AddFile(infile.Data());
    temp = TFile::Open(infile.Data(), "read");
    Counter_cent->Add((TH1D *) temp->Get("Counter_cent"));
    temp->Close();
    return kTRUE;
  } else {
    return kFALSE;
  }
}

void addMultFiles(TString commaSeparatedList, TChain * fChain, TH1D * Counter_cent)
{
  TObjArray * arr = commaSeparatedList.Tokenize(",");
  Int_t n         = arr->GetEntries();

  for (Int_t i = 0; i < n; i++) {
    TString name = ((TObjString *) arr->At(i))->GetString();
    if (!addFile(name, fChain, Counter_cent)) {
      delete arr;
    }
  }
  delete arr;
}

typedef struct {
  Float_t newEvent, multTofRpc, multEmc, multPhoton, centBin,
    targetBin, beamTime, RPlanePhi, impPar, emc_mq, qualityDTime, betaCand,
    sector, cell, size, beta, energy, time, xpos, ypos, zpos, track,
    trackParent, mom_kine, idParent;
} PID;


Int_t getBeamBin(Float_t & beamtime)
{
  if (beamtime < 74.4) return 0;
  else if (beamtime < 82.5) return 1;
  else return 2;
}

TH1D * hEnergyCorrection;

void correctNonlinearity(Float_t & energy)
{
  Double_t corr;

  if (energy < 110) corr = hEnergyCorrection->Interpolate(110);
  if (energy > 640) corr = hEnergyCorrection->Interpolate(640);
  else corr = hEnergyCorrection->Interpolate(energy);
  energy *= corr;
}

Int_t photons_to_histograms(TString inFile, TString outputFile, Int_t nevt)
{
  TFile * fCorr = new TFile(
    "/lustre/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/NonlinearityPi0PeakEnergyCorrection.root", "READ");

  hEnergyCorrection = (TH1D *) fCorr->Get("hEnergyCorrection");


  TFile * out = new TFile(outputFile, "RECREATE");
  out->cd();

  TH1D * Counter_cent = new TH1D("Counter_cent", "Counter_cent", 5, 0, 5);

  Counter_cent->GetXaxis()->SetBinLabel(5, "all_cent");
  Counter_cent->GetXaxis()->SetBinLabel(1, "cent0_10");
  Counter_cent->GetXaxis()->SetBinLabel(2, "cent10_20");
  Counter_cent->GetXaxis()->SetBinLabel(3, "cent20_30");
  Counter_cent->GetXaxis()->SetBinLabel(4, "cent30_40");

  TH1D * nEventsBeamCounter;
  TH1D * nEventsBeam;
  TH2F * betaBeam[6][163];
  TH2F * betaEnergy[3][6][163];

  TH2D * betaCell;
  TH2D * betaCellLow;

  const Int_t nBeamBins = 200;
  Float_t xBeamBins[nBeamBins + 1] = { 0., 74.4, 82.5, 100 };
  const Int_t nBetaBins = 100;

  //nEventsBeamCounter = new TH1D("nEventsBeamCounter", "nEventsBeam;beam time, days;nEvents", nBeamBins, xBeamBins);
  //nEventsBeam        = new TH1D("nEventsBeam", "nEventsBeam;beam time, days;nEvents", nBeamBins, xBeamBins);
  nEventsBeamCounter = new TH1D("nEventsBeamCounter", "nEventsBeam;beam time, days;nEvents", nBeamBins, 68,88);
  nEventsBeam        = new TH1D("nEventsBeam", "nEventsBeam;beam time, days;nEvents", nBeamBins, 68,88);

  betaCell    = new TH2D("betaCell", "betaCell;sector*200+cell;Beta", 1200, 0, 1200, 500, 0.8, 1.5);
  betaCellLow = new TH2D("betaCellLow", "betaCellLow;sector*200+cell;Beta", 1200, 0, 1200, 500, 0.8, 1.5);


  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      TString histName = Form("sector%i_cell%i_beam_vs_beta", sector, cell);
      betaBeam[sector][cell] = new TH2F(histName, histName + ";beam time, days;Beta", nBeamBins, 68, 88, nBetaBins, 0.8,
          1.3);
    //   for (Int_t beamBin = 0; beamBin < nBeamBins; beamBin++)
    //     betaEnergy[beamBin][sector][cell] = new TH2F(Form("sector%icell%ibeam%ibetaEnergy", sector, cell, beamBin),
    //         ";Energy, MeV ;Beta", 5, 50, 800, nBetaBins, 0.8, 1.3);
     }
  }


  TFile * betaCuts;

  TH2F * betaEnergyMean[nBeamBins];
  TH2F * betaEnergySigma[nBeamBins];
  // if (isSim) betaCuts = new TFile("/lustre/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/betaBeamMeanSigmaSim.root", "Read");
  // else
  betaCuts = new TFile(
    "/lustre/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/betaBeamMeanSigma.root", "Read");
  for (Int_t iBeamBin = 0; iBeamBin < 3; iBeamBin++) {
    betaEnergyMean[iBeamBin]  = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergyMean", iBeamBin));
    betaEnergySigma[iBeamBin] = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergySigma", iBeamBin));
    if (betaEnergyMean[iBeamBin] == NULL ||
      betaEnergySigma[iBeamBin] == NULL) std::cout << "No beta corrections loaded" << std::endl;
  }


  TChain * nt = new TChain("nt");
  addMultFiles(inFile, nt, Counter_cent);
  PID pid;
  nt->SetBranchAddress("pid", &pid.newEvent);
  Long64_t nentries = nt->GetEntries();
  nt->GetEntry(0); // get first number of photons for beam start time

  nEventsBeamCounter->Fill(pid.beamTime, Counter_cent->GetBinContent(5));
  Bool_t isSim = (pid.mom_kine > 0);
  for (Long64_t i = 1; i < nentries;) {
    printProgress(i, nentries, 1, "Analyze pairs :");
		if (!(isGoodBeamTime(pid.beamTime) || isSim) ) { nt->GetEntry(++i); continue; }
    if (pid.newEvent == 1) Counter_cent->Fill(pid.centBin);

    Int_t nPhot = pid.multPhoton; // number of photons per new event
    if (nPhot == 0)  { nt->GetEntry(++i); continue; }

    for (int k = 0; k < nPhot; k++) {
			if (!isSim) correctNonlinearity(pid.energy);
      if (pid.size <= 1 &&
        pid.energy > 150. &&
        pid.energy < 2000. &&
        (pid.qualityDTime > 3 ||
        pid.qualityDTime < -3)
        &&isGoodCell(pid.sector,pid.cell)
        &&
        pid.emc_mq > 10 &&
        (isSim || isGoodBeamTime(pid.beamTime))
      )
      {
        Int_t sectorCell = (int) pid.sector * 200 + (int) pid.cell;
        Int_t beamBin    = getBeamBin(pid.beamTime);
        // betaEnergy[beamBin][(int)pid.sector][(int)pid.cell]->Fill(pid.energy,pid.beta);


        Int_t energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(pid.energy);
        if (pid.energy < 200) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(250);
        if (pid.energy >= 800) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(750);
        // cout<<"energyBin="<<energyBin<<endl;
        Float_t beta_mean = betaEnergyMean[beamBin]->GetBinContent(energyBin,
            (int) pid.sector * 200 + (int) pid.cell + 1);

        Float_t beta_sigma = betaEnergySigma[beamBin]->GetBinContent(energyBin,
            (int) pid.sector * 200 + (int) pid.cell + 1);

        Float_t beta_cut_down = beta_mean - 1. * beta_sigma;
        Float_t beta_cut_up   = beta_mean + 1. * beta_sigma;

        if (pid.beta > beta_cut_down && pid.beta < beta_cut_up) {
        //  betaEnergy[beamBin][(int) pid.sector][(int) pid.cell]->Fill(pid.energy, pid.beta);
          if (pid.energy > 300) {
            betaBeam[(int) pid.sector][(int) pid.cell]->Fill(pid.beamTime, pid.beta);
            betaCell->Fill(sectorCell, pid.beta);
          } else {
            betaCellLow->Fill(sectorCell, pid.beta);
          }
        }
        nt->GetEntry(++i); // first increase counter than get next photon in te same event
      } else { nt->GetEntry(++i); }
    }
    // =====================================================================================//
  }

  out->cd();
  out->GetList()->Write();

  out->Save();
  out->Close();

  return 1;
} /* photons_to_histograms */
