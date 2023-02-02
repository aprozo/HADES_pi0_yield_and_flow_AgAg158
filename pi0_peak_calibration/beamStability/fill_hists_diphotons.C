#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "htool.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"

#include <iostream>
#include <stdio.h>

TH1D * hEnergyCorrection;
#include "addFiles.h"    // library for adding multiple ttrees separated by comma
#include "heventmixer.h" // Szymon's event mixer, a bit improved

#include "EcalCell.h"
#include "exp_setup_cells.h"


using namespace std;

Int_t getBeamBin(const Float_t & beamtime)
{
  if (beamtime < 74.4) return 0;
  else if (beamtime < 82.5) return 1;
  else return 2;
}

void correctNonlinearity(Float_t & energy)
{
  Double_t corr;

  if (energy < 110) corr = hEnergyCorrection->Interpolate(110);
  if (energy > 640) corr = hEnergyCorrection->Interpolate(640);
  else corr = hEnergyCorrection->Interpolate(energy);
  energy *= corr;
}

void fill_hists_diphoton(TString inFile, TString outputFile, Int_t iteration = 0)
{
  TH1::SetDefaultSumw2();

  TChain * nt = new TChain("nt");
  addMultFiles(inFile, nt);
  PID pid;
  nt->SetBranchAddress("pid", &pid.newEvent);
  Long64_t nentries = nt->GetEntries();
  TFile * out       = new TFile(outputFile, "RECREATE");

  vector < EcalCell > ecalCells;
  EcalCell allCells(0, 0);
  allCells.InitializeAllMix();


  EcalCell sector2(2, 0);
  EcalCell sector4(4, 0);
  EcalCell sector5(5, 0);

  sector2.InitializeAllMix();
  sector4.InitializeAllMix();
  sector5.InitializeAllMix();


  HGenericEventMixerObj < Photon > eventmixer;
  eventmixer.setPIDs(1, 1, 7);
  eventmixer.setBuffSize(30);

  TFile * fCorr = new TFile(
    "/lustre/hebe/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/NonlinearityPi0PeakEnergyCorrection.root", "READ");
  hEnergyCorrection = (TH1D *) fCorr->Get("hEnergyCorrection");

  TFile * betaCuts;
  const Int_t nBeamBins = 3;
  TH2F * betaEnergyMean[nBeamBins];
  TH2F * betaEnergySigma[nBeamBins];

  betaCuts = new TFile(
    "/lustre/hebe/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/betaBeamMeanSigma.root", "Read");
  for (Int_t iBeamBin = 0; iBeamBin < nBeamBins; iBeamBin++) {
    betaEnergyMean[iBeamBin]  = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergyMean", iBeamBin));
    betaEnergySigma[iBeamBin] = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergySigma", iBeamBin));
    if (betaEnergyMean[iBeamBin] == NULL ||
      betaEnergySigma[iBeamBin] == NULL) std::cout << "No beta corrections loaded" << std::endl;
  }


  nt->GetEntry(0); // get first number of photons for beam start time
  Float_t beamTimeStart = pid.beamTime;
  Float_t beamTime      = beamTimeStart;
  for (Long64_t i = 0; i < nentries;) { // Entries loop start -> counter is increased at filling photons
    printProgress(i, nentries, 1, "Analyze pairs :");
    eventmixer.nextEvent(1000 * pid.centBin + 20 * pid.multPhoton + (pid.targetBin + 1));
    
    // ==============Reading and filling one event with photons=========================================//
    vector < Photon > photons;
    //  photons.reserve(nPhot);
    Int_t nPhot = pid.multPhoton; // number of photons per new event
    if (nPhot == 0)  { nt->GetEntry(++i); continue; }
    for (int k = 0; k < nPhot; k++) {
      correctNonlinearity(pid.energy);
      if (pid.energy >= 100 && pid.energy < 1200 &&
        isGoodCell((Int_t) pid.sector, (Int_t) pid.cell) &&
        pid.centBin >= 0 && abs(pid.emc_mq) > 10 && abs(pid.qualityDTime) > 10
      )
      {

        Int_t beamBin   = getBeamBin(pid.beamTime); // cout<<"beamBin="<<beamBin<<endl;
        Int_t energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(pid.energy);
        if (pid.energy < 100) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(150);
        if (pid.energy >= 800) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(750);
        // cout<<"energyBin="<<energyBin<<endl;

        Float_t BETA_MEAN = betaEnergyMean[beamBin]->GetBinContent(energyBin,
            (int) pid.sector * 200 + (int) pid.cell + 1);
        if (pid.energy < 200) BETA_MEAN = betaEnergyMean[beamBin]->GetBinContent(energyBin + 1,
              (int) pid.sector * 200 + (int) pid.cell + 1);

        Float_t BETA_SIGMA = betaEnergySigma[beamBin]->GetBinContent(energyBin,
            (int) pid.sector * 200 + (int) pid.cell + 1);
        Float_t BETA_CUT_DOWN = BETA_MEAN - BETA_SIGMA;
        Float_t BETA_CUT_UP   = BETA_MEAN + BETA_SIGMA;

        if (pid.beta > BETA_CUT_DOWN && pid.beta < BETA_CUT_UP) {
          photons.push_back({ pid });
        }
        nt->GetEntry(++i); // first increase counter than get the next photon in the same event
      } else { nt->GetEntry(++i); }
    }

    // =============================================================================================//
    Float_t beamTimeEnd = pid.beamTime;
    beamTime = (beamTimeEnd + beamTimeStart) / 2;


    if (photons.empty()) continue;


    // ==============One event combinations========================================================//
    for (size_t e1 = 0; e1 < (photons.size() - 1); e1++) { // first gamma
      Photon photon1 = photons.at(e1);
      for (size_t e2 = e1 + 1; e2 < photons.size(); e2++) { // second gamma
        Photon photon2 = photons.at(e2);
        Diphoton diphot(photon1, photon2);
        if ((diphot.opangle * 180. / 3.14) < 8) continue;
        // allCells.Fill(diphot, "all");
        if (diphot.sector == 2) sector2.Fill(diphot, "all");
        if (diphot.sector == 4) sector4.Fill(diphot, "all");
        if (diphot.sector == 5) sector5.Fill(diphot, "all");
      }
    }
    // =============================================================================================//

    // ====================Mixed-event combinations=================================================//
    eventmixer.addVector(photons, 1);
    vector < pair < Photon, Photon >>& pairsVec = eventmixer.getMixedVector();
    Int_t size = pairsVec.size();

    for (Int_t j = 0; j < size; j++) {
      pair < Photon, Photon >& pair = pairsVec[j];
      Photon gammaMix1 = pair.first;
      Photon gammaMix2 = pair.second;
      Diphoton diphot(gammaMix1, gammaMix2);
      if ((diphot.opangle * 180. / 3.14) < 6) continue;
      // allCells.Fill(diphot, "cb");
      if (diphot.sector == 2) sector2.Fill(diphot, "cb");
      if (diphot.sector == 4) sector4.Fill(diphot, "cb");
      if (diphot.sector == 5) sector5.Fill(diphot, "cb");
    }
  }
  out->cd();
  // =============================================================================================//
  // allCells.Evaluate(beamTime);
  // allCells.Write();
  sector2.Evaluate(beamTime);
  sector2.Write();
  sector4.Evaluate(beamTime);
  sector4.Write();
  sector5.Evaluate(beamTime);
  sector5.Write();


  out->Save();
  out->Close();
} /* fill_hists_diphoton */
