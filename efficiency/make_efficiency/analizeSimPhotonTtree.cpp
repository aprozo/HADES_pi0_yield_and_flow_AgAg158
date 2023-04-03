#include <iostream>
#include <map>
#include <stdio.h>
#include "TH3F.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPad.h"
#include "TTree.h"
#include "TChain.h"
#include "smoothMatrix.C"
#include "exp_setup_cells.h"
using namespace std;


typedef struct {
  Int_t   cent_bin;
  Int_t   sector;
  Int_t   emc_cell;
  Int_t   rowEMC;
  Int_t   columnEMC;

  Int_t   id_geant;
  Float_t theta_kine;
  Float_t phi_kine;

  Float_t x_kine;
  Float_t y_kine;
  Float_t mom_kine;
  Int_t   matchedCells;
  Int_t   matchedPart;
  Float_t emc_energy;
  Float_t within_resolution;
  Float_t emc_beta;
  Int_t   emc_size;
  Float_t xEMC;
  Float_t yEMC;
} TREE_VARIABLES;


void analizeSimPhotonTtree(void)
{
  TFile * fileIn = new TFile("../sim_photon_ttree_gen5a.root", "Read");
  TChain * nt    = (TChain *) fileIn->Get("nt");
  TREE_VARIABLES vars;

  nt->SetBranchAddress("cent_bin", &vars.cent_bin);
  nt->SetBranchAddress("sector", &vars.sector);
  nt->SetBranchAddress("emc_cell", &vars.emc_cell);
  nt->SetBranchAddress("rowEMC", &vars.rowEMC);
  nt->SetBranchAddress("columnEMC", &vars.columnEMC);
  nt->SetBranchAddress("id_geant", &vars.id_geant);
  nt->SetBranchAddress("theta_kine", &vars.theta_kine);
  nt->SetBranchAddress("phi_kine", &vars.phi_kine);
  nt->SetBranchAddress("x_kine", &vars.x_kine);
  nt->SetBranchAddress("y_kine", &vars.y_kine);
  nt->SetBranchAddress("mom_kine", &vars.mom_kine);
  nt->SetBranchAddress("matchedCells", &vars.matchedCells);
  nt->SetBranchAddress("matchedPart", &vars.matchedPart);
  nt->SetBranchAddress("emc_energy", &vars.emc_energy);
  nt->SetBranchAddress("within_resolution", &vars.within_resolution);
  nt->SetBranchAddress("emc_beta", &vars.emc_beta);
  nt->SetBranchAddress("emc_size", &vars.emc_size);
  nt->SetBranchAddress("xEMC", &vars.xEMC);
  nt->SetBranchAddress("yEMC", &vars.yEMC);

  TFile * in         = new TFile("AngleToCell.root", "read");
  TFile * acceptance = new TFile("acceptancePhiThetaMar19.root", "read");
  // TH2D * AngleToCell = (TH2D *) in->Get("AngleToCell");
  TH2D * AngleToRow    = (TH2D *) in->Get("AngleToRow");
  TH2D * AngleToColumn = (TH2D *) in->Get("AngleToColumn");

  TH2D * acceptancePhiTheta = (TH2D *) acceptance->Get("acceptancePhiThetaMar19");
  TH2D * numeratorCells[3];
  TH2D * numeratorCells2[3];
  TH2D * denominatorCells[3];

  for (int cent = 0; cent < 3; cent++) {
    TString centName = Form("%i0_%i0", cent, cent + 1);
    numeratorCells[cent] =
      new TH2D("cell_mom_efficiency" + centName,
        "Efficiency" + centName + ";sector*200+cell;mom", 1200, 0, 1200, 12, 100, 700);

    numeratorCells2[cent] = new TH2D("cell_mom_efficiency2" + centName,
        "Efficiency" + centName + ";sector*200+cell;mom", 1200, 0, 1200, 12, 100,
        700);
    numeratorCells[cent]->Sumw2();

    denominatorCells[cent] = new TH2D("cell_mom_denominator" + centName,
        "Efficiency" + centName + ";sector*200+cell;mom", 1200, 0, 1200, 12, 100,
        700);
    denominatorCells[cent]->Sumw2();
  }

  TFile * out = new TFile("efficiencyPhotons.root", "RECREATE");

  for (Long64_t i = 0; i < nt->GetEntries(); i++) { // Entries loop start -> counter is increased at filling photons
    nt->GetEntry(i);
    if (vars.cent_bin > 3) continue;
    // if (vars.mom_kine <100) continue;
    Int_t binAngle = acceptancePhiTheta->FindBin(vars.phi_kine, vars.theta_kine);
    if (acceptancePhiTheta->GetBinContent(binAngle) == 0) continue;

    Int_t rowEMC    = AngleToRow->GetBinContent(AngleToRow->FindBin(vars.phi_kine, vars.theta_kine));
    Int_t columnEMC = AngleToColumn->GetBinContent(AngleToColumn->FindBin(vars.phi_kine, vars.theta_kine));
    Int_t cell      = systemMap[rowEMC * 17 + columnEMC];
    if (cell < 10) continue;

    cell += vars.sector * 200;
    denominatorCells[(int) vars.cent_bin - 1]->Fill(cell, vars.mom_kine);

    if (vars.matchedCells == 0 &&
      vars.matchedPart == 0
      vars.emc_beta >= 0.95 && vars.emc_beta <= 1.05 &&
 // vars.within_resolution < 5.5 //(3 sigma)
    )
    {
      numeratorCells[(int) vars.cent_bin - 1 ]->Fill(cell, vars.emc_energy);
      numeratorCells2[(int) vars.cent_bin - 1]->Fill(vars.emc_cell + vars.sector * 200, vars.emc_energy);
    }
  }
  for (int cent = 0; cent < 3; cent++) {
    numeratorCells[cent]->Divide(denominatorCells[cent]);
    numeratorCells2[cent]->Divide(denominatorCells[cent]);
    out->cd();
    numeratorCells[cent]->Write();
    numeratorCells2[cent]->Write();
  }
} // MakeAngleCellMap
