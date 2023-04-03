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

using namespace std;


typedef struct {
	Int_t cent_bin;
	Int_t sector;
	Int_t emc_cell;
	Int_t rowEMC;
	Int_t columnEMC;

	Int_t id_geant;
	Float_t theta_kine;
	Float_t phi_kine;

	Float_t x_kine;
	Float_t y_kine;
	Float_t mom_kine;
	Int_t matchedCells;
	Int_t matchedPart;
	Float_t emc_energy;
	Float_t within_resolution;
	Float_t emc_beta;
	Int_t emc_size;
	Float_t xEMC;
	Float_t yEMC;
} TREE_VARIABLES;


void MakeAngleCellMap(void)
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

	TFile * out        = new TFile("AngleToCell.root", "recreate");
	TH2D * AngleToCell = new TH2D("AngleToCell", "AngleToCell", 360 * 2, 0, 360, 40 * 2, 10, 50);

	TH2D * AngleToRow = new TH2D("AngleToRow", "AngleToCell", 360 * 2, 0, 360, 40 * 2, 10, 50);
	TH2D * AngleToColumn = new TH2D("AngleToColumn", "AngleToCell", 360 * 2, 0, 360, 40 * 2, 10, 50);

	AngleToCell->Sumw2();
	TH2D * Counter = new TH2D("acceptancePhiTheta", "acceptancePhiTheta", 360 * 2, 0, 360, 40 * 2, 10, 50);
	Counter->Sumw2();


	for (Long64_t i = 0; i < nt->GetEntries(); i++) { // Entries loop start -> counter is increased at filling photons
		nt->GetEntry(i);

		if (
			vars.emc_size == 1 &&
			vars.within_resolution < 6 &&
			vars.mom_kine > 300 &&
			vars.id_geant == 1)
		{
			AngleToRow->Fill(vars.phi_kine, vars.theta_kine, vars.rowEMC);
			AngleToColumn->Fill(vars.phi_kine, vars.theta_kine, vars.columnEMC);

			AngleToCell->Fill(vars.phi_kine, vars.theta_kine, vars.emc_cell);
			Counter->Fill(vars.phi_kine, vars.theta_kine, 1);
		}
	}

	AngleToRow->Divide(Counter);
	AngleToColumn->Divide(Counter);
	AngleToCell->Divide(Counter);
	out->cd();
	AngleToCell->Write();
	AngleToRow->Write();
	AngleToColumn->Write();
	//Counter->Write();
} // MakeAngleCellMap
