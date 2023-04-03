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
  Float_t mult_tofrpc, imppar, cent_bin, charge, sector,
    emc_cell, theta_rpc, phi_rpc, mom, beta, mom_kine,
    emc_size, emc_energy, emc_time, emc_mq, qualityDTime,
    emc_path, mm_qa, geantid, within_resolution,
    true_emc_hit_track, true_rpc_hit_track, cflag, beamtime;
}
PID_exp_lepton;

typedef struct {
  Float_t imppar, cent_bin, charge, emc_match, sector,
    emc_cell, theta_rpc, phi_rpc, mom, beta, mom_kine,
    emc_size, emc_energy, emc_time, emc_energy_corr, emc_mq, qualityDTime,
    emc_path, mm_qa, geantid, within_resolution, emc_energy_centr,
    true_emc_hit_track, true_rpc_hit_track, cflag;
}
PID_sim_lepton;


void analizeSimPhotonTtree(void)
{
  fileIn[sim_lepton] = new TFile("../sim_photon_ttree_gen5a.root", "Read");
  fileIn[exp_lepton] = new TFile("../exp_lepton_ttree_gen5a.root", "Read");

  TChain * nt = (TChain *) fileInSim->Get("nt");
  TChain * nt = (TChain *) fileInSim->Get("nt");
  PID_exp_lepton expPid;
  PID_sim_lepton simPid;

  nt->SetBranchAddress("pid", &expPid.mult_tofrpc);
  nt->SetBranchAddress("pid", &simPid.imppar);


  TFile * out = new TFile("efficiencyLeptons.root", "RECREATE");

  for (Long64_t i = 0; i < nt->GetEntries(); i++) { // Entries loop start -> counter is increased at filling photons
    nt->GetEntry(i);

    if (vars.matchedCells == 0 &&
      vars.matchedPart == 0 &&
      vars.emc_beta > 0.95 &&
      vars.emc_beta < 1.05 &&
      vars.within_resolution < 5.5)
    {
      numeratorCells->Fill(cell, vars.emc_energy);
      numeratorCells2->Fill(vars.emc_cell + vars.sector * 200, vars.emc_energy);
    }
  }
  out->cd();
} // MakeAngleCellMap
