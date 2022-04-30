#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "exp_setup_cells.h"

using namespace std;



typedef struct {
	Float_t sector, cell,  tof, emc_mq, mom, multEmc,
	        beta, charge, energy, time, emc_path,distRpcEmc,timeCorr,
	        beamtime, TOT, qualityDThDPh, qualityDTime, raw_time, sys;
}
PID;
void MakeCellTtrees()
{

	TString treeLocation ="/lustre/hebe/hades/user/prozorov/mar19/ChargedParticlesTime/";
	TChain * nt = new TChain("nt");
	nt->Add(treeLocation+"chargedParticles.root");



	PID pid;
	PID newPid;
	nt->SetBranchAddress("pid", &pid);

	TTree * tree[6][163];
	for (Int_t iSector = 0; iSector < 6; ++iSector) {
		for (Int_t iCell = 0; iCell < 163; ++iCell) {
			tree[iSector][iCell] = new TTree(Form("Sector%iCell%i", iSector, iCell), Form("Sector%iCell%i", iSector, iCell));
			tree[iSector][iCell]->Branch("pid", &newPid,"sector:cell:tof:emc_mq:mom:multEmc:beta:charge:energy:time:emc_path:"
			          "distRpcEmc:timeCorr:beamtime:TOT:qualityDThDPh:qualityDTime:raw_time:sys");
	}}

	// ==========================================
	// Loop over ttree
	// ==========================================
	for (Long_t iEntry = 0; iEntry < nt->GetEntries(); ++iEntry) {
		nt->GetEntry(iEntry);
		// counter
		if (iEntry % 5000000 == 0) cout << ((float) iEntry / nt->GetEntries() * 100) << "%" << endl;

		Int_t sector = pid.sector;
		Int_t cell   = pid.cell;
		if (cell < 0) continue;
		if (sector > 6 || sector < 0 || cell < 0 || cell > 162) {
			std::cout << "sector" << sector << " cell" << cell << "ERROR" << std::endl;
			continue;
		}
		newPid      = pid;
		tree[sector][cell]->Fill();
	}


	for (Int_t iSector = 0; iSector < 6; ++iSector) {
		for (Int_t iCell = 0; iCell < 163; ++iCell) {
      if(!isGoodCell(iSector,iCell))continue;
			TFile * outfile =
				new TFile(Form("/lustre/hebe/hades/user/prozorov/mar19/ChargedParticlesTime/splitted/tree_Sector%i_Cell%i.root",
				               iSector, iCell), "RECREATE");
			outfile->cd();
			tree[iSector][iCell]->Write();
			outfile->Save();
			outfile->Close();
			delete outfile;
			delete tree[iSector][iCell];
		}
	}
	delete nt;
} // MakeCellTtrees
