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


#include "addFiles.h"    // library for adding multiple ttrees separated by comma
#include "heventmixer.h" // Szymon's event mixer, a bit improved
#include "CorrFactors.h"
#include "EcalCell.h"
#include "exp_setup_cells.h"
#include "BetaCut.h"

using namespace std;

//


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
	for (Int_t iSector = 0; iSector < 6; iSector++) {
		for (Int_t iCell = 0; iCell < 163; iCell++) {
			EcalCell a(iSector, iCell);
			if (isGoodCell(iSector, iCell)) a.InitializeAllMix();
			ecalCells.push_back(a);
		}
	}


	CorrFactors corrFactors[6][163];
  readCorrFactors("AllCorrFactors.root",corrFactors);

	HGenericEventMixerObj < Photon > eventmixer;
	eventmixer.setPIDs(1, 1, 7);
	eventmixer.setBuffSize(20);
	BetaCut betaCut("cellBordersBeta.root");

	nt->GetEntry(0);                      // get first number of photons for beam start time

	for (Long64_t i = 0; i < nentries;) {     // Entries loop start -> counter is increased at filling photons
		printProgress(i, nentries, 1, "Analyze pairs :");
		eventmixer.nextEvent(1000*pid.centBin+20*pid.multPhoton+(pid.targetBin+1));
		Int_t nPhot = pid.multPhoton;                          // number of photons per new event
		// ==============Reading and filling one event with photons=========================================//
		vector < Photon > photons;
		photons.reserve(nPhot);

		for (int k = 0; k < nPhot; k++) {
			if (pid.energy >= minEnergy && pid.energy < maxEnergy &&
			    pid.centBin>=0&& abs(pid.emc_mq)>20&& abs(pid.qualityDTime)>10 ) {
				Int_t energyBin = getEnergyBin(pid.energy);
				pid.mom_kine = energyBin;
			//	if( corrFactors[(int) pid.sector][(int) pid.cell].val[energyBin]<0.5) corrFactors[(int) pid.sector][(int) pid.cell].val[energyBin]=1.0;
				pid.energy  *= corrFactors[(int) pid.sector][(int) pid.cell].val[energyBin];

				photons.push_back({ pid });
				nt->GetEntry(++i);     // first increase counter than get the next photon in the same event
			} else { nt->GetEntry(++i); }
		}
		// =============================================================================================//

		if (photons.size() == 0) continue;


		// ==============One event combinations========================================================//
		for (size_t e1 = 0; e1 < (photons.size() - 1); e1++) {     // first gamma
			Photon photon1 = photons.at(e1);
			for (size_t e2 = e1 + 1; e2 < photons.size(); e2++) {     // second gamma
				Photon photon2 = photons.at(e2);
				Diphoton diphot(photon1, photon2);
				if ((diphot.opangle * 180. / 3.14) < 8) continue;

				allCells.Fill(diphot, "all, cells");
				ecalCells[photon1.cell].Fill(diphot, "all");
				ecalCells[photon2.cell].Fill(diphot, "all");
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
			if ((diphot.opangle * 180. / 3.14) < 8) continue;

			allCells.Fill(diphot, "cb, cells");
			ecalCells[gammaMix1.cell].Fill(diphot, "cb");
			ecalCells[gammaMix2.cell].Fill(diphot, "cb");
		}
	}
	out->cd();
	// =============================================================================================//
	for (Int_t iSector = 0; iSector < 6; iSector++) {
		for (Int_t iCell = 0; iCell < 163; iCell++) {
			if (!isGoodCell(iSector, iCell)) continue;
			ecalCells[iSector * 163 + iCell].WriteHists();
		}
	}
	allCells.WriteHists();
	out->Save();
	out->Close();
} /* fill_hists_diphoton */
