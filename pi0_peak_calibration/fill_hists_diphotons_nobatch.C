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



const Float_t pi0TrueMass = 134.977;
const Int_t nEnergyBins   = 1;
const Float_t n = 1.6;
const Float_t nInterations = 4;
const Float_t minEnergy = 100.;
const Float_t maxEnergy = 1500.;

//const Double_t energyBinAxis[nEnergyBins+1] = {minEnergy,150.,200.,250.,300.,375.,500.,maxEnergy};
const Double_t energyBinAxis[nEnergyBins+1] = {minEnergy,maxEnergy};


#include "CorrFactors.h"
#include "EcalCell.h"
#include "exp_setup_cells.h"
#include "BetaCut.h"

using namespace std;
//


void fill_hists_diphoton(TString inFile, TString outputFile, Int_t nev = 1000.)
{
	TH1::SetDefaultSumw2();

	TChain * nt = new TChain("nt");
	addMultFiles(inFile, nt);
	PID pid;
	nt->SetBranchAddress("pid", &pid.newEvent);
	Long64_t nentries = nt->GetEntries();
	TFile * out       = new TFile(outputFile, "RECREATE");


	vector < EcalCell > ecalCells;
	EcalCell allCells(1000, 1000);

	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			EcalCell a(sector, cell);
			ecalCells.push_back(a);
		}
	}


	CorrFactors corrFactors[6][163];
	//readCorrFactors("corrFactors.txt",corrFactors);


	for (Int_t iter = 0; iter < nInterations; iter++) {
		allCells.Reset();

		for (Int_t sector = 0; sector < 6; sector++) {
			for (Int_t cell = 0; cell < 163; cell++) {
				ecalCells[sector * 163 + cell].Reset();
			}
		}

		HGenericEventMixerObj < Photon > eventmixer;
		eventmixer.setPIDs(1, 1, 7);
		eventmixer.setBuffSize(10);
		BetaCut betaCut("cellBordersBeta.root");

		cout << "Iteration " << iter << endl;
		nt->GetEntry(0);                  // get first number of photons for beam start time
		//nentries/=1000;

		for (Long64_t i = 0; i < nentries;) { // Entries loop start -> counter is increased at filling photons
			printProgress(i, nentries, 1, "Analyze pairs :");
			eventmixer.nextEvent(1000*pid.centBin+20*pid.multPhoton+(pid.targetBin+1));
			Int_t nPhot = pid.multPhoton;                      // number of photons per new event
			// ==============Reading and filling one event with photons=========================================//
			vector < Photon > photons;
			photons.reserve(nPhot);

			for (int k = 0; k < nPhot; k++) {
				if (pid.energy >= minEnergy && pid.energy < maxEnergy && (betaCut.isGoodBeta(pid))&&
				    pid.centBin>2&& abs(pid.emc_mq)>10&& abs(pid.qualityDTime)>10) {
					Int_t energyBin = getEnergyBin(pid.energy);
					pid.mom_kine = energyBin;
					pid.energy  *= corrFactors[(int) pid.sector][(int) pid.cell].val[energyBin];

					photons.push_back({ pid });
					nt->GetEntry(++i); // first increase counter than get the next photon in the same event
				} else { nt->GetEntry(++i); }
			}
			// =============================================================================================//

			if (photons.size() == 0) continue;

			// ==============One event combinations========================================================//
			for (size_t e1 = 0; e1 < (photons.size() - 1); e1++) { // first gamma
				Photon photon1 = photons.at(e1);
				for (size_t e2 = e1 + 1; e2 < photons.size(); e2++) { // second gamma
					Photon photon2 = photons.at(e2);
					Diphoton diphot(photon1, photon2);
					if (diphot.opangle * 180. / 3.14 < 8) continue;

					allCells.Fill(diphot, "all");

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
				if (diphot.opangle * 180. / 3.14 < 8) continue;

				allCells.Fill(diphot, "cb");
				ecalCells[gammaMix1.cell].Fill(diphot, "cb");
				ecalCells[gammaMix2.cell].Fill(diphot, "cb");
			}
		}
		// =============================================================================================//
		allCells.Evaluate(iter, corrFactors[0][0]);


		for (Int_t sector = 0; sector < 6; sector++) {
			for (Int_t cell = 0; cell < 163; cell++) {
				if (!isGoodCell(sector, cell)) continue;

				ecalCells[sector * 163 + cell].Evaluate(iter, corrFactors[sector][cell]);
			}
		}
		if (iter == 0)   {
			out->cd(); allCells.Write();
			for (Int_t iSector = 0; iSector < 6; iSector++) {
				for (Int_t iCell = 0; iCell < 163; iCell++) {
					if (!isGoodCell(iSector, iCell)) continue;
					ecalCells[iSector * 163 + iCell].Write();
				}
			}
		}

	}

	out->cd();

	for (Int_t iSector = 0; iSector < 6; iSector++) {
		for (Int_t iCell = 0; iCell < 163; iCell++) {
			if (!isGoodCell(iSector, iCell)) continue;
			//	cout << "Write Sector" << iSector << " Cell " << iCell << endl;
			ecalCells[iSector * 163 + iCell].SetName("after");
			ecalCells[iSector * 163 + iCell].Write();
			cout<<"Sector "<<iSector<<" Cell "<<iCell<<" corr = "<< corrFactors[iSector][iCell].val[0]<<endl;


		}
	}

	allCells.SetName("after");
	allCells.Write();
	writeCorrFactors("corrFactors.txt",corrFactors);
	out->Save();
	out->Close();
} /* fill_hists_diphoton */
