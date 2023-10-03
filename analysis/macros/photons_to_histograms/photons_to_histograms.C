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

#include <iostream>
#include <map>
#include <stdio.h>
#include <iomanip>
#include "../GlobalParameters.h"


#include "heventmixer.h" // Szymon's event mixer, a bit improved
#include "exp_setup_cells.h"
using namespace std;
TH1D * hEnergyCorrection;
Bool_t isSim;
#include "HistogramManager.h" // classes for PID Photon Diphoton and HistogramManager

void correctNonlinearity(Float_t & energy)
{
	Double_t corr;

	if (energy < 110) corr = hEnergyCorrection->Interpolate(110);
	if (energy > 640) corr = hEnergyCorrection->Interpolate(640);
	else corr = hEnergyCorrection->Interpolate(energy);
	energy *= corr;
}

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
		Counter_cent->Add((TH1D *) temp->Get("Counter_cent_cut"));
		if((TH1D *)temp->Get("Counter_cent_cut")==NULL) Counter_cent->Add((TH1D *) temp->Get("Counter_cent"));
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

Int_t getBeamBin(const Float_t & beamtime)
{
	if (beamtime < 74.4) return 0;
	else if (beamtime < 82.5) return 1;
	else return 2;
}

Int_t photons_to_histograms(TString inFile, TString outputFile, Bool_t isExperiment = kTRUE)
{
	// =====================================================================================
	//   R E A D I N G   P A R A M E T E R S
	// =====================================================================================

	// -------------------------------------------

	GlobalParameters par("../parameters.conf");
	// -------------------------------------------
	TFile * fCorr = new TFile(Form("%s/NonlinearityPi0PeakEnergyCorrection.root", par.PATH_TO_CORRECTION_DATA), "READ");
	hEnergyCorrection = (TH1D *) fCorr->Get("hEnergyCorrection");

	TFile * out = new TFile(outputFile, "RECREATE");
	if (VERBOSITY > 0) cout << "Output file created\n";

	TH1::SetDefaultSumw2();

	enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
	TH1D * Counter_cent = new TH1D("Counter_cent_cut", "Counter_cent_cut", cNumCounters_cent, 0, cNumCounters_cent);

	Counter_cent->GetXaxis()->SetBinLabel(5, "all_cent");
	Counter_cent->GetXaxis()->SetBinLabel(1, "cent0_10");
	Counter_cent->GetXaxis()->SetBinLabel(2, "cent10_20");
	Counter_cent->GetXaxis()->SetBinLabel(3, "cent20_30");
	Counter_cent->GetXaxis()->SetBinLabel(4, "cent30_40");


	TH1D * nEvents = new TH1D("nEvents", "nEvents", cNumCounters_cent, 0, cNumCounters_cent);
	nEvents->GetXaxis()->SetBinLabel(5, "all_cent");
	nEvents->GetXaxis()->SetBinLabel(1, "cent0_10");
	nEvents->GetXaxis()->SetBinLabel(2, "cent10_20");
	nEvents->GetXaxis()->SetBinLabel(3, "cent20_30");
	nEvents->GetXaxis()->SetBinLabel(4, "cent30_40");


	TDirectory * subD1 = out->mkdir("All");
	TDirectory * subD2 = out->mkdir("Mix");
	out->cd();

	if (VERBOSITY > 0) cout << "Directories created\n";

	// HistogramManager HistMan(subD1, par.N_PT_BINS, par.N_Y_BINS, par.N_MASS_BINS);
	// HistogramManager HistManMix(subD2, par.N_PT_BINS, par.N_Y_BINS, par.N_MASS_BINS);

	HistogramManager HistMan(subD1, par );
	HistogramManager HistManMix(subD2, par );


	if (VERBOSITY > 0) cout << "HistogramManagers created\n";


	TFile * FileEff = new TFile(Form("%s/embeddingAccEffMatrices.root", par.PATH_TO_CORRECTION_DATA), "Read");
	photon_eff[0] = (TH1F *) FileEff->Get("cent0_correctionMomEmbedding");
	photon_eff[1] = (TH1F *) FileEff->Get("cent1_correctionMomEmbedding");
	photon_eff[2] = (TH1F *) FileEff->Get("cent2_correctionMomEmbedding");

	TFile * fileBetaCutEff = new TFile(Form("%s/betCutExpSimRatio.root", par.PATH_TO_CORRECTION_DATA), "Read");
	betaCutEff[0] = (TH1F *) fileBetaCutEff->Get("BetaCutioCent0");
	betaCutEff[1] = (TH1F *) fileBetaCutEff->Get("BetaCutioCent1");
	betaCutEff[2] = (TH1F *) fileBetaCutEff->Get("BetaCutioCent2");

  TFile * fileAccEffEmbedding= new TFile(Form("%s/accEffMatrices.root", par.PATH_TO_CORRECTION_DATA), "Read");
	accEffEmbedding[0] = (TH1F *) fileAccEffEmbedding->Get("cent0_acceptanceEff");
	accEffEmbedding[1] = (TH1F *) fileAccEffEmbedding->Get("cent1_acceptanceEff");
	accEffEmbedding[2] = (TH1F *) fileAccEffEmbedding->Get("cent2_acceptanceEff");



	if (VERBOSITY > 0 && photon_eff[0] != NULL) cout << "Read embedding succefully!\n";

	TChain * nt = new TChain("nt");
	if (VERBOSITY > 0) cout << "Chain init\n";
	addMultFiles(inFile, nt, Counter_cent);

	if (VERBOSITY > 0) cout << "Files added\n";
	PID pid;
	nt->SetBranchAddress("pid", &pid.newEvent);
	if (VERBOSITY > 0) cout << "Photon init\n";
	Long64_t nentries = nt->GetEntries();

	HGenericEventMixerObj < Photon > eventmixer;
	eventmixer.setPIDs(1, 1, 7);
	eventmixer.setBuffSize(30);

	if (VERBOSITY > 0) cout << "EventMixer set\n";

	nt->GetEntry(0); // get first number of photons for beam start time
	nt->Show();
	isSim = (pid.mom_kine > 0);


	TFile * betaCuts;
	const Int_t nBeamBins = 3;
	TH2F * betaEnergyMean[nBeamBins];
	TH2F * betaEnergySigma[nBeamBins];
	if (isSim) betaCuts = new TFile(Form("%s/betaBeamMeanSigmaSim.root", par.PATH_TO_CORRECTION_DATA), "Read");
	else betaCuts = new TFile(Form("%s/betaBeamMeanSigma.root", par.PATH_TO_CORRECTION_DATA), "Read");
	for (Int_t iBeamBin = 0; iBeamBin < nBeamBins; iBeamBin++) {
		betaEnergyMean[iBeamBin]  = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergyMean", iBeamBin));
		betaEnergySigma[iBeamBin] = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergySigma", iBeamBin));
		if (betaEnergyMean[iBeamBin] == NULL ||
		    betaEnergySigma[iBeamBin] == NULL) std::cout << "No beta corrections loaded" << std::endl;
	}


	// Float_t BETA_MIDDLE   = 0.97;
	// if(isSim) BETA_MIDDLE=1.0;
	// Float_t BETA_CUT_DOWN = BETA_MIDDLE-TMath::Abs(par.BETA_CUT-BETA_MIDDLE);
	// Float_t BETA_CUT_UP   = BETA_MIDDLE+TMath::Abs(par.BETA_CUT-BETA_MIDDLE);

	Float_t RPC_TIME_MIDDLE = -0.30;
	if (isSim) RPC_TIME_MIDDLE = 0.0;

	Float_t RPC_TIME_UP   = RPC_TIME_MIDDLE + TMath::Abs(par.RPC_TIME_CUT - RPC_TIME_MIDDLE);
	Float_t RPC_TIME_DOWN = RPC_TIME_MIDDLE - TMath::Abs(par.RPC_TIME_CUT - RPC_TIME_MIDDLE);


	if (VERBOSITY > 0) cout << "0 event is taken\n";
	
	for (Long64_t i = 1; i < nentries;) {
		printProgress(i, nentries, 1, "Analyze pairs :");

		if (!(isGoodBeamTime(pid.beamTime) || isSim)  ) { nt->GetEntry(++i); continue; }
		if (pid.newEvent == 1) nEvents->Fill(pid.centBin);


		Int_t nPhot = pid.multPhoton; // number of photons per new event
		if (nPhot == 0)  { nt->GetEntry(++i); continue; }

		eventmixer.nextEvent(1000 * pid.centBin + 20 * pid.multPhoton + (pid.targetBin + 1)); // event mixer classes -> same centrality(10% classes), photon multiplicity and target bin
		if (VERBOSITY > 1) cout << "event " << i << ": " << nPhot << " photons" << endl;
		// ======Reading and filling one event with photons=========================================//
		vector < Photon > photons;
		photons.reserve(nPhot);
		for (int k = 0; k < nPhot; k++) {
			if (!isSim) correctNonlinearity(pid.energy);
			if (
				isGoodCell((Int_t) pid.sector, (Int_t) pid.cell) &&
				pid.size >= 1 &&
				pid.size < par.CLUSTER_SIZE_CUT &&
				pid.energy > par.ENERGY_CUT &&
				pid.energy < 2000. &&
				pid.centBin >= 0 &&
				pid.centBin < 3 &&

				((pid.sector == 0 && par.SECTOR_0) ||
				 (pid.sector == 1 && par.SECTOR_1) ||
				 (pid.sector == 2 && par.SECTOR_2) ||
				 (pid.sector == 3 && par.SECTOR_3) ||
				 (pid.sector == 4 && par.SECTOR_4) ||
				 (pid.sector == 5 && par.SECTOR_5))
				&&
				pid.cell >= 10 &&
				(pid.qualityDTime > RPC_TIME_UP ||
				 pid.qualityDTime < RPC_TIME_DOWN)
				&&
				pid.emc_mq > par.RPC_XY_CUT
				)
			{
				Int_t beamBin   = getBeamBin(pid.beamTime);// cout<<"beamBin="<<beamBin<<endl;
				Int_t energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(pid.energy);
				if (pid.energy < 100) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(150);
				if (pid.energy >= 800) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(750);
				// cout<<"energyBin="<<energyBin<<endl;

				Float_t BETA_MEAN = betaEnergyMean[beamBin]->GetBinContent(energyBin,
				                                                           (int) pid.sector * 200 + (int) pid.cell + 1);
				if (pid.energy < 200) BETA_MEAN = betaEnergyMean[beamBin]->GetBinContent(energyBin+1,
					                                                                     (int) pid.sector * 200 + (int) pid.cell + 1);


				Float_t BETA_SIGMA = betaEnergySigma[beamBin]->GetBinContent(energyBin,
				                                                             (int) pid.sector * 200 + (int) pid.cell + 1);
				Float_t BETA_CUT_DOWN = BETA_MEAN - par.BETA_CUT * BETA_SIGMA;
				Float_t BETA_CUT_UP   = BETA_MEAN + par.BETA_CUT * BETA_SIGMA;

				if (pid.beta > BETA_CUT_DOWN && pid.beta < BETA_CUT_UP) {
					photons.push_back({ pid });
				}
				nt->GetEntry(++i); // first increase counter than get next photon in te same event
			} else { nt->GetEntry(++i); }
		}
		// =====================================================================================//
		if (photons.empty()) continue;
		// ==============One event combinations========================================================//
		for (size_t e1 = 0; e1 < (photons.size() - 1); e1++) { // first gamma
			Photon photon1 = photons.at(e1);
			for (size_t e2 = e1 + 1; e2 < photons.size(); e2++) { // second gamma
				Photon photon2 = photons.at(e2);
				Diphoton diphot(photon1, photon2);			
				HistMan.Hfill(diphot);
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
			HistManMix.Hfill(diphot);
		}
	}
	// =============================================================================================//
	if (VERBOSITY > 0) cout << "Write results to the file\n";
	out->cd();
	Counter_cent->Write();
	nEvents->Write();
	HistMan.Write(subD1);
	HistManMix.WriteMix(subD2);
	out->Save();
	out->Close();

	return 1;
} /* photons_to_histograms */
