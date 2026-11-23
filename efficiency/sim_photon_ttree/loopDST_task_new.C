#include "hades.h"
#include "hloop.h"
#include "htaskset.h"
#include "hdst.h"
#include "htime.h"
#include "hcategory.h"
#include "hmatrixcategory.h"
#include "hcategorymanager.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticletracksorter.h"
#include "hparticlebooker.h"
#include "hparticletool.h"
#include "hparticledef.h"
#include "rpcdef.h"
#include "hparticleevtinfo.h"
#include "hrpchit.h"
#include "hrpccluster.h"
#include "hlocation.h"
#include "henergylosscorrpar.h"
#include "hrpcgeompar.h"
#include "hphysicsconstants.h"

#include "exp_setup_cells.h"
#include "lepton_efficiency.h"


#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TMath.h"
#include "hmetamatcher.h"
#include "TVector3.h"
#include "Includes.h"
#include <iostream>
#include <vector>
using namespace std;
TH2D * acceptance_mar19;

Bool_t isInAcc_mar19(const Float_t & theta, const Float_t & phi)
{
	UInt_t bin = 0;
	bin = acceptance_mar19->FindBin(phi, theta);
	Float_t w = acceptance_mar19->GetBinContent(bin);
	if (w ==1) return kTRUE;
	else return kFALSE;
}

Float_t impact_bins[] = { 0., 3.8, 5.4, 6.6, 7.6 }; // latest version from WIKI


Int_t getCentBin(const Double_t & b)
{
	Int_t impBin = 5; // Hades binning

	if (b >= impact_bins[0] && b < impact_bins[1]) impBin = 1;
	if (b >= impact_bins[1] && b < impact_bins[2]) impBin = 2;
	if (b >= impact_bins[2] && b < impact_bins[3]) impBin = 3;
	if (b >= impact_bins[3] && b < impact_bins[4]) impBin = 4;
	return impBin;
}

typedef struct {
	Int_t cent_bin;
	Int_t sector;
	Int_t emc_cell;
	Int_t rowEMC;
	Int_t columnEMC;

	Int_t id_geant;
	Float_t theta_kine;
	Float_t phi_kine;

	Float_t theta_ecal;
	Float_t phi_ecal;

	Float_t x_kine;
	Float_t y_kine;
	Float_t mom_kine;
	Int_t matchedCells;
	Int_t matchedPart;
	Float_t emc_energy;
	Float_t emc_max_energy;
	Float_t within_resolution;
	Float_t emc_beta;
	Int_t emc_size;
	Float_t xEMC;
	Float_t yEMC;

	void reset()
	{
		// cent_bin           = -1;
		sector     = -1;
		emc_cell   = -1;
		theta_kine = -1.;
		phi_kine   = -1.;

		theta_ecal        = -1.;
		phi_ecal          = -1.;
		mom_kine          = -1.;
		id_geant          = -1;
		x_kine            = -99999.;
		y_kine            = -99999.;
		mom_kine          = 0.;
		matchedCells      = -1;
		matchedPart       = -1;
		rowEMC            = -1;
		columnEMC         = -1;
		emc_energy        = -1.;
		emc_max_energy    = -1.;
		within_resolution = 100.;
		emc_beta          = -1.;
		emc_size          = -1;
		xEMC = -99999.;
		yEMC = -99999.;
	}
} TREE_VARIABLES;

TTree * makeTree(TREE_VARIABLES& vars)
{
	TTree * nt = new TTree("nt", "Simulated tree for prediction of EMC cell");

	nt->Branch("cent_bin", &vars.cent_bin, "cent_bin/I");
	nt->Branch("sector", &vars.sector, "sector/I");
	nt->Branch("emc_cell", &vars.emc_cell, "emc_cell/I");
	nt->Branch("rowEMC", &vars.rowEMC, "rowEMC/I");
	nt->Branch("columnEMC", &vars.columnEMC, "columnEMC/I");
	nt->Branch("id_geant", &vars.id_geant, "id_geant/I");
	nt->Branch("theta_kine", &vars.theta_kine, "theta_kine/F");
	nt->Branch("phi_kine", &vars.phi_kine, "phi_kine/F");
	nt->Branch("theta_ecal", &vars.theta_ecal, "theta_ecal/F");
	nt->Branch("phi_ecal", &vars.phi_ecal, "phi_ecal/F");
	nt->Branch("x_kine", &vars.x_kine, "x_kine/F");
	nt->Branch("y_kine", &vars.y_kine, "y_kine/F");
	nt->Branch("mom_kine", &vars.mom_kine, "mom_kine/F");
	nt->Branch("matchedCells", &vars.matchedCells, "matchedCells/I");
	nt->Branch("matchedPart", &vars.matchedPart, "matchedPart/I");
	nt->Branch("emc_energy", &vars.emc_energy, "emc_energy/F");
	nt->Branch("emc_max_energy", &vars.emc_max_energy, "emc_max_energy/F");
	nt->Branch("within_resolution", &vars.within_resolution, "within_resolution/F");
	nt->Branch("emc_beta", &vars.emc_beta, "emc_beta/F");
	nt->Branch("emc_size", &vars.emc_size, "emc_size/I");
	nt->Branch("xEMC", &vars.xEMC, "xEMC/F");
	nt->Branch("yEMC", &vars.yEMC, "yEMC/F");

	return nt;
}

void fillEntry(HEmcClusterSim * cls, TREE_VARIABLES& vars)
{
	vars.emc_cell   = getSystemMap(cls->getCell());
	vars.columnEMC  = cls->getCell() % 17;
	vars.rowEMC     = floor(cls->getCell() / 17);
	vars.sector     = cls->getSector();
	vars.emc_size   = cls->getNCells();
	vars.xEMC       = cls->getXMod();
	vars.yEMC       = cls->getYMod();
	vars.phi_ecal   = cls->getPhi();
	vars.theta_ecal = cls->getTheta();


	Float_t clustEnergy = correctEMCenergy(cls->getEnergy(), cls->getCell());
	vars.emc_max_energy = correctEMCenergy(cls->getMaxEnergy(), cls->getCell());
	vars.emc_energy     = clustEnergy;

	Float_t energydiff = (100 * abs(clustEnergy - vars.mom_kine) / clustEnergy);
	vars.within_resolution = energydiff * sqrt(vars.mom_kine / 1000.) / 3;
}

Int_t loopDST_task(TString infileList, TString outfile, Int_t nEvents)
{
	//  infileList : comma seprated file list "file1.root,file2.root" or "something*.root"
	//  outfile    : optional (not used here) , used to store hists in root file
	//  nEvents    : number of events to processed. if  nEvents < entries or < 0 the chain will be processed
	Bool_t isSimulation = kFALSE;

	// -------------------------------------------------
	// create loop obejct and hades
	HLoop loop(kTRUE);
	gHades->setBeamTimeID(HADES::kMar19);
	// -------------------------------------------------
	// -------------------------------------------------
	// setup the spectrometer for filter tasks
	// this part is needed if your task needs
	// detectors or parameter container
	if (1) {
		TString beamtime    = "mar19";
		Int_t mdcMods[6][4] =
		{ { 1, 1, 1, 1 },
		  { 1, 1, 1, 1 },
		  { 1, 1, 1, 1 },
		  { 1, 1, 1, 1 },
		  { 1, 1, 1, 1 },
		  { 1, 1, 1, 1 } };
		TString asciiParFile = "";
		TString rootParFile  = "/cvmfs/hadessoft.gsi.de/param/sim/mar19/allParam_mar19_sim_run_16000_gen4_25112020.root";
		TString paramSource  = "root";       // root, ascii, oracle
		TString paramrelease = "MAR19_dst_gen5"; // now, APR12_gen2_dst APR12_gen5_dst
		HDst::setupSpectrometer(beamtime, mdcMods, "rich,mdc,tof,rpc,emc,wall,start,tbox");
		HDst::setupParameterSources(paramSource, asciiParFile, rootParFile, paramrelease); // now, APR12_gen2_dst
	}
	// -------------------------------------------------

	// -------------------------------------------------
	// list of all files with working sectors
	// if(!isSimulation) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);
	// -------------------------------------------------

	// -------------------------------------------------
	// reading input files and decalring containers
	Bool_t ret = kFALSE;
	if (infileList.Contains(",")) {
		ret = loop.addMultFiles(infileList); // file1,file2,file3
	} else {
		ret = loop.addFiles(infileList); // myroot*.root
	}
	if (ret == 0) {
		cout << "READBACK: ERROR : cannot find inputfiles : " << infileList.Data() << endl;
		return 1;
	}
	TFile * out = new TFile(outfile, "RECREATE");
	out->cd();
	// ==================================================================================================================//
	// Output                                                                                                           //
	// ==================================================================================================================//
	//
	TREE_VARIABLES vars;       //
	TTree * nt = makeTree(vars); //
	// ------------------------------------------------------------------------------------------------------------------//


	if (!loop.setInput("-*,+HParticleCandSim,+HEmcClusterSim,+HGeantKine")) {
		cout << "READBACK: ERROR : cannot read input !" << endl;
		exit(1);
	} // read all categories
	loop.printCategories();
	loop.printChain();
	// -------------------------------------------------

	// -------------------------------------------------
	HParticleTrackSorter sorter;
	sorter.init(); // get catgegory pointers etc...
	// -------------------------------------------------

	// -------------------------------------------------
	// input data
	HEventHeader * evHeader   = loop.getEventHeader();
	HCategory * candCat       = (HCategory *) HCategoryManager::getCategory(catParticleCand);
	HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
	HCategory * geantKineCat  = (HCategory *) HCategoryManager::getCategory(catGeantKine);
	// -------------------------------------------------
	// -------------------------------------------------
	// parameters
	TFile * acc = new TFile("acceptancePhiThetaMar19.root", "Read");
	acceptance_mar19 = (TH2D *) acc->Get("acceptancePhiThetaMar19");
	setupSimEffCorr("corretion_energy_ag1650ag.root");
	HPartialEvent * fSimul    = NULL;
	HGeantHeader * fSubHeader = NULL;
	// -------------------------------------------------
	// -------------------------------------------------
	// added for the task
	HTaskSet * masterTaskSet = gHades->getTaskSet("all");
	HMetaMatcher * matcher   = new HMetaMatcher();
	// matcher->setDebug();
	masterTaskSet->add(matcher);
	// -------------------------------------------------
	// -------------------------------------------------
	// centrality histogramm
	enum Counters_cent { all_cents, cent0_10, cent10_20, cent20_30, cent30_40, cNumCounters_cent };
	TH1D * Counter_cent = new TH1D("EventNum", "EventNum", cNumCounters_cent, 0, cNumCounters_cent);

	Counter_cent->GetXaxis()->SetBinLabel(1, "all_cent");
	Counter_cent->GetXaxis()->SetBinLabel(2, "cent0_10");
	Counter_cent->GetXaxis()->SetBinLabel(3, "cent10_20");
	Counter_cent->GetXaxis()->SetBinLabel(4, "cent20_30");
	Counter_cent->GetXaxis()->SetBinLabel(5, "cent30_40");

	TH1D * hMom_kine    = new TH1D("hMom_kine", "mom_kine", 200, 0, 1000);
	TH1D * hMom_reco    = new TH1D("hMom_reco", "hMom_reco", 200, 0, 1000);
	TH1D * hEnergy_reco = new TH1D("hEnergy_reco", "hEnergy_reco", 200, 0, 1000);
	TH1D * cells        = new TH1D("cells", "cells", 163, 0, 163);

	// -------------------------------------------------
	// -------------------------------------------------
	// event loop starts here
	Int_t entries = loop.getEntries();
	if (nEvents < entries && nEvents >= 0) entries = nEvents;
	for (Int_t i = 0; i < entries; i++) {
		Int_t nbytes = loop.nextEvent(i); // get next event. categories will be cleared before
		if (nbytes <= 0) { cout << nbytes << endl; break; } // last event reached
		// if(i%1000 == 0) cout<<"event "<<i<<endl;
		HGeomVector vertex(evHeader->getVertexCluster().getPos()); // Event vertex position
		fSimul     = ( (HRecEvent *) gHades->getCurrentEvent())->getPartialEvent(catSimul);
		fSubHeader = (HGeantHeader *) (fSimul->getSubHeader());
		Float_t impactParameter = -1.;
		if (fSubHeader) impactParameter = fSubHeader->getImpactParameter();
		UInt_t cent_bin = 0;
		if (impactParameter < 3.8) cent_bin = cent0_10;
		if (impactParameter >= 3.8 && impactParameter < 5.4) cent_bin = cent10_20;
		if (impactParameter >= 5.4 && impactParameter < 6.6) cent_bin = cent20_30;
		if (impactParameter >= 6.6) continue;
		vars.cent_bin = cent_bin;

		Counter_cent->Fill(vars.cent_bin);
		Counter_cent->Fill(all_cents);

		Int_t emc_ModMult[6][163] = { 0 };
		sorter.cleanUp();
		sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
		Int_t nCandLep     = sorter.fill(HParticleTrackSorter::selectLeptons);// fill only good leptons
		Int_t nCandLepBest = sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
		Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);// fill only good hadrons (already marked good leptons will be skipped)
		Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);
		HParticleCand * cand;

		for (Int_t j = 0; j < candCat->getEntries(); j++) { // loop over Particlecand
			cand = HCategoryManager::getObject(cand, candCat, j);
			if (!cand->isFlagBit(Particle::kIsUsed)) continue;
			Int_t emc_index = cand->getEmcInd();
			Double_t emc_mq = cand->getMetaMatchQualityEmc();
			if (emc_index >= 0) {
				HEmcClusterSim * emcClst = (HEmcClusterSim *) emcClusterCat->getObject(emc_index);
				Int_t sec  = emcClst->getSector();
				Int_t cell = getSystemMap(emcClst->getCell());
				emc_ModMult[sec][cell]++;
			}
		}

		// -------------------------------------------------
		// geant kine loop starts here
		for (Int_t j = 0; j < geantKineCat->getEntries(); j++) {
			HGeantKine * photon_kine = HCategoryManager::getObject(photon_kine, geantKineCat, j);
			if (photon_kine->getID() != 1) continue;

			HGeantKine * mother = photon_kine->getParent(photon_kine->getTrack());
			if (!(mother && mother->getID() == 7 && mother->isPrimary())) continue;
			Float_t mom_kine = photon_kine->getTotalMomentum();
			Float_t phi_kine, theta_kine;
			photon_kine->getPhiThetaDeg(theta_kine, phi_kine);
			if (!isInAcc_mar19(theta_kine, phi_kine)) continue;

			hMom_kine->Fill(mom_kine);
			// theta_kine*=TMath::DegToRad();
			// phi_kine=TMath::DegToRad()*(fmod(phi_kine,60.f)+60.);
			vars.reset();
			vars.id_geant   = 1;
			vars.theta_kine = theta_kine;
			vars.phi_kine   = phi_kine;
			vars.mom_kine   = mom_kine;
			vars.x_kine     = sin(theta_kine*TMath::DegToRad()) * sin(phi_kine*TMath::DegToRad());
			vars.y_kine     = sin(theta_kine*TMath::DegToRad()) * cos(phi_kine*TMath::DegToRad());
			vars.sector     = photon_kine->getSector();

			Int_t photon_kinetrack = photon_kine->getTrack();
			Int_t geantid = 0;
			HEmcClusterSim * cls;
			Float_t length_to_cell = 0;
			Double_t Beta = 0;
			Int_t matchedPart = -1, matchedCells = -1;


			Bool_t fromparent = kFALSE; // flag that emc saw a hit from parent photon
			Bool_t directhit  = kFALSE;

			for (Int_t k = 0; k < emcClusterCat->getEntries(); k++) {
				cls = HCategoryManager::getObject(cls, emcClusterCat, k); // get an cluster emc object
				HGeomVector cellVector(cls->getXLab(), cls->getYLab(), cls->getZLab());
				Int_t cls_track = cls->getTrack(); // get the track of cluster emc object
				Int_t emc_cell  = getSystemMap(cls->getCell());// get cell of cluster emc object;
				Int_t sector    = cls->getSector();// get cell of cluster emc object;

				Double_t trLength = cellVector.length();
				Beta         = (trLength / 1000.) / (cls->getTime() * 1.e-9) / TMath::C();
				matchedPart  = emc_ModMult[sector][emc_cell];
				matchedCells = cls->getNMatchedCells();


				HGeantKine * kineEcal = HCategoryManager::getObject(kineEcal, geantKineCat, cls_track - 1); // find a kine track corresponding to that cluster emc object
				Int_t kineEcal_track  = kineEcal->getTrack();

				if (kineEcal_track == photon_kinetrack) {
					directhit = kTRUE;
					fillEntry(cls, vars);
					vars.matchedPart  = matchedPart;
					vars.matchedCells = matchedCells;
					vars.emc_beta     = Beta;
					break;
				}

				HGeantKine * prntKineEcal = kineEcal->getParent(kineEcal_track);
				Int_t trackParentKineEcal = prntKineEcal != NULL ? prntKineEcal->getTrack() : 0; // get a track  of parnet of geant emc object track
				if (trackParentKineEcal && trackParentKineEcal == photon_kinetrack) {    // parent track == photon track  - flag that emc saw a hit from parent photon
					HGeomVector prntKineEcalVertex;
					kineEcal->getVertex(prntKineEcalVertex);
					prntKineEcalVertex -= cellVector;
					length_to_cell      = prntKineEcalVertex.length();
					if (length_to_cell < 330) {
						fromparent    = kTRUE;
						vars.id_geant = kineEcal->getID();
						fillEntry(cls, vars);
						vars.matchedPart  = matchedPart;
						vars.matchedCells = matchedCells;
						vars.emc_beta     = Beta;
						break;
					}
				}
			}

			nt->Fill();

			// cell_mom_kine->Fill(pid.emc_cell,mom_kine);

			// end of loop over Cluster category
		} // end of loop over Geant Kine category
		  // -------------------------------------------------
	} // end eventloop
	  // -------------------------------------------------
	  // added for the task
	// if(gHades)gHades->finalizeTasks();
	// -------------------------------------------------
	out->cd();
	Counter_cent->Write();

	/*
	 *      cell_mom_kine->Write();
	 *      cell_energy_reco->Write();
	 *      cell_mom_reco->Write();
	 */
	// hMom_kine->Write();
	// hEnergy_reco->Write();
	// hMom_reco->Write();
	cells->Write();

	out->Write();
	delete gHades;
	return 0;
} /* loopDST_task */
