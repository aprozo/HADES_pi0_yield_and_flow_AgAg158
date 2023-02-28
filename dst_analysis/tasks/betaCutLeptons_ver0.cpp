#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/Includes.h"
#include "hparticleevtchara.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/EcalCellPar.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/EcalMap.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/CutValues.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/exp_setup_cells.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/lepton_efficiency.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/hmetamatcher.h"


TH2D * acceptance_mar19;
Bool_t isAcepted(const Float_t & phi, const Float_t & theta)
{
	UInt_t bin = 0;
	bin = acceptance_mar19->FindBin(phi, theta);
	Float_t w = acceptance_mar19->GetBinContent(bin);
	if (theta < 16.5 || theta > 50.) return kFALSE;
	if (w > 0.) return kTRUE;
	return kFALSE;
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

EcalCellPar ecalCellPars[6][163];
void recalibrate(HEmcCluster * cls, HCategory * emcCalCat, Float_t beamTime = 80.1568327436)
{
	Float_t energySum = 0.;
	Int_t sector      = cls->getSector();
	Int_t clusterSize = cls->getNCells();

	for (int i = 0; i < clusterSize; i++) {
		UChar_t c1 = cls->getCell(i);
		Int_t cell = getSystemMap(c1);
		HLocation loc;
		loc.set(2, 0, 0); // set location of cluster
		loc[0] = sector;
		loc[1] = c1;
		HEmcCal * cal     = (HEmcCal *) emcCalCat->getObject(loc);
		Float_t TOT       = cal->getWidth();
		Params newParams  = ecalCellPars[sector][cell].newParams;
		Float_t newEnergy = newParams.adc0 + exp(newParams.adc1 + newParams.adc2 * TOT); // recalibration of a single module energy within cluster
		energySum += newEnergy;
	}
	energySum *= h_energy_correction_exp->GetBinContent(sector * 163 + getSystemMap(cls->getCell()) + 1);
	correctNonlinearity(energySum);
	cls->setEnergy(energySum);
}


Int_t getBeamBin(const Float_t & beamTime)
{
	if (beamTime < 74.4) return 0;
	else if (beamTime < 82.5) return 1;
	else return 2;
}



Int_t ana(TString infileList, TString outputFile, Int_t nev = 1000000.)
{

	Bool_t isSim = kFALSE;
// -------------------------------------------------

// -------------------------------------------------
// list of all files with working sectors
// if(!isSim) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);
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

	if (!loop.setInput("-*,+:HParticleCand,+HEmcCluster,+HParticleEvtInfo,+HParticleCandSim,+HEmcCal,+HEmcClusterSim,+HGeantKine"))
	{
		// if(!loop.setInput("-*,+HParticleCand,+HEmcCluster,+HParticleEvtInfo")) {
		//	if(!loop.setInput("*")) {
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
	HCategory * emcCalCat     = (HCategory *) HCategoryManager::getCategory(catEmcCal);
	HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
	HCategory * evtInfoCat    = (HCategory *) HCategoryManager::getCategory(catParticleEvtInfo);
	HCategory * geantKineCat  = (HCategory *) HCategoryManager::getCategory(catGeantKine);
// -------------------------------------------------

	if (geantKineCat) {
		isSim = kTRUE;
		Info("init()", "GeantKine found - is Simulation");
	} else {
		isSim = kFALSE;
		Info("init()", "GeantKine not found - is not Simulation");
	}

	TH1F*hMomBetaCutEfficiency1SigmaCut  = new TH1F("hMomBetaCutEfficiency1SigmaCut","; mom [MeV/c]; beta eff",100,100,1000);
	TH1F*hMomBetaCutEfficiencyNoCut  = new TH1F("hMomBetaCutEfficiencyNoCut","; mom [MeV/c]; beta eff",100,100,1000);

	TH1F*hEnergyBetaCutEfficiency1SigmaCut  = new TH1F("hEnergyBetaCutEfficiency1SigmaCut","; emc energy [MeV/c]; beta eff",100,100,1000);
	TH1F*hEnergyBetaCutEfficiencyNoCut  = new TH1F("hEnergyBetaCutEfficiencyNoCut",";  emc energy  [MeV/c]; beta eff",100,100,1000);

	TH1F*hMomBetaCutEfficiency1SigmaCutWtihinEnergyResoltion  = new TH1F("hMomBetaCutEfficiency1SigmaCutWtihinEnergyResoltion","; mom [MeV/c]; beta eff",100,100,1000);
	TH1F*hMomBetaCutEfficiencyNoCutWtihinEnergyResoltion  = new TH1F("hMomBetaCutEfficiencyNoCutWtihinEnergyResoltion",";  mom  [MeV/c]; beta eff",100,100,1000);

	// TH2F*hThetaMomElectronsWithEcalCluster  = new TH2F("hThetaMomElectronsWithEcalCluster","hThetaMomElectronsWithEcalCluster; Momentum x Polarity [MeV/c*q];beta",1000,-1000,1000,100,10,50);
	// TH2F*hThetaMomElectronsWithEcalCluster_aCut  = new TH2F("hThetaMomElectronsWithEcalCluster_aCut","hThetaMomElectronsWithEcalCluster_aCut; Momentum x Polarity [MeV/c*q];beta",1000,-1000,1000,100,10,50);
	//





//==============================CORRECTIONS====================================================//
//=============================================================================================//
//=============================================================================================//
//=============================================================================================//

	TString path2corrections = "/lustre/hebe/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/";
	if (readParams(path2corrections + "mar19EmcCalParNewEnergy.txt", ecalCellPars) > 0) return 3;

	if (setupEnergyCorrectionExp(path2corrections + "CorrFactorsPi0Peak.root") > 0) return 4;

	TFile * fCorr = new TFile(path2corrections + "NonlinearityPi0PeakEnergyCorrection.root", "READ");
	hEnergyCorrection = (TH1D *) fCorr->Get("hEnergyCorrection");

	TFile * acc = new TFile(path2corrections + "acceptancePhiThetaMar19.root", "Read");
	acceptance_mar19 = (TH2D *) acc->Get("acceptancePhiThetaMar19");
	setupMassCuts(path2corrections + "m2_cut.root");
	setupSimEffCorr(path2corrections + "simCorrEnergyEcal.root");

	TFile * betaCuts;
	const Int_t nBeamBins = 3;
	TH2F * betaEnergyMean[nBeamBins];
	TH2F * betaEnergySigma[nBeamBins];
	if (isSim) betaCuts = new TFile(path2corrections + "betaBeamMeanSigmaSim.root", "Read");
	else betaCuts = new TFile(path2corrections + "betaBeamMeanSigma.root", "Read");

	for (Int_t iBeamBin = 0; iBeamBin < nBeamBins; iBeamBin++) {
		betaEnergyMean[iBeamBin]  = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergyMean", iBeamBin));
		betaEnergySigma[iBeamBin] = (TH2F *) betaCuts->Get(Form("beam%ibetaEnergySigma", iBeamBin));
		if (betaEnergyMean[iBeamBin] == NULL ||
		    betaEnergySigma[iBeamBin] == NULL) std::cout << "No beta corrections loaded" << std::endl;
	}


	//=============================================================================================//
	//============================CORRECTIONS================================================//
	//=============================================================================================//
	//=============================================================================================//




// -------------------------------------------------


// -------------------------------------------------
// event loop starts here
	Int_t entries = loop.getEntries();
	if (nEvents < entries && nEvents >= 0) entries = nEvents;

	enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
	TH1D * Counter_cent = new TH1D("Counter_cent", "Counter_cent", cNumCounters_cent, 0, cNumCounters_cent);

	Counter_cent->GetXaxis()->SetBinLabel(5, "all_cent");
	Counter_cent->GetXaxis()->SetBinLabel(1, "cent0_10");
	Counter_cent->GetXaxis()->SetBinLabel(2, "cent10_20");
	Counter_cent->GetXaxis()->SetBinLabel(3, "cent20_30");
	Counter_cent->GetXaxis()->SetBinLabel(4, "cent30_40");

	TString tempFileName, currentFileName;
	Int_t dayOfYear, hour, min, second;


	HParticleEvtChara evtChara;
	if (!evtChara.setParameterFile(
			"/lustre/hebe/hades/user/bkardan/param/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root"))
	{
		cout << "HParticleEvtChara: file not found!!! " << endl;
	}
	if (!evtChara.init()) {
		cout << "HParticleEvtChara not init!!! " << endl;
	}
	evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);


// ==========================================================================================//
// ====================================EVENT LOOP============================================//
// ==========================================================================================//


	for (Int_t i = 0; i < entries; i++) {
		HTool::printProgress(i, entries, 1, "Analyze pairs :");
		Int_t nbytes = loop.nextEvent(i); // get next event. categories will be cleared before

		if (nbytes <= 0) { cout << nbytes << endl; break; } // last event reached
		// if(i%1000 == 0) cout<<"event "<<i<<endl;

		if (loop.isNewFile(tempFileName)) {
			TString type;
			Int_t year;
			Int_t eb;
			currentFileName = tempFileName;
			HTime::splitFileName(HTime::stripFileName(tempFileName), type, year, dayOfYear, hour, min, second, eb);
		}
		Float_t valueMinute = dayOfYear + hour / 24. + min / 24. / 60.;
		Float_t valueSecond = valueMinute + second / 24. / 60. / 60.;
		Float_t beamTime    = valueSecond;
		if (isSim) beamTime = 0.;
		if (!isGoodBeamTime(beamTime)) continue;
		// -------------------------------------------------
		// summary event info object

		HGeomVector vertex(evHeader->getVertexCluster().getPos()); // Event vertex position
		HParticleEvtInfo * evtInfo = HCategoryManager::getObject(evtInfo, catParticleEvtInfo, 0);
		Int_t mult_tofrpc = evtInfo->getSumTofMultCut() + evtInfo->getSumRpcMultHitCut();

		Int_t cent_bin = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);
		Float_t ImpactParameter = 0;


		if (cent_bin > 3) continue;
		cent_bin = cent_bin - 1; // starting from 0, not 1

		Counter_cent->Fill(cent_bin);
		Counter_cent->Fill(4);

		// -------------------------------------------------
		// loop over particle candidates in event
		if (candCat) {
			// ------------------------------------------------------------------------
			// clean vectors and index arrays
			sorter.cleanUp();
			// if you do not want to modify the flags
			// for the following tasks etc (needs call of
			// restoreFlags() at the end of the eventloop)
			// sorter.backupFlags();
			// ------------------------------------------------------------------------

			sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
			sorter.fill(selectLeptonsBetaLoose);     // fill only good loose leptons
			Int_t nCandLepBest        = sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
			Int_t size                = candCat->getEntries();
			HParticleCandClass * cand = 0;

			for (Int_t j = 0; j < size; j++) {
				// ------------------------------------------------------------------------
				cand = HCategoryManager::getObject(cand, candCat, j);
				if (!leptonParticlecandCut(cand)) continue;     // see  #include "lepton_efficiency.h"
				if (!cand->isFlagBit(kIsUsed)) continue;     // skip particles rejected by HParticleTrackSorter
				Int_t emc_index = cand->getEmcInd();
				if (emc_index<0) continue;
				HEmcClusterClass * cls = (HEmcClusterClass *) emcClusterCat->getObject(emc_index);
				Int_t sector         = cls->getSector();
				Int_t cell = getSystemMap(cls->getCell());
				if (!isGoodCell(sector, cell)) continue;


				Float_t mom = cand->getMomentum();
				Float_t beta   = cand->getBeta();
				Float_t charge = cand->getCharge();
				Float_t mm_qa  = cand->getMetaMatchQuality();

				Float_t time      = cls->getTime();
				HGeomVector trackVec(cls->getXLab()-vertex.X(), cls->getYLab()-vertex.Y(), cls->getZLab() - vertex.Z());
				Double_t trackLength = trackVec.length();
				Double_t emc_beta        = (trackLength / 1000.) / (time * 1.e-9) / TMath::C();
				Int_t size = cls->getNCells();

				if (!isSim) recalibrate(cls, emcCalCat,beamTime);
				Float_t emc_energy = cls->getEnergy();
				if (isSim) emc_energy = correctEMCenergy(energy, sector, cell);
				if(emc_energy<100) continue;



				Float_t emc_mq     = cand->getMetaMatchQualityEmc();
				Float_t emc_path   = cand->getEmcCellPath();
				Float_t emc_time   = cand->getEmcTime();

				if(emc_mq>4) continue;

				Int_t beamBin   = getBeamBin(beamTime);            // cout<<"beamBin="<<beamBin<<endl;
				Int_t energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(emc_energy);
				if (emc_energy >= 800) energyBin = betaEnergyMean[beamBin]->GetXaxis()->FindBin(750);
				Float_t betaMean = betaEnergyMean[beamBin]->GetBinContent(energyBin,sector * 200 + cell + 1);
				if (emc_energy < 200) betaMean = betaEnergyMean[beamBin]->GetBinContent(energyBin + 1,sector * 200 + cell + 1);
				Float_t betaSigma = betaEnergySigma[beamBin]->GetBinContent(energyBin,sector * 200 + cell + 1);
				Float_t betCutDown = betaMean - 1 * betaSigma;
				Float_t betaCutUp  = betaMean + 1 * betaSigma;

				Float_t withinSigma = abs(emc_beta-betaMean)/betaSigma;
				Float_t energyDiff = (100 * abs(emc_energy - mom) / emc_energy);
				Float_t withinResolution = energyDiff * sqrt(mom / 1000.) / 3;


				hEnergyBetaCutEfficiencyNoCut->Fill(emc_energy);
				hMomBetaCutEfficiencyNoCut->Fill(mom);
				if(withinResolution<5.9) hMomBetaCutEfficiencyNoCutWtihinEnergyResoltion->Fill(mom);

				if(withinSigma<=1.) {
					hMomBetaCutEfficiency1SigmaCut->Fill(mom);
					hEnergyBetaCutEfficiency1SigmaCut->Fill(emc_energy);
					if(withinResolution<5.9) hMomBetaCutEfficiency1SigmaCutWtihinEnergyResoltion->Fill(mom);
				}

				// -------------------------------------------------

			} // end cand loop
		} // end cand cat
		  // -------------------------------------------------
	} // end eventloop


// -------------------------------------------------
// added for the task
// if(gHades)gHades->finalizeTasks();
// -------------------------------------------------
	out->cd();
	Counter_cent->Write();
	out->Write();
	delete gHades;
	return 0;
} /* fill_ttree_diphoton */
