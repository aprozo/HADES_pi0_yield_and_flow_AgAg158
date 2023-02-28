#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/Includes.h"
#include "hparticleevtchara.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/EcalMap.h"
#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/CutValues.h"

TH2D * acceptance_mar19;
Bool_t isAcepted(const Float_t & phi, const Float_t & theta)
{
	UInt_t bin = 0;

	bin = acceptance_mar19->FindBin(phi, theta);
	Float_t w = acceptance_mar19->GetBinContent(bin);
	if (theta < 16.5||theta>50.) return kFALSE;

	if (w > 0.) return kTRUE;

	return kFALSE;
}

Bool_t isKineMatched(HEmcClusterSim * cls, Float_t trackKine, HCategory * geantKineCat)
{
	for (Short_t i = 0; i < 1; i++) {
		Int_t clusterTrack       = cls->getTrack(i) - 1;
		HGeantKine * clusterKine = (HGeantKine *) geantKineCat->getObject(clusterTrack);
		if (clusterKine == NULL) continue;
		if (clusterKine->getTrack() == trackKine || (clusterKine->getParentTrack() == trackKine)) {
			return kTRUE;
		}
	}
	return kFALSE;
}

struct Photon {
	Photon(HGeantKine * g)
	{
		isAcc         = kFALSE;
		isRec         = kFALSE;
		isRecWithCuts = kFALSE;
		setKine(g);
	}

	void setKine(HGeantKine * g)
	{
		kine = g;
		Float_t px1, py1, pz1, thetaKine, phiKine;
		kine->getMomentum(px1, py1, pz1);
		tlvKine.SetXYZM(px1, py1, pz1, kine->getM());
		kine->getPhiThetaDeg(thetaKine, phiKine);
		isAcc = isAcepted(phiKine, thetaKine);
	}

	void setCluster(HEmcClusterSim * cls, const EcalMap & ecalMap, const Float_t & VertexZ)
	{
		Float_t time = cls->getTime();
		HGeomVector trackVec(cls->getXLab(), cls->getYLab(), cls->getZLab() - VertexZ);
		Double_t trackLength = trackVec.length();

		clusterBeta = (trackLength / 1000.) / (time * 1.e-9) / TMath::C();
		Float_t energy = cls->getEnergy();
		clusterTheta  = cls->getTheta();
		clusterEnergy = ecalMap.getCorrectedEnergy(cls);
		trackVec     *= (clusterEnergy / trackLength);
		tlvCluster.SetXYZM(trackVec.getX(), trackVec.getY(), trackVec.getZ(), 0);
		isRec = kTRUE;
		Int_t matchedCells = cls->getNMatchedCells();
		Int_t matchedPart  = ecalMap.getEcalCandMultiplicity(cls);
		isMatched = (matchedPart > 0 || matchedCells > 0);
		if (!isMatched && clusterBeta > 0.85 && clusterBeta < 1.3 && energy > 100 &&
		    clusterTheta > 16.5)
			isRecWithCuts = kTRUE;
	}

	inline Float_t getClusterEnergy(){
		return clusterEnergy;
	}

	inline Float_t getClusterTheta(){
		return clusterTheta;
	}

	inline Float_t getClusterBeta(){
		return clusterBeta;
	}

	inline HGeantKine * getKine(){
		return kine;
	}

	HGeantKine *   kine;
	TLorentzVector tlvCluster, tlvKine;
	Bool_t isAcc, isRec, isRecWithCuts, isMatched;
	Float_t clusterEnergy, clusterTheta, clusterBeta;
};


Int_t ana(TString infileList, TString outputFile, Int_t nev = 1000000.)
{
	HLoop loop(kTRUE); // create HADES (needed if one wants to use gHades)
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

	if (!loop.setInput("-*,+HGeantKine,+HParticleCand,+HEmcCal,+HEmcCluster,+HParticleEvtInfo")) {
		cout << "READBACK: ERROR : cannot read input !" << endl;
		exit(1);
	}
	loop.printCategories();
	loop.printChain();

	CutValues iniVal;
	EcalMap ecalMap;

	TFile * out = new TFile(outputFile, "RECREATE");
	out->cd();


	Int_t nYBins  = iniVal.getNYBins();
	Int_t nPtBins = iniVal.getNPtBins();

	Double_t yMin        = iniVal.getMinY();
	Double_t yMax        = iniVal.getMaxY();
	Double_t ptMin       = iniVal.getMinPt();
	Double_t ptMax       = iniVal.getMaxPt();
	Double_t midrapidity = iniVal.getMidrapidity();


	enum centClass { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCentClass };
	TH1D * centCounter = new TH1D("centCounter", "centCounter", cNumCentClass, 0, cNumCentClass);

	centCounter->GetXaxis()->SetBinLabel(5, "00 - 40 %");
	centCounter->GetXaxis()->SetBinLabel(1, "00 - 10 %");
	centCounter->GetXaxis()->SetBinLabel(2, "10 - 20 %");
	centCounter->GetXaxis()->SetBinLabel(3, "20 - 30 %");
	centCounter->GetXaxis()->SetBinLabel(4, "30 - 40 %");


	TH2F * phaseSpace[cNumCentClass];
	TH2F * phaseSpaceFull[cNumCentClass];

	TH2F * phaseSpaceAccepted[cNumCentClass];
	TH2F * phaseSpaceAcceptedFull[cNumCentClass];

	TH2F * phaseSpaceReconstructed[cNumCentClass];
	TH2F * phaseSpaceReconstructedFull[cNumCentClass];

	TH2F * phaseSpaceReconstructedWithCuts[cNumCentClass];
	TH2F * phaseSpaceReconstructedFullWithCuts[cNumCentClass];

	TH1F * singlePhotonEffMom[cNumCentClass];
	TH1F * singlePhotonEffTheta[cNumCentClass];

	TH1F * singlePhotonAccMom[cNumCentClass];
	TH1F * singlePhotonAccTheta[cNumCentClass];

	TH2F * singlePhotonAccMomCells[cNumCentClass];
	TH2F * singlePhotonEffMomCells[cNumCentClass];
	TH1F * ecalCells = new TH1F("ecalCells", "", 1200, 0, 1200);

	for (Int_t iCent = 0; iCent < cNumCentClass; iCent++) {
		TDirectory * sub;
		if (iCent != cNumCentClass - 1) sub = out->mkdir(Form("cent%i", iCent));
		else sub = out->mkdir("all_cents");
		sub->cd();
		phaseSpace[iCent] =
			new TH2F("phaseSpace", "; y_{cm} ;p_{t} [MeV/c]", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
		phaseSpaceFull[iCent] = new TH2F("phaseSpaceFull", "; y_{cm};p_{t} [MeV/c]", 400, -3., 3., 400, 0, 1500);

		phaseSpaceAccepted[iCent] = new TH2F("phaseSpaceAccepted", "; y_{cm} ;p_{t} [MeV/c]", nYBins, yMin, yMax,
		                                     nPtBins, ptMin, ptMax);
		phaseSpaceAcceptedFull[iCent] = new TH2F("phaseSpaceAcceptedFull", "; y_{cm} ;p_{t} [MeV/c]", 400, -3., 3., 400, 0,
		                                         1500);

		phaseSpaceReconstructed[iCent] = new TH2F("phaseSpaceReconstructed", "; y_{cm} ;p_{t} [MeV/c]", nYBins, yMin,
		                                          yMax, nPtBins, ptMin, ptMax);
		phaseSpaceReconstructedFull[iCent] = new TH2F("phaseSpaceReconstructedFull", "; y_{cm} ;p_{t} [MeV/c]", 400, -3.,
		                                              3., 400, 0, 1500);

		phaseSpaceReconstructedWithCuts[iCent] = new TH2F("phaseSpaceReconstructedWithCuts", "; y_{cm} ;p_{t} [MeV/c]",
		                                                  nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
		phaseSpaceReconstructedFullWithCuts[iCent] = new TH2F("phaseSpaceReconstructedFullWithCuts",
		                                                      "; y_{cm} ;p_{t} [MeV/c]", 400, -3., 3., 400, 0, 1500);

		singlePhotonAccMom[iCent] = new TH1F("singlePhotonAccMom", "; E [MeV];#varepsilon", 200, 0, 1000);
		singlePhotonEffMom[iCent] = new TH1F("singlePhotonEffMom", "; E [MeV];#varepsilon", 200, 0, 1000);

		singlePhotonAccTheta[iCent] = new TH1F("singlePhotonAccTheta", "; #theta [deg];#varepsilon", 90, 5, 50);
		singlePhotonEffTheta[iCent] = new TH1F("singlePhotonEffTheta", "; #theta [deg];#varepsilon", 90, 5, 50);

		singlePhotonAccMomCells[iCent] = new TH2F("singlePhotonAccMomCells", ";sector*200+cell; E [MeV];#varepsilon", 1200,
		                                          0, 1200, 200, 0, 1000);
		singlePhotonEffMomCells[iCent] = new TH2F("singlePhotonEffMomCells", ";sector*200+cell; E [MeV];#varepsilon", 1200,
		                                          0, 1200, 200, 0, 1000);


		out->cd();
	}


	// --------------------------CONFIGURATION---------------------------------------------------
	// At begin of the program (outside the event loop)
	HCategory * candCat       = (HCategory *) HCategoryManager::getCategory(catParticleCand);
	HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
	HCategory * geantKineCat  = (HCategory *) HCategoryManager::getCategory(catGeantKine);
	HCategory * emcCalCat     = (HCategory *) HCategoryManager::getCategory(catEmcCal);


	HParticleTrackSorter sorter;
	sorter.setIgnoreInnerMDC(); // do not reject Double_t inner MDC hits
	sorter.init();            // get catgegory pointers etc...
	// =======================================================================
	Int_t nEntries = loop.getEntries();
	printf("nEntries = %i\n", nEntries);
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
	if (nEntries > nev && nev > 0) nEntries = nev;

	TString tempFileName, currentFileName;
	Int_t dayOfYear, hour, min, second;


	HGeantHeader * fSubHeader = NULL;
	HPartialEvent * fSimul    = NULL;

	Bool_t isSimulation = kFALSE;
	if (geantKineCat) { isSimulation = kTRUE;  Info("init()", "GeantKine found - is Simulation"); } else {
		isSimulation = kFALSE;
		Info("init()", "GeantKine not found - is not Simulation");
	}

	TFile * GammaAcc = new TFile("/lustre/hebe/hades/user/prozorov/include/acceptancePhiThetaMar19.root", "read");
	acceptance_mar19 = (TH2D *) GammaAcc->Get("acceptancePhiThetaMar19");

	TFile * in           = new TFile("/lustre/hebe/hades/user/prozorov/include/AngleToCell.root", "read");
	TH2D * AngleToRow    = (TH2D *) in->Get("AngleToRow");
	TH2D * AngleToColumn = (TH2D *) in->Get("AngleToColumn");

	// --------------------------EVENT LOOP---------------------------------------------------//
	for (Int_t i = 0; i < nEntries; i++) {
		Int_t nbytes = loop.nextEvent(i); // get next event. categories will be cleared before
		if (nbytes <= 0) { cout << nbytes << endl; break; } // last event reached
		if (loop.isNewFile(tempFileName)) {
			TString type;
			Int_t year;
			Int_t eb;
			currentFileName = tempFileName;
			HTime::splitFileName(HTime::stripFileName(tempFileName), type, year, dayOfYear, hour, min, second, eb);
		}
		Float_t valueMinute = dayOfYear + hour / 24. + min / 24. / 60.;
		Float_t beamTime    = valueMinute + second / 24. / 60. / 60.;


		HTool::printProgress(i, nEntries, 1, "Analyze pairs :");

		// -------------------------------------------------
		// summary event info object

		Float_t VertexZ = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();
		if (VertexZ < -70) continue;
		Int_t cent_bin = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10) - 1;
		if (cent_bin >= 4 || cent_bin < 0) continue;

		centCounter->Fill(cent_bin);
		centCounter->Fill(all_cents);

		std::vector<Photon> daughters;
		daughters.clear();

		for (Int_t i = 0; i < geantKineCat->getEntries(); i++) {
			HGeantKine * kine = (HGeantKine *) geantKineCat->getObject(i);
			if (kine->isPrimary() == false) continue;
			if (kine->getGeneratorInfo() < 500) continue;
			if (kine->getID() != 1) continue;
			daughters.push_back({ kine });
		}
		// ========================================================================================//

		TLorentzVector pi0 = daughters.at(0).tlvKine + daughters.at(1).tlvKine;
		Float_t pt         = pi0.Pt();
		Float_t y = pi0.Rapidity() - midrapidity;
		phaseSpace[cent_bin]->Fill(y, pt);
		phaseSpaceFull[cent_bin]->Fill(y, pt);
		phaseSpace[all_cents]->Fill(y, pt);
		phaseSpaceFull[all_cents]->Fill(y, pt);

		// ===================Accepted=======================================================//
		if (daughters.at(0).isAcc && daughters.at(1).isAcc) {
			phaseSpaceAccepted[cent_bin]->Fill(y, pt);
			phaseSpaceAcceptedFull[cent_bin]->Fill(y, pt);
			phaseSpaceAccepted[all_cents]->Fill(y, pt);
			phaseSpaceAcceptedFull[all_cents]->Fill(y, pt);
		}
		// ===================Reconstructed=======================================================//

		sorter.cleanUp();
		sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
		sorter.fill(HParticleTrackSorter::selectLeptons); // fill only good leptons
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
		sorter.fill(HParticleTrackSorter::selectHadrons); // fill only good hadrons (already marked good leptons will be skipped)
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);

		ecalMap.InitializeEcalCandMatching(candCat, emcClusterCat);
		candCat->Clear();


		for (Int_t e = 0; e < emcClusterCat->getEntries(); e++) { // loop over all clusters in emc for finding photons
			HEmcClusterSim * cls = (HEmcClusterSim *) emcClusterCat->getObject(e);
			Int_t track = cls->getTrack(0);
			if (track <= 0) continue;
			HGeantKine * clusterKine0 = (HGeantKine *) geantKineCat->getObject(track - 1);
			if (clusterKine0 == NULL) continue;
			if (!ecalMap.isGoodCell(cls)) continue;

			if (clusterKine0->getTrack() == daughters.at(0).kine->getTrack() ||
			    clusterKine0->getParentTrack() == daughters.at(0).kine->getTrack()
			    )
				daughters.at(0).setCluster(cls, ecalMap, VertexZ);

			if (clusterKine0->getTrack() == daughters.at(1).kine->getTrack() ||
			    clusterKine0->getParentTrack() == daughters.at(1).kine->getTrack()
			    )
				daughters.at(1).setCluster(cls, ecalMap, VertexZ);
		} // loop over emc


		if (daughters.at(0).isRec && daughters.at(1).isRec) {
			TLorentzVector pi0Rec;
			Float_t ptRec, yRec;
			pi0Rec = daughters.at(0).tlvCluster + daughters.at(1).tlvCluster;
			ptRec  = pi0Rec.Pt();
			yRec   = pi0Rec.Rapidity() - midrapidity;

			phaseSpaceReconstructed[cent_bin]->Fill(y, pt);
			phaseSpaceReconstructedFull[cent_bin]->Fill(y, pt);
			phaseSpaceReconstructed[all_cents]->Fill(y, pt);
			phaseSpaceReconstructedFull[all_cents]->Fill(y, pt);


			if (daughters.at(0).isRecWithCuts && daughters.at(1).isRecWithCuts) {
				phaseSpaceReconstructedWithCuts[cent_bin]->Fill(yRec, ptRec);
				phaseSpaceReconstructedFullWithCuts[cent_bin]->Fill(yRec, ptRec);
				phaseSpaceReconstructedWithCuts[all_cents]->Fill(yRec, ptRec);
				phaseSpaceReconstructedFullWithCuts[all_cents]->Fill(yRec, ptRec);
			}
		}

		// ========================================================================================//


		// ===============Single photon efficiency====================================//


		sorter.cleanUp();
		sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
		sorter.fill(HParticleTrackSorter::selectLeptons); // fill only good leptons
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
		sorter.fill(HParticleTrackSorter::selectHadrons); // fill only good hadrons (already marked good leptons will be skipped)
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);

		ecalMap.InitializeEcalCandMatching(candCat, emcClusterCat);
		candCat->Clear();


		for (Int_t i = 0; i < geantKineCat->getEntries(); i++) {
			HGeantKine * kine = (HGeantKine *) geantKineCat->getObject(i);
			if (kine->isPrimary() == false) continue;
			if (kine->getID() != 1) continue;
			Float_t momKine = kine->getTotalMomentum();
			Float_t thetaKine, phiKine;
			kine->getPhiThetaDeg(thetaKine, phiKine);
			Photon photon(kine);
			if (!photon.isAcc) continue;
			Int_t sector = kine->getSector();
			singlePhotonAccMom[cent_bin]->Fill(momKine);
			singlePhotonAccTheta[cent_bin]->Fill(thetaKine);

			Int_t rowEMC    = round(AngleToRow->GetBinContent(AngleToRow->FindBin(phiKine, thetaKine)));
			Int_t columnEMC = round(AngleToColumn->GetBinContent(AngleToColumn->FindBin(phiKine, thetaKine)));
			Int_t cellKine  = sector * 200 + ecalMap.getSystemCell(rowEMC * 17 + columnEMC);

			singlePhotonAccMomCells[cent_bin]->Fill(cellKine, momKine);


			for (Int_t e = 0; e < emcClusterCat->getEntries(); e++) { // loop over all clusters in emc for finding photons
				HEmcClusterSim * cls = (HEmcClusterSim *) emcClusterCat->getObject(e);
				if (!ecalMap.isGoodCell(cls)) continue;
				ecalCells->Fill(ecalMap.getCell(cls));

				if (isKineMatched(cls, kine->getTrack(), geantKineCat)) {
					photon.setCluster(cls, ecalMap, VertexZ);
					// cout << "Kine CLuster " << cellKine << "  " << ecalMap.getCell(cls) << endl;
					// cout << "Energy Kine Cluster " << momKine << "  " << photon.getClusterEnergy() << endl;
					Float_t energyDiff       = (100 * abs(photon.getClusterEnergy() - momKine) / photon.getClusterEnergy());
					Float_t withinResolution = energyDiff * sqrt(momKine / 1000.) / 3;
					if (!photon.isMatched &&
					    photon.getClusterBeta() > 0.8 &&
					    photon.getClusterBeta() < 1.3 &&
					    withinResolution < 5.9)
					{
						singlePhotonEffMomCells[cent_bin]->Fill(cellKine, momKine);
						if (photon.getClusterTheta() > 16.5) singlePhotonEffMom[cent_bin]->Fill(momKine);
						if (photon.getClusterEnergy() > 100) singlePhotonEffTheta[cent_bin]->Fill(thetaKine);
					}
				}// found photon in cluster!
			} // loop over emc
		} // loop over geant kine
	} // event loop
	out->cd();
	out->GetList()->Write();

	out->Save();
	out->Close();
	return 0;
} /* fill_ttree_diphoton */
