// #ifndef __CINT__
#include "Includes.h"
#include "exp_setup_cells.h"

#include "hades.h"
#include "hloop.h"
#include "hcategory.h"
#include "hcategorymanager.h"

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


#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "hparticleevtchara.h"
#include "TMath.h"
#include "hmetamatcher.h"
#include "TVector3.h"
#include "Includes.h"
#include <iostream>
#include <vector>


#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "hemcraw.h"
#include "hemccal.h"
#include "hemccluster.h"
#include "TMath.h"

#include "TVector3.h"
#include "Includes.h"
#include <iostream>
#include <vector>
#include "hmetamatcher.h"



typedef struct {
	Float_t sector, cell,  tof, emc_mq, mom, multEmc,
	        beta, charge, energy, time, emc_path,distRpcEmc,timeCorr,
	        beamtime, TOT, qualityDThDPh, qualityDTime, raw_time, sys;
}
PID;

void analysis(TString inFile, TString outputFile, Int_t nev = 1000000.)
{
	Bool_t createHades = kTRUE; // kFALSE;  // kTRUE = create HADES new
	HLoop loop(createHades);  // kTRUE : create Hades  (needed to work with standard eventstructure)

	// -------------------------------------------------
	// setup the spectrometer for filter tasks
	// this part is needed if your task needs
	// detectors or parameter container
	if (1) {
		TString beamtime = "mar19";

		Int_t mdcMods[6][4] =
		{ { 1, 1, 1, 1 },
			{ 1, 1, 1, 1 },
			{ 1, 1, 1, 1 },
			{ 1, 1, 1, 1 },
			{ 1, 1, 1, 1 },
			{ 1, 1, 1, 1 } };
		TString asciiParFile = "";
		TString rootParFile  = "/cvmfs/hades.gsi.de/param/real/mar19/allParam_Mar19_gen5_13042021.root";
		TString paramSource  = "root";           // root, ascii, oracle
		TString paramrelease = "MAR19_dst_gen5"; // now, APR12_gen2_dst APR12_gen5_dst
		HDst::setupSpectrometer(beamtime, mdcMods, "rich,mdc,tof,rpc,emc,wall,start,tbox");
		HDst::setupParameterSources(paramSource, asciiParFile, rootParFile, paramrelease); // now, APR12_gen2_dst
	}
	// -------------------------------------------------



	Bool_t ret = kFALSE;

	if (inFile.Contains(",")) {
		ret = loop.addMultFiles(inFile); // file1,file2,file3
	} else {
		ret = loop.addFiles(inFile); // myroot*.root
	}
	if (ret == 0) {
		cout << "READBACK: ERROR : cannot find inputfiles : " << inFile.Data() << endl;
		return;
	}
	if (!loop.setInput("-*,+HParticleCand,+HEmcCluster,+HEmcCal,+HParticleEvtInfo,+HStart2Hit")) {
		cout << "READBACK: ERROR : cannot read input !" << endl;
		exit(1);
	}
	loop.printCategories();
	loop.printChain();


	HCategory * candCat       = (HCategory *) HCategoryManager::getCategory(catParticleCand);
	HCategory * emcCalCat     = (HCategory *) HCategoryManager::getCategory(catEmcCal);
	HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
	HCategory * startCat      = (HCategory *) HCategoryManager::getCategory(catStart2Hit);
	HIterator * iterStart     = (HIterator *) startCat->MakeIterator("native");

	// -------------------------------------------------
	// added for the task
	HTaskSet * masterTaskSet = gHades->getTaskSet("all");
	HMetaMatcher * matcher   = new HMetaMatcher();
	// matcher->setDebug();
	masterTaskSet->add(matcher);
	// -------------------------------------------------



	TFile * out = new TFile(outputFile, "RECREATE");
	out->cd();


	PID pid;
	TTree nt("nt", "pid");



	nt.Branch("pid",&pid,"sector:cell:tof:emc_mq:mom:multEmc:beta:charge:energy:time:emc_path:"
	          "distRpcEmc:timeCorr:beamtime:TOT:qualityDThDPh:qualityDTime:raw_time:sys");


	// --------------------------CONFIGURATION---------------------------------------------------
	// At begin of the program (outside the event loop)
	HParticleTrackSorter sorter;
	sorter.setIgnoreInnerMDC(); // do not reject Double_t inner MDC hits
	sorter.init();            // get catgegory pointers etc...
	// =======================================================================
	Int_t nEntries = loop.getEntries();
	printf("nEntries = %i\n", nEntries);

	TString tempFileName, currentFileName;
	Int_t dayOfYear, hour, min, second;

	for (Int_t i = 0; i < nEntries; i++) {
		Int_t nbytes = loop.nextEvent(i); // get next event. categories will be cleared before
		if (nbytes <= 0) { cout << nbytes << endl; break; } // last event reached
		HParticleEvtInfo * evtinfo  = HCategoryManager::getObject(evtinfo, catParticleEvtInfo, 0);
		HEventHeader * event_header = NULL;
		if (!(event_header = loop.getEventHeader())) continue;


		Float_t start_time = 0.;
		HStart2Hit * start = 0;
		iterStart->Reset();
		while ((start = (HStart2Hit *) iterStart->Next()) != 0) {
			start_time = start->getTime();
		}
		startCat->Clear();



		if (loop.isNewFile(tempFileName)) {
			TString type;
			Int_t year;
			Int_t eb;
			currentFileName = tempFileName;
			HTime::splitFileName(HTime::stripFileName(tempFileName), type, year, dayOfYear, hour, min, second, eb);
		}
		Float_t valueMinute = dayOfYear + hour / 24. + min / 24. / 60.;
		Float_t valueSecond = valueMinute + second / 24. / 60. / 60.;
		pid.beamtime = valueSecond;


		HTool::printProgress(i, nEntries, 1, "Analyze pairs :");
		Float_t VertexZ = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();
		sorter.cleanUp();
		sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
		sorter.fill(HParticleTrackSorter::selectHadrons); // fill only good hadrons (already marked good leptons will be skipped)
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);


		Int_t nEmcCluster = emcClusterCat->getEntries();
		pid.multEmc = nEmcCluster;

		vector<PID> photons;
		for (Int_t e = 0; e < nEmcCluster; e++) { // loop over all clusters in emc for finding photons
			HEmcCluster * cls = (HEmcCluster *) emcClusterCat->getObject(e);
			Float_t energy    = cls->getEnergy();
			Int_t sector      = cls->getSector();
			Int_t cell        = getSystemMap(cls->getCell());
			Float_t time      = cls->getTime();
			Int_t size = cls->getNCells();
			if (size != 1) continue;
			if (cls->getNMatchedCells() == 0) continue;

			pid.emc_mq   = 1000;
			for (Int_t j = 0; j < candCat->getEntries(); j++) { // loop over Particlecand to see if there is a match
				HParticleCand * cand = HCategoryManager::getObject(cand, candCat, j);
				if (!cand->isFlagBit(kIsUsed)) continue;


				Int_t emc_index = cand->getEmcInd();
				pid.emc_mq   = cand->getMetaMatchQualityEmc();
				if (emc_index == cls->getIndex() && pid.emc_mq<2.5) {


					Float_t timeCorr   = 0.;
					Float_t distRpcEmc = -100;

					if (cand->isTofClstUsed()) {
						HTofCluster * tofCluster = matcher->matchWithTof(cand);
						Float_t xrl, yrl, zrl;
						tofCluster->getXYZLab(xrl, yrl, zrl);
						Float_t deltaX = cls->getXLab() - xrl;
						Float_t deltaY = cls->getYLab() - yrl;
						Float_t deltaZ = cls->getZLab() - zrl;
						distRpcEmc = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
						timeCorr   = distRpcEmc / TMath::Sqrt(xrl * xrl + yrl * yrl + zrl * zrl) * tofCluster->getTof();
					}



					if (cand->isRpcClstUsed()) {
						HRpcCluster * fakerpccluster = matcher->matchWithRpc(cand);
						Float_t xrl, yrl, zrl;
						fakerpccluster->getXYZLab(xrl, yrl, zrl);
						Float_t deltaX = cls->getXLab() - xrl;
						Float_t deltaY = cls->getYLab() - yrl;
						Float_t deltaZ = cls->getZLab() - zrl;
						distRpcEmc = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
						timeCorr   = distRpcEmc / TMath::Sqrt(xrl * xrl + yrl * yrl + zrl * zrl) * fakerpccluster->getTof();
					}



					pid.distRpcEmc = distRpcEmc;
					pid.timeCorr   = timeCorr;


					pid.emc_path = cand->getEmcCellPath();
					pid.beta     = cand->getBeta();
					pid.charge   = cand->getCharge();
					pid.tof      = cand->getTof();
					pid.sys      = cand->getSystemUsed();
					pid.mom      = cand->getMomentum();


					UChar_t c1 = cls->getCell(0);
					HLocation loc;// location of cluster
					loc.set(2, 0, 0);
					loc[0] = sector; // set location of cluster
					loc[1] = c1; // set location of cluster
					HEmcCal * cal = (HEmcCal *) emcCalCat->getObject(loc);
					pid.TOT      = cal->getWidth();
					pid.raw_time = cal->getLead() - start_time;
					pid.sector = sector;
					pid.cell   = cell;
					pid.time   = time;
					pid.energy = energy;



					pid.qualityDThDPh = cls->getQualDThDPh();
					pid.qualityDTime  = cls->getQualDTime();
					nt.Fill();

				}
			}
			candCat->Clear();


		}
	} // event loop
	out->cd();
	out->GetList()->Write();
	out->Save();
	out->Close();
} // end of macro
