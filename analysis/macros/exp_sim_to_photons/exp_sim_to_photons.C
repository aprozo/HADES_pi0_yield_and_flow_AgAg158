
#include "hades.h"
#include "htool.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "hgeantkine.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "hemcclustersim.h"
#include "TLorentzVector.h"
#include "hcategorymanager.h"
#include "TMacro.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hparticleevtinfo.h"

#include "hstart2cal.h"
#include "hstart2hit.h"

#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "hstartdef.h"
#include "richdef.h"
#include "rpcdef.h"
#include "showerdef.h"
#include "hrpccluster.h"
#include "tofdef.h"
#include "walldef.h"
#include "hemccluster.h"
#include "hemccal.h"
#include "hemccalibrater.h"
#include "emcdef.h"

#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLine.h"
#include "hparticlebooker.h"
#include "htime.h"
#include "hparticleanglecor.h"

#include "hparticleevtchara.h"

#include "exp_setup_cells.h"
#include "eventclassifier.h"
#include "sim_correction.h"
#include "EcalCellPar.h"
#include "../GlobalParameters.h"

EcalCellPar ecalCellPars[6][163];

TH1D* h_energy_correction_exp;
Int_t setupEnergyCorrectionExp(TString pathEnegryCorrectionExp="")
{
	TFile *fileEnergyCorrectionExp= new TFile(pathEnegryCorrectionExp);
	if (fileEnergyCorrectionExp)
	{
		h_energy_correction_exp = (TH1D*) fileEnergyCorrectionExp->Get("AllCorrFactors1D");;
		if (!h_energy_correction_exp)
		{
			std::cout << "No histogram with exp energy correction found in file. Exitting..." << std::endl;
			return 1;
		}
	}
	else
	{
		std::cout << "No file with with exp energy correction found. Exitting..." << std::endl;
		return 2;
	}
	return 0;
}


void correctNonlinearity(Float_t & energy){
	energy*= 0.99367-1.42255/sqrt(energy) + 23.8587/energy -110/pow(energy,1.5)+1571.92/pow(energy,2);
}
void correctBeamDependence(Float_t & energy, const Float_t & beamTime){
	energy*= (155.932-0.261425*beamTime)/134.977;
}

void recalibrate(HEmcCluster * cls, HCategory * emcCalCat,Float_t beamTime=80.1568327436)
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
		HEmcCal * cal = (HEmcCal *) emcCalCat->getObject(loc);
		Float_t TOT   = cal->getWidth();
		Params newParams = ecalCellPars[sector][cell].newParams;
		Float_t newEnergy = newParams.adc0 + exp(newParams.adc1 + newParams.adc2 * TOT); // recalibration of a single module energy within cluster
		energySum += newEnergy;
	}
	energySum*=h_energy_correction_exp->GetBinContent(sector*163+getSystemMap(cls->getCell())+1);

	//correctNonlinearity(energySum);
	//correctBeamDependence(energySum,beamTime);
	cls->setEnergy(energySum);
}

TFile pi0PeakCorrFactorsFile("pi0PeakCorrFactors.root","read");
TH1D * pi0PeakCorrFactors=(TH1D *)pi0PeakCorrFactorsFile->Get("pi0PeakCorrFactors");
...
Double_t energy = emcCluster->GetEnergy();
energy *=pi0PeakCorrFactors->GetBinContent(sector*163+cell+1); //cell is from 0 to 162 
emcCluster->setEnergy(energy);
}



typedef struct {
	Float_t newEvent, multTofRpc, multEmc, multPhoton, centBin,
	        targetBin, beamTime, RPlanePhi, impPar, emc_mq, qualityDTime, betaCand,
	        sector, cell, size, beta, energy, time, xpos, ypos, zpos, track, trackParent, mom_kine, idParent;
} PID;
struct Candidate {
	Candidate(const Int_t & _emc_index, const Float_t & _emc_mq, const Float_t & _beta)
		: emc_index(_emc_index), emc_mq(_emc_mq), beta(_beta){
	};
	Int_t emc_index;
	Float_t emc_mq, beta;
};
const Float_t impact_bins[] = { 0., 3.8, 5.4, 6.6, 7.6 }; // latest version from WIKI

Int_t getCentBin(const Double_t & b)
{
	Int_t impBin = 5; // Hades binning

	if (b >= impact_bins[0] && b < impact_bins[1]) impBin = 1;
	if (b >= impact_bins[1] && b < impact_bins[2]) impBin = 2;
	if (b >= impact_bins[2] && b < impact_bins[3]) impBin = 3;
	if (b >= impact_bins[3] && b < impact_bins[4]) impBin = 4;
	return impBin;
}


// exclude days
Bool_t isGoodBeamTime(const Float_t & beamtime){
	if(beamtime >=68.1 && beamtime <68.2) return kFALSE;
	if(beamtime >=69.1 && beamtime <69.6) return kFALSE;
	if(beamtime >=70.3 && beamtime <70.4) return kFALSE;
	if(beamtime >=70.6 && beamtime <72.6) return kFALSE;
	if(beamtime >=73.5 && beamtime <73.6) return kFALSE;
	if(beamtime >=80.1 && beamtime <80.5) return kFALSE;
	if(beamtime >=82.5 && beamtime <85.6) return kFALSE;
	if(beamtime >=85.8 && beamtime <86.4) return kFALSE;
	return kTRUE;
}


Int_t exp_sim_to_photons(TString infileList, TString outputFile, Bool_t isExperiment=kTRUE, Int_t nEventsToAnalyze=-1)
{

	GlobalParameters par("../parameters.conf");
	HLoop  loop(kTRUE);// create HADES (needed if one wants to use gHades)
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

	if (!loop.setInput("-*,+HGeantKine,+HEmcCal,+HEmcCluster,+HParticleEvtInfo, +HParticleCand, +HWallHit")) {
		cout << "READBACK: ERROR : cannot read input !" << endl;
		exit(1);
	}
	loop.printCategories();
	loop.printChain();

	HCategory * candCat       = (HCategory *) HCategoryManager::getCategory(catParticleCand);
	HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
	HCategory * geantKineCat  = (HCategory *) HCategoryManager::getCategory(catGeantKine);
	HCategory * emcCalCat     = (HCategory *) HCategoryManager::getCategory(catEmcCal);
	HGeantHeader * fSubHeader = NULL;
	HPartialEvent * fSimul    = NULL;

	Bool_t isSimulation = !isExperiment;
	if (geantKineCat) { isSimulation = kTRUE;  Info("init()", "GeantKine found - is Simulation"); } else {
		isSimulation = kFALSE;
		Info("init()", "GeantKine not found - is not Simulation");
	}

	setupSimEffCorr(Form("%s/simCorrEnergyEcal.root",par.PATH_TO_CORRECTION_DATA));
	if(!isSimulation) {
		if(readParams(Form("%s/mar19EmcCalParNewEnergy.txt",par.PATH_TO_CORRECTION_DATA),ecalCellPars)>0) return 3;
		if(setupEnergyCorrectionExp(Form("%s/CorrFactorsPi0Peak.root",par.PATH_TO_CORRECTION_DATA))>0) return 4;
	}

	TFile * out = new TFile(outputFile, "RECREATE");
	out->cd();


	PID pid;
	TTree nt("nt", "pid");

	nt.Branch("pid", &pid, "newEvent:multTofRpc:multEmc:multPhoton:centBin:targetBin:beamTime:RPlanePhi:impPar:emc_mq:qualityDTime:betaCand:"
	          "sector:cell:size:beta:energy:time:xpos:ypos:zpos:track:trackParent:mom_kine:idParent");

	enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
	TH1D * Counter_cent = new TH1D("Counter_cent", "Counter_cent", cNumCounters_cent, 0, cNumCounters_cent);

	Counter_cent->GetXaxis()->SetBinLabel(5, "all_cent");
	Counter_cent->GetXaxis()->SetBinLabel(1, "cent0_10");
	Counter_cent->GetXaxis()->SetBinLabel(2, "cent10_20");
	Counter_cent->GetXaxis()->SetBinLabel(3, "cent20_30");
	Counter_cent->GetXaxis()->SetBinLabel(4, "cent30_40");

	TH1D * Counter_cent_cut = new TH1D("Counter_cent_cut", "Counter_cent_cut", cNumCounters_cent, 0, cNumCounters_cent);
	Counter_cent_cut->GetXaxis()->SetBinLabel(5, "all_cent");
	Counter_cent_cut->GetXaxis()->SetBinLabel(1, "cent0_10");
	Counter_cent_cut->GetXaxis()->SetBinLabel(2, "cent10_20");
	Counter_cent_cut->GetXaxis()->SetBinLabel(3, "cent20_30");
	Counter_cent_cut->GetXaxis()->SetBinLabel(4, "cent30_40");


	TH2D * Counter_cent_beamtime = new TH2D("Counter_cent_beamtime", "Counter_cent_beamtime", 200, 68,88, cNumCounters_cent, 0, cNumCounters_cent);
	Counter_cent_beamtime->GetYaxis()->SetBinLabel(5, "all_cent");
	Counter_cent_beamtime->GetYaxis()->SetBinLabel(1, "cent0_10");
	Counter_cent_beamtime->GetYaxis()->SetBinLabel(2, "cent10_20");
	Counter_cent_beamtime->GetYaxis()->SetBinLabel(3, "cent20_30");
	Counter_cent_beamtime->GetYaxis()->SetBinLabel(4, "cent30_40");



	TH1D * hDeltaTimeEcalRpc = new TH1D("hDeltaTimeEcalRpc", "", 200, -4, 4);
	TH2D * hDeltaTimeEcalRpcCells = new TH2D("hDeltaTimeEcalRpcCells", "",1200,0,1200, 200, -4, 4);

	TH1D * hEcalMatchQuality = new TH1D("hEcalMatchQuality", "", 200, 0, 10);
	TH2D * hEcalMatchQualityCells = new TH2D("hEcalMatchQualityCells", "", 1200,0,1200,200, 0, 10);

	TH1D * hBetaEcal = new TH1D("hBetaEcal", "", 200, 0.8, 1.5);
	TH2D * hBetaEcalCells = new TH2D("hBetaEcalCells", "",1200,0,1200, 200, 0.8, 1.5);

	TH1D * hEnergyEcal = new TH1D("hEnergyEcal", "", 300, 0, 1500);
	TH2D * hEnergyEcalCells = new TH2D("hEnergyEcalCells", "",1200,0,1200, 300, 0, 1500);

	TH1D * hSizeEcal = new TH1D("hSizeEcal", "", 10, 0, 10);
	TH2D * hSizeEcalCells = new TH2D("hSizeEcalCells", "",1200,0,1200, 10, 0, 10);

	TH2D * hSizeEnergyEcal = new TH2D("hSizeEnergyEcal", "",300, 0, 1500, 10, 0, 10);


	TH1F *hEventPlaneCharaAll = new TH1F("hEventPlaneCharaAll", ";#Psi_{RP}[#circ];Counts", 200, -200, 200);
	TH1F *hEventPlaneInfo = new TH1F("hEventPlaneInfo", ";#Psi_{RP}[#circ];Counts", 200, -200, 200 );

	const Int_t nCentralityClasses=4;
	TH1F *hDeltaABcent10[nCentralityClasses];
	TH1F *hDeltaABcent5[nCentralityClasses*2];
	TH1F *hEventPlaneChara[nCentralityClasses];

	for (Int_t k = 0; k < nCentralityClasses; k++)
	{
		hDeltaABcent10[k] = new TH1F(Form("hDeltaABcent10_%d", k), ";#Delta#Psi=#Psi_{A}-#Psi_{B},[#circ];Counts", 180, 0, 180);
		hEventPlaneChara[k] = new TH1F(Form("hEventPlaneChara%d", k), ";#Psi_{RP}[#circ];Counts", 200, -200, 200);

	}
	for (Int_t k = 0; k < nCentralityClasses*2; k++)
	{   hDeltaABcent5[k] = new TH1F(Form("hDeltaABcent5_%d", k), ";#Delta#Psi=#Psi_{A}-#Psi_{B},[#circ];Counts", 180, 0, 180);}


	// --------------------------CONFIGURATION---------------------------------------------------
	// At begin of the program (outside the event loop)
	HParticleTrackSorter sorter;
	sorter.setIgnoreInnerMDC(); // do not reject Double_t inner MDC hits
	sorter.init();            // get catgegory pointers etc...
	// =======================================================================

	Int_t nEntriesInFile = loop.getEntries();
	Int_t nEvents = (nEventsToAnalyze>nEntriesInFile ? nEntriesInFile : nEventsToAnalyze);
	if(nEventsToAnalyze<0) nEvents=nEntriesInFile;

	printf("nEventsToAnalyze = %i\n",nEventsToAnalyze);
	printf("nEntriesInFile = %i\n",nEntriesInFile);
	printf("nEvents = %i\n",nEvents);
	printf(isExperiment ? "Experiment\n" : "Simulation\n");

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

	TString tempFileName, currentFileName;
	Int_t dayOfYear, hour, min, second;



	for (Int_t i = 0; i < nEvents; i++) {
		Int_t nbytes = loop.nextEvent(i); // get next event. categories will be cleared before
		pid.newEvent = 1;
		if (nbytes <= 0) { cout << nbytes << endl; break; } // last event reached
		if (loop.isNewFile(tempFileName)) {
			TString type;
			Int_t year;
			Int_t eb;
			currentFileName = tempFileName;
			HTime::splitFileName(HTime::stripFileName(tempFileName), type, year, dayOfYear, hour, min, second, eb);
		}
		Float_t valueMinute = dayOfYear + hour / 24. + min / 24. / 60.;
		Float_t beamTime = valueMinute + second / 24. / 60. / 60.;
		pid.beamTime = beamTime;
		if (isSimulation) pid.beamTime = 0;

		HTool::printProgress(i, nEvents, 1, "Analyze pairs :");
		Float_t VertexZ = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();
		if (VertexZ < -70) continue;

		HParticleEvtInfo * evtInfo = HCategoryManager::getObject(evtInfo, catParticleEvtInfo, 0);
		Int_t mult_tofrpc = evtInfo->getSumTofMultCut() + evtInfo->getSumRpcMultHitCut();
		Int_t cent_bin    = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);
		Int_t cent_bin5    = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5);
		Int_t target_bin  = getTargetBin();
		if (cent_bin == 0) continue; // exclude underflow bin


		// Extracting event plane
		Float_t eventPlaneInfo = evtInfo->getRPlanePhi();
		hEventPlaneInfo->Fill(eventPlaneInfo);

		pid.impPar    = 0;

		if (isSimulation) {
			fSimul     = ( (HRecEvent *) gHades->getCurrentEvent())->getPartialEvent(catSimul);
			fSubHeader = (HGeantHeader *) (fSimul->getSubHeader());
			if (fSubHeader) {
				pid.RPlanePhi = fSubHeader->getEventPlane();
				pid.impPar    = fSubHeader->getImpactParameter();
				//cent_bin      = getCentBin(pid.impPar);
			}
		}

		Float_t eventPlaneChara = evtChara.getEventPlane(HParticleEvtChara::kDefault)*TMath::RadToDeg();
		Float_t EventPlaneA =evtChara.getEventPlane(HParticleEvtChara::kDefault,1)*TMath::RadToDeg();
		Float_t EventPlaneB = evtChara.getEventPlane(HParticleEvtChara::kDefault,2)*TMath::RadToDeg();


		if(eventPlaneChara > 180.) eventPlaneChara -= 360.;
		pid.RPlanePhi=eventPlaneChara;

		if (evtChara.getEventPlane(HParticleEvtChara::kDefault) == -1) {pid.RPlanePhi =-1000;eventPlaneChara=-1000;}
		if (EventPlaneA*TMath::DegToRad() == -1 || EventPlaneB*TMath::DegToRad() == -1) pid.RPlanePhi =-1000;
		if (cent_bin > 4) continue;
		cent_bin = cent_bin - 1; // starting from 0, not 1

		Float_t deltaPhiAB = (EventPlaneA - EventPlaneB)*TMath::RadToDeg();
		if( EventPlaneA >= 0. && EventPlaneB >= 0. ) {
			deltaPhiAB = TMath::Abs( EventPlaneA - EventPlaneB );
			if(deltaPhiAB > 180.) deltaPhiAB = 360. - deltaPhiAB; //this makes sure that phiAB is in (0,180)
			hDeltaABcent10[cent_bin]->Fill(deltaPhiAB);
			if(cent_bin5<nCentralityClasses*2) hDeltaABcent5[cent_bin5]->Fill(deltaPhiAB);
		}

		hEventPlaneCharaAll->Fill(eventPlaneChara);
		hEventPlaneChara[cent_bin]->Fill(eventPlaneChara);

		if(isGoodBeamTime(pid.beamTime))
		{
			Counter_cent_cut->Fill(cent_bin);
			Counter_cent_cut->Fill(4);
		}

		Counter_cent->Fill(cent_bin);
		Counter_cent->Fill(4);

		Counter_cent_beamtime->Fill(pid.beamTime,cent_bin);
		Counter_cent_beamtime->Fill(pid.beamTime,4);

		Int_t emc_ModMult[6][163] = { 0 };
		vector < Candidate > candidates;

		sorter.cleanUp();
		sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
		sorter.fill(HParticleTrackSorter::selectLeptons); // fill only good leptons
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
		sorter.fill(HParticleTrackSorter::selectHadrons); // fill only good hadrons (already marked good leptons will be skipped)
		sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);

		for (Int_t j = 0; j < candCat->getEntries(); j++) { // loop over Particlecand to see if there is a match
			HParticleCand * cand = HCategoryManager::getObject(cand, candCat, j);
			if (!cand->isFlagBit(kIsUsed)) continue;
			Int_t emc_index = cand->getEmcInd();
			Float_t beta    = cand->getBeta();
			if (emc_index >= 0) {
				HEmcCluster * cls    = (HEmcCluster *) emcClusterCat->getObject(emc_index);
				Int_t sector         = cls->getSector();
				Int_t cell           = getSystemMap(cls->getCell());
				Float_t qualityDTime = cls->getQualDTime();
				Float_t emc_mq       = cand->getMetaMatchQualityEmc();
				Int_t cell1200 = sector*200+cell;

				hDeltaTimeEcalRpc->Fill(qualityDTime);
				hDeltaTimeEcalRpcCells->Fill(cell1200,qualityDTime);

				hEcalMatchQuality->Fill(emc_mq);
				hEcalMatchQualityCells->Fill(cell1200,emc_mq);


				if ( emc_mq < 2.5) emc_ModMult[sector][cell]++;
				else candidates.push_back({ emc_index, emc_mq, beta });
			}
		}
		candCat->Clear();

		Int_t nEmcCluster = emcClusterCat->getEntries();
		pid.multTofRpc = mult_tofrpc;
		pid.multEmc    = nEmcCluster;
		pid.centBin    = cent_bin;
		pid.targetBin  = target_bin;


		vector < PID > photons;
		for (Int_t e = 0; e < nEmcCluster; e++) { // loop over all clusters in emc for finding photons
			HEmcCluster * cls       = (HEmcCluster *) emcClusterCat->getObject(e);
			HEmcClusterSim * clsSim = dynamic_cast < HEmcClusterSim * > (cls);

			Float_t matchedCells = cls->getNMatchedCells();
			Float_t qualityDTime = cls->getQualDTime();
			Int_t sector         = cls->getSector();
			Int_t cell = getSystemMap(cls->getCell());
			Int_t cell1200 = sector*200+cell;
			if (!isGoodCell(sector, cell)) continue;

			Int_t matchedPart = emc_ModMult[sector][cell];
			Float_t time      = cls->getTime();
			HGeomVector trackVec(cls->getXLab(), cls->getYLab(), cls->getZLab() - VertexZ);
			Double_t trackLength = trackVec.length();
			Double_t beta        = (trackLength / 1000.) / (time * 1.e-9) / TMath::C();
			Int_t size = cls->getNCells();


			if (matchedPart > 0 || beta < 0.8 || beta > 1.5) continue;
			if (!isSimulation) recalibrate(cls, emcCalCat,beamTime);
			Float_t energy = cls->getEnergy();
			if (isSimulation) energy = correctEMCenergy(energy, sector, cell);

			pid.emc_mq       = 1000;
			pid.betaCand     = 0;
			pid.qualityDTime = qualityDTime;

			for (Int_t i = 0; i < candidates.size(); i++) {
				if (candidates[i].emc_index == e && pid.emc_mq > candidates[i].emc_mq) {
					pid.emc_mq   = candidates[i].emc_mq;
					pid.betaCand = candidates[i].beta;
				}
			}

			if(pid.emc_mq>20&&abs(pid.qualityDTime)>10) {
				hBetaEcal->Fill(beta);
				hBetaEcalCells->Fill(cell1200,beta);

				hEnergyEcal->Fill(energy);
				hEnergyEcalCells->Fill(cell1200,energy);

				hSizeEcal->Fill(size);
				hSizeEcalCells->Fill(cell1200,size);

				hSizeEnergyEcal->Fill(energy,size);
			}
			if (energy < 50) continue;

			pid.sector      = sector;
			pid.cell        = cell;
			pid.size        = size;
			pid.time        = time;
			pid.energy      = energy;
			pid.beta        = beta;
			pid.xpos        = cls->getXLab();
			pid.ypos        = cls->getYLab();
			pid.zpos        = cls->getZLab() - VertexZ;
			pid.track       = 0;
			pid.trackParent = 0;
			pid.mom_kine    = 0;
			pid.idParent    = 0;

			if (isSimulation) {
				Int_t track        = clsSim->getTrack();
				HGeantKine * pKine = (HGeantKine *) geantKineCat->getObject(track - 1);
				Float_t mom_kine   = pKine->getTotalMomentum();
				Int_t id = pKine->getID();
				HGeantKine * prnt = pKine->getParent(track);
				Int_t idp         = prnt != NULL ? prnt->getID() : 0;
				Int_t trackParent = prnt != NULL ? prnt->getTrack() : 0;

				pid.track       = id == 1 ? track : (-100);
				pid.idParent    = idp;
				pid.trackParent = (idp == 7 || idp == 17) && prnt->isPrimary() ? trackParent : (-100);
				pid.mom_kine    = mom_kine;

				HGeantKine * grandprnt = prnt->getParent(trackParent);
				Int_t idgrandp         = grandprnt != NULL ? grandprnt->getID() : 0;

				if ((id == 2 || id == 3) && (idp == 1) && (idgrandp == 7 || idgrandp == 17)) {
					pid.track       = prnt->getTrack();
					pid.trackParent = grandprnt->getTrack();
					pid.idParent    = idgrandp;
				}
			}
			photons.push_back(pid);
		} // loop over emc

		Bool_t isFirstEvent = kTRUE;
		for (PID p:photons) {
			pid = p;
			if (isFirstEvent) { pid.newEvent = 1; isFirstEvent = kFALSE; } else { pid.newEvent = 0; }
			pid.multPhoton = photons.size();
			nt.Fill();
		}

		if(photons.size()==0) {
			pid.newEvent    = 1;
			pid.multPhoton  = 0;
			pid.emc_mq      = 0;
			pid.betaCand    = 0;
			pid.qualityDTime= 0;
			pid.sector      = 0;
			pid.cell        = 0;
			pid.size        = 0;
			pid.time        = 0;
			pid.energy      = 0;
			pid.beta        = 0;
			pid.xpos        = 0;
			pid.ypos        = 0;
			pid.zpos        = 0;
			pid.track       = 0;
			pid.trackParent = 0;
			pid.mom_kine    = 0;
			pid.idParent    = 0;
			nt.Fill();
		}



	} // event loop
	out->cd();
	out->GetList()->Write();
	TMacro * m = new TMacro("exp_sim_to_photons.C");
	TTree macros("macros", "macros");
	macros.GetUserInfo()->Add(m);
	macros.Write();





	out->Save();
	out->Close();
	return 0;
} /* fill_ttree_diphoton */
