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
TH2D *acceptance_mar19;

Bool_t isInAcc_mar19(const Float_t & theta,  const Float_t & phi) {
	UInt_t bin = 0;
	bin = acceptance_mar19->FindBin(phi,theta);
	Float_t w = acceptance_mar19->GetBinContent(bin);
	if(w>0.9) return kTRUE;
	return kFALSE;
}

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


Int_t loopDST_task( TString infileList,TString outfile,Int_t nEvents)
{
	//  infileList : comma seprated file list "file1.root,file2.root" or "something*.root"
	//  outfile    : optional (not used here) , used to store hists in root file
	//  nEvents    : number of events to processed. if  nEvents < entries or < 0 the chain will be processed

	Bool_t isSimulation = kFALSE;

	//-------------------------------------------------
	// create loop obejct and hades
	HLoop loop(kTRUE);
	//-------------------------------------------------

	//-------------------------------------------------
	// setup the spectrometer for filter tasks
	// this part is needed if your task needs
	// detectors or parameter container
	if(1) {
		TString beamtime="mar19";

		Int_t mdcMods[6][4]=
		{ {1,1,1,1},
		  {1,1,1,1},
		  {1,1,1,1},
		  {1,1,1,1},
		  {1,1,1,1},
		  {1,1,1,1} };
		TString asciiParFile     = "";
		TString rootParFile      = "/cvmfs/hades.gsi.de/param/sim/mar19/allParam_mar19_sim_run_16000_gen4_25112020.root";
		TString paramSource      = "root";// root, ascii, oracle
		TString paramrelease     = "MAR19_dst_gen4";// now, APR12_gen2_dst APR12_gen5_dst
		HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,emc,wall,start,tbox");
		HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease); // now, APR12_gen2_dst
	}
	//-------------------------------------------------

	//-------------------------------------------------
	// list of all files with working sectors
	//if(!isSimulation) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);
	//-------------------------------------------------


	//-------------------------------------------------
	// reading input files and decalring containers
	Bool_t ret =kFALSE;
	if(infileList.Contains(",")) {
		ret = loop.addMultFiles(infileList);  // file1,file2,file3
	} else{
		ret = loop.addFiles(infileList); // myroot*.root
	}

	if(ret == 0) {
		cout<<"READBACK: ERROR : cannot find inputfiles : "<<infileList.Data()<<endl;
		return 1;
	}

	TFile* out = new TFile(outfile,"RECREATE");
	out->cd();
	typedef struct {
		Float_t imppar,cent_bin,charge,emc_match,sector,
		        emc_cell,theta_rpc,phi_rpc,mom,beta,mom_kine,
		        emc_size,emc_energy,emc_time,emc_energy_corr,emc_mq,qualityDTime,
		        emc_path,mm_qa,geantid,within_resolution,emc_energy_centr,
		        true_emc_hit_track,true_rpc_hit_track, cflag;
	} PID;
	static PID pid;

	TTree *nt = new TTree("nt","pid");
	nt->Branch("pid",&pid, "imppar:cent_bin:charge:emc_match:sector:"
	           "emc_cell:theta_rpc:phi_rpc:mom:beta:mom_kine:"
	           "emc_size:emc_energy:emc_time:emc_energy_corr:emc_mq:qualityDTime:"
	           "emc_path:mm_qa:geantid:within_resolution:emc_energy_centr:"
	           "true_emc_hit_track:true_rpc_hit_track:cflag");

	if(!loop.setInput("-*,+HParticleCandSim,+HEmcClusterSim,+HGeantKine")) {
		cout<<"READBACK: ERROR : cannot read input !"<<endl;
		exit(1);
	} // read all categories
	loop.printCategories();
	loop.printChain();

	//-------------------------------------------------

	//-------------------------------------------------
	HParticleTrackSorter sorter;
	sorter.init();                                                  // get catgegory pointers etc...
	//-------------------------------------------------

	//-------------------------------------------------
	// input data
	HEventHeader *evHeader  = loop.getEventHeader();
	HCategory* candCat      = (HCategory*)HCategoryManager::getCategory(catParticleCand);
	HCategory* emcClusterCat= (HCategory*)HCategoryManager::getCategory(catEmcCluster);
	HCategory* geantKineCat = (HCategory*)HCategoryManager::getCategory(catGeantKine);
	//-------------------------------------------------


	//-------------------------------------------------
	//parameters
	//HEnergyLossCorrPar dEdxCorr;
	//dEdxCorr.setDefaultPar("mar19");
	TFile* acc       = new TFile("acceptancePhiThetaMar19.root","Read");
	//acceptance_mar19 =(TH2F*)acc->Get("Acceptance2D");
	acceptance_mar19 =(TH2D*)acc->Get("acceptancePhiThetaMar19");
	setupMassCuts("m2_cut.root");
	setupSimEffCorr("corretion_energy_ag1650ag.root");
	HPartialEvent* fSimul = NULL;
	HGeantHeader* fSubHeader = NULL;
	//-------------------------------------------------

	//-------------------------------------------------
	// added for the task
	HTaskSet *masterTaskSet = gHades->getTaskSet("all");
	HMetaMatcher* matcher = new HMetaMatcher();
	//matcher->setDebug();
	masterTaskSet->add(matcher);
	//-------------------------------------------------

	//-------------------------------------------------
	// event loop starts here
	Int_t entries = loop.getEntries();
	if(nEvents < entries && nEvents >= 0 ) entries = nEvents;




enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
  TH1D * Counter_cent = new TH1D("Counter_cent", "Counter_cent", cNumCounters_cent, 0, cNumCounters_cent);

  Counter_cent->GetXaxis()->SetBinLabel(5, "all_cent");
  Counter_cent->GetXaxis()->SetBinLabel(1, "cent0_10");
  Counter_cent->GetXaxis()->SetBinLabel(2, "cent10_20");
  Counter_cent->GetXaxis()->SetBinLabel(3, "cent20_30");
  Counter_cent->GetXaxis()->SetBinLabel(4, "cent30_40");

	HLocation loc;
	loc.set(3,0,0,0);

	for (Int_t i = 0; i < entries; i++) {
		Int_t nbytes =  loop.nextEvent(i);         // get next event. categories will be cleared before
		if(nbytes <= 0) { cout<<nbytes<<endl; break; } // last event reached
		//if(i%1000 == 0) cout<<"event "<<i<<endl;

		/*
		   //-------------------------------------------------
		   // summary event info object
		   HParticleEvtInfo* evtInfo=0;
		   evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

		   if(evtInfo&&!evtInfo->isGoodEvent(Particle::kGoodTRIGGER|          // standard event selection apr12
		                  Particle::kGoodVertexClust|
		                  Particle::kGoodVertexCand|
		                  Particle::kGoodSTART|
		                  Particle::kNoPileUpSTART|
		                  Particle::kNoVETO|
		                  Particle::kGoodSTARTVETO|
		                  Particle::kGoodSTARTMETA
		                 )) continue;
		   //-------------------------------------------------
		 */
		HGeomVector vertex(evHeader->getVertexCluster().getPos()); // Event vertex position
		fSimul     = ( (HRecEvent*)gHades->getCurrentEvent())->getPartialEvent(catSimul);
		fSubHeader = (HGeantHeader*)(fSimul->getSubHeader());
		Float_t impactParameter = -1.;
		if (fSubHeader) impactParameter = fSubHeader->getImpactParameter();

		UInt_t cent_bin      = getCentBin(impactParameter)-1;
		Counter_cent->Fill(cent_bin);
		Counter_cent->Fill(4);

		pid.imppar=impactParameter;
		pid.cent_bin=cent_bin;

		//-------------------------------------------------
		// loop over particle candidates in event
		if(candCat) {
			//------------------------------------------------------------------------
			// clean vectors and index arrays
			sorter.cleanUp();
			// if you do not want to modify the flags
			// for the following tasks etc (needs call of
			// restoreFlags() at the end of the eventloop)
			//sorter.backupFlags();
			//------------------------------------------------------------------------

			sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
			sorter.fill(selectLeptonsBetaLoose); // fill only good loose leptons
			Int_t nCandLepBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsLeptonSorter);
			Int_t size = candCat->getEntries();
			HParticleCandSim* cand=0;

			for(Int_t j = 0; j < size; j++) {

				//------------------------------------------------------------------------
				cand = HCategoryManager::getObject(cand,candCat,j);
				if(cand) {
					if(!leptonParticlecandCut(cand)) continue;    // see  #include "lepton_efficiency.h"
					if(!cand->isFlagBit(kIsUsed)) { continue; }   // skip particles rejected by HParticleTrackSorter
					if(cand->isRpcClstUsed()) {

						HRpcCluster* fakerpccluster = matcher->matchWithRpc(cand);
						Float_t mom                 = cand->getMomentum();
						Float_t emc_energy          = cand->getEmcEnergy();
						Int_t emc_index           = cand->getEmcInd();
						Float_t beta               = cand->getBeta();
						pid.charge                  = cand->getCharge();
						pid.mm_qa                   = cand->getMetaMatchQuality();
						pid.mom_kine                = cand->getGeantTotalMom();
						pid.cflag                   = cand->getGeantCorrTrackIds();
						pid.geantid                 = cand->getGeantPID();
						pid.theta_rpc               = fakerpccluster->getTheta();
						pid.phi_rpc                 = fakerpccluster->getPhi();

						UInt_t cflag                = cand->getGeantCorrTrackIds();
						pid.true_rpc_hit_track      = (cflag&kIsInRPC)==kIsInRPC;
						pid.true_emc_hit_track      = (cflag&kIsInRPC)==kIsInRPC &&
						                              (cflag&kIsInEMC)==kIsInEMC;
						Int_t sec                 = cand->getSector();
						/*
						   for(Int_t sector = 0; sector < 6; sector++) {
						   Int_t   sec                 = sector;//cand->getSector();
						   if (sector<5){
						   while(pid.phi_rpc<=(sector+1)*60)pid.phi_rpc+=60;
						   while(pid.phi_rpc>(sector+2)*60) pid.phi_rpc-=60;}

						   else
						   { while(pid.phi_rpc>60) pid.phi_rpc-=60;
						   while(pid.phi_rpc<=0) pid.phi_rpc+=60;
						   }
						 */


						if(!isInAcc_mar19(pid.theta_rpc,pid.phi_rpc)) continue;

						pid.phi_rpc                 = fakerpccluster->getPhi();
						pid.mom                     = mom;
						pid.sector                  = sec;
						pid.beta                    = beta;
						pid.emc_cell                =-100;
						pid.emc_size                = cand->getEmcNCells();
						pid.emc_energy              = emc_energy;
						pid.emc_energy_corr         = emc_energy;
						pid.emc_energy_centr        = emc_energy;
						pid.emc_mq                  = cand->getMetaMatchQualityEmc();
						pid.emc_path                = cand->getEmcCellPath();
						pid.emc_time                = cand->getEmcTime();
						pid.emc_match               = 0;
						pid.qualityDTime            = 1000;
						pid.within_resolution       = 100000;


						if(emc_index>=0) {
							HEmcClusterSim* cls  = (HEmcClusterSim*)emcClusterCat->getObject(emc_index);
							Int_t cell           =cls->getCell();

							if(isGoodCell(sec,getSystemMap(cell))) {

								if (pid.true_emc_hit_track==0&&(cflag&kIsInRPC)==kIsInRPC) { //if there is a hit for only RPC and no EMC hit
									Int_t cls_track     = cls->getTrack();
									HGeantKine * kine   = HCategoryManager::getObject(kine,geantKineCat,cls_track-1);
									HGeantKine*  prnt   = kine->getParent(kine->getTrack());
									Int_t trkParent     = prnt!=NULL ? prnt->getTrack() : 0;
									if(cand->getGeantTrackMeta(0)==trkParent) //if parent of EMC cluster is our particlecand in RPC
										pid.true_emc_hit_track=1;
								}
								pid.qualityDTime = cls->getQualDTime();
								pid.emc_match           = 1;
								pid.emc_cell            = getSystemMap(cell);
								pid.emc_energy_corr     = correctEMCenergy(emc_energy,cell);

								Float_t correction      = -47+exp(3.8-0.0123*mom);
								Float_t emc_energy_centr= pid.emc_energy_corr-correction;
								pid.emc_energy_centr=emc_energy_centr;
								Float_t energydiff      = (100*abs(emc_energy_centr-mom)/emc_energy_centr);
								pid.within_resolution   = energydiff*sqrt(mom/1000.)/3;
							}




							nt->Fill();
						}
					}
					//-------------------------------------------------
				}
			} // end cand loop
		} // end cand cat
		  //-------------------------------------------------

	} // end eventloop


	//-------------------------------------------------
	// added for the task
	//if(gHades)gHades->finalizeTasks();
	//-------------------------------------------------
	out->cd();
	Counter_cent->Write();
	out->Write();
	delete gHades;
	return 0;
}
