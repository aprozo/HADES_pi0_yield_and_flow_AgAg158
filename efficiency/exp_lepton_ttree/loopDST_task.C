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
#include "hparticleevtchara.h"
#include "TMath.h"
#include "hmetamatcher.h"
#include "TVector3.h"
#include "Includes.h"
#include <iostream>
#include <vector>
using namespace std;


const Bool_t isSim = kTRUE;

typedef conditional<isSim, HParticleCandSim, HParticleCand>    ::type HParticleCandClass;
typedef conditional<isSim, HEmcClusterSim, HEmcCluster>    ::type HEmcClusterClass;
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
	#if isSim
		TString rootParFile      = "/cvmfs/hadessoft.gsi.de/param/sim/mar19//cvmfs/hadessoft.gsi.de/param/sim/mar19/allParam_mar19_sim_run_16000_gen4_25112020.root";
		TString paramrelease     = "MAR19_dst_gen4";// now, APR12_gen2_dst APR12_gen5_dst
#else
		//TString rootParFile      = "/cvmfs/hadessoft.gsi.de/param/real/mar19/allParam_Mar19_gen5_13042021.root";
		TString rootParFile      = "/cvmfs/hades.gsi.de/param/real/mar19/allParam_Mar19_gen5_13042021.root";
		TString paramrelease     = "MAR19_dst_gen5";    // now, APR12_gen2_dst APR12_gen5_dst
  #endif
		TString paramSource      = "root";// root, ascii, oracle


		HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,emc,wall,start,tbox");
		HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease); // now, APR12_gen2_dst
	}




	//-------------------------------------------------

	//-------------------------------------------------
	// list of all files with working sectors
	//if(!isSim) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);
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
		Float_t mult_tofrpc,imppar,cent_bin,charge,sector,
		        emc_cell,theta_rpc,phi_rpc,mom,beta,mom_kine,
		        emc_size,emc_energy,emc_time,emc_mq,qualityDTime,
		        emc_path,mm_qa,geantid,within_resolution,
		        true_emc_hit_track,true_rpc_hit_track, cflag,beamtime;
		void reset(){
			qualityDTime = 1000;
			mult_tofrpc=-1000;
			imppar=-1000;
			cent_bin=-1000;
			charge=-1000;
			sector=-1000;
			emc_cell=-1000;
			theta_rpc=-1000;
			phi_rpc=-1000;
			mom=-1000;
			beta=-1000;
			mom_kine=-1000;
			emc_size=-1000;
			emc_energy=-1000;
			emc_time=-1000;
			emc_mq=-1000;
			emc_path=-1000;
			mm_qa=-1000;
			geantid=-1000;
			within_resolution=-1000;
			true_emc_hit_track=-1000;
			true_rpc_hit_track=-1000;
			cflag=-1000;
			beamtime=-1000;
		};
	} PID;


	static PID pid;

	TTree *nt = new TTree("nt","pid");
	nt->Branch("pid",&pid, "mult_tofrpc:imppar:cent_bin:charge:sector:"
	           "emc_cell:theta_rpc:phi_rpc:mom:beta:mom_kine:"
	           "emc_size:emc_energy:emc_time:emc_mq:qualityDTime:"
	           "emc_path:mm_qa:geantid:within_resolution:"
	           "true_emc_hit_track:true_rpc_hit_track: cflag:beamtime");

	if(!loop.setInput("-*,+HParticleCand,+HEmcCluster,+HParticleEvtInfo,+HParticleCandSim,+HEmcClusterSim,+HGeantKine")) {
		//if(!loop.setInput("-*,+HParticleCand,+HEmcCluster,+HParticleEvtInfo")) {
//	if(!loop.setInput("*")) {
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
	HCategory* evtInfoCat   = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
	HCategory* geantKineCat = (HCategory*)HCategoryManager::getCategory(catGeantKine);
	//-------------------------------------------------


	//-------------------------------------------------
	//parameters
	//HEnergyLossCorrPar dEdxCorr;
	//dEdxCorr.setDefaultPar("mar19");
	TFile* acc       = new TFile("acceptancePhiThetaMar19.root","Read");
	acceptance_mar19 =(TH2D*)acc->Get("acceptancePhiThetaMar19");
	setupMassCuts("m2_cut.root");

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

	TString tempFileName,currentFileName;
	Int_t dayOfYear, hour,min,second;


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

	for (Int_t i = 0; i < entries; i++) {
		HTool::printProgress(i,entries,1,"Analyze pairs :");
		Int_t nbytes =  loop.nextEvent(i);         // get next event. categories will be cleared before

		if(nbytes <= 0) { cout<<nbytes<<endl; break; } // last event reached
		//if(i%1000 == 0) cout<<"event "<<i<<endl;

		if (loop.isNewFile(tempFileName)) {
			TString type;
			Int_t year;
			Int_t eb;
			currentFileName = tempFileName;
			HTime::splitFileName(HTime::stripFileName(tempFileName),type,year,dayOfYear,hour,min,second,eb);
		}
		Float_t valueMinute = dayOfYear + hour/24. + min/24./60.;
		Float_t valueSecond = valueMinute + second/24./60./60.;
		Float_t BeamTime=valueSecond;

		//-------------------------------------------------
		// summary event info object

		HGeomVector vertex(evHeader->getVertexCluster().getPos()); // Event vertex position
		HParticleEvtInfo* evtInfo = HCategoryManager::getObject(evtInfo,catParticleEvtInfo,0);
		Int_t mult_tofrpc         = evtInfo->getSumTofMultCut() + evtInfo->getSumRpcMultHitCut();

		Int_t cent_bin            = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);
		Float_t ImpactParameter = 0;

		if(isSim) {
			fSimul = ( (HRecEvent*)gHades->getCurrentEvent())->getPartialEvent(catSimul);
			fSubHeader = (HGeantHeader*)(fSimul->getSubHeader());
			if (fSubHeader) {ImpactParameter = fSubHeader->getImpactParameter();
				             cent_bin = getCentBin(ImpactParameter);}
		}

		if(cent_bin==0) continue;
		cent_bin      -= 1;
	  if (cent_bin > 4) continue;
		Counter_cent->Fill(cent_bin);
		Counter_cent->Fill(4);

		pid.imppar=ImpactParameter;
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
			HParticleCandClass* cand=0;

			for(Int_t j = 0; j < size; j++) {

				//------------------------------------------------------------------------
				cand = HCategoryManager::getObject(cand,candCat,j);
				if(cand) {
					if(!leptonParticlecandCut(cand)) continue;    // see  #include "lepton_efficiency.h"
					if(!cand->isFlagBit(kIsUsed)) { continue; }   // skip particles rejected by HParticleTrackSorter
					if(cand->isRpcClstUsed()) {
						pid.reset();
						HRpcCluster* fakerpccluster = matcher->matchWithRpc(cand);
						pid.theta_rpc               = fakerpccluster->getTheta();
						pid.phi_rpc                 = fakerpccluster->getPhi();
						if(!isInAcc_mar19(pid.theta_rpc,pid.phi_rpc)) continue;


						pid.mult_tofrpc=mult_tofrpc;
						pid.cent_bin=cent_bin;
						pid.imppar=ImpactParameter;
						pid.beamtime=BeamTime;

						Float_t mom                 = cand->getMomentum();

						Int_t emc_index             = cand->getEmcInd();
						Float_t beta                = cand->getBeta();
						pid.charge                  = cand->getCharge();
						pid.mm_qa                   = cand->getMetaMatchQuality();
						pid.mom                     = mom;
						Int_t sec                   = cand->getSector();
						pid.sector                  = sec;
						pid.beta                    = beta;


						#if isSim
						UInt_t cflag = cand->getGeantCorrTrackIds();
						pid.true_rpc_hit_track      = (cflag&kIsInRPC)==kIsInRPC;
						pid.true_emc_hit_track      = (cflag&kIsInRPC)==kIsInRPC &&
						                              (cflag&kIsInEMC)==kIsInEMC;
						pid.mom_kine                = cand->getGeantTotalMom();
						pid.cflag                   = cand->getGeantCorrTrackIds();
						pid.geantid                 = cand->getGeantPID();
						#endif

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

						if(emc_index>=0) {
							HEmcClusterClass* cls  = (HEmcClusterClass*)emcClusterCat->getObject(emc_index);

							Float_t emc_energy          = cand->getEmcEnergy();
							pid.emc_size                = cand->getEmcNCells();
							pid.emc_energy              = emc_energy;
							pid.emc_mq                  = cand->getMetaMatchQualityEmc();
							pid.emc_path                = cand->getEmcCellPath();
							pid.emc_time                = cand->getEmcTime();
							#if isSim
							if (pid.true_emc_hit_track==0&&(cflag&kIsInRPC)==kIsInRPC) {     //if there is a hit for only RPC and no EMC hit
								Int_t cls_track        = cls->getTrack();
								HGeantKine * kine   = HCategoryManager::getObject(kine,geantKineCat,cls_track-1);
								HGeantKine*  prnt   = kine->getParent(kine->getTrack());
								Int_t trkParent     = prnt!=NULL ? prnt->getTrack() : 0;
								if(cand->getGeantTrackMeta(0)==trkParent)     //if parent of EMC cluster is our particlecand in RPC
									pid.true_emc_hit_track=1;
							}
							pid.emc_energy    = correctEMCenergy(emc_energy,cell);
							Float_t correction      = -47+exp(3.8-0.0123*mom);
							Float_t emc_energy_centr= pid.emc_energy-correction;
							emc_energy=emc_energy_centr;
							#endif
							pid.qualityDTime = cls->getQualDTime();
							pid.emc_cell            = getSystemMap(cls->getCell());
							Float_t energydiff      = (100*abs(emc_energy-mom)/emc_energy);
							pid.within_resolution   = energydiff*sqrt(mom/1000.)/3;
						}
						nt->Fill();
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
