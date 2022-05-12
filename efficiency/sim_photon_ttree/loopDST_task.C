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
TH2F *acceptance_mar19;

Bool_t isInAcc_mar19(const Float_t & theta,  const Float_t & phi) {
	UInt_t bin = 0;
	bin = acceptance_mar19->FindBin(phi,theta);
	Float_t w = acceptance_mar19->GetBinContent(bin);
	//if(theta<=15) return kFALSE;
	if(w>0.9) return kTRUE;
	//if(theta>=12.5&&theta<=17&&w>0.085) return kTRUE; ///for electrons
	return kFALSE;
}

HLoop* CreateLoop(TString inputfile)
{
	HLoop* loop = new HLoop(kTRUE);
	Bool_t ret =kFALSE;
	if(inputfile.Contains(",")) {
		ret = loop->addMultFiles(inputfile);         // file1,file2,file3
	}
	else{    ret = loop->addFiles(inputfile);    // myroot*.root
	}
	if(ret == 0) {
		cout<<"READBACK: ERROR : cannot find inputfiles : "<<inputfile.Data()<<endl;
		return NULL;
	}

	if(!loop->setInput()) return NULL;
	loop->printCategories();
	loop->printChain();
	return loop;
}

HCategory* GetCategory(const TString categoryName, HLoop* loop)
{
	HCategory* category = loop->getCategory(categoryName);
	if(!category)
	{
		cout<<"Category " << categoryName << " is NULL!" << endl;
		return 0;
	}
	return category;
}


typedef struct
{
	Int_t cent_bin;
	Int_t sector;
	Int_t emc_cell;
	Int_t rowEMC;
	Int_t columnEMC;

	Int_t id_geant;
	Float_t theta_kine;
	Float_t phi_kine;

	Float_t x_kine;
	Float_t y_kine;
	Float_t mom_kine;
	Int_t matchedCells;
	Int_t matchedPart;
	Float_t emc_energy;
	Float_t within_resolution;
	Float_t emc_beta;
	Int_t emc_size;
	Float_t xEMC;
	Float_t yEMC;

	void reset()
	{
		//cent_bin           = -1;
		sector             = -1;
		emc_cell           = -1;
		theta_kine         = -1.;
		phi_kine           = -1.;
		mom_kine           = -1.;
		id_geant           =  0;
		x_kine             =  -99999.;
		y_kine             =  -99999.;
		mom_kine           =  0.;
		matchedCells       =  0;
		matchedPart        =  0;
		rowEMC             = -1;
		columnEMC          = -1;
		emc_energy         = -1.;
		within_resolution  =100.;
		emc_beta           =  0.;
		emc_size           =  0;
		xEMC               =-99999.;
		yEMC               =-99999.;
	}
}TREE_VARIABLES;

TTree* makeTree(TREE_VARIABLES& vars)
{
	TTree* nt = new TTree("nt","Simulated tree for prediction of EMC cell");
	nt->Branch("cent_bin",&vars.cent_bin,"cent_bin/I"         );
	nt->Branch("sector",&vars.sector,"sector/I"         );
	nt->Branch("emc_cell",&vars.emc_cell,"emc_cell/I"          );
	nt->Branch("rowEMC",&vars.rowEMC,"rowEMC/I"         );
	nt->Branch("columnEMC",&vars.columnEMC,"columnEMC/I"         );
	nt->Branch("id_geant",&vars.id_geant,"id_geant/I"         );
	nt->Branch("theta_kine",&vars.theta_kine,"theta_kine/F"          );
	nt->Branch("phi_kine",&vars.phi_kine,"phi_kine/F"          );
	nt->Branch("x_kine",&vars.x_kine,"x_kine/F"          );
	nt->Branch("y_kine",&vars.y_kine,"y_kine/F"          );
	nt->Branch("mom_kine",&vars.mom_kine,"mom_kine/F"         );
	nt->Branch("matchedCells",&vars.matchedCells,"matchedCells/I"         );
	nt->Branch("matchedPart",&vars.matchedPart,"matchedPart/I"         );
	nt->Branch("emc_energy",&vars.emc_energy,"emc_energy/F"         );
	nt->Branch("within_resolution",&vars.within_resolution,"within_resolution/F"         );
	nt->Branch("emc_beta",&vars.emc_beta,"emc_beta/F"         );
	nt->Branch("emc_size",&vars.emc_size,"emc_size/I"         );
	nt->Branch("xEMC",&vars.xEMC,"xEMC/F"      );
	nt->Branch("yEMC",&vars.yEMC,"yEMC/F"      );

	return nt;
}



//=====================================================================================================================================================//
//____________________________________________   Main function   ______________________________________________________________________________________//

Int_t loopDST_task( TString infileList,TString outfile,Int_t nEvents)
{
	//=================================================================================================================================================//
	// Initialization of gHades                                                                                                                        //
	//=================================================================================================================================================//
	HLoop* loop = CreateLoop(infileList);

	//HCategory* emcClusterCat  = GetCategory("HEmcClusterSim",   loop);//HEmcClusterSim*   emcCluster = 0; if(!emcClusterCat ){delete gHades; return 1;} //
	//HCategory* geantKineCat   = GetCategory("HGeantKine",       loop);
	//HCategory* candCat    = GetCategory("HParticleCandSim", loop);

	HCategory* candCat      = (HCategory*)HCategoryManager::getCategory(catParticleCand);
	HCategory* emcClusterCat= (HCategory*)HCategoryManager::getCategory(catEmcCluster);
	HCategory* geantKineCat      = (HCategory*)HCategoryManager::getCategory(catGeantKine);

	HGeantKine*       geantKine  = 0; if(!geantKineCat  ) {delete gHades; return 1;}     //
	HParticleCand*    cand       = 0; if(!candCat   ) {delete gHades; return 1;}     //
	//-------------------------------------------------


	//==================================================================================================================//
	// Simulation parameter
	//==================================================================================================================//
	{  Bool_t isSimulation = kFALSE;
	   TString beamtime="mar19";
	   Int_t mdcMods[6][4]=
	   { {1,1,1,1},
		 {1,1,1,1},
		 {1,1,1,1},
		 {1,1,1,1},
		 {1,1,1,1},
		 {1,1,1,1} };
	   TString asciiParFile     = "";
	   TString rootParFile      = "/cvmfs/hadessoft.gsi.de/param/sim/mar19/allParam_mar19_sim_run_16000_gen4_25112020.root";
	   TString paramSource      = "root";         // root, ascii, oracle
	   TString paramrelease     = "MAR19_dst_gen4";         // now, APR12_gen2_dst APR12_gen5_dst
	   HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,emc,wall,start,tbox");
	   HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease);          // now, APR12_gen2_dst
	}
	//-------------------------------------------------

	//==================================================================================================================//
	// Output                                                                                                           //
	//==================================================================================================================//
	TFile* outputFile = new TFile(outfile, "RECREATE");                                                              //
	TREE_VARIABLES vars;                                                                                               //
	TTree* nt = makeTree(vars);                                                                                        //
	//------------------------------------------------------------------------------------------------------------------//



	//parameters
	TFile* acc       = new TFile("acceptance_gen4.root","Read");
	acceptance_mar19 =(TH2F*)acc->Get("Acceptance2D_expsetup");
	setupSimEffCorr("corretion_energy_ag1650ag.root");
	HPartialEvent* fSimul = NULL;
	HGeantHeader* fSubHeader = NULL;
	//------------------------------------------------------------------------------------------------------------------//


	//==================================================================================================================//
	// Recalculation of RpcCluster needs to run HParticleMetaMatcher first                                              //
	//==================================================================================================================//

	HTaskSet *masterTaskSet = gHades->getTaskSet("all");
	HMetaMatcher* matcher = new HMetaMatcher();
	//matcher->setDebug();
	masterTaskSet->add(matcher);

	HParticleTrackSorter sorter;
	sorter.init();
	//-------------------------------------------------


	//==================================================================================================================//
	// nEvents -centrality histogram                                                                                                         //
	//==================================================================================================================//
	enum Counters_cent {
		cent0_10,
		cent10_20,
		cent20_30,
		cent30_40,
		all_cents,
		cNumCounters_cent
	};

	TH1D* Counter_cent = new TH1D( "Counter_cent", "Counter_cent", cNumCounters_cent, 0, cNumCounters_cent);

	for (Int_t c = 0; c < cNumCounters_cent; c++) {
		Counter_cent->GetXaxis()->SetBinLabel( 1, "cent0_10" );
		Counter_cent->GetXaxis()->SetBinLabel( 2, "cent10_20" );
		Counter_cent->GetXaxis()->SetBinLabel( 3, "cent20_30" );
		Counter_cent->GetXaxis()->SetBinLabel( 4, "cent30_40" );
		Counter_cent->GetXaxis()->SetBinLabel( 5, "all_cent" );
	}
	//------------------------------------------------------------------------------------------------------------------//


	//******************************************************************************************************************//
	// Event Loop                                                                                                       //
	//******************************************************************************************************************//

	Int_t entries = loop->getEntries();
	HEventHeader *evHeader  = loop->getEventHeader();
	if(nEvents < entries && nEvents >= 0 ) entries = nEvents;
	for (Int_t i = 0; i < entries; i++) {
		//======================================================================//
		// Selection of events                                                  //
		//======================================================================//
		Int_t nbytes =  loop->nextEvent(i);        // get next event. categories will be cleared before
		if(nbytes <= 0) { cout<<nbytes<<endl; break; }
		// last event reached
		//if(i%1000 == 0) cout<<"event "<<i<<endl;
		HGeomVector vertex(evHeader->getVertexCluster().getPos());         // Event vertex position
		fSimul     = ( (HRecEvent*)gHades->getCurrentEvent())->getPartialEvent(catSimul);
		fSubHeader = (HGeantHeader*)(fSimul->getSubHeader());
		Float_t impactParameter = -1.;
		if (fSubHeader) impactParameter = fSubHeader->getImpactParameter();
		UInt_t cent_bin=0;
		if(impactParameter<3.8)                    {cent_bin=cent0_10;}
		if(impactParameter>=3.8&&impactParameter<5.4) {cent_bin=cent10_20;}
		if(impactParameter>=5.4&&impactParameter<6.6) {cent_bin=cent20_30;}
		//if(impactParameter>=6.6&&impactParameter<7.6) {cent_bin=cent30_40;}
		if(impactParameter>=6.6) continue;         //centrality 22-44% (5.5<b<8)
		Counter_cent->Fill(cent_bin);
		Counter_cent->Fill(all_cents);

		vars.cent_bin=cent_bin;

		//=======================================================================================//
		// particle candidate  loop + sorter                                                      //
		//=======================================================================================//
		Int_t emc_ModMult[6][163]={0};
		sorter.cleanUp();
		sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);         // reset all flags for flags (0-28) ,reject,used,lepton
		Int_t nCandLep     = sorter.fill(HParticleTrackSorter::selectLeptons);        // fill only good leptons
		Int_t nCandLepBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsLeptonSorter);
		Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);        // fill only good hadrons (already marked good leptons will be skipped)
		Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsHadronSorter);


		for(Int_t j = 0; j < candCat->getEntries(); j++) {         // loop over Particlecand
			cand = HCategoryManager::getObject(cand,candCat,j);
			if(!cand->isFlagBit(Particle::kIsUsed)) continue;
			Int_t emc_index = cand->getEmcInd();
			//Double_t emc_mq = cand->getMetaMatchQualityEmc();
			if(emc_index >=0) {
				HEmcClusterSim* emcClst = (HEmcClusterSim*)emcClusterCat->getObject(emc_index);
				Int_t sec = emcClst->getSector();
				Int_t cell= getSystemMap(emcClst->getCell());
				emc_ModMult[sec][cell]++;
			}
		}

		//=======================================================================================//
		// HGeantKine loop                                                                       //
		//=======================================================================================//
		for(Int_t j = 0; j < geantKineCat->getEntries(); j++) {
			vars.reset();
			geantKine = HCategoryManager::getObject(geantKine,geantKineCat,j);
			if(geantKine->getID()!=1) continue;
			vars.id_geant=1;
			HGeantKine * mother  = geantKine->getParent(geantKine->getTrack());
			if(!(mother && mother->getID()==7 && mother->isPrimary())) continue;
			Float_t mom_kine = geantKine->getTotalMomentum();
			Float_t phi_kine, theta_kine;
			Bool_t fromparent=kFALSE;             //flag that emc saw a hit from parent photon
			Bool_t directhit =kFALSE;

			geantKine->getPhiThetaDeg(theta_kine,phi_kine);
			if(!isInAcc_mar19(theta_kine,phi_kine)) continue;
			if(mom_kine<20) continue;


			vars.theta_kine=theta_kine;
			vars.phi_kine  =phi_kine;
			vars.mom_kine  =mom_kine;
			vars.x_kine    =sin(theta_kine)*sin(phi_kine);
			vars.y_kine    =sin(theta_kine)*cos(phi_kine);
			vars.sector    =geantKine->getSector();

			Int_t geantKineTrack = geantKine->getTrack();

			Float_t length_to_cell=0;

			Int_t matchedPart=-1,matchedCells=-1;
			//===============================================================================//
			// EMC Cluster                                                                   //
			//===============================================================================//
			Double_t beta;
			HEmcClusterSim *emcCluster;
			for(Int_t k = 0; k < emcClusterCat->getEntries(); k++)   {

				emcCluster =HCategoryManager::getObject(emcCluster,emcClusterCat,k);                 // get an cluster emc object
				HGeomVector cellVector(emcCluster->getXLab(),emcCluster->getYLab(),emcCluster->getZLab());
				Int_t emcCluster_track     =emcCluster->getTrack();                // get the track of cluster emc object
				Int_t emc_cell      =getSystemMap(emcCluster->getCell());                // get cell of cluster emc object;
				Int_t sector        =emcCluster->getSector();                // get cell of cluster emc object;
				if(!isGoodCell(sector,emc_cell)) continue;
				Double_t trLength    = cellVector.length();
				beta        = (trLength/1000.)/(emcCluster->getTime()*1.e-9) / TMath::C();
				matchedPart = emc_ModMult[sector][emc_cell];
				matchedCells= emcCluster->getNMatchedCells();

				HGeantKine * kine2  =HCategoryManager::getObject(kine2,geantKineCat,emcCluster_track-1);                // find a kine track corresponding to that cluster emc object
				Int_t kine2_track=kine2->getTrack();

				if( kine2_track==geantKineTrack)  {directhit=kTRUE; break;}

				HGeantKine* prnt    = kine2->getParent(kine2_track);
				Int_t trkMother = prnt!=NULL ? prnt->getTrack() : 0;                 // get a track  of parnet of geant emc object track
				if (trkMother==geantKineTrack&&trkMother!=0) {                 //parent track == photon track  - flag that emc saw a hit from parent photon
					HGeomVector parent_photon_scat_vertex;
					kine2->getVertex(parent_photon_scat_vertex);
					parent_photon_scat_vertex-=cellVector;
					length_to_cell=parent_photon_scat_vertex.length();
					if(length_to_cell<330) {fromparent=kTRUE;
						                    vars.id_geant=kine2->getID();}
					break;
				}

				//===============================================================================//
				// if EMC Cluster   WAS FOUND                                                    //
				//===============================================================================//

				if(fromparent||directhit) {
					vars.emc_cell  =getSystemMap(emcCluster->getCell());
					vars.columnEMC = emcCluster->getCell()%17;
					vars.rowEMC    = floor(emcCluster->getCell()/17);
					vars.sector    = emcCluster->getSector();
					vars.emc_size  = emcCluster->getNCells();
					vars.xEMC      = emcCluster->getXMod();
					vars.yEMC      = emcCluster->getYMod();

					Float_t clustEnergy=correctEMCenergy(emcCluster->getEnergy(),emcCluster->getCell());
					vars.emc_energy=clustEnergy;
					vars.matchedPart=matchedPart;
					vars.matchedCells=matchedCells;
					vars.emc_beta=beta;

					//Float_t resolution = 3*5.5/sqrt(mom_kine/1000.);
					Float_t energydiff = (100*abs(clustEnergy-mom_kine)/clustEnergy);
					vars.within_resolution=energydiff*sqrt(mom_kine/1000.)/3;
				}
				//=======================================================================================//
			}            // end of EMC Cluster Loop
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
			nt->Fill();
		}         //end of HGeantKine Loop
		          //-------------------------------------------------------------------------------------------------------------//

	}         // end of Event Loop                                                                                           //
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

	//========================//
	// Write tree to the file //
	//========================//
	outputFile->cd();
	Counter_cent->Write();
	outputFile->cd();             //
	nt->Write();                  //
	outputFile->Close();          //
	//------------------------//

	//------------------------//
	//========================//
	// close program          //
	//========================//
	delete gHades;
	return 0;
	//------------------------//
}
