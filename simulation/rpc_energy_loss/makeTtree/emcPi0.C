#include "hades.h"
#include "htool.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "hgeantkine.h"
#include "TH1.h"
#include "TH2.h"

#include "hemcclustersim.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hparticleevtinfo.h"


#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"

#include "TMath.h"
#include "TFile.h"


#include <iostream>
#include <map>
#include <stdio.h>

//#include "particle.h"      // definition of object
//#include "particleKine.h"  // definition of object


Int_t emcPi0(TString inFile, TString outputFile, Int_t nev=10000000.) {


    TROOT dst_production( "DstProductionApp", "Produce a DST File" );
    Bool_t createHades = kTRUE; //kFALSE;  // kTRUE = create HADES new
    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)

    loop->addMultFiles(inFile);
    if(!loop->setInput("-*,+HEmcCluster,+HGeantKine")) { exit(1); }
    loop->printCategories();

    //------------------------------------------------------------------
    // setting up the output

    TFile* out = new TFile(outputFile,"RECREATE");
    out->cd();
   
    TH2F* devMomThetaVsPhi      = new TH2F("devMomThetaVsPhi",     "devMomThetaVsPhi",     60, 0, 60, 60, 0, 60);
    TH2F* devMomThetaVsPhiN     = new TH2F("devMomThetaVsPhiN",    "devMomThetaVsPhiN",    60, 0, 60, 60, 0, 60);
    
    TH2F* devMomCellVsEne_noABS = new TH2F("devMomCellVsEne_noABS","devMomCellVsEne_noABS",255, 0, 255, 200, 0, 2000);
    TH2F* devMomCellVsEneN      = new TH2F("devMomCellVsEneN",     "devMomCellVsEneN",     255, 0, 255, 200, 0, 2000);

    TH2F* devMomVsEne = new TH2F("devMomVsEne","devMomVsEne",200, 0, 2000, 200, 0, 2000);
    
    TH2F* energy_vs_correction = new TH2F("energy_vs_correction","energy_vs_correction",200,-1000,1000,1000, -500, 500);
    
    TH1F* devMomCell            = new TH1F("devMomCell",            "devMomCell",          255, 0, 255);
    TH1F* devMomCell_noABS      = new TH1F("devMomCell_noABS",      "devMomCell_noABS",    255, 0, 255);
    TH1F* devMomCellN           = new TH1F("devMomCellN",           "devMomCellN",         255, 0, 255);

    TH2F* clEnVsMomGamma        = new TH2F("clEnVsMomGamma",        "#gamma;p, MeV/c;Cluster Energy, MeV",200,0.,2000., 100,0.,2000.);
    TH2F* clEnVsMomGammaCorr    = new TH2F("clEnVsMomGammaCorr",    "#gammaCorr;p, MeV/c;Cluster Energy, MeV",200,0.,2000., 100,0.,2000.);
    
    HEventHeader *evHeader = loop->getEventHeader();
    HGeantHeader* fSubHeader = NULL;
    HPartialEvent* fSimul = NULL;
    HCategory *catEmcClus = loop->getCategory("HEmcCluster");
    HCategory *catGeaKine = loop->getCategory("HGeantKine");

    Int_t nEntries = loop->getEntries();
    printf("nEntries = %i\n",nEntries);
    if(nEntries>nev && nev>0) nEntries = nev;
    for (Int_t i=0; i<nEntries ; i++) if(loop->nextEvent(i) > 0)  {

	HGeomVector vertex(evHeader->getVertexCluster().getPos());   // Event vertex position
	HEventHeader* header = loop->getEventHeader();
	fSimul = ( (HRecEvent*)gHades->getCurrentEvent())->getPartialEvent(catSimul);
	fSubHeader = (HGeantHeader*)(fSimul->getSubHeader());
	Int_t nEmcCluster = catEmcClus->getEntries();

	for(Int_t e=0; e<nEmcCluster; e++)
	{
	    HEmcClusterSim* cls = (HEmcClusterSim*)catEmcClus->getObject(e);

	    // Calculte Beta:
	    HGeomVector trackVec(cls->getXLab(),cls->getYLab(),cls->getZLab()+30.);

	    Double_t    trLength= trackVec.length();
	    Double_t    Beta    = (trLength/1000.)/(cls->getTime()*1.e-9) / TMath::C();
	    Int_t       track   = cls->getTrack();
	    HGeantKine* pKine   = (HGeantKine*)catGeaKine->getObject(track-1);
	    Float_t     mom     = pKine->getTotalMomentum();
	    Int_t       id      = pKine->getID();
	    Int_t       cell    = cls->getCell();	    
	    Float_t clustEnergy = cls->getEnergy();

	      
	    if(Beta > 0.95&&cell!=0) {
	      if(id==0){

		
		    devMomCell_noABS     ->Fill(cell,(mom - clustEnergy)/clustEnergy);
		    devMomCellN          ->Fill(cell);
		    
		    devMomCellVsEne_noABS->Fill(cell,mom,(mom - clustEnergy)/clustEnergy);
		    devMomCellVsEneN     ->Fill(cell,mom);
		    
		    energy_vs_correction ->Fill(mom,(mom - clustEnergy));
		    
		    devMomVsEne->Fill(mom,(mom - clustEnergy)/clustEnergy);
	      }
		       

		} // beta

	} // loop over cluster
	
    }// loop over events

	out->cd();
	out->GetList()->Write();
	out->Save();
	out->Close();

	return 0;

}


