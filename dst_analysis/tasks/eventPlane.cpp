#include "/u/prozorov/mar19/diphotons/gen5a/dstAnalysis/include/Includes.h"
#include "hparticleevtchara.h"


Int_t ana(TString infileList, TString outputFile, Int_t nEvents = 10000000.)
{
  // -------------------------------------------------
  TROOT dst_analysis("DstAnalysisMacro", "Simple DST analysis Macro");
  HLoop * loop = new HLoop(kTRUE);
  gHades->setBeamTimeID(HADES::kMar19);


  // -------------------------------------------------
  // reading input files and decalring containers
  Bool_t ret = kFALSE;

  if (infileList.Contains(",")) {
    ret = loop->addMultFiles(infileList); // file1,file2,file3
  } else {
    ret = loop->addFiles(infileList); // myroot*.root
  }

  if (ret == 0) {
    cout << "READBACK: ERROR : cannot find inputfiles : " << infileList.Data() << endl;
    return 1;
  }

  TFile * out = new TFile(outputFile, "RECREATE");
  out->cd();

  if (!loop->setInput("-*,+HParticleEvtInfo, +HParticleCand, +HWallHit")) {
    // if(!loop->setInput("-*,+HParticleCand,+HEmcCluster,+HParticleEvtInfo")) {
    //	if(!loop->setInput("*")) {
    cout << "READBACK: ERROR : cannot read input !" << endl;
    exit(1);
  } // read all categories
  // --------------------------------------------------------------------------------
  // Setting the cache size of the HLoop internal TChain to read data from the DSTs
  // Improves performance of the lustre storage by decreasing load on lustre META servers
  // --------------------------------------------------------------------------------
  loop->getChain()->SetCacheSize(8e6); // 8Mb
  loop->getChain()->AddBranchToCache("*", kTRUE);
  loop->getChain()->StopCacheLearningPhase();
  loop->printCategories(); // Just for informative purposes
  // -------------------------------------------------

  // -------------------------------------------------
  HParticleTrackSorter sorter;
  sorter.init(); // get catgegory pointers etc...
  // -------------------------------------------------

  // -------------------------------------------------
  // input data
  HCategory * evtInfoCat = (HCategory *) HCategoryManager::getCategory(catParticleEvtInfo);
  HParticleEvtChara evtChara;


  TString ParameterfileCVMFS =
    "/lustre/hades/user/bkardan/param/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root";
  if (!evtChara.setParameterFile(ParameterfileCVMFS)) {
    cout << "Parameterfile not found !!! " << endl;
    return kFALSE;
  }
  evtChara.setPostFix("");
  if (!evtChara.init()) {
    cout << "HParticleEvtChara not init!!! " << endl;
    return kFALSE;
  }

  evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5);
  evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);


  // -------------------------------------------------

  TH1F * hEventPlaneCharaAll = new TH1F("hEventPlaneCharaAll", ";#Psi_{RP}[deg];Counts", 200, -200, 200);
  TH1F * hEventPlaneInfo     = new TH1F("hEventPlaneInfo", ";#Psi_{RP}[deg];Counts", 360, -180, 180);

  const Int_t nCentralityClasses = 4;
  TH1F * hDeltaABcent10[nCentralityClasses];
  TH1F * hDeltaABcent5[nCentralityClasses * 2];
  TH1F * hEventPlaneChara[nCentralityClasses];

  for (Int_t k = 0; k < nCentralityClasses; k++) {
    hDeltaABcent10[k] = new TH1F(Form("hDeltaABcent10_%d",
        k), ";#Delta#Psi=#Psi_{A}-#Psi_{B},[#circ];Counts", 360, -180, 180);
    hEventPlaneChara[k] = new TH1F(Form("hEventPlaneChara%d", k), ";#Psi_{RP}[#circ];Counts", 200, -200, 200);
  }
  for (Int_t k = 0; k < nCentralityClasses * 2; k++) hDeltaABcent5[k] = new TH1F(Form("hDeltaABcent5_%d",
        k), ";#Delta#Psi=#Psi_{A}-#Psi_{B},[#circ];Counts", 360, -180, 180);

  // ==========================================================================================//

  // ==========================================================================================//
  // ====================================EVENT LOOP============================================//
  // ==========================================================================================//

  // -------------------------------------------------
  // event loop starts here
  Int_t entries = loop->getEntries();
  if (nEvents < entries && nEvents >= 0) entries = nEvents;

  for (Int_t i = 0; i < entries; i++) {
    HTool::printProgress(i, entries, 1, "Analyze pairs :");
    Int_t nbytes = loop->nextEvent(i); // get next event. categories will be cleared before
    if (nbytes <= 0) { cout << nbytes << endl; break; } // last event reached

    // -------------------------------------------------
    // summary event info object
    HParticleEvtInfo * evtInfo = HCategoryManager::getObject(evtInfo, catParticleEvtInfo, 0);
    // Extracting event plane
    Float_t eventPlaneInfo = evtInfo->getRPlanePhi();
    hEventPlaneInfo->Fill(eventPlaneInfo);

    // Event Plane resolution:

    Int_t cent_bin  = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10) - 1;
    Int_t cent_bin5 = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5) - 1;

    if (cent_bin < 0 || cent_bin > 3) continue;


    Float_t eventPlaneChara = evtChara.getEventPlane(HParticleEvtChara::kDefault) * TMath::RadToDeg();
    Float_t EventPlaneA     = evtChara.getEventPlane(HParticleEvtChara::kDefault, 1) * TMath::RadToDeg();
    Float_t EventPlaneB     = evtChara.getEventPlane(HParticleEvtChara::kDefault, 2) * TMath::RadToDeg();


    if (eventPlaneChara * TMath::DegToRad() == -1)
      continue;
    if (EventPlaneA * TMath::DegToRad() == -1 || EventPlaneB * TMath::DegToRad() == -1)
      continue;

    if (evtChara.getEventPlane(HParticleEvtChara::kDefault) == -1) eventPlaneChara = -1000;
    if (eventPlaneChara > 180.) eventPlaneChara -= 360.;

    Float_t deltaPhiAB = (EventPlaneA - EventPlaneB) * TMath::RadToDeg();
    if (EventPlaneA >= 0. && EventPlaneB >= 0.) {
      deltaPhiAB = ( EventPlaneA - EventPlaneB );
      //	if(deltaPhiAB > 180.) deltaPhiAB = 360. - deltaPhiAB; //this makes sure that phiAB is in (0,180)
      hDeltaABcent10[cent_bin]->Fill(deltaPhiAB);
      if (cent_bin5 < nCentralityClasses * 2) hDeltaABcent5[cent_bin5]->Fill(deltaPhiAB);
    }

    hEventPlaneCharaAll->Fill(eventPlaneChara);
    hEventPlaneChara[cent_bin]->Fill(eventPlaneChara);


    // -------------------------------------------------
  } // end eventloop


  // -------------------------------------------------
  out->cd();
  //	out->GetList()->Write();
  hEventPlaneCharaAll->Write();
  hEventPlaneInfo->Write();
  for (Int_t k = 0; k < nCentralityClasses; k++) {
    hDeltaABcent10[k]->Write();
    hEventPlaneChara[k]->Write();
  }
  for (Int_t k = 0; k < nCentralityClasses * 2; k++) {
    hDeltaABcent5[k]->Write();
  }
  out->Write();
  return 0;
} /* fill_ttree_diphoton */
