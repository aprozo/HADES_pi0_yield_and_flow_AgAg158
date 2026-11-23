#include "Includes.h"

#include "exp_setup_cells.h"



Int_t analysis(TString infileList, TString outputFile, Int_t nev = 1000000.)
{
  TROOT dst_production("DstProductionApp", "Produce a DST File");
  Bool_t createHades = kTRUE;
  HLoop  loop(createHades); // create HADES (needed if one wants to use gHades)
  gHades->setBeamTimeID(HADES::kMar19);
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

  if (!loop.setInput("-*,+HGeantKine,+HParticleCand,+HEmcCluster,+HParticleEvtInfo")) {
    cout << "READBACK: ERROR : cannot read input !" << endl;
    exit(1);
  }
  loop.printCategories();
  loop.printChain();

  HCategory * candCat       = (HCategory *) HCategoryManager::getCategory(catParticleCand);
  HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
  HCategory * geantKineCat  = (HCategory *) HCategoryManager::getCategory(catGeantKine);

  HGeantHeader * fSubHeader = NULL;
  HPartialEvent * fSimul    = NULL;

  Bool_t isSimulation = kFALSE;
  if (geantKineCat) { isSimulation = kTRUE;  Info("init()", "GeantKine found - is Simulation"); } else {
    isSimulation = kFALSE;
    Info("init()", "GeantKine not found - is not Simulation");
  }
  setupSimEffCorr("simCorrEnergyEcal.root");
  // setupSimEffCorr("corretion_energy_ag1650ag.root");
  TFile * out = new TFile(outputFile, "RECREATE");
  out->cd();


  PID pid;
  TTree nt("nt", "pid");

  nt.Branch("pid", &pid, "newEvent:multTofRpc:multEmc:multPhoton:centBin:targetBin:beamTime:RPlanePhi:impPar:emc_mq:qualityDTime:"
    "sector:cell:size:beta:energy:time:xpos:ypos:zpos:track:trackParent:mom_kine:idParent");

  enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
  TH1D * Counter_cent = new TH1D("Counter_cent", "Counter_cent", cNumCounters_cent, 0, cNumCounters_cent);

  Counter_cent->GetXaxis()->SetBinLabel(5, "all_cent");
  Counter_cent->GetXaxis()->SetBinLabel(1, "cent0_10");
  Counter_cent->GetXaxis()->SetBinLabel(2, "cent10_20");
  Counter_cent->GetXaxis()->SetBinLabel(3, "cent20_30");
  Counter_cent->GetXaxis()->SetBinLabel(4, "cent30_40");

  // --------------------------CONFIGURATION---------------------------------------------------
  // At begin of the program (outside the event loop)
  HParticleTrackSorter sorter;
  sorter.setIgnoreInnerMDC(); // do not reject Double_t inner MDC hits
  sorter.init();              // get catgegory pointers etc...
  // =======================================================================
  Int_t nEntries = loop.getEntries();
  printf("nEntries = %i\n", nEntries);
  HParticleEvtChara evtChara;
  if (!evtChara.setParameterFile(
      "/lustre/hades/user/bkardan/param/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root"))
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

  for (Int_t i = 0; i < nEntries; i++) {
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
    Float_t valueSecond = valueMinute + second / 24. / 60. / 60.;
    pid.beamTime = valueSecond;
    if (isSimulation) pid.beamTime = 0;

    HTool::printProgress(i, nEntries, 1, "Analyze pairs :");
    Float_t VertexZ = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();
    if (VertexZ < -70) continue;

    HParticleEvtInfo * evtInfo = HCategoryManager::getObject(evtInfo, catParticleEvtInfo, 0);
    Int_t mult_tofrpc = evtInfo->getSumTofMultCut() + evtInfo->getSumRpcMultHitCut();
    Int_t cent_bin    = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);
    Int_t target_bin  = getTargetBin();
    if (cent_bin == 0) continue;  // exclude underflow bin


    pid.RPlanePhi = evtInfo->getRPlanePhi();
    pid.impPar    = 0;

    if (isSimulation) {
      fSimul     = ( (HRecEvent *) gHades->getCurrentEvent())->getPartialEvent(catSimul);
      fSubHeader = (HGeantHeader *) (fSimul->getSubHeader());
      if (fSubHeader) {
        pid.RPlanePhi = fSubHeader->getEventPlane();
        pid.impPar    = fSubHeader->getImpactParameter();
        cent_bin      = getCentBin(pid.impPar);
      }
    }


    if (cent_bin > 4) continue;
    cent_bin = cent_bin - 1; // starting from 0, not 1

    Counter_cent->Fill(cent_bin);
    Counter_cent->Fill(4);

    Int_t emc_ModMult[6][163] = { 0 };
    vector < pair < Int_t, Float_t > > candidates;
    sorter.cleanUp();
    sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE);    // reset all flags for flags (0-28) ,reject,used,lepton
    sorter.fill(HParticleTrackSorter::selectLeptons); // fill only good leptons
    sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
    sorter.fill(HParticleTrackSorter::selectHadrons); // fill only good hadrons (already marked good leptons will be skipped)
    sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);

    for (Int_t j = 0; j < candCat->getEntries(); j++) { // loop over Particlecand to see if there is a match
      HParticleCand * cand = HCategoryManager::getObject(cand, candCat, j);
			//HParticleCandSim * candSim = dynamic_cast < HParticleCandSim * > (cand);
      if (!cand->isFlagBit(kIsUsed)) continue;
      Int_t emc_index = cand->getEmcInd();
      Float_t emc_mq  = cand->getMetaMatchQualityEmc();
      if (emc_index >= 0) {
          HEmcCluster * cls = (HEmcCluster *) emcClusterCat->getObject(emc_index);
        Int_t sector      = cls->getSector();
        Int_t cell        = getSystemMap(cls->getCell());

        if (emc_mq < 2) emc_ModMult[sector][cell]++;
        else candidates.push_back(std::make_pair(emc_index, emc_mq));
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
      Float_t energy       = cls->getEnergy();
      Float_t matchedCells = cls->getNMatchedCells();
      Float_t qualityDTime = cls->getQualDTime();
      Int_t sector         = cls->getSector();
      Int_t cell = getSystemMap(cls->getCell());
      if (!isGoodCell(sector, cell)) continue;
      if (isSimulation) energy = correctEMCenergy(energy, sector, cell);
      Int_t matchedPart = emc_ModMult[sector][cell];
      Float_t time      = cls->getTime();
      HGeomVector trackVec(cls->getXLab(), cls->getYLab(), cls->getZLab() - VertexZ);
      Double_t trackLength = trackVec.length();
      Double_t beta        = (trackLength / 1000.) / (time * 1.e-9) / TMath::C();
      Int_t size = cls->getNCells();

      pid.qualityDTime = qualityDTime;

      if (energy < 100 || matchedPart > 0 || beta < 0.9 || beta > 1.5) continue;
			pid.emc_mq = 1000;
      // for (auto& [emc_index, emc_mq] : candidates) {
      //   if (emc_index == e) pid.emc_mq = emc_mq;
      // }

			for (Int_t i =0;i<candidates.size();i++) {
			  if (candidates[i].first == e) pid.emc_mq = candidates[i].second;
			}



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
  } // event loop
  out->cd();
  out->GetList()->Write();
  TMacro * m = new TMacro("fill_ttree_photon.C");
  TTree macros("macros", "macros");
  macros.GetUserInfo()->Add(m);
  macros.Write();
  out->Save();
  out->Close();
  return 0;
} /* fill_ttree_diphoton */
