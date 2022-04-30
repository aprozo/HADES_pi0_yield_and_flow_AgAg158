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
#include "Useful_Stuff.h"


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

using namespace std;

static TH1F * hcutMassSys0 = NULL;
static TH1F * hcutMassSys1 = NULL;

static Bool_t selectLeptonsBetaLoose(HParticleCand * pcand)
{
  //  selection function for lepton candidates.

  Bool_t selectEpEm = kFALSE;


  if (pcand->isFlagAND(5,
    Particle::kIsAcceptedHitRICH,
    Particle::kIsAcceptedHitInnerMDC,
    Particle::kIsAcceptedHitOuterMDC,
    Particle::kIsAcceptedHitMETA,
    Particle::kIsAcceptedRK)
    &&
    pcand->getInnerSegmentChi2() > 0
    &&
    pcand->getOuterSegmentChi2() > 0
    &&
    pcand->getMetaMatchQuality() < 4
    &&
    pcand->getMetaMatchQualityEmc() < 4
    &&
    HParticleTool::isGoodMetaCell(pcand, 4, kTRUE)
    &&
    pcand->getChi2() < 1000
    &&
    pcand->getBeta() > 0.9
    &&
    pcand->getMomentum() > 0.
    &&
    pcand->getMomentum() < 1500.) selectEpEm = kTRUE;

  return selectEpEm;
} // selectLeptonsBetaLoose

void setupMassCuts(TString filename)
{
  TFile * fileCutMass = new TFile(filename);

  if (fileCutMass) {
    hcutMassSys0 = (TH1F *) fileCutMass->Get("m2_cut_sys0");
    hcutMassSys1 = (TH1F *) fileCutMass->Get("m2_cut_sys1");
    if (!hcutMassSys0 || !hcutMassSys1) {
      cout << "No histogram with mass cuts found in file. Exitting..." << endl;
      exit(-1);
    }
  } else {
    cout << "No file with mass cuts found. Exitting..." << endl;
    exit(-1);
  }
}

Int_t analysis(TString infileList, TString outfile, Int_t nEvents = 1000000)
{
  Bool_t lowfield = kFALSE;// kFALSE;


  //  infileList : comma seprated file list "file1.root,file2.root" or "something*.root"
  //  outfile    : optional (not used here) , used to store hists in root file
  //  nEvents    : number of events to processed. if  nEvents < entries or < 0 the chain will be processed

  // -------------------------------------------------
  // create loop obejct and hades
  HLoop loop(kTRUE);

  // -------------------------------------------------

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

  // -------------------------------------------------
  // list of all files with working sectors
  // if(!isSimulation) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);
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


  if (!loop.setInput()) {
    cout << "READBACK: ERROR : cannot read input !" << endl;
    exit(1);
  } // read all categories
  loop.printCategories();
  loop.printChain();

  // -------------------------------------------------

  // -------------------------------------------------
  HParticleTrackSorter sorter;
  sorter.setIgnoreInnerMDC(); // do not reject Double_t inner MDC hits
  sorter.init();              // get catgegory pointers etc...
  // -------------------------------------------------

  // -------------------------------------------------
  // input data
  HEventHeader * evHeader   = loop.getEventHeader();
  HCategory * candCat       = (HCategory *) HCategoryManager::getCategory(catParticleCand);
  HCategory * emcClusterCat = (HCategory *) HCategoryManager::getCategory(catEmcCluster);
  HCategory * emcCalCat     = (HCategory *) HCategoryManager::getCategory(catEmcCal);
  HCategory * evtInfoCat    = (HCategory *) HCategoryManager::getCategory(catParticleEvtInfo);
  HCategory * start2HitCat  = (HCategory *) HCategoryManager::getCategory(catStart2Hit);
  HCategory * richHitCat    = (HCategory *) HCategoryManager::getCategory(catRichHit);
  HCategory * richCalCat    = (HCategory *) HCategoryManager::getCategory(catRichCal);

  HIterator * iterEmcCluster = (HIterator *) emcClusterCat->MakeIterator("native");
  HIterator * iterEmcCal     = (HIterator *) emcCalCat->MakeIterator("native");
  // HIterator *iterEmcRaw = (HIterator*)emcRawCat->MakeIterator("native");
  HIterator * iterStart = (HIterator *) start2HitCat->MakeIterator("native");


  // -------------------------------------------------
  // added for the task
  HTaskSet * masterTaskSet = gHades->getTaskSet("all");
  HMetaMatcher * matcher   = new HMetaMatcher();
  // matcher->setDebug();
  masterTaskSet->add(matcher);
  // -------------------------------------------------


  // -------------------------------------------------
  setupMassCuts("m2_cut.root");
  // setting up the output
  TFile * out = new TFile(outfile, "RECREATE");
  out->cd();


  typedef struct {
    Float_t sector, cell, cell2, tof, theta, phi, richQa,
      emc_mq, mom, beta, charge, mm_qa, energy, maxenergy, time, time2, size, emc_theta,
      emc_phi, emc_path, distRpcEmc, timeCorr, beamtime, TOT, TOT2,
      dthetadphi_qa, dtime, chi2, raw_time, raw_time2, goodrich, sys;
  } PID;
  PID pid;

  TTree * nt = new TTree("nt", "pid");
  nt->Branch("pid", &pid, "sector:cell:cell2:tof:"
    "theta:phi:richQa:emc_mq:mom:"
    "beta:charge:mm_qa:energy:maxenergy:"
    "time:time2:size:emc_theta:emc_phi:emc_path:distRpcEmc:"
    "timeCorr:beamtime:TOT:TOT2:dthetadphi_qa:dtime:"
    "chi2:raw_time:raw_time2:goodrich:sys");


  HLocation loc; // location of cluster
  loc.set(2, 0, 0);

  Int_t nEntries = loop.getEntries();
  printf("nEntries = %i\n", nEntries);
  if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

  TString tempFileName, currentFileName;
  Int_t dayOfYear, hour, min, second;

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
    Float_t valueSecond = valueMinute + second / 24. / 60. / 60.;
    pid.beamtime = valueSecond;

    // -------------------------------------------------
    // summary event info object
    HParticleEvtInfo * evtInfo = 0;
    evtInfo = HCategoryManager::getObject(evtInfo, evtInfoCat, 0);

    // clean vectors and index arrays
    sorter.cleanUp();
    // if you do not want to modify the flags
    // for the following tasks etc (needs call of
    // restoreFlags() at the end of the eventloop)
    // sorter.backupFlags();
    // ------------------------------------------------------------------------
    sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
    sorter.fill(selectLeptonsBetaLoose);           // fill only good loose leptons

    Int_t nCandLepBest = sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);


    Float_t VertexZ = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();
    if (VertexZ < -70) continue;
    Float_t start_time = 0.0;
    HStart2Hit * start = 0;
    iterStart->Reset();
    while ((start = (HStart2Hit *) iterStart->Next()) != 0) {
      start_time = start->getTime();
    }

    Int_t size = candCat->getEntries();
    HParticleCand * cand = 0;
    for (Int_t j = 0; j < size; j++) { // loop over Particlecand
      cand = HCategoryManager::getObject(cand, candCat, j);

      if (!cand->isFlagBit(Particle::kIsUsed)) continue;
      if (!cand->isFlagBit(Particle::kIsLepton)) continue;
      if (cand->getEmcInd() < 0) continue;
      // if (cand->isTofClstUsed() ==0) continue;

      Float_t theta  = cand->getTheta();
      Float_t phi    = cand->getPhi();
      Float_t beta   = cand->getBeta();
      Float_t mom    = cand->getMomentum();
      Float_t m2     = mom * mom * (1 - beta * beta) / (beta * beta);
      Float_t charge = cand->getCharge();
      Int_t sys      = cand->getSystemUsed();
      Int_t sec      = cand->getSector();
      Float_t richQa = cand->getRichMatchingQuality();
      Float_t richR  = cand->getRingRadius();

      Float_t emc_energy = cand->getEmcEnergy();
      Float_t emc_time   = cand->getEmcTime();

      Float_t emc_mq          = cand->getMetaMatchQualityEmc();
      Float_t MetaMatchQA     = cand->getMetaMatchQuality();
      Float_t candpart_length = cand->getDistanceToMetaHit();
      Float_t tof      = cand->getTof();
      Float_t emc_path = cand->getEmcCellPath();

      Bool_t oaRingFlag     = kFALSE;
      Bool_t conversionFlag = kFALSE;
      if (isOARing(cand, richHitCat)) oaRingFlag = kTRUE;
      if (isConversion(cand, richHitCat, richCalCat)) conversionFlag = kTRUE;
      pid.goodrich = (!oaRingFlag && !conversionFlag) ? 1 : 0;
      // Ring matching
      // rich condition for leptons
      if (richQa < 0. || richQa > 1.5 || richR < 17) continue;
      // beta vs. m2 condition for leptons
      if (sys == 1 && !(m2 < hcutMassSys1->Interpolate(mom * charge))) continue;
      if (sys == 0 && !(m2 < hcutMassSys0->Interpolate(mom * charge))) continue;
      pid.sys      = sys;
      pid.chi2     = cand->getChi2();
      pid.richQa   = richQa;
      pid.sector   = sec;
      pid.tof      = tof;
      pid.theta    = theta;
      pid.phi      = phi;
      pid.mom      = mom;
      pid.beta     = beta;
      pid.charge   = charge;
      pid.mm_qa    = MetaMatchQA;
      pid.emc_path = emc_path;
      pid.TOT      = -1;
      pid.TOT2     = -1;


      HEmcCluster * emcClst = (HEmcCluster *) emcClusterCat->getObject(cand->getEmcInd());
      Int_t sector      = emcClst->getSector();
      Int_t clusterSize = emcClst->getNCells();
      if (clusterSize > 2 || emc_mq > 2.5) continue;
      UChar_t c1 = emcClst->getCell(0);
      loc[0] = sector; // set location of cluster
      loc[1] = c1;     // set location of cluster
      HEmcCal * cal1 = (HEmcCal *) emcCalCat->getObject(loc);
      pid.TOT      = cal1->getWidth();
      pid.raw_time = cal1->getLead() - start_time;

      Float_t timeCorr   = 0.;
      Float_t distRpcEmc = -100;

      if (sys == 1 && cand->isTofClstUsed()) {
        HTofCluster * tofCluster = matcher->matchWithTof(cand);
        Float_t xrl, yrl, zrl;
        tofCluster->getXYZLab(xrl, yrl, zrl);
        Float_t deltaX = emcClst->getXLab() - xrl;
        Float_t deltaY = emcClst->getYLab() - yrl;
        Float_t deltaZ = emcClst->getZLab() - zrl;
        distRpcEmc = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        timeCorr   = distRpcEmc / TMath::Sqrt(xrl * xrl + yrl * yrl + zrl * zrl) * tofCluster->getTof();
      }


      if (sys == 0 && cand->isRpcClstUsed()) {
        HRpcCluster * fakerpccluster = matcher->matchWithRpc(cand);
        Float_t xrl, yrl, zrl;
        fakerpccluster->getXYZLab(xrl, yrl, zrl);
        Float_t deltaX = emcClst->getXLab() - xrl;
        Float_t deltaY = emcClst->getYLab() - yrl;
        Float_t deltaZ = emcClst->getZLab() - zrl;
        distRpcEmc = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        timeCorr   = distRpcEmc / TMath::Sqrt(xrl * xrl + yrl * yrl + zrl * zrl) * fakerpccluster->getTof();
      }


      pid.distRpcEmc = distRpcEmc;
      pid.timeCorr   = timeCorr;
      pid.cell2      = -1;
      pid.time2      = -100;
      pid.raw_time2  = -100;
      if (clusterSize == 2) { // do the same for second cal in cluster
        UChar_t c2 = emcClst->getCell(1);
        loc[1] = c2;
        HEmcCal * cal2 = (HEmcCal *) emcCalCat->getObject(loc);
        pid.TOT2      = cal2->getWidth();
        pid.raw_time2 = cal2->getLead() - start_time;
        pid.cell2     = getSystemMap(emcClst->getCell(1));
        pid.time2     = cal2->getTime();
      }

      pid.dtime         = emcClst->getQualDTime();
      pid.dthetadphi_qa = emcClst->getQualDThDPh();
      pid.cell = getSystemMap(emcClst->getCell(0));

      pid.emc_theta = emcClst->getTheta();
      pid.emc_phi   = emcClst->getPhi();
      pid.emc_mq    = emc_mq;
      pid.energy    = emc_energy;
      pid.maxenergy = emcClst->getMaxEnergy();
      pid.time      = emcClst->getTime();
      pid.size      = clusterSize;

      nt->Fill();
    } // end of loop over ParticleCandidate*/
  }   // loop over events

  out->cd();
  out->GetList()->Write();

  // TMacro * m = new TMacro("lepton_ttree.C");
  // TTree macros("lepton_ttree", "lepton_ttree");
  // macros.GetUserInfo()->Add(m);
  // macros.Write();
  out->Save();
  out->Close();

  return 0;
} // emc_modern
