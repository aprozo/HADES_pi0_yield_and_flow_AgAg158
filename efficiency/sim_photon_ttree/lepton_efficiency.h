static TH1F  *hcutMassSys0 = NULL;
static TH1F  *hcutMassSys1 = NULL;
static TH1F *corr1D_energy=  NULL;


static Bool_t selectLeptonsBetaLoose(HParticleCand* pcand){
    //  selection function for lepton candidates.
  Bool_t selectEpEm = kFALSE;
  if(pcand->isFlagAND(5,
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
     HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
     &&
     pcand->getChi2() < 1500
     &&
     pcand->getBeta() > 0.9
     &&
     pcand->getMomentum() > 0.
     &&
     pcand->getMomentum() < 3000.      ) selectEpEm = kTRUE;
  return selectEpEm;
}





void setupMassCuts(TString filename) {
    TFile *fileCutMass = new TFile(filename);
    if (fileCutMass) {
        hcutMassSys0 = (TH1F*)fileCutMass->Get("m2_cut_sys0");
        hcutMassSys1 = (TH1F*)fileCutMass->Get("m2_cut_sys1");
        if (!hcutMassSys0 || !hcutMassSys1) {
            cout << "No histogram with mass cuts found in file. Exitting..." << endl;
            exit(-1);
        }
    }
    else {
        cout << "No file with mass cuts found. Exitting..." << endl;
        exit(-1);
    }
}

void setupSimEffCorr(TString filename) {
    TFile *fileSimEffCorr= new TFile(filename);
    if (fileSimEffCorr) {
        corr1D_energy = (TH1F*) fileSimEffCorr->Get("energy_correction1D_Cell");;
        if (!corr1D_energy) {
            cout << "No histogram with sim eff correction found in file. Exitting..." << endl;
            exit(-1);
        }
    }
    else {
        cout << "No file with with sim eff correction found. Exitting..." << endl;
        exit(-1);
    }
}

Float_t correctEMCenergy(const Float_t & energy, const Int_t & cell) {  
     Int_t binX = corr1D_energy->GetXaxis()->FindBin(cell);
     Float_t corr = corr1D_energy->GetBinContent(binX);
     if(!TMath::Finite(corr)) corr=1.;
     Float_t emcEnCorr = (1+corr)*energy;
     return emcEnCorr;
}

Int_t leptonParticlecandCut( HParticleCandSim * cand1) {
  Int_t   sys         = cand1->getSystemUsed();
  Float_t richQa      = cand1->getRichMatchingQuality();
  Float_t beta        = cand1->getBeta();
  Float_t mom         = cand1->getMomentum();
  Float_t m2          = mom*mom*(1-beta*beta)/(beta*beta);
  Float_t charge      = cand1->getCharge();	    
  if(sys==1)return 0;
  if(!cand1->isFlagBit(Particle::kIsLepton))  return 0;
  if (richQa < 0.|| richQa > 1.5 ||cand1->getRingRadius() < 17) return 0;
  if (beta>0.9&&sys==0 && !(m2 < hcutMassSys0->Interpolate(mom*charge))) return 0;
  return 1;
}

