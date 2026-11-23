#ifndef __HISTOGRAM_MANAGER__
#define __HISTOGRAM_MANAGER__

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TH2.h"
#include "TBits.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "hades.h"
#include "TVector3.h"
#include "TProfile.h"
#include "../GlobalParameters.h"

typedef struct {
  Float_t newEvent, multTofRpc, multEmc, multPhoton, centBin,
    targetBin, beamTime, RPlanePhi, impPar, emc_mq, qualityDTime, betaCand,
    sector, cell, size, beta, energy, time, xpos, ypos, zpos, track, trackParent, mom_kine, idParent;
} PID;

// =======================================================================//
TH1F * photon_eff[3];
Float_t getCorrectionEmbedding(const PID & pid)
{
  Float_t energy = pid.energy;

  if (energy >= 790) energy = 790;
  Float_t corr = photon_eff[(int) pid.centBin]->GetBinContent(photon_eff[(int) pid.centBin]->FindBin(energy));
  // Float_t corrError = photon_eff[(int) pid.centBin]->GetBinError(photon_eff[(int) pid.centBin]->FindBin(energy));
  return (corr > 0 && corr < 1) ? corr : 1;
}

Float_t getCorrectionErrorEmbedding(const PID & pid)
{
  Float_t energy = pid.energy;

  if (energy >= 790) energy = 790;
  Float_t corr      = photon_eff[(int) pid.centBin]->GetBinContent(photon_eff[(int) pid.centBin]->FindBin(energy));
  Float_t corrError = photon_eff[(int) pid.centBin]->GetBinError(photon_eff[(int) pid.centBin]->FindBin(energy));
  return (corr > 0 && corr < 1) ? corrError : 0;
}

// =======================================================================//
TH1F * betaCutEff[3];
Float_t getBetaCutEff(const PID & pid)
{
  Float_t energy = pid.energy;

  if (energy >= 990) energy = 990;
  Float_t corr = betaCutEff[(int) pid.centBin]->GetBinContent(betaCutEff[(int) pid.centBin]->FindBin(energy));
  return (corr > 0 && corr < 1) ? corr : 1;
}

// =======================================================================//

TH1F * accEffEmbedding[3];
Float_t getAcceptanceEfficiency(const Float_t & _y, const Float_t & _pt, const Int_t & centBin)
{
  Float_t y = _y;

  if (y <= 0.15) y = 0.151;
  if (y >= 0.95) y = 0.949;
  Float_t pt = _pt;
  if (pt >= 1200) pt = 1199;
  Float_t corr = accEffEmbedding[centBin]->Interpolate(y, pt);
  return (corr > 0 && corr < 1) ? corr : 1;
}

struct Photon {
  Photon(){ };
  Photon(const PID & pid) : track(0), trackParent(0){
    TVector3 trackVec(pid.xpos, pid.ypos, pid.zpos);
    Double_t trackLength = trackVec.Mag();
    Float_t energy       = pid.energy;
    RPlanePhi = pid.RPlanePhi;
    trackVec *= (energy / trackLength);
    lorVec.SetXYZM(trackVec.x(), trackVec.y(), trackVec.z(), 0);
    centBin = (Int_t) pid.centBin;
    //if (!isSim) {
  //    eff        = getCorrectionEmbedding(pid) * getBetaCutEff(pid);
  //    weight     = (1. / eff);
  //    weight_err = getCorrectionErrorEmbedding(pid) / (eff * eff);
    //} else {
       weight = 1.0; weight_err = 0.0;
     //}

    if (isSim) { // for SIM
      track       = pid.track;
      trackParent = pid.trackParent;
      mom_kine    = pid.mom_kine;
      idParent    = pid.idParent;
    }
  }
  TLorentzVector lorVec;
  Int_t          centBin;
  Float_t        RPlanePhi;
  Int_t          track, trackParent, idParent;
  Float_t        eff, mom_kine, weight, weight_err;
};


struct Diphoton {
  Diphoton(const Photon & gamma1, const Photon & gamma2) : fromPi0(0)
  {
    TLorentzVector pairMix = gamma1.lorVec + gamma2.lorVec;

    opangle = gamma1.lorVec.Angle(gamma2.lorVec.Vect()) * TMath::RadToDeg();
    mass    = pairMix.M();
    mom     = pairMix.P();
    theta   = pairMix.Theta() * TMath::RadToDeg();
    phi     = pairMix.Phi() * TMath::RadToDeg() - 0.5 * (gamma2.RPlanePhi + gamma1.RPlanePhi);
    while (phi > 180.) phi -= 360.;
    while (phi <= -180.) phi += 360.;
    y       = pairMix.Rapidity() - 0.822;// correct to midrapidity
    pt      = pairMix.Pt();
    centBin = gamma1.centBin;
    weight  = gamma1.weight * gamma2.weight;
    accEffCorrection = 1.;// / getAcceptanceEfficiency(y, pt, centBin);
    weight_err2      = TMath::Power(gamma1.weight * gamma2.weight_err, 2) + TMath::Power(
      gamma1.weight_err * gamma2.weight, 2);


    fromPi0 = gamma1.track > 0 && gamma2.track > 0 && (gamma1.track != gamma2.track) &&
      (gamma1.trackParent == gamma2.trackParent) && gamma1.trackParent > 0 && gamma1.idParent == 7;
  }

  Float_t mass, opangle, pt, y, theta, phi, mom, weight, weight_err2, accEffCorrection;
  Int_t   centBin;
  Bool_t  fromPi0;
};


class HistogramManager {
private:

  Int_t ptBins;
  Int_t ptMin;
  Int_t ptMax;
  Float_t ptBinSize;

  Int_t yBins;
  Float_t yMin;
  Float_t yMax;
  Float_t yBinSize;

  Int_t ptBinsFlow;
  Int_t ptMinFlow;
  Int_t ptMaxFlow;
  Float_t ptBinSizeFlow;

  Int_t yBinsFlow;
  Float_t yMinFlow;
  Float_t yMaxFlow;
  Float_t yBinSizeFlow;

  Int_t MassBins;
  Float_t MassMin;
  Float_t MassMax;

  Int_t centMin;
  Int_t centMax;
  Int_t centBins;


  const Int_t phiMin = -180;
  const Int_t phiMax = 180;
  Int_t phiBins;
  Float_t phiBinSize;

  const Int_t quadPhiMin       = -45;
  const Int_t quadPhiMax       = 315;
  const Int_t quadPhiBins      = 4;
  const Float_t quadPhiBinSize = 90;

  TH1F ***** massPhi;// [phiBins][ptBins][yBins][centBins];
  TH1F ***** massPhiPi0;// [phiBins][ptBins][yBins][centBins];

  TH1F ***** massQuadPhi;// [phiBins][ptBins][yBins][centBins];

  TProfile **** profCosDphi;// [ptBins][yBins][centBins];
  TProfile **** profCos2Dphi;// [ptBins][yBins][centBins];


  TH1F **** massSum;

  Float_t startBeamTime;
  Float_t endBeamTime;

  TH1D * massIntegrated[4];


public:
  HistogramManager(TDirectory * dir, const GlobalParameters & par);
  // virtual ~HistogramManager();
  void Hfill(const Diphoton & diphot);
  void HFill(const Photon & photon);
  void Write(TDirectory * dir);
  void WriteMix(TDirectory * dir);
  // ClassDef(HistogramManager,1) //Manages all histograms
};

// ______________________________________________________________________________
HistogramManager::HistogramManager(TDirectory * dir, const GlobalParameters & par)
{
  // Create histogram manager object. Histograms will be created
  // in the "dir" directory.

  // Save current directory and cd to "dir".

  TDirectory * saved = gDirectory;

  dir->cd();
  MassBins = par.N_MASS_BINS;
  MassMin  = 0;
  MassMax  = 600;

  centMin  = 0;
  centMax  = 30;
  centBins = 3;

  ptBins    = par.N_PT_BINS;
  ptMin     = 0;
  ptMax     = 1200;
  ptBinSize = (ptMax - ptMin) / Float_t(ptBins);

  yBins    = par.N_Y_BINS;
  yMin     = 0.05;
  yMax     = 1.05;
  yBinSize = (yMax - yMin) / Float_t(yBins);


  phiBins    = par.N_PHI_BINS;
  phiBinSize = 360 / Float_t(phiBins);

  ptBinsFlow    = par.N_PT_FLOW_BINS;
  ptMinFlow     = 100;
  ptMaxFlow     = 900;
  ptBinSizeFlow = (ptMaxFlow - ptMinFlow) / Float_t(ptBinsFlow);

  yBinsFlow    = par.N_Y_FLOW_BINS;
  yMinFlow     = 0.1;
  yMaxFlow     = 0.9;
  yBinSizeFlow = (yMaxFlow - yMinFlow) / Float_t(yBinsFlow);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==================================== DIFFERENTIAL HISTOGRAMS======================================================//
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TString xTitle = ";M_{#gamma#gamma} [MeV/c^2]", yTitle = ";Counts";


  for (Int_t cent = 0; cent < centBins; cent++) {
    TString centName = Form("%i0%%-%i0%% ", cent + centMin, (cent + 1) + centMin);
    massIntegrated[cent] = new TH1D(Form("cent%i", cent), "Integrated " + centName + xTitle + yTitle, MassBins, MassMin,
        MassMax);
  }

  // ================================================================================================================//
  // ====================================================FLOW========================================================//
  // ================================================================================================================//

  massPhi      = new TH1F * ** *[centBins];// [phiBins][ptBins][yBins];
  massPhiPi0   = new TH1F * ** *[centBins];// [phiBins][ptBins][yBins];
  massSum      = new TH1F * * *[centBins];// [ptBins][yBins];
  massQuadPhi  = new TH1F * ** *[centBins];// [phiBins][ptBins][yBins];
  profCosDphi  = new TProfile * * *[centBins];// [ptBins][yBins];
  profCos2Dphi = new TProfile * * *[centBins]; // [ptBins][yBins];

  for (Int_t cent = 0; cent < centBins; cent++) {
    massPhi[cent]      = new TH1F * * *[ptBinsFlow];// [phiBins][yBins];
    massSum[cent]      = new TH1F * *[ptBinsFlow];// [yBins];
    massPhiPi0[cent]   = new TH1F * * *[ptBinsFlow];// [phiBins]yBins];
    massQuadPhi[cent]  = new TH1F * * *[ptBinsFlow];// [phiBins][yBins];
    profCosDphi[cent]  = new TProfile * *[ptBinsFlow];// yBins];
    profCos2Dphi[cent] = new TProfile * *[ptBinsFlow]; // yBins];

    for (Int_t pt = 0; pt < ptBinsFlow; pt++) {
      massPhi[cent][pt]      = new TH1F * *[yBinsFlow];// [phiBins];
      massPhiPi0[cent][pt]   = new TH1F * *[yBinsFlow];// [phiBins];
      massSum[cent][pt]      = new TH1F *[yBinsFlow];// [yBins];
      massQuadPhi[cent][pt]  = new TH1F * *[yBinsFlow];// [phiBins];
      profCosDphi[cent][pt]  = new TProfile *[yBinsFlow];
      profCos2Dphi[cent][pt] = new TProfile *[yBinsFlow];

      for (Int_t y = 0; y < yBinsFlow; y++) {
        TString centName = Form("%i0%%-%i0%% ", cent + centMin, (cent + 1) + centMin);
        TString yName    =
          Form("| %.2f < y_cm < %.2f |", (y * yBinSizeFlow + yMinFlow), ((y + 1) * yBinSizeFlow + yMinFlow));
        TString ptName =
          Form("| %.0f < pt < %.0f |", (pt * ptBinSizeFlow + ptMinFlow), ((pt + 1) * ptBinSizeFlow + ptMinFlow));
        TString hName = Form("pt%iy%icent%i", pt, y, cent);

        massPhi[cent][pt][y]     = new TH1F *[phiBins];
        massPhiPi0[cent][pt][y]  = new TH1F *[phiBins];
        massQuadPhi[cent][pt][y] = new TH1F *[4];


        massSum[cent][pt][y] = new TH1F(hName, ptName + yName + centName + xTitle + yTitle, MassBins, MassMin, MassMax);

        for (Int_t phi = 0; phi < phiBins; phi++) {
          TString phiName = Form("| %.0f < phi < %.0f |", phi * phiBinSize + phiMin, (phi + 1) * phiBinSize + phiMin);
          massPhi[cent][pt][y][phi] = new TH1F(Form("phi%i", phi) + hName,
              phiName + ptName + yName + centName + xTitle + yTitle, MassBins,
              MassMin, MassMax);
          massPhiPi0[cent][pt][y][phi] = new TH1F(Form("phi%i", phi) + hName + "geant",
              phiName + ptName + yName + centName + xTitle + yTitle, MassBins,
              MassMin, MassMax);
        }
        // ==================================== QUADRANT METHOD======================================================//
        for (Int_t quadPhi = 0; quadPhi < 4; quadPhi++) {
          TString quadPhiName = Form("| %i < quad phi < %i |", (quadPhi * 90 - 45), (quadPhi + 1) * 90 - 45);
          massQuadPhi[cent][pt][y][quadPhi] = new TH1F(Form("qphi%i", quadPhi) + hName,
              quadPhiName + ptName + yName + centName + xTitle + yTitle,
              MassBins, MassMin, MassMax);
        }

        // ==================================== TPROFILE METHOD======================================================//
        profCosDphi[cent][pt][y] = new TProfile("cosDphi" + hName,
            "cosDphi " + ptName + yName + centName + xTitle + ";cos(#Delta#phi)",
            MassBins, MassMin, MassMax, -1, 1);
        profCos2Dphi[cent][pt][y] = new TProfile("cos2Dphi" + hName,
            "cos2Dphi " + ptName + yName + centName + xTitle + ";cos(2#Delta#phi)",
            MassBins, MassMin, MassMax, -1, 1);
      }
    }
  }

  // ==================================================================================================================================//
  // ==================================================================================================================================//


  // cd back to original directory
  saved->cd();
}

void HistogramManager::Hfill(const Diphoton & diphot)
{ // Float_t  assymetry = TMath::Abs((pid.energy1-pid.energy2))/(pid.energy1+pid.energy2);
  Float_t invmass     = diphot.mass;
  Float_t weight      = diphot.weight;
  Float_t weight_err2 = diphot.weight_err2;
  Float_t error       = 0;
  Int_t centBin       = diphot.centBin;


  // ==================================================================//


  Int_t ptBin = floor((diphot.pt - ptMinFlow) / ptBinSizeFlow);
  Int_t yBin  = floor((diphot.y + 0.0000000001 - yMinFlow) / yBinSizeFlow);

  // mapping iQuarter number to correct bin number:
  // iQuarter = 0:  -45 <= Dphi < +45                         ... iBinQ = 2
  // iQuarter = 1: -135 <= Dphi < -45                         ... iBinQ = 1 (bin number starting from 1, zero is for underflow)
  // iQuarter = 2: -180 <= Dphi < -135 and 135 <= Dphi < 180  ... iBinQ = 4 (due to selected range from -135 to 225)
  // iQuarter = 3:  +45 <= Dphi < +135                        ... iBinQ = 3
  Int_t quadPhiBin = 2;
  Double_t phi     = diphot.phi;

  if (-45 <= phi && phi < 45) quadPhiBin = 0;
  if (45 <= phi && phi < 135) quadPhiBin = 1;
  if (-135 <= phi && phi < -45) quadPhiBin = 3;


  Int_t phiBin     = floor((diphot.phi - phiMin) / phiBinSize);
  Double_t cosPhi  = TMath::Cos(phi * TMath::DegToRad());
  Double_t cos2Phi = TMath::Cos(2 * phi * TMath::DegToRad());


  if (centBin >= 0 && centBin < centBins) {
    Int_t binGGmass = massIntegrated[centBin]->FindBin(invmass);
    if (invmass > MassMin && invmass < MassMax) error = massIntegrated[centBin]->GetBinError(binGGmass);
    error *= error;
    massIntegrated[centBin]->Fill(invmass, weight);
    massIntegrated[centBin]->SetBinError(binGGmass, TMath::Sqrt(error + weight_err2));
  }


  if (centBin >= 0 && centBin < centBins &&
    yBin >= 0 && yBin < yBinsFlow &&
    ptBin >= 0 && ptBin < ptBinsFlow)
  {
    if (phiBin >= 0 && phiBin < phiBins) {
      massPhi[centBin][ptBin][yBin][phiBin]->Fill(invmass, diphot.accEffCorrection);
      if (diphot.fromPi0) massPhiPi0[centBin][ptBin][yBin][phiBin]->Fill(invmass, diphot.accEffCorrection);
    }

    massSum[centBin][ptBin][yBin]->Fill(invmass, diphot.accEffCorrection);

    massQuadPhi[centBin][ptBin][yBin][quadPhiBin]->Fill(invmass, diphot.accEffCorrection);

    profCosDphi[centBin][ptBin][yBin]->Fill(invmass, cosPhi, diphot.accEffCorrection);
    profCos2Dphi[centBin][ptBin][yBin]->Fill(invmass, cos2Phi, diphot.accEffCorrection);
  }

  // ==================================================================//
} /* HistogramManager::Hfill */

void HistogramManager::Write(TDirectory * dir)
{
  TDirectory * saved = gDirectory;

  dir->cd();
  for (Int_t cent = 0; cent < centBins; cent++) {
    massIntegrated[cent]->Write();
  }

  // ==================================================================//
  for (Int_t cent = 0; cent < centBins; cent++) {
    for (Int_t y = 0; y < yBinsFlow; y++) {
      for (Int_t pt = 0; pt < ptBinsFlow; pt++) {
        massSum[cent][pt][y]->Write();
        for (Int_t phi = 0; phi < phiBins; phi++) {
          massPhi[cent][pt][y][phi]->Write();
          if (isSim) massPhiPi0[cent][pt][y][phi]->Write();
        }
        for (Int_t quadPhi = 0; quadPhi < 4; quadPhi++) {
          massQuadPhi[cent][pt][y][quadPhi]->Write();
        }
        profCosDphi[cent][pt][y]->Write();
        profCos2Dphi[cent][pt][y]->Write();
      }
    }
  }

  // ==================================================================//


  saved->cd();
} /* HistogramManager::Write */

void HistogramManager::WriteMix(TDirectory * dir)
{
  TDirectory * saved = gDirectory;

  dir->cd();
  for (Int_t cent = 0; cent < centBins; cent++) {
    massIntegrated[cent]->Write();
  }


  // ==================================================================//

  for (Int_t cent = 0; cent < centBins; cent++) {
    for (Int_t y = 0; y < yBinsFlow; y++) {
      for (Int_t pt = 0; pt < ptBinsFlow; pt++) {
        massSum[cent][pt][y]->Write();
        for (Int_t phi = 0; phi < phiBins; phi++) {
          massPhi[cent][pt][y][phi]->Write();
        }
        for (Int_t quadPhi = 0; quadPhi < 4; quadPhi++) {
          massQuadPhi[cent][pt][y][quadPhi]->Write();
        }
        // profCosDphi[cent][pt][y]->Write();
        // profCos2Dphi[cent][pt][y]->Write();
      }
    }
  }

  // ==================================================================//


  saved->cd();
} /* HistogramManager::Write */

#endif /* ifndef __HISTOGRAM_MANAGER__ */
