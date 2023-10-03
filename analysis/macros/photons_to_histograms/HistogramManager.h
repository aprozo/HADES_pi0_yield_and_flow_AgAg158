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
	Photon(){
	};
	Photon(const PID & pid) : track(0), trackParent(0){
		TVector3 trackVec(pid.xpos, pid.ypos, pid.zpos);
		Double_t trackLength = trackVec.Mag();
		Float_t energy       = pid.energy;
		RPlanePhi = pid.RPlanePhi;
		trackVec *= (energy / trackLength);
		lorVec.SetXYZM(trackVec.x(), trackVec.y(), trackVec.z(), 0);
		centBin = (Int_t) pid.centBin;
		if (!isSim) {
			eff        = getCorrectionEmbedding(pid) * getBetaCutEff(pid);
			weight     = (1. / eff);
			weight_err = getCorrectionErrorEmbedding(pid) / (eff * eff);
		} else { weight = 1.0; weight_err = 0.0; }

		if (isSim) { // for SIM
			track       = pid.track;
			trackParent = pid.trackParent;
			mom_kine    = pid.mom_kine;
			idParent    = pid.idParent;
		}
		sector=pid.sector;
	}
	TLorentzVector lorVec;
	Int_t centBin;
	Int_t sector;
	Float_t RPlanePhi;
	Int_t track, trackParent, idParent;
	Float_t eff, mom_kine, weight, weight_err;
};


struct Diphoton {
	Diphoton(const Photon & gamma1, const Photon & gamma2) : isPi0(0),isEta(0)
	{
		TLorentzVector pairMix = gamma1.lorVec + gamma2.lorVec;

		opangle = gamma1.lorVec.Angle(gamma2.lorVec.Vect()) * TMath::RadToDeg();
		mass    = pairMix.M();
		mom     = pairMix.P();
		theta   = pairMix.Theta() * TMath::RadToDeg();
		phi     = pairMix.Phi() * TMath::RadToDeg() - gamma2.RPlanePhi;
		while (phi > 180.) phi -= 360.;
		while (phi <= -180.) phi += 360.;
		y       = pairMix.Rapidity() - 0.822;// correct to midrapidity
		pt      = pairMix.Pt();
		centBin = gamma1.centBin;
		weight  = gamma1.weight * gamma2.weight;
		accEffCorrection = 1. / getAcceptanceEfficiency(y, pt, centBin);
		weight_err2      = TMath::Power(gamma1.weight * gamma2.weight_err, 2) + TMath::Power(
			gamma1.weight_err * gamma2.weight, 2);


		isPi0 = gamma1.track > 0 && gamma2.track > 0 && (gamma1.track != gamma2.track) &&
		        (gamma1.trackParent == gamma2.trackParent) && gamma1.trackParent > 0 && gamma1.idParent == 7;

		isEta = gamma1.track > 0 && gamma2.track > 0 && (gamma1.track != gamma2.track) &&
		        (gamma1.trackParent == gamma2.trackParent) && gamma1.trackParent > 0 && gamma1.idParent == 17;

		isSameSector=(gamma1.sector==gamma2.sector);
	}

	Float_t mass, opangle, pt, y, theta, phi, mom, weight, weight_err2, accEffCorrection;
	Int_t centBin;
	Bool_t isPi0,isEta;
	Bool_t isSameSector;
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


	Int_t MassBins;
	Float_t MassMin;
	Float_t MassMax;

	Int_t centMin;
	Int_t centMax;
	Int_t centBins;

	Float_t opAngleCut;

	// =========================================================//
	TH1F **** mass;// [ptBins][yBins][centBins];
	TH1F **** massPi0;// [ptBins][yBins][centBins];


	// quality HISTOGRAMS
	TH1F * qualityOpAngleGeant;
	TH1F * qualityOpAngleGeantEta;
	TH1F * qualityMassGeant;
	TH1F * qualityMass;
	TH2F * qualityPtYSpace;
	TH2F * qualityPtMass;
	TH2F * qualityYMass;
	TH2F * qualityOpAngleMass;
	TH2F * qualityOpAngleMassGeant;
	TH2F * qualityOpAnglePt;
	TH2F * qualityOpAnglePtEta;

	TH2F * qualityOpAnglePtOneSector;
	TH2F * qualityOpAnglePtEtaOneSector;

	TH2F * qualityOpAnglePtDiffSector;
	TH2F * qualityOpAnglePtEtaDiffSector;



	TH2F * qualityPtMassOneSector;
	TH2F * qualityYMassOneSector;

	TH2F * qualityPtMassDiffSector;
	TH2F * qualityYMassDiffSector;


	// =============================================================//

	Float_t startBeamTime;
	Float_t endBeamTime;

	TH1D * massIntegrated[4];


public:
	HistogramManager(TDirectory * dir, const GlobalParameters & par);
	virtual ~HistogramManager();
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


	MassBins = par.N_MASS_BINS;
	MassMin  = 0;
	MassMax  = 800;

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

	opAngleCut = par.ANGLE_CUT;


	TDirectory * saved = gDirectory;
	dir->cd();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ==================================== DIFFERENTIAL HISTOGRAMS======================================================//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	TString massTitle = ";M_{#gamma#gamma} [MeV/c^2]";
	TString ptTitle = ";p_{t} [MeV/c]", rapTitle = ";y_{cm}", opAngleTitle = ";#Theta_{#gamma#gamma},[#circ]";
	qualityOpAngleGeant = new TH1F("qualityOpAngleGeant", opAngleTitle + ";Counts", 240, 0, 120);
	qualityOpAngleGeantEta = new TH1F("qualityOpAngleGeantEta", opAngleTitle + ";Counts", 240, 0, 120);
	qualityOpAnglePt = new TH2F("qualityOpAnglePt", ptTitle + opAngleTitle + ";Counts",  300, 0, 1200, 120, 0, 120);
	qualityOpAnglePtEta = new TH2F("qualityOpAnglePtEta", ptTitle + opAngleTitle + ";Counts", 300, 0, 1200, 120, 0, 120);

	qualityOpAnglePtOneSector = new TH2F("qualityOpAnglePtOneSector", ptTitle + opAngleTitle + ";Counts",  300, 0, 1200, 120, 0, 120);
	qualityOpAnglePtEtaOneSector = new TH2F("qualityOpAnglePtEtaOneSector", ptTitle + opAngleTitle + ";Counts", 300, 0, 1200, 120, 0, 120);

	qualityOpAnglePtDiffSector = new TH2F("qualityOpAnglePtDiffSector", ptTitle + opAngleTitle + ";Counts",  300, 0, 1200, 120, 0, 120);
	qualityOpAnglePtEtaDiffSector = new TH2F("qualityOpAnglePtEtaDiffSector", ptTitle + opAngleTitle + ";Counts", 300, 0, 1200, 120, 0, 120);



	qualityPtYSpace     = new TH2F("qualityPtYSpace", rapTitle + ptTitle + ";Counts", 200, 0, 1.2, 300, 0, 1200);
	qualityPtMass       = new TH2F("qualityPtMass", ptTitle + massTitle + ";Counts", 300, 0, 1200, MassBins, MassMin,
	                               500);
	qualityYMass       = new TH2F("qualityYMass", rapTitle + massTitle + ";Counts", 200, 0, 1.2, MassBins, MassMin, 600);
	qualityOpAngleMass = new TH2F("qualityOpAngleMass", opAngleTitle + massTitle + ";Counts", 120, 0, 60, MassBins,MassMin, 600);
	qualityOpAngleMassGeant = new TH2F("qualityOpAngleMassGeant", opAngleTitle + massTitle + ";Counts", 120, 0, 60,
	                                   MassBins, MassMin, 500);
	qualityMassGeant = new TH1F("qualityMassGeant", massTitle + ";Counts", MassBins, MassMin, 500);
	qualityMass      = new TH1F("qualityMass", massTitle + ";Counts", MassBins, MassMin, 500);

	qualityPtMassOneSector       = new TH2F("qualityPtMassOneSector", ptTitle + massTitle + ";Counts", 300, 0, 1200, MassBins, MassMin, 600);
	qualityYMassOneSector       = new TH2F("qualityYMassOneSector", rapTitle + massTitle + ";Counts", 200, 0, 1.2, MassBins, MassMin, 600);


	qualityPtMassDiffSector       = new TH2F("qualityPtMassDiffSector", ptTitle + massTitle + ";Counts", 300, 0, 1200, MassBins, MassMin, 600);
	qualityYMassDiffSector       = new TH2F("qualityYMassDiffSector", rapTitle + massTitle + ";Counts", 200, 0, 1.2, MassBins, MassMin, 600);


	mass    = new TH1F * * *[centBins];
	massPi0 = new TH1F * * *[centBins];

	for (Int_t cent = 0; cent < centBins; cent++) {
		TString centName = Form("%i0%%-%i0%% ", cent + centMin, (cent + 1) + centMin);
		massIntegrated[cent] = new TH1D(Form("cent%i", cent), "Integrated " + centName + massTitle + ";Counts", MassBins,
		                                MassMin,
		                                MassMax);
		mass[cent]    = new TH1F * *[ptBins];
		massPi0[cent] = new TH1F * *[ptBins];
		for (Int_t pt = 0; pt < ptBins; pt++) {
			mass[cent][pt]    = new TH1F *[yBins];
			massPi0[cent][pt] = new TH1F *[yBins];
			for (Int_t y = 0; y < yBins; y++) {
				TString yName  = Form("| %.2f < y_cm < %.2f |", (y * yBinSize + yMin), ((y + 1) * yBinSize + yMin));
				TString ptName = Form("| %.0f < pt < %.0f |", (pt * ptBinSize + ptMin), ((pt + 1) * ptBinSize + ptMin));
				TString hName  = Form("pt%iy%icent%i", pt, y, cent);
				massPi0[cent][pt][y] = new TH1F(hName + "geant", yName + ptName + centName + massTitle + ";Counts", MassBins,
				                                MassMin,
				                                MassMax);
				mass[cent][pt][y] = new TH1F(hName, yName + ptName + centName + massTitle + ";Counts", MassBins, MassMin,
				                             MassMax);
			}
		}
	}


	// cd back to original directory
	saved->cd();
}

// ______________________________________________________________________________
HistogramManager::~HistogramManager()
{
	// Clean up all histograms.
	for (Int_t cent = 0; cent < centBins; cent++) {
		// Clean up all histograms.

		for (Int_t pt = 0; pt < ptBins; pt++) {
			delete[] mass[cent][pt];
		}
		delete[] mass[cent];
	}
	delete[] mass;
}

void HistogramManager::Hfill(const Diphoton & diphot)
{ // Float_t  assymetry = TMath::Abs((pid.energy1-pid.energy2))/(pid.energy1+pid.energy2);
	Float_t invmass     = diphot.mass;
	Float_t weight      = diphot.weight;
	Float_t weight_err2 = diphot.weight_err2;
	Float_t error       = 0;
	Int_t ptBin         = floor((diphot.pt - ptMin) / ptBinSize);
	Int_t yBin    = floor((diphot.y + 0.0000000001 - yMin) / yBinSize);
	Int_t centBin = diphot.centBin;


	if (diphot.isPi0) {
		qualityOpAngleGeant->Fill(diphot.opangle);
		qualityOpAngleMassGeant->Fill(diphot.opangle, invmass);
		qualityOpAnglePt->Fill(diphot.pt,diphot.opangle);


		if(diphot.isSameSector) 	qualityOpAnglePtOneSector->Fill(diphot.pt,diphot.opangle);
		else 	qualityOpAnglePtDiffSector->Fill(diphot.pt,diphot.opangle);

	}

	if (diphot.isEta) {
		qualityOpAngleGeantEta->Fill(diphot.opangle);
		qualityOpAnglePtEta->Fill(diphot.pt,diphot.opangle);
		if(diphot.isSameSector) 	qualityOpAnglePtEtaOneSector->Fill(diphot.pt,diphot.opangle);
		else 	qualityOpAnglePtEtaDiffSector->Fill(diphot.pt,diphot.opangle);
	}

	if (centBin == 2) {
		qualityPtYSpace->Fill(diphot.y, diphot.pt);
		qualityPtMass->Fill(diphot.pt, invmass);
		qualityYMass->Fill(diphot.y, invmass);
		qualityOpAngleMass->Fill(diphot.opangle, invmass);
		if (diphot.isPi0) qualityMassGeant->Fill(invmass);
		qualityMass->Fill(invmass);

		if(diphot.isSameSector) {
			qualityPtMassOneSector->Fill(diphot.pt, invmass);
			qualityYMassOneSector->Fill(diphot.y, invmass);
		}
		else {
			qualityPtMassDiffSector->Fill(diphot.pt, invmass);
			qualityYMassDiffSector->Fill(diphot.y, invmass);
		}


	}

	if (diphot.opangle < opAngleCut) return;

	if (centBin >= 0 && centBin < centBins) {
		Int_t binGGmass = massIntegrated[centBin]->FindBin(invmass);
		if (invmass > MassMin && invmass < MassMax) error = massIntegrated[centBin]->GetBinError(binGGmass);
		error *= error;
		massIntegrated[centBin]->Fill(invmass, weight);
		massIntegrated[centBin]->SetBinError(invmass, TMath::Sqrt(error + weight_err2));
	}

	if (centBin >= 0 && centBin < centBins &&
	    yBin >= 0 && yBin < yBins &&
	    ptBin >= 0 && ptBin < ptBins)
	{
		Int_t massBin = mass[centBin][ptBin][yBin]->FindBin(invmass);
		error  = mass[centBin][ptBin][yBin]->GetBinError(massBin);
		error *= error;
		mass[centBin][ptBin][yBin]->Fill(invmass, weight);
		mass[centBin][ptBin][yBin]->SetBinError(invmass, TMath::Sqrt(error + weight_err2));
		if (diphot.isPi0) {
			massPi0[centBin][ptBin][yBin]->Fill(invmass, weight);
		}
	}
} /* HistogramManager::Hfill */

void HistogramManager::Write(TDirectory * dir)
{
	qualityPtYSpace->Write();
	qualityPtMass->Write();
	qualityYMass->Write();
	qualityPtMassOneSector->Write();
	qualityYMassOneSector->Write();
	qualityPtMassDiffSector->Write();
	qualityYMassDiffSector->Write();



	qualityOpAngleMass->Write();
	if (isSim) {
		qualityOpAngleGeantEta->Write();
		qualityOpAngleGeant->Write();
		qualityOpAngleMassGeant->Write();
		qualityOpAnglePt->Write();
		qualityOpAnglePtEta->Write();


		qualityOpAnglePtOneSector->Write();
		qualityOpAnglePtEtaOneSector->Write();

		qualityOpAnglePtDiffSector->Write();
		qualityOpAnglePtEtaDiffSector->Write();

	}

	TDirectory * saved = gDirectory;

	dir->cd();
	qualityMass->Write();
	if (isSim) qualityMassGeant->Write();
	for (Int_t cent = 0; cent < centBins; cent++) {
		for (Int_t y = 0; y < yBins; y++) {
			for (Int_t pt = 0; pt < ptBins; pt++) {
				mass[cent][pt][y]->Write();
				if (isSim) massPi0[cent][pt][y]->Write();
			}
		}
		massIntegrated[cent]->Write();
	}


	saved->cd();
} /* HistogramManager::Write */

void HistogramManager::WriteMix(TDirectory * dir)
{
	TDirectory * saved = gDirectory;

	dir->cd();
	qualityMass->Write();
	for (Int_t cent = 0; cent < centBins; cent++) {
		for (Int_t y = 0; y < yBins; y++) {
			for (Int_t pt = 0; pt < ptBins; pt++) {
				mass[cent][pt][y]->Write();
			}
		}
		massIntegrated[cent]->Write();
	}


	saved->cd();
} /* HistogramManager::Write */

#endif /* ifndef __HISTOGRAM_MANAGER__ */
