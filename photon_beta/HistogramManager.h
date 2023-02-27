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


typedef struct {
	Float_t newEvent, multTofRpc, multEmc, multPhoton, centBin,
	        targetBin, beamTime, RPlanePhi, impPar, emc_mq, qualityDTime, betaCand,
	        sector, cell, size, beta, energy, time, xpos, ypos, zpos, track, trackParent, mom_kine, idParent;
} PID;


TH1F * photon_eff[3];

Float_t getCorrectionEmbedding(const PID & pid)
{
	Float_t energy = pid.energy;
	if (energy >= 790) energy = 790;
	Float_t corr = photon_eff[(int) pid.centBin]->GetBinContent(photon_eff[(int) pid.centBin]->FindBin(energy));
	//Float_t corrError = photon_eff[(int) pid.centBin]->GetBinError(photon_eff[(int) pid.centBin]->FindBin(energy));
	return (corr > 0 && corr < 1) ? corr : 1;
}

Float_t getCorrectionErrorEmbedding(const PID & pid)
{
	Float_t energy = pid.energy;
	if (energy >= 790) energy = 790;
 	Float_t corr = photon_eff[(int) pid.centBin]->GetBinContent(photon_eff[(int) pid.centBin]->FindBin(energy));
	Float_t corrError = photon_eff[(int) pid.centBin]->GetBinError(photon_eff[(int) pid.centBin]->FindBin(energy));
	return (corr > 0 && corr < 1) ? corrError : 0;
}



struct Photon {
	Photon(){
	};
	Photon(const PID & pid) : track(0), trackParent(0){
		TVector3 trackVec(pid.xpos, pid.ypos, pid.zpos);
		Double_t trackLength = trackVec.Mag();
		Float_t energy       =  pid.energy;
		trackVec *= (energy / trackLength);
		lorVec.SetXYZM(trackVec.x(), trackVec.y(), trackVec.z(), 0);
		centBin = (Int_t)pid.centBin;
		if(!isSim) {eff     = getCorrectionEmbedding(pid);
			        weight     = 1./eff;
			        weight_err = getCorrectionErrorEmbedding(pid)/(eff*eff);}
		else {weight=1.0; weight_err=0.0;}

		if (isSim) {     // for SIM
			track       = pid.track;
			trackParent = pid.trackParent;
			mom_kine    = pid.mom_kine;
			idParent    = pid.idParent;
		}
	}
	TLorentzVector lorVec;
	Int_t centBin;
	Int_t track, trackParent,idParent;
	Float_t eff, mom_kine,weight, weight_err;
};


struct Diphoton {
	Diphoton(const Photon & gamma1, const Photon & gamma2) : fromPi0(0)
	{
		TLorentzVector pairMix = gamma1.lorVec + gamma2.lorVec;

		opangle = gamma1.lorVec.Angle(gamma2.lorVec.Vect()) * TMath::RadToDeg();
		mass    = pairMix.M();
		mom     = pairMix.P();
		theta   = pairMix.Theta() * TMath::RadToDeg();
		phi     = pairMix.Phi() * TMath::RadToDeg();
		y       = pairMix.Rapidity() - 0.822; //correct to midrapidity
		pt      = pairMix.Pt();
		centBin = gamma1.centBin;
		weight  =  gamma1.weight*gamma2.weight;
		weight_err2 = TMath::Power(gamma1.weight*gamma2.weight_err, 2)+TMath::Power(gamma1.weight_err*gamma2.weight, 2);


		fromPi0 = gamma1.track > 0 && gamma2.track > 0 && (gamma1.track != gamma2.track) &&
		          (gamma1.trackParent == gamma2.trackParent) && gamma1.trackParent > 0 &&gamma1.idParent ==7;
	}

	Float_t mass, opangle, pt, y, theta, phi, mom,weight,weight_err2;
	Int_t centBin;
	Bool_t fromPi0;
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

	TH1F**** mass;//[16][5];
	TH1F**** massPi0;//[16][5];

	Float_t startBeamTime;
	Float_t endBeamTime;

	TH1D * massIntegrated[4];


public:
	HistogramManager(TDirectory *dir, Int_t _ptBins, Int_t _yBins, Int_t _MassBins);
	virtual ~HistogramManager();
	void Hfill(const Diphoton & diphot);
	void Write(TDirectory * dir);
	void WriteMix(TDirectory * dir);
	// ClassDef(HistogramManager,1) //Manages all histograms
};

// ______________________________________________________________________________
HistogramManager::HistogramManager(TDirectory *dir, Int_t _ptBins, Int_t _yBins, Int_t _MassBins)
{

	// Create histogram manager object. Histograms will be created
	// in the "dir" directory.

	// Save current directory and cd to "dir".

	TDirectory *saved = gDirectory;
	dir->cd();
	MassBins = _MassBins;
	MassMin  = 0;
	MassMax  = 800;
	ptBins=_ptBins;
	ptMin=0;
	ptMax=1200;

	ptBinSize=(ptMax-ptMin)/Float_t(ptBins);
	yBins=_yBins;
	yMin=0.05;
	yMax=1.05;
	yBinSize=(yMax-yMin)/Float_t(yBins);

	centMin  = 0;
	centMax  = 30;
	centBins = 3;




	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ==================================== DIFFERENTIAL HISTOGRAMS======================================================//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TString xTitle = ";M_{#gamma#gamma} [MeV/c^2]", yTitle = ";Counts";

	mass = new TH1F***[centBins];
	massPi0 = new TH1F***[centBins];
	for (Int_t cent = 0; cent < centBins; cent++) {
		TString centName = Form("%i0%%-%i0%% ", cent  + centMin, (cent + 1) + centMin);

		massIntegrated[cent]=new TH1D(Form("cent%i", cent), "Integrated "+ centName + xTitle + yTitle, MassBins, MassMin, MassMax);
		mass[cent] = new TH1F**[ptBins];
		massPi0[cent] = new TH1F**[ptBins];
		for (Int_t pt = 0; pt < ptBins; pt++) {
			mass[cent][pt] = new TH1F*[yBins];
			massPi0[cent][pt] = new TH1F*[yBins];
			for (Int_t y = 0; y < yBins; y++) {
				TString yName = Form("%.2f < y_cm < %.2f ", (y * yBinSize + yMin), ((y + 1) * yBinSize + yMin));
				TString ptName = Form("%f < pt < %f ", (pt * ptBinSize + ptMin), ((pt + 1) * ptBinSize + ptMin));
				TString hName  = Form("pt%iy%icent%i", pt, y, cent);
				massPi0[cent][pt][y] = new TH1F(hName + "pi0", yName + ptName + centName + xTitle + yTitle, MassBins, MassMin,
				                                MassMax);
				mass[cent][pt][y] = new TH1F(hName, yName + ptName + centName + xTitle + yTitle, MassBins, MassMin, MassMax);
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
		for (Int_t pt = 0; pt < ptBins; pt++)
		{
		//	std::cout << "delete ";
			for (Int_t y = 0; y < yBins; y++)
			{

			}
			//std::cout << "delete[] mass[pt]" << std::endl;
			delete[] mass[cent][pt];
		}
		delete[] mass[cent];
	}
	delete[] mass;
	//std::cout << "mass[][] deleted\n";
}

void HistogramManager::Hfill(const Diphoton & diphot)
{     // Float_t  assymetry = TMath::Abs((pid.energy1-pid.energy2))/(pid.energy1+pid.energy2);
	Float_t invmass = diphot.mass;
	Float_t weight  = diphot.weight;
	Float_t weight_err2= diphot.weight_err2;
  Float_t error=0;
	Int_t ptBin   = floor((diphot.pt - ptMin) / ptBinSize);
	Int_t yBin    = floor((diphot.y + 0.00000001 - yMin) / yBinSize);
	Int_t centBin = diphot.centBin;

	if (centBin >= 0 && centBin < centBins) {
		Int_t binGGmass = massIntegrated[centBin]->FindBin(invmass);
		if(invmass>MassMin && invmass<MassMax) error = massIntegrated[centBin]->GetBinError(binGGmass);
		error*=error;
		massIntegrated[centBin]->Fill(invmass,weight);
		massIntegrated[centBin]->SetBinError(invmass,TMath::Sqrt(error+weight_err2));
	}

	if (centBin >= 0 && centBin < centBins &&
	    yBin >= 0 && yBin < yBins &&
	    ptBin >= 0 && ptBin < ptBins)
	{
		Int_t massBin = mass[centBin][ptBin][yBin]->FindBin(invmass);
		error = mass[centBin][ptBin][yBin]->GetBinError(massBin);
		error*=error;
		mass[centBin][ptBin][yBin]->Fill(invmass,weight);
		mass[centBin][ptBin][yBin]->SetBinError(invmass,TMath::Sqrt(error+weight_err2));
		if (diphot.fromPi0) massPi0[centBin][ptBin][yBin]->Fill(invmass, weight);
	}
}

void HistogramManager::Write(TDirectory * dir)
{
	TDirectory * saved = gDirectory;

	dir->cd();
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
}     /* HistogramManager::Write */

void HistogramManager::WriteMix(TDirectory * dir)
{
	TDirectory * saved = gDirectory;

	dir->cd();
	for (Int_t cent = 0; cent < centBins; cent++) {

		for (Int_t y = 0; y < yBins; y++) {
			for (Int_t pt = 0; pt < ptBins; pt++) {
				mass[cent][pt][y]->Write();
			}
		}
		massIntegrated[cent]->Write();
	}
	saved->cd();
}     /* HistogramManager::Write */

#endif /* ifndef __HISTOGRAM_MANAGER__ */
