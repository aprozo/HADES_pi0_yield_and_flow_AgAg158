#ifndef __ECAL_CELL__
#define __ECAL_CELL__

#include "TRef.h"
#include "TH1.h"
#include "TH2.h"
#include "TBits.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

#include "CorrFactors.h" //
#include "TVector3.h"
#include "exp_setup_cells.h"

typedef struct {
	Float_t newEvent, multTofRpc, multEmc, multPhoton, centBin,
	        targetBin, beamTime, RPlanePhi, impPar, emc_mq, qualityDTime, betaCand,
	        sector, cell, size, beta, energy, time, xpos, ypos, zpos, track, trackParent, mom_kine, idParent;
} PID;


Int_t getEnergyBin(Float_t energy){
	for (size_t i = 0; i < nEnergyBins; i++) {
		if (energy >=energyBinAxis[i] && energy<energyBinAxis[i+1]) return i;
	}
	return -1;
}


struct Photon {
	Photon(){
	};
	Photon(const PID & pid){
		TVector3 trackVec(pid.xpos, pid.ypos, pid.zpos);
		Double_t trackLength = trackVec.Mag();
		trackVec *= (pid.energy / trackLength);
		lorVec.SetXYZM(trackVec.x(), trackVec.y(), trackVec.z(), 0);
		cell      = (int) pid.sector * 163 + (int) pid.cell;
		energyBin = (int) pid.mom_kine;
		energy    = pid.energy;
	}
	TLorentzVector lorVec;
	Int_t cell, energyBin;
	Float_t energy;
};

struct Diphoton {
	Diphoton(const Photon & gamma1, const Photon & gamma2){
		TLorentzVector pairMix = gamma1.lorVec + gamma2.lorVec;

		opangle    = gamma1.lorVec.Angle(gamma2.lorVec.Vect());
		cell1      = gamma1.cell;
		cell2      = gamma2.cell;
		energy1    = gamma1.energy;
		energy2    = gamma2.energy;
		energyBin1 = gamma1.energyBin;
		energyBin2 = gamma2.energyBin;
	}
	Float_t opangle;
	Int_t cell1, cell2;
	Int_t energyBin1, energyBin2;
	Float_t energy1, energy2;
};


TH1D * getSig(TH1D * All, TH1D * CB)
{
	Int_t below_min = 80;
	Int_t below_max = 90;
	Int_t above_min = 210;
	Int_t above_max = 250;

	Int_t BinBelowMin = All->FindBin(below_min);
	Int_t BinBelowMax = All->FindBin(below_max);
	Int_t BinAboveMin = All->FindBin(above_min);
	Int_t BinAboveMax = All->FindBin(above_max);

	Float_t int_below_mix = CB->Integral(BinBelowMin, BinBelowMax);
	Float_t int_below     = All->Integral(BinBelowMin, BinBelowMax);
	Float_t int_above_mix = CB->Integral(BinAboveMin, BinAboveMax);
	Float_t int_above     = All->Integral(BinAboveMin, BinAboveMax);
	Float_t nbCB = 1.;

	if (int_above_mix != 0) nbCB = int_above / int_above_mix;
	if (nbCB == 0) nbCB = 1.;
	CB->Scale(nbCB);
	TString histname = All->GetName();
	TH1D * Sig       = (TH1D *) All->Clone("sig" + histname);
	Sig->SetLineColor(kRed);
	Sig->Add(CB, -1);

	return Sig;
} // getSig

class EcalCell {
public:

	EcalCell(const Int_t & _sector, const Int_t & _cell) : cell(_sector * 163 + _cell){
		// ======================================//

		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			TString additionName = Form("Energy%.0f_%.0fMeV_Sector%iCell%i", energyBinAxis[iEnergyBin], energyBinAxis[iEnergyBin+1], _sector, _cell);
			hMassAll[iEnergyBin] = new TH1D("hMassAll" + additionName,  additionName, 250, 0, 400);
			hMassAll[iEnergyBin]->SetLineColor(kBlack);
			hMassCB[iEnergyBin] = new TH1D("hMassCB" + additionName,  additionName, 250, 0, 400);
			hMassCB[iEnergyBin]->SetLineColor(kBlue);
			fSigFit[iEnergyBin] = new TF1("gausfit" + additionName, "gaus", 100, 250);
			fSigFit[iEnergyBin]->SetParLimits(1,100,250);
			fSigFit[iEnergyBin]->SetLineColor(kGreen);
			can[iEnergyBin] = new TCanvas("can" + additionName,  additionName, 1000, 1000);
		}


		hMeanValues = new TH2D(Form("hMeanSector%iCell%i", _sector, _cell), Form("pi0 mean sector %i cell %i  ; Iteration ; Old energy, MeV", _sector,
		                                                                         _cell), nInterations, 0, nInterations, nEnergyBins, energyBinAxis);
		hSigmaValues = new TH2D(Form("hSigmaSector%iCell%i", _sector, _cell), Form("pi0 sigma sector %i cell %i  ; Iteration ; Old energy, MeV", _sector,
		                                                                           _cell), nInterations, 0, nInterations, nEnergyBins, energyBinAxis);
		hCorrectionFactor = new TH2D(Form("hCorrectionFactorSector%iCell%i", _sector, _cell),
		                             Form("Correction factor sector %i cell %i; Iteration ; Old energy, MeV", _sector, _cell), nInterations, 0, nInterations,
		                             nEnergyBins, energyBinAxis);

		hCorrectionFactor->GetXaxis()->SetNdivisions(5);
		hMeanValues->GetXaxis()->SetNdivisions(5);
		hSigmaValues->GetXaxis()->SetNdivisions(5);

	}

private:


	TH1D * hMassAll[nEnergyBins];
	TH1D * hMassCB[nEnergyBins];
	TH1D * hMassSig[nEnergyBins];
	TF1 * fSigFit[nEnergyBins];

	TCanvas * can[nEnergyBins];

	TH2D * hMeanValues;
	TH2D * hSigmaValues;
	TH2D * hCorrectionFactor;
	Int_t cell;

public:
	void Fill(const Diphoton & di, Option_t * option)
	{
		TString opt = option;

		opt.ToLower();

		Float_t energy1 = di.energy1, energy2 = di.energy2;
		Int_t energyBin = di.energyBin1;

		if (di.cell2 == cell) { energy2 = di.energy1; energy1 = di.energy2; energyBin = di.energyBin2; }
		Double_t mass = sqrt(2 * energy1 * energy2 * (1 - cos(di.opangle)));

		if (opt.Contains("cb") )
			hMassCB[energyBin]->Fill(mass);
		else hMassAll[energyBin]->Fill(mass);
	}

	void Evaluate(const Int_t & iteration, CorrFactors & corrFactors)
	{
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			if (hMassAll[iEnergyBin]->GetEntries() == 0 || hMassCB[iEnergyBin]->GetEntries() == 0) continue;
			hMassSig[iEnergyBin] = getSig(hMassAll[iEnergyBin], hMassCB[iEnergyBin]);
			TF1 * gausfit = fSigFit[iEnergyBin];
			gausfit->SetParameter(1, 150);
			gausfit->SetParameter(2, 20);
			gausfit->SetParameter(0, hMassSig[iEnergyBin]->GetMaximum());
			hMassSig[iEnergyBin]->Fit(gausfit, "RQN");
			Float_t mean       = gausfit->GetParameter(1);
			Float_t sigma      = gausfit->GetParameter(2);
			Float_t correction =    corrFactors.val[iEnergyBin] * pow( (pi0TrueMass / mean ), n);

			//if(correction>1.8||correction<0.7) {cout<<"bad calibration  sector "<<(int)(cell/163.)<<"  cell "<<cell%163<<"   pi0 mean = "<<mean<<endl;
				//                                correction = 0.9;}
			corrFactors.val[iEnergyBin] = correction;
			fSigFit[iEnergyBin]          = move(gausfit);
			hMeanValues->SetBinContent(iteration + 1, iEnergyBin + 1, mean);
			hSigmaValues->SetBinContent(iteration + 1, iEnergyBin + 1, sigma);
			hCorrectionFactor->SetBinContent(iteration + 1, iEnergyBin + 1, corrFactors.val[iEnergyBin]);
		}
	}

	void SetName(TString name)
	{
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			can[iEnergyBin]->SetName((TString) can[iEnergyBin]->GetName() + name);
			hMassAll[iEnergyBin]->SetTitle((TString) hMassAll[iEnergyBin]->GetTitle() + name);

		}
	}


	void Reset()
	{
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			hMassCB[iEnergyBin]->Reset("ICESM");
			hMassAll[iEnergyBin]->Reset("ICESM");

			//	hMassSig[iEnergyBin]->Reset("ICESM");
		}
	}

	void Write()
	{
		hMeanValues->Write();
		hSigmaValues->Write();
		hCorrectionFactor->Write();
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			// hMassCB[iEnergyBin]->Write();
			// hMassAll[iEnergyBin]->Write();
			if (hMassAll[iEnergyBin]->GetEntries() == 0 || hMassCB[iEnergyBin]->GetEntries() == 0) continue;
			can[iEnergyBin]->cd();
			hMassAll[iEnergyBin]->Draw();
			hMassCB[iEnergyBin]->Draw("same");
			if(hMassSig[iEnergyBin]) hMassSig[iEnergyBin]->Draw("same");
			fSigFit[iEnergyBin]->Draw("same");
			can[iEnergyBin]->Write();
		}
	}
};

#endif /* ifndef __HISTOGRAM_MANAGER__ */
