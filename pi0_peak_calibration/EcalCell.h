#ifndef __ECAL_CELL__
#define __ECAL_CELL__

#include "TRef.h"
#include "TH1.h"
#include "TH2.h"
#include "TBits.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "CorrFactors.h" //
#include "TVector3.h"
#include "exp_setup_cells.h"

const Float_t pi0TrueMass = 134.977;
const Float_t nInterations = 15;

#define TCOND(a, b, c) ((a) ? (b) : (c))
typedef struct {
	Float_t newEvent, multTofRpc, multEmc, multPhoton, centBin,
	        targetBin, beamTime, RPlanePhi, impPar, emc_mq, qualityDTime, betaCand,
	        sector, cell, size, beta, energy, time, xpos, ypos, zpos, track, trackParent, mom_kine, idParent;
} PID;




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
	Int_t above_min = 200;
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
	}

private:


	TH1D * hMassAll[nEnergyBins];
	TH1D * hMassCb[nEnergyBins];
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
			hMassCb[energyBin]->Fill(mass);
		else hMassAll[energyBin]->Fill(mass);

		if (opt.Contains("cells") ) {
			Int_t energyBin = di.energyBin2;
			if (opt.Contains("cb") )
				hMassCb[energyBin]->Fill(mass);
			else hMassAll[energyBin]->Fill(mass);

		}



	}


	void InitializeAllMix()
	{
		Int_t _sector=(int)(cell/163.);
		Int_t _cell = cell%163;

		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			TString additionName = Form("Energy%.0f_%.0fMeV_Sector%iCell%i", energyBinAxis[iEnergyBin], energyBinAxis[iEnergyBin+1], _sector, _cell);
			hMassAll[iEnergyBin] = new TH1D("hMassAll" + additionName,  additionName, 150, 0, 300);
			hMassAll[iEnergyBin]->SetLineColor(kBlack);
			hMassCb[iEnergyBin] = new TH1D("hMassCb" + additionName,  additionName, 150, 0, 300);
			hMassCb[iEnergyBin]->SetLineColor(kBlue);
		}
	}




	void Initialize()
	{
		Int_t _sector=(int)(cell/163.);
		Int_t _cell = cell%163;

		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			TString additionName = Form("Energy%.0f_%.0fMeV_Sector%iCell%i", energyBinAxis[iEnergyBin], energyBinAxis[iEnergyBin+1], _sector, _cell);

			fSigFit[iEnergyBin] = new TF1("fit" + additionName, "crystalball", 90, 180);
			fSigFit[iEnergyBin]->SetParLimits(1,110,160);
			fSigFit[iEnergyBin]->SetLineColor(kGreen);
			can[iEnergyBin] = new TCanvas("can" + additionName,  additionName, 800, 600);
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



	void WriteHists()
	{

		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			hMassCb[iEnergyBin]->Write();
			hMassAll[iEnergyBin]->Write();
		}
	}

	void ReadHists(TFile * in)
	{
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {

			TString additionName = Form("Energy%.0f_%.0fMeV_Sector%iCell%i", energyBinAxis[iEnergyBin], energyBinAxis[iEnergyBin+1], (int)(cell/163.), cell%163);
			hMassAll[iEnergyBin]=(TH1D*)in->Get("hMassAll" + additionName);
			hMassCb[iEnergyBin]=(TH1D*)in->Get("hMassCb" + additionName);
		}
	}

	void Evaluate(const Int_t & iteration, CorrFactors & corrFactors, Float_t & parN)
	{
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			if (hMassAll[iEnergyBin]->GetEntries() == 0 || hMassCb[iEnergyBin]->GetEntries() == 0) continue;
			hMassSig[iEnergyBin] = getSig(hMassAll[iEnergyBin], hMassCb[iEnergyBin]);
			TF1 * fit = fSigFit[iEnergyBin];
			fit->SetParameters(hMassSig[iEnergyBin]->GetMaximum(), 140, 18, 1, 1);

			hMassSig[iEnergyBin]->Fit(fit, "RQNEM");
			Float_t mean       = fit->GetParameter(1);
			Float_t sigma      = fit->GetParameter(2);

			Int_t minBin = hMassSig[iEnergyBin]->FindBin(mean-2*sigma);
			Int_t maxBin = hMassSig[iEnergyBin]->FindBin(mean+2*sigma);
			//	Double_t SpBErr, BErr;
			Double_t SpB    = hMassAll[iEnergyBin]->Integral(minBin, maxBin);
			Double_t S      = hMassSig[iEnergyBin]->Integral(minBin, maxBin);
			Double_t Sig    = S / TMath::Sqrt(SpB);

			Float_t diff = pow( (pi0TrueMass / mean ), parN);
			Float_t correction =    corrFactors.val[iEnergyBin] * diff;


		//	if(diff<1.2&&diff>0.8&&S>300 && correction <1.2 && correction>0.75) corrFactors.val[iEnergyBin] = correction;
			if(diff<1.5&&diff>0.5&&S>300 && correction <1.5 && correction>0.5) corrFactors.val[iEnergyBin] = correction;
			else corrFactors.val[iEnergyBin]=1.0;

			//else if(iEnergyBin==0) corrFactors.val[iEnergyBin]=1.0;
			//else corrFactors.val[iEnergyBin]=corrFactors.val[iEnergyBin-1];

			fSigFit[iEnergyBin]          = move(fit);
			hMeanValues->SetBinContent(iteration + 1, iEnergyBin + 1, mean);
			hSigmaValues->SetBinContent(iteration + 1, iEnergyBin + 1, sigma);
			hCorrectionFactor->SetBinContent(iteration + 1, iEnergyBin + 1, corrFactors.val[iEnergyBin]);
		}
	}

	void Write()
	{
		hMeanValues->Write();
		hSigmaValues->Write();
		hCorrectionFactor->Write();
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			if (hMassAll[iEnergyBin]->GetEntries() == 0 || hMassCb[iEnergyBin]->GetEntries() == 0) continue;
			can[iEnergyBin]->cd();


			Float_t ymax = hMassSig[iEnergyBin]->GetMaximum();
			hMassAll[iEnergyBin]->GetXaxis()->SetRangeUser(10, 300);
			hMassAll[iEnergyBin]->GetYaxis()->SetRangeUser(-0.01 * ymax, 1.2 * ymax);


			hMassAll[iEnergyBin]->Draw();
			hMassCb[iEnergyBin]->Draw("same");
			hMassSig[iEnergyBin]->Draw("same");
			fSigFit[iEnergyBin]->Draw("same");

			Float_t mean       = fSigFit[iEnergyBin]->GetParameter(1);
			Float_t sigma      = fSigFit[iEnergyBin]->GetParameter(2);

			Int_t minBin = hMassSig[iEnergyBin]->FindBin(mean-2*sigma);
			Int_t maxBin = hMassSig[iEnergyBin]->FindBin(mean+2*sigma);

			Double_t SpBErr, BErr;
			Double_t SpB    = hMassAll[iEnergyBin]->IntegralAndError(minBin, maxBin, SpBErr, "");
			Double_t S      = hMassSig[iEnergyBin]->Integral(minBin, maxBin);
			Double_t B      = hMassCb[iEnergyBin]->IntegralAndError(minBin, maxBin, BErr, "");
			Double_t SErr   = sqrt(SpBErr * SpBErr + BErr * BErr);


			TLatex * latex = new TLatex();
			latex->SetTextSize(0.035);
			latex->SetTextFont(42);
			latex->SetTextAlign(12);

			Float_t labelOffsetX = 0.53; // 0.43 Lambda, K0s 0.53 Hypertriton, XiM


			TString label1 =
				"#splitline{Signal (#mu #pm 2#sigma) = %s}"
				"{#splitline{Signal / Background = %.2f}{Significance = %.1f}}";
			TString Signal =
				TCOND(S >= 10000.,
				      Form("(%.2f #pm %.2f) x 10^{%d}", S / TMath::Power(10., TMath::Floor(TMath::Log10(S))),
				           TMath::Ceil(SErr / TMath::Power(10., TMath::Floor(TMath::Log10(S))) * 100.) / 100.,
				           TMath::FloorNint(TMath::Log10(S))), Form("%.0f #pm %.0f", S, SErr));

			TString label2 =
				"#splitline{Mean = %.1f}"
				"{Sigma = %.1f}";

			Double_t correction = hCorrectionFactor->ProjectionY("")->GetBinContent(iEnergyBin + 1);

			latex->DrawLatexNDC(labelOffsetX, 0.70, Form(label1, Signal.Data(), S / B, S / TMath::Sqrt(SpB)));
			latex->DrawLatexNDC(labelOffsetX, 0.45, Form(label2, mean, sigma));
			latex->DrawLatexNDC(labelOffsetX, 0.35, Form("Correction = %.2f", correction));

			can[iEnergyBin]->Write();
		}
	}


};

#endif /* ifndef __HISTOGRAM_MANAGER__ */
