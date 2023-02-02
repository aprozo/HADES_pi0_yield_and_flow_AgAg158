#ifndef __ECAL_CELL__
#define __ECAL_CELL__

#include "TRef.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TBits.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "CorrFactors.h" //
#include "TVector3.h"


const Float_t pi0TrueMass = 134.977;

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
		cell   = (int) pid.sector * 163 + (int) pid.cell;
    sector = (int) pid.sector;
		energy = pid.energy;
	}
	TLorentzVector lorVec;
	Int_t cell;
	Int_t sector;
	Float_t energy;
};

struct Diphoton {
	Diphoton(const Photon & gamma1, const Photon & gamma2){
		TLorentzVector pairMix = gamma1.lorVec + gamma2.lorVec;

		opangle = gamma1.lorVec.Angle(gamma2.lorVec.Vect());
		cell1   = gamma1.cell;
		cell2   = gamma2.cell;
		energy1 = gamma1.energy;
		energy2 = gamma2.energy;
    if     (gamma1.sector==2||gamma1.sector==2) sector=2;
    else if(gamma1.sector==4||gamma2.sector==4) sector=4;
		else if(gamma1.sector==5||gamma2.sector==5) sector=5;
    else sector=0;


	}
	Float_t opangle;
	Int_t cell1, cell2;
  Int_t sector;

	Float_t energy1, energy2;
};


TH1D * getSig(TH1D * All, TH1D * CB)
{
	Int_t below_min = 80;
	Int_t below_max = 90;
	Int_t above_min = 180;
	Int_t above_max = 220;

	Int_t BinBelowMin = All->FindBin(below_min);
	Int_t BinBelowMax = All->FindBin(below_max);
	Int_t BinAboveMin = All->FindBin(above_min);
	Int_t BinAboveMax = All->FindBin(above_max);

	Float_t int_below_mix = CB->Integral(BinBelowMin, BinBelowMax);
	Float_t int_below     = All->Integral(BinBelowMin, BinBelowMax);
	Float_t int_above_mix = CB->Integral(BinAboveMin, BinAboveMax);
	Float_t int_above     = All->Integral(BinAboveMin, BinAboveMax);
	Float_t nbCB = 1.;

	if (int_above_mix!=0 && int_above >100) nbCB = int_above / int_above_mix;
	else if (int_below > 100) nbCB = int_below / int_below_mix;

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


	TH1D * hMassAll;
	TH1D * hMassCb;
	TH1D * hMassSig;
	TF1 * fit;

	TProfile * hMeanValues;
	TProfile * hSigmaValues;
	TH2D * hSigProj;
	Int_t cell;

public:

	void Fill(const Diphoton & di, Option_t * option)
	{
		TString opt = option;

		opt.ToLower();
		Float_t energy1 = di.energy1, energy2 = di.energy2;
		Double_t mass = sqrt(2 * energy1 * energy2 * (1 - cos(di.opangle)));
		if (opt.Contains("cb") ) hMassCb->Fill(mass);
		else hMassAll->Fill(mass);
	}

	void InitializeAllMix()
	{
		Int_t _sector = (int) (cell / 163.);
		Int_t _cell   = cell % 163;

		TString additionName = Form("Sector%iCell%i", _sector, _cell);

		hMassAll = new TH1D("hMassAll" + additionName, additionName, 150, 0, 300);
		hMassAll->SetLineColor(kBlack);
		hMassCb = new TH1D("hMassCb" + additionName, additionName, 150, 0, 300);
		hMassCb->SetLineColor(kBlue);


		fit = new TF1("fit" + additionName, "gaus", 90, 180);
		fit->SetParLimits(1, 0, 160);
		fit->SetLineColor(kGreen);
		// can = new TCanvas("can" + additionName,  additionName, 800, 600);
		hMeanValues = new TProfile("hMean" +additionName,
		                           Form("pi0 mean sector %i cell %i  ; Beam time ; #pi^{0} mean, MeV", _sector,
		                                _cell), 500, 67, 88,100,200); hMeanValues->SetMarkerStyle(20);
		hSigmaValues = new TProfile("hSigma" +additionName,
		                            Form("pi0 sigma sector %i cell %i  ; Beam time ;#pi^{0} sigma, MeV", _sector,
		                                 _cell), 500, 67, 88,100,200); hSigmaValues->SetMarkerStyle(20);

		hSigProj = new TH2D("hSigProj"+additionName,
		                    Form("Signal sector %i cell %i  ; Beam time ; M_{#gamma #gamma}, MeV", _sector,
		                         _cell), 500, 67, 88, 150, 0, 300);


	}

	void WriteHists()
	{
		hMassCb->Write();
		hMassAll->Write();
	}

	void ReadHists(TFile * in)
	{
		TString additionName = Form("Sector%iCell%i", (int) (cell / 163.), cell % 163);

		hMassAll = (TH1D *) in->Get("hMassAll" + additionName);
		hMassCb  = (TH1D *) in->Get("hMassCb" + additionName);
	}

	void Evaluate(const Float_t & beamTime)
	{
		if (hMassAll->GetEntries() == 0 || hMassCb->GetEntries() == 0) return;

		hMassSig = getSig(hMassAll, hMassCb);
		if(hMassSig->GetEntries()<100) return;
		fit->SetRange(110, 170);
		fit->SetParameters(hMassSig->GetMaximum(), 135, 15);

		hMassSig->Fit(fit, "RQNEM");
		Float_t mean       = fit->GetParameter(1);
		Float_t meanErr     = fit->GetParError(1);
		Float_t sigma      = fit->GetParameter(2);
		Float_t sigmaErr    = fit->GetParError(2);
		Int_t beamBin = hMeanValues->GetXaxis()->FindBin(beamTime);

		Int_t minBin = hMassSig->FindBin(mean-1*sigma);
		Int_t maxBin = hMassSig->FindBin(mean+1*sigma);

		Double_t S      = hMassSig->Integral(minBin, maxBin);


		hMeanValues->Fill(beamTime,mean,S);
		hSigmaValues->Fill(beamTime,sigma,S);

		for (size_t binY = 1; binY < hSigProj->GetNbinsY(); binY++) {
      Int_t bin = hSigProj->GetBin(beamBin,binY);
			hSigProj->AddBinContent(bin,hMassSig->GetBinContent(binY));
		}

	}

	void Write()
	{
		hMeanValues->Write();
		hSigmaValues->Write();
		hSigProj->Write();

		//
		// if (hMassAll->GetEntries() == 0 || hMassCb->GetEntries() == 0) continue;
		// can->cd();
		//
		//
		// Float_t ymax = hMassAll->GetMaximum();
		// hMassAll->GetXaxis()->SetRangeUser(10, 300);
		// hMassAll->GetYaxis()->SetRangeUser(-0.01 * ymax, 1.2 * ymax);
		//
		//
		// hMassAll->Draw();
		// hMassCb->Draw("same");
		// hMassSig->Draw("same");
		// fit->Draw("same");
		//
		// Float_t mean       = fit->GetParameter(1);
		// Float_t sigma      = fit->GetParameter(2);
		//
		// Int_t minBin = hMassSig->FindBin(mean-2*sigma);
		// Int_t maxBin = hMassSig->FindBin(mean+2*sigma);
		//
		// Double_t SpBErr, BErr;
		// Double_t SpB    = hMassAll->IntegralAndError(minBin, maxBin, SpBErr, "");
		// Double_t S      = hMassSig->Integral(minBin, maxBin);
		// Double_t B      = hMassCb->IntegralAndError(minBin, maxBin, BErr, "");
		// Double_t SErr   = sqrt(SpBErr * SpBErr + BErr * BErr);
		//
		//
		// TLatex * latex = new TLatex();
		// latex->SetTextSize(0.035);
		// latex->SetTextFont(42);
		// latex->SetTextAlign(12);
		//
		// Float_t labelOffsetX = 0.53;     // 0.43 Lambda, K0s 0.53 Hypertriton, XiM
		//
		//
		// TString label1 =
		//  "#splitline{Signal (#mu #pm 2#sigma) = %s}"
		//  "{#splitline{Signal / Background = %.2f}{Significance = %.1f}}";
		// TString Signal =
		//  TCOND(S >= 10000.,
		//        Form("(%.2f #pm %.2f) x 10^{%d}", S / TMath::Power(10., TMath::Floor(TMath::Log10(S))),
		//             TMath::Ceil(SErr / TMath::Power(10., TMath::Floor(TMath::Log10(S))) * 100.) / 100.,
		//             TMath::FloorNint(TMath::Log10(S))), Form("%.0f #pm %.0f", S, SErr));
		//
		// TString label2 =
		//  "#splitline{Mean = %.1f}"
		//  "{Sigma = %.1f}";
		//
		// Double_t correction = hCorrectionFactor->ProjectionY("")->GetBinContent(iBeamBin + 1);
		//
		// latex->DrawLatexNDC(labelOffsetX, 0.70, Form(label1, Signal.Data(), S / B, S / TMath::Sqrt(SpB)));
		// latex->DrawLatexNDC(labelOffsetX, 0.45, Form(label2, mean, sigma));
		// latex->DrawLatexNDC(labelOffsetX, 0.35, Form("Correction = %.2f", correction));
		//
		// can->Write();
	} /* Write */
};

#endif /* ifndef __HISTOGRAM_MANAGER__ */
