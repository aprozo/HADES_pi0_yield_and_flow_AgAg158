#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TChain.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TNtuple.h"
#include <iostream>
#include <vector>


typedef struct {
	Float_t cent_bin, charge, mom, beta, emc_size, emc_energy, emc_mq, emc_beta,
	        emc_path, withinResolution, withinBetaSigma, sector, cell;
} PID;

TH1F * photon_eff[3];
Float_t getCorrectionEmbedding(const PID & pid)
{
	Float_t energy = pid.emc_energy;
	if (energy >= 790) energy = 790;
	Float_t corr = photon_eff[(int) pid.cent_bin]->GetBinContent(photon_eff[(int) pid.cent_bin]->FindBin(energy));
	//Float_t corrError = photon_eff[(int) pid.centBin]->GetBinError(photon_eff[(int) pid.centBin]->FindBin(energy));
	return (corr > 0 && corr < 1) ? corr : 1;
}



Int_t betaCutLeptons(void)
{
	TFile * out = new TFile("betCutExpSimRatio.root", "RECREATE");
	out->cd();

	TString path2corrections = "/lustre/hebe/hades/user/prozorov/mar19/systematic_uncertainty/correction_data/";
	TFile * FileEff = new TFile(Form(path2corrections+"embeddingAccEffMatrices.root"), "Read");
	photon_eff[0] = (TH1F *) FileEff->Get("cent0_correctionMomEmbedding");
	photon_eff[1] = (TH1F *) FileEff->Get("cent1_correctionMomEmbedding");
	photon_eff[2] = (TH1F *) FileEff->Get("cent2_correctionMomEmbedding");




	TH1F * hMom1Sigma[3][2];
	TH1F * hMomNoCut[3][2];

	TH1F * hEnergy1Sigma[3][2];
	TH1F * hEnergyNoCut[3][2];

	TH1F * hMom1SigmaWithinEnergyRes[3][2];
	TH1F * hMomNoCutWithinEnergyRes [3][2];

	TH1F * hEnergy1SigmaWithinEnergyRes[3][2];
	TH1F * hEnergyNoCutWithinEnergyRes [3][2];

	for (Int_t SimOrExp = 0; SimOrExp < 2; SimOrExp++) {
		for (Int_t iCent = 0; iCent < 3; iCent++) {
			hMom1Sigma[iCent][SimOrExp] = new TH1F(Form("hMom1SigmaCent%i%i", iCent,
			                                            SimOrExp), "; mom [MeV/c]; beta eff", 100, 0, 1000);
			hMomNoCut[iCent][SimOrExp] = new TH1F(Form("hMomNoCutCent%i%i", iCent,
			                                           SimOrExp), "; mom [MeV/c]; beta eff", 100, 0, 1000);

			hEnergy1Sigma[iCent][SimOrExp] =
				new TH1F(Form("hEnergy1SigmaCent%i%i", iCent,
				              SimOrExp), "; emc energy [MeV/c]; beta eff (1 sigma cut)", 100, 0, 1000);
			hEnergyNoCut[iCent][SimOrExp] =
				new TH1F(Form("hEnergyNoCutCent%i%i", iCent,
				              SimOrExp), ";  emc energy  [MeV/c]; beta eff (1 sigma cut)", 100, 0, 1000);

			hMom1SigmaWithinEnergyRes[iCent][SimOrExp] = new TH1F(Form("hMom1SigmaWithinEnergyResCent%i%i",
			                                                           iCent, SimOrExp), "; mom [MeV/c]; beta eff", 100, 0, 1000);
			hMomNoCutWithinEnergyRes[iCent][SimOrExp] = new TH1F(Form("hMomNoCutWithinEnergyResCent%i%i",
			                                                          iCent, SimOrExp), ";  mom  [MeV/c]; beta eff", 100, 0, 1000);

			hEnergy1SigmaWithinEnergyRes[iCent][SimOrExp] = new TH1F(Form("hEnergy1SigmaWithinEnergyResCent%i%i",
			                                                              iCent, SimOrExp), "; mom [MeV/c]; beta eff", 100, 0, 1000);
			hEnergyNoCutWithinEnergyRes[iCent][SimOrExp] = new TH1F(Form("hEnergyNoCutWithinEnergyResCent%i%i",
			                                                             iCent, SimOrExp), ";  mom  [MeV/c]; beta eff", 100, 0, 1000);
		}
	}


	PID pid;
	TChain * ntExp = new TChain("nt");
	ntExp->Add("/lustre/hebe/hades/user/prozorov/mar19/gen5/leptonsBetaCut/betaCutLeptonsExp.root");
	ntExp->SetBranchAddress("pid", &pid);

	for (Int_t iEntry = 0; iEntry < ntExp->GetEntries(); ++iEntry) {
		ntExp->GetEntry(iEntry);
		if (pid.withinResolution > 5.9 && pid.charge != 1) continue;
		hEnergyNoCut[(Int_t) pid.cent_bin][0]->Fill(pid.emc_energy);
	//	  hEnergyNoCut[(Int_t) pid.cent_bin][0]->Fill(pid.emc_energy,1./getCorrectionEmbedding(pid));
		if (pid.withinBetaSigma < 1) hEnergy1Sigma[(Int_t) pid.cent_bin][0]->Fill(pid.emc_energy);
	}


	TChain * ntSim = new TChain("nt");
	ntSim->Add("/lustre/hebe/hades/user/prozorov/mar19/gen5/leptonsBetaCut/betaCutLeptonsSim.root");
	ntSim->SetBranchAddress("pid", &pid);

	for (Int_t iEntry = 0; iEntry < ntSim->GetEntries(); ++iEntry) {
		ntSim->GetEntry(iEntry);
		if ((pid.withinResolution > 5.5 ||(pid.withinResolution > 5.8&&pid.sector<5)) && pid.charge != 1) continue;
		hEnergyNoCut[(Int_t) pid.cent_bin][1]->Fill(pid.emc_energy);
		if (pid.withinBetaSigma < 1) hEnergy1Sigma[(Int_t) pid.cent_bin][1]->Fill(pid.emc_energy);
	}
	gStyle->SetOptStat(0);

	TCanvas * can  = new TCanvas("can", "can", 800, 600);
	TLatex * latex = new TLatex();
	latex->SetTextSize(0.05);
	latex->SetTextFont(42);
	can->SaveAs("BetaCutEfficiencyNoEmbeddingEff.pdf[");
	for (Int_t iCent = 0; iCent < 3; iCent++) {
		can->cd();
		TLegend * leg = new TLegend(0.7, 0.2, 0.9, 0.3);
		for (Int_t SimOrExp = 0; SimOrExp < 2; SimOrExp++) {
			hEnergy1Sigma[iCent][SimOrExp]->Divide(hEnergyNoCut[iCent][SimOrExp]);
			hEnergy1Sigma[iCent][SimOrExp]->SetLineColor(2 + SimOrExp);
			hEnergy1Sigma[iCent][SimOrExp]->SetLineWidth(2);
			hEnergy1Sigma[iCent][SimOrExp]->GetYaxis()->SetRangeUser(0, 0.8);
			hEnergy1Sigma[iCent][SimOrExp]->Draw(SimOrExp == 0 ? "hist" : "same hist");
			latex->DrawLatexNDC(0.5, 0.3, Form("%i0-%i0 %%", iCent, iCent + 1));
			latex->DrawLatexNDC(0.5, 0.4, "within 3 sigma Energy");
			latex->DrawLatexNDC(0.5, 0.5, "positrons");
			leg->AddEntry(hEnergy1Sigma[iCent][SimOrExp], SimOrExp == 0 ? "exp" : "sim", "l");
		}
		leg->Draw("same");
		can->SaveAs("BetaCutEfficiencyNoEmbeddingEff.pdf");
		can->Clear();
	}
	can->SaveAs("BetaCutEfficiencyNoEmbeddingEff.pdf]");



	out->cd();
	for (Int_t iCent = 0; iCent < 3; iCent++) {
		hEnergy1Sigma[iCent][0]->Divide(hEnergy1Sigma[iCent][1]);

		// for (Int_t binX = 1; binX <= hEnergy1Sigma[iCent][0]->GetNbinsX(); binX++) {
		// 	Float_t energy = hEnergy1Sigma[iCent][0]->GetXaxis()->GetBinCenter(binX);
		// 	Float_t correction=photon_eff[iCent]->GetBinContent(photon_eff[iCent]->FindBin(energy));
		// 	hEnergy1Sigma[iCent][0]->SetBinContent(binX,hEnergy1Sigma[iCent][0]->GetBinContent(binX)*correction);
		// }

		hEnergy1Sigma[iCent][0]->GetYaxis()->SetRangeUser(0.5, 1.3);
    hEnergy1Sigma[iCent][0]->SetName(Form("BetaCutioCent%i",iCent));
		hEnergy1Sigma[iCent][0]->Write();
	}





	// out->cd();
	//
	// Counter_cent->Write();
	// for (Int_t iCent = 0; iCent < 3; iCent++) {
	//   TDirectory * dir = out->mkdir(Form("cent%i0_%i0", iCent, iCent + 1));
	//   dir->cd();
	//   hEnergyNoCut[iCent]->Write();
	//   hMomNoCut[iCent]->Write();
	//   hMomNoCutWithinEnergyRes[iCent]->Write();
	//   hEnergy1Sigma[iCent]->Write();
	//   hMom1Sigma[iCent]->Write();
	//   hMom1SigmaWithinEnergyRes[iCent]->Write();
	//   hEnergyNoCutWithinEnergyRes[iCent]->Write();
	//   hEnergy1SigmaWithinEnergyRes[iCent]->Write();
	//   out->cd();
	// }

	return 0;
} /* fill_ttree_diphoton */
