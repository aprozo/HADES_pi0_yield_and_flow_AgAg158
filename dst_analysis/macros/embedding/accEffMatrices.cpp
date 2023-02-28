#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TStyle.h"
#include "TRatioPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
void accEffMatrices(void)
{
	enum centClass { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCentClass };
	enum mode { sim, exp };
	TFile * in[2];
	TString day = "";

	in[sim] = new TFile("embeddingSim.root", "read");
	in[exp] = new TFile("embeddingExp"+day+".root", "read");

	TH1::SetDefaultSumw2();

	TFile * out = new TFile("accEffMatrices"+day+".root", "recreate");
	TH2F * phaseSpace[2][cNumCentClass];
	TH2F * phaseSpaceReconstructedWithCuts[2][cNumCentClass];
	TH2F * acceptanceEff[2][cNumCentClass];
	out->cd();
	for (Int_t iCent = 0; iCent < (cNumCentClass - 2); iCent++) {

		TString directory;
		if (iCent != cNumCentClass - 1) directory = Form("cent%i", iCent);
		else directory = "all_cents";
		for (size_t SimOrExp = 1; SimOrExp < 2; SimOrExp++) {
			TString stringSimOrExp = SimOrExp == 0 ? "_sim" : "_exp";
			phaseSpace[SimOrExp][iCent]         = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpace");
			phaseSpaceReconstructedWithCuts[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceReconstructed");
			acceptanceEff[SimOrExp][iCent]     = (TH2F *) phaseSpaceReconstructedWithCuts[SimOrExp][iCent]->Clone(directory+"_acceptanceEff");
			acceptanceEff[SimOrExp][iCent]->Divide(phaseSpace[SimOrExp][iCent]);
			acceptanceEff[SimOrExp][iCent]->Write();
		}

	}


	out->Save();
	out->Close();
} // embeddingCreateAccEffMatrices
