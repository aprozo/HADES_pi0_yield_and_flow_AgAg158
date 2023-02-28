#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TStyle.h"
#include "TRatioPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
void embeddingCreateAccEffMatrices(void)
{
	enum centClass { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCentClass };
	enum mode { sim, exp };
	TFile * in[2];
	TString day = "073";

	in[sim] = new TFile("embeddingSim.root", "read");
	in[exp] = new TFile("embeddingExp"+day+".root", "read");

	TH1::SetDefaultSumw2();

	TFile * out = new TFile("embeddingAccEffMatrices"+day+".root", "recreate");
	TH2F * phaseSpace[2][cNumCentClass];
	TH2F * phaseSpaceFull[2][cNumCentClass];
	TH2F * phaseSpaceAccepted[2][cNumCentClass];
	TH2F * phaseSpaceAcceptedFull[2][cNumCentClass];
	TH2F * phaseSpaceReconstructed[2][cNumCentClass];
	TH2F * phaseSpaceReconstructedFull[2][cNumCentClass];
	TH2F * phaseSpaceReconstructedWithCuts[2][cNumCentClass];
	TH2F * phaseSpaceReconstructedFullWithCuts[2][cNumCentClass];
	TH1F * singlePhotonEffMom[2][cNumCentClass];
	TH1F * singlePhotonEffTheta[2][cNumCentClass];
	TH1F * singlePhotonAccMom[2][cNumCentClass];
	TH1F * singlePhotonAccTheta[2][cNumCentClass];

	TH2F * singlePhotonAccMomCells[2][cNumCentClass];
	TH2F * singlePhotonEffMomCells[2][cNumCentClass];

	TH2F * acceptance[2][cNumCentClass];
	TH2F * acceptanceFull[2][cNumCentClass];
	TH2F * reconstructed[2][cNumCentClass];
	TH2F * reconstructedFull[2][cNumCentClass];
	TH2F * reconstructedWithCuts[2][cNumCentClass];
	TH2F * reconstructedWithCutsFull[2][cNumCentClass];
	TH1F * effMom[2][cNumCentClass];
	TH1F * effTheta[2][cNumCentClass];
	TH2F * effMomCells[2][cNumCentClass];
	TH2F * correctionMomCells[cNumCentClass];
	TH1F * correctionMom[cNumCentClass];
	TH1F * correctionTheta[cNumCentClass];

	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("1.2f");
	TRatioPlot *rp[cNumCentClass];
	TCanvas *can[cNumCentClass];

	Int_t rebinFactor = 20;


	for (Int_t iCent = 0; iCent < (cNumCentClass - 2); iCent++) {

		TString directory;
		if (iCent != cNumCentClass - 1) directory = Form("cent%i", iCent);
		else directory = "all_cents";

		//TDirectory * sub;
		//if (iCent != cNumCentClass - 1) sub = out->mkdir(Form("cent%i", iCent));
		//else sub = out->mkdir("all_cents");
		//sub->cd();
		TCanvas * can = new TCanvas("can","",800,600);

		for (size_t SimOrExp = 0; SimOrExp < 2; SimOrExp++) {
			TString stringSimOrExp = SimOrExp == 0 ? "_sim" : "_exp";
			//TDirectory * ssub      = sub->mkdir(stringSimOrExp);
			//ssub->cd();


			phaseSpace[SimOrExp][iCent]         = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpace");
			phaseSpaceFull[SimOrExp][iCent]     = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceFull");
			phaseSpaceAccepted[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceAccepted");
			phaseSpaceAcceptedFull[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceAcceptedFull");
			phaseSpaceReconstructed[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceReconstructed");
			phaseSpaceReconstructedFull[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceReconstructedFull");
			phaseSpaceReconstructedWithCuts[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceReconstructedWithCuts");
			phaseSpaceReconstructedFullWithCuts[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/phaseSpaceReconstructedFullWithCuts");

			acceptance[SimOrExp][iCent]     = (TH2F *) phaseSpaceAccepted[SimOrExp][iCent]->Clone(directory+stringSimOrExp+directory+stringSimOrExp+"acceptance");
			acceptanceFull[SimOrExp][iCent] = (TH2F *) phaseSpaceAcceptedFull[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"acceptanceFull");
			acceptance[SimOrExp][iCent]->GetZaxis()->SetTitle("#varepsilon_{acc}");
			acceptanceFull[SimOrExp][iCent]->GetZaxis()->SetTitle("#varepsilon_{acc}");

			reconstructed[SimOrExp][iCent]     = (TH2F *) phaseSpaceReconstructed[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"reconstructed");
			reconstructedFull[SimOrExp][iCent] = (TH2F *) phaseSpaceReconstructedFull[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"reconstructedFull");

			reconstructedWithCuts[SimOrExp][iCent]     = (TH2F *) phaseSpaceReconstructedWithCuts[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"reconstructedWithCuts");
			reconstructedWithCutsFull[SimOrExp][iCent] = (TH2F *) phaseSpaceReconstructedFullWithCuts[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"reconstructedWithCutsFull");

			reconstructed[SimOrExp][iCent]->GetZaxis()->SetTitle("#varepsilon_{reco}");
			reconstructedFull[SimOrExp][iCent]->GetZaxis()->SetTitle("#varepsilon_{reco}");

			acceptance[SimOrExp][iCent]->Divide(phaseSpace[SimOrExp][iCent]);
			acceptanceFull[SimOrExp][iCent]->Divide(phaseSpaceFull[SimOrExp][iCent]);
			reconstructed[SimOrExp][iCent]->Divide(phaseSpaceAccepted[SimOrExp][iCent]);
			reconstructedFull[SimOrExp][iCent]->Divide(phaseSpaceAcceptedFull[SimOrExp][iCent]);

			reconstructedWithCuts[SimOrExp][iCent]->Divide(phaseSpaceReconstructed[SimOrExp][iCent]);
			reconstructedWithCutsFull[SimOrExp][iCent]->Divide(phaseSpaceReconstructedFull[SimOrExp][iCent]);

			singlePhotonAccMomCells[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/singlePhotonAccMomCells");
			singlePhotonAccMomCells[SimOrExp][iCent]->RebinY(200); singlePhotonAccMomCells[SimOrExp][iCent]->Sumw2();
			singlePhotonEffMomCells[SimOrExp][iCent] = (TH2F *) in[SimOrExp]->Get(directory + "/singlePhotonEffMomCells");
			singlePhotonEffMomCells[SimOrExp][iCent]->RebinY(200); singlePhotonAccMomCells[SimOrExp][iCent]->Sumw2();

			singlePhotonAccMom[SimOrExp][iCent] = (TH1F *) in[SimOrExp]->Get(directory + "/singlePhotonAccMom");
			singlePhotonAccMom[SimOrExp][iCent]->Rebin(rebinFactor); singlePhotonAccMom[SimOrExp][iCent]->Sumw2();
			singlePhotonEffMom[SimOrExp][iCent] = (TH1F *) in[SimOrExp]->Get(directory + "/singlePhotonEffMom");
			singlePhotonEffMom[SimOrExp][iCent]->Rebin(rebinFactor); singlePhotonEffMom[SimOrExp][iCent]->Sumw2();

			singlePhotonAccTheta[SimOrExp][iCent] = (TH1F *) in[SimOrExp]->Get(directory + "/singlePhotonAccTheta");
			singlePhotonAccTheta[SimOrExp][iCent]->Rebin(2);
			singlePhotonEffTheta[SimOrExp][iCent] = (TH1F *) in[SimOrExp]->Get(directory + "/singlePhotonEffTheta");
			singlePhotonEffTheta[SimOrExp][iCent]->Rebin(2);

			effMomCells[SimOrExp][iCent] = (TH2F *) singlePhotonEffMomCells[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"_effMomCells");
			effMom[SimOrExp][iCent]      = (TH1F *) singlePhotonEffMom[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"_effMom");
			effMom[SimOrExp][iCent]->GetYaxis()->SetTitle("#varepsilon");
			effTheta[SimOrExp][iCent]    = (TH1F *) singlePhotonEffTheta[SimOrExp][iCent]->Clone(directory+stringSimOrExp+"_effTheta");

			effMomCells[SimOrExp][iCent]->Divide(singlePhotonAccMomCells[SimOrExp][iCent]);
			effTheta[SimOrExp][iCent]->Divide(singlePhotonAccTheta[SimOrExp][iCent]);
			effMom[SimOrExp][iCent]->Divide(singlePhotonAccMom[SimOrExp][iCent]);

			can->cd();

			phaseSpaceReconstructed[SimOrExp][iCent]->Draw("colz");
			can->SaveAs("../plots/phaseSpaceReconstructed"+stringSimOrExp+directory+day+".pdf");
			can->SaveAs("../plots/phaseSpaceReconstructed"+stringSimOrExp+directory+day+".C");
			can->SaveAs("../plots/phaseSpaceReconstructed"+stringSimOrExp+directory+day+".png");


			phaseSpaceAcceptedFull[SimOrExp][iCent]->Draw("colz");
			can->SaveAs("../plots/phaseSpaceAcceptedFull"+stringSimOrExp+directory+day+".pdf");
			can->SaveAs("../plots/phaseSpaceAcceptedFull"+stringSimOrExp+directory+day+".C");
			can->SaveAs("../plots/phaseSpaceAcceptedFull"+stringSimOrExp+directory+day+".png");


			phaseSpaceFull[SimOrExp][iCent]->Draw("colz");
			can->SaveAs("../plots/phaseSpaceFull"+stringSimOrExp+directory+day+".pdf");
			can->SaveAs("../plots/phaseSpaceFull"+stringSimOrExp+directory+day+".C");
			can->SaveAs("../plots/phaseSpaceFull"+stringSimOrExp+directory+day+".png");

			acceptance[SimOrExp][iCent]->Draw("colz");
			can->SaveAs("../plots/acceptance"+stringSimOrExp+directory+day+".pdf");
			can->SaveAs("../plots/acceptance"+stringSimOrExp+directory+day+".C");
			can->SaveAs("../plots/acceptance"+stringSimOrExp+directory+day+".png");

			reconstructed[SimOrExp][iCent]->Draw("colz");
			can->SaveAs("../plots/reconstructed"+stringSimOrExp+directory+day+".pdf");
			can->SaveAs("../plots/reconstructed"+stringSimOrExp+directory+day+".C");
			can->SaveAs("../plots/reconstructed"+stringSimOrExp+directory+day+".png");



		}

		//  sub->cd();
		correctionMomCells[iCent] = (TH2F *) effMomCells[exp][iCent]->Clone(directory+"_correctionMomCellsEmbedding");
		correctionMom[iCent]      = (TH1F *) effMom[exp][iCent]->Clone(directory+"_correctionMomEmbedding");
		correctionTheta[iCent]    = (TH1F *) effTheta[exp][iCent]->Clone(directory+"_correctionThetaEmbedding");
		correctionMomCells[iCent]->Divide(effMomCells[sim][iCent]);
		correctionMom[iCent]->Divide(effMom[sim][iCent]);
		correctionTheta[iCent]->Divide(effTheta[sim][iCent]);


		// can[iCent] = new TCanvas("can"+directory,"",800,1000);
		// can[iCent]->cd();
		// gStyle->SetOptStat(0);
		// rp[iCent] = new TRatioPlot(effMom[0][iCent],effMom[1][iCent]);
		// rp[iCent]->Draw("");
		// can[iCent]->Update();
		// can[iCent]->SaveAs("../plots/EmbeddingEffMom"+directory+".pdf");
		// can[iCent]->SaveAs("../plots/EmbeddingEffMom"+directory+".C");
		// can[iCent]->SaveAs("../plots/EmbeddingEffMom"+directory+".png");


		// Define two gaussian histograms. Note the X and Y title are defined
		// at booking time using the convention "Hist_title ; X_title ; Y_title"

		// TH1D *h1 = effMomCells[exp][iCent]->ProjectionX("pxExp",1,1);
		// h1->GetYaxis()->SetTitle("#varepsilon");
		// TH1D *h2 = effMomCells[sim][iCent]->ProjectionX("pxSim",1,1);
		// h1->GetXaxis()->SetRangeUser(810,970);
		// h2->GetXaxis()->SetRangeUser(810,970);
		TH1D *h1 = (TH1D *)effMom[exp][iCent];
		TH1D *h2 = (TH1D *)effMom[sim][iCent];

		Float_t yMin=h1->GetMinimum()*0.95;
		Float_t yMax=h2->GetMaximum()*1.05;
		h1->SetTitle(Form("%i-%i %%",iCent*10,(iCent+1)*10));
		// Define the Canvas
		TCanvas *c = new TCanvas("c", "canvas", 800, 600);
		// Upper plot will be in pad1
		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(0); // Upper and lower plot are joined
		pad1->SetGridx();        // Vertical grid
		pad1->Draw();            // Draw the upper pad: pad1
		pad1->cd();              // pad1 becomes the current pad
		h1->SetStats(0);         // No statistics on upper plot

		h1->SetMinimum(0);
		h1->SetMaximum(0.85);
		// h1 settings
		h1->SetLineColor(kBlue+1);
		h1->SetLineWidth(2);
		// Y axis h1 plot settings
		h1->GetYaxis()->SetTitleSize(25);
		h1->GetYaxis()->SetTitleFont(43);
		h1->GetYaxis()->SetTitleOffset(1.3);

		// h2 settings
		h2->SetLineColor(kRed);
		h2->SetLineWidth(2);

		//h1->GetXaxis()->SetRangeUser(409,570);
		h1->Draw("hist");              // Draw h1
		h2->Draw("same hist");        // Draw h2 on top of h1

	#if ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)
		// Avoid the first label (0) to be clipped.
		TAxis *axis = h1->GetYaxis();
		axis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		axis->SetLabelSize(20);
	#else
		// Do not draw the Y axis label on the upper plot and redraw a small
		// axis instead, in order to avoid the first label (0) to be clipped.
		h1->GetYaxis()->SetLabelSize(0.);
		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		axis->SetLabelSize(20);
		axis->Draw();
	#endif


		auto legend = new TLegend(0.6,0.3,0.75,0.45);
		legend->AddEntry(h1,"Exp","l");
		legend->AddEntry(h2,"Sim","l");

		legend->Draw();

		// lower plot will be in pad
		c->cd();         // Go back to the main canvas before defining pad2
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0);
		pad2->SetBottomMargin(0.2);
		pad2->SetGridx(); // vertical grid
		pad2->Draw();
		pad2->cd();      // pad2 becomes the current pad
		// Define the ratio plot
		TH1D *h3 = (TH1D*)h1->Clone("h3");
		h3->GetYaxis()->SetRangeUser(0.48,1.2);
		h3->SetLineColor(kBlack);
		h3->Sumw2();
		h3->SetStats(0);     // No statistics on lower plot
		h3->Divide(h2);
		h3->SetMarkerStyle(21);

		// Ratio plot (h3) settings
		h3->SetTitle(""); // Remove the ratio title
		// Y axis ratio plot settings
		h3->GetYaxis()->SetTitle("ratio exp/sim ");
		h3->GetYaxis()->SetNdivisions(505);
		h3->GetYaxis()->SetTitleSize(25);
		h3->GetYaxis()->SetTitleFont(43);
		h3->GetYaxis()->SetTitleOffset(1.3);
		h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetYaxis()->SetLabelSize(20);

		h3->Draw("hist");      // Draw the ratio plot


		c->SaveAs("../plots/EmbeddingEffMom"+directory+day+".pdf");
		c->SaveAs("../plots/EmbeddingEffMom"+directory+day+".C");
		c->SaveAs("../plots/EmbeddingEffMom"+directory+day+".png");





	}


	out->cd();
	out->GetList()->Write();
	out->Save();
	//out->Close();
} // embeddingCreateAccEffMatrices
