#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TStyle.h"
#include "TRatioPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
void drawEfficiency(void) {

	TFile * inFile = new TFile("effLeptons.root", "read");

	gStyle->SetOptStat(0);

	TCanvas * can = new TCanvas("can","",800,600);

	TH1D *h1 = (TH1D *)inFile->Get("exp_leptonmom_efficiency_cent2");
	TH1D *h2 = (TH1D *)inFile->Get("sim_leptonmom_efficiency_cent2");

	Float_t yMin=h1->GetMinimum()*0.95;
	Float_t yMax=h2->GetMaximum()*1.05;
	h1->SetTitle("0-30 %; E[MeV]; #varepsilon");
  h2->SetTitle("0-30 %; E[MeV]; #varepsilon");
	// Define the Canvas
	TCanvas *c = new TCanvas("c", "canvas", 800, 600);
	// Upper plot will be in pad1
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0);     // Upper and lower plot are joined
	pad1->SetGridx();            // Vertical grid
	pad1->Draw();                // Draw the upper pad: pad1
	pad1->cd();                  // pad1 becomes the current pad
	h1->SetStats(0);             // No statistics on upper plot

	h1->SetMinimum(0);
	h1->SetMaximum(0.85);
	// h1 settings
	h1->SetLineColor(kBlue+1);
	h1->SetLineWidth(2);
	// Y axis h1 plot settings
	h1->GetYaxis()->SetTitleSize(25);
	h1->GetYaxis()->SetTitleFont(43);
	h1->GetYaxis()->SetTitleOffset(1.3);

  // X axis h1 plot settings


	// h2 settings
	h2->SetLineColor(kRed);
	h2->SetLineWidth(2);

	h1->Draw("hist");                  // Draw h1
	h2->Draw("same hist");            // Draw h2 on top of h1

	#if ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)
	// Avoid the first label (0) to be clipped.
	TAxis *axis = h1->GetYaxis();
	axis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
	axis->SetLabelFont(43);     // Absolute font size in pixel (precision 3)
	axis->SetLabelSize(20);
	#else
	// Do not draw the Y axis label on the upper plot and redraw a small
	// axis instead, in order to avoid the first label (0) to be clipped.
	h1->GetYaxis()->SetLabelSize(0.);
	TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
	axis->SetLabelFont(43);     // Absolute font size in pixel (precision 3)
	axis->SetLabelSize(20);
	axis->Draw();
	#endif


	auto legend = new TLegend(0.6,0.3,0.75,0.45);
	legend->AddEntry(h1,"Exp","l");
	legend->AddEntry(h2,"Sim","l");
	legend->Draw();

	// lower plot will be in pad
	c->cd();             // Go back to the main canvas before defining pad2
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();     // vertical grid
	pad2->Draw();
	pad2->cd();          // pad2 becomes the current pad
	// Define the ratio plot
	TH1D *h3 = (TH1D*)h1->Clone("h3");
	h3->GetYaxis()->SetRangeUser(0.48,1.2);
	h3->SetLineColor(kBlack);
	h3->Sumw2();
	h3->SetStats(0);         // No statistics on lower plot
	h3->Divide(h2);
	h3->SetMarkerStyle(21);

	// Ratio plot (h3) settings
	h3->SetTitle("");     // Remove the ratio title
	// Y axis ratio plot settings
	h3->GetYaxis()->SetTitle("ratio exp/sim ");
	h3->GetYaxis()->SetNdivisions(505);
	h3->GetYaxis()->SetTitleSize(25);
	h3->GetYaxis()->SetTitleFont(43);
	h3->GetYaxis()->SetTitleOffset(1.1);
	h3->GetYaxis()->SetLabelFont(43);     // Absolute font size in pixel (precision 3)
	h3->GetYaxis()->SetLabelSize(20);

  h3->GetXaxis()->SetLabelSize(0.1);
  h3->GetXaxis()->SetTitleSize(0.1);

	h3->Draw("hist");          // Draw the ratio plot


	c->SaveAs("efficiency.pdf");

	//out->Close();


}
