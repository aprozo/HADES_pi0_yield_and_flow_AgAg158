#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>
#include <stdio.h>


const Int_t nIterations = 14;
const Int_t nParameters = 7;
using namespace std;


void set_nice_color()
{
  new TColor(2000, (80. / 255.), (120. / 255.), (150. / 255.) ); // magenta-ish
  new TColor(2001, (66. / 255.), (98. / 255.), (255. / 255.) );  // blue-ish
  new TColor(2002, (255. / 255.), (89. / 255.), (74. / 255.));   // red-ish
  new TColor(2003, (25. / 255.), (170. / 255.), (25. / 255.));   // green-ish
  new TColor(2004, (255. / 255.), (166. / 255.), (33. / 255.));  // yellow-ish
  new TColor(2005, (0. / 255.), (170. / 255.), (255. / 255.) );  // azur-ish
  new TColor(2006, (153. / 255.), (0. / 255.), (153. / 255.) );  // magenta-ish
}

void parNCompare()
{
  set_nice_color();
  Int_t model_colors[8] = { 2000, 2001, 2002 };
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(55);
  std::vector<Float_t> parN = { 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 };
  TFile * in[nParameters][nIterations];
  TH1D * mean[nParameters];
  TH1D * sigma[nParameters];
  TH1D * diffMean[nParameters];

  TH2D * mean2d     = new TH2D("mean2d", "mean2d", 15, -0.5, 15 - 0.5, 7, 0.7, 2.1);
  TH2D * sigma2d    = new TH2D("sigma2d", "sigma2d", 15, -0.5, 15 - 0.5, 7, 0.7, 2.1);
  TH2D * diffMean2d = new TH2D("diffMean2d", "diffMean2d", 15, -0.5, 15 - 0.5, 7, 0.7, 2.1);


  TCanvas * can  = new TCanvas("can", "can", 800, 600);
  TCanvas * can2 = new TCanvas("can2", "can2", 800, 600);
  TLegend * leg  = new TLegend(0.82, 0.50, 1.0, 0.92);
  leg->SetBorderSize(0);
  TLegend * leg2 = new TLegend(0.79, 0.15, 0.87, 0.66);
  leg2->SetBorderSize(0);

  for (Int_t iPar = 0; iPar < parN.size(); iPar++) {
    mean[iPar]     = new TH1D(Form("mean%.1f", parN[iPar]), "mean", 15, 0, 15);
    sigma[iPar]    = new TH1D(Form("sigma%.1f", parN[iPar]), "sigma", 15, 0, 15);
    diffMean[iPar] = new TH1D(Form("diffMean%.1f", parN[iPar]), "diffMean", 15, -0.5, 15 - 0.5);

    for (Int_t iIter = 0; iIter < nIterations; iIter++) {
      in[iPar][iIter] = new TFile(Form("checkHistogramsIteration%iparN%.1f.root", iIter, parN[iPar]), "Read");
      TH2D * meanTemp  = (TH2D *) in[iPar][iIter]->Get("hMeanSector0Cell0");
      TH2D * sigmaTemp = (TH2D *) in[iPar][iIter]->Get("hSigmaSector0Cell0");
      TH1D * temp      = meanTemp->ProjectionX("pxMean", 1, 1);

      Float_t value = temp->GetBinContent(temp->FindBin(iIter));
      mean[iPar]->SetBinContent(iIter + 1, value);
      Float_t bin = mean2d->FindBin(iIter, parN[iPar]);
      mean2d->SetBinContent(bin, value);
      diffMean2d->SetBinContent(bin, abs(134.977 - value));

      diffMean[iPar]->SetBinContent(iIter + 1, abs(134.977 - value));
    }


    can->cd();
    gPad->SetRightMargin(0.2);
    leg->AddEntry(diffMean[iPar], Form("n = %.1f", parN[iPar]), "lep");
    diffMean[iPar]->SetLineWidth(2);
    diffMean[iPar]->SetMarkerStyle(20);
    diffMean[iPar]->SetMarkerColor(2000 + iPar);
    diffMean[iPar]->SetLineColor(2000 + iPar);

    diffMean[iPar]->GetYaxis()->SetTitle("|m_{sum} - m_{PDG}|, MeV");
    //	diffMean[iPar]->GetYaxis()->SetTitle("#sigma_{sum}, MeV");

    diffMean[iPar]->GetXaxis()->SetTitle("Iteration");
    diffMean[iPar]->Draw("same LP");
  }
  leg->Draw("same");
  gPad->SetLogy();
  can->SaveAs("calibrationConvergenceParNSigma.pdf");


  can2->cd();
  gPad->SetLogy();
  for (Int_t iIter = 0; iIter < nIterations; iIter++) {
    TH1D * proj = diffMean2d->ProjectionY(Form("iter%i", iIter), iIter + 1, iIter + 1);
    proj->GetYaxis()->SetRangeUser(10E-4, 15);
    leg2->AddEntry(proj, Form("Iter %i", iIter), "lep");
    proj->SetLineWidth(2);
    proj->SetMarkerStyle(20);
    proj->SetMarkerColor(kOrange - 9 + iIter);
    proj->SetLineColor(kOrange - 9 + iIter);
    proj->GetYaxis()->SetTitle("|m_{sum} - m_{PDG}|, MeV");
    // proj->GetYaxis()->SetTitle("#sigma_{sum}, MeV");

    proj->GetXaxis()->SetTitle("n");
    proj->Draw("same PL");
  }
  leg2->Draw("same");
  //	gPad->SetLogy();
  can2->SaveAs("calibrationConvergenceParNIterationsSigma.pdf");
} /* fill_hists_diphoton */
