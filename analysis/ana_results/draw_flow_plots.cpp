#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include "TF1.h"
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TROOT.h>
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TGaxis.h"


#define KPLUS      0
#define KZERO      1
#define HigherHarm 0


using namespace std;

Double_t funcFourier2(Double_t * x, Double_t * par)
{
  return par[0] / (2 * TMath::Pi())
         * (1 + 2 * par[1] * TMath::Cos(x[0] * TMath::DegToRad()) + 2 * par[2] * TMath::Cos(
           2 * x[0] * TMath::DegToRad()));
}

Double_t funcFourier4(Double_t * x, Double_t * par)
{
  return par[0] / (2 * TMath::Pi()) * (1
         + 2 * par[1] * TMath::Cos(1. * x[0] * TMath::DegToRad())
         + 2 * par[2] * TMath::Cos(2. * x[0] * TMath::DegToRad())
         + 2 * par[3] * TMath::Cos(3. * x[0] * TMath::DegToRad())
         + 2 * par[4] * TMath::Cos(4. * x[0] * TMath::DegToRad()));
}

const int NYbins = 4;
const int NPTbins = 4;
const int NX = NYbins, NY = NPTbins, NT = NX * NY;
double XLabelSize[NT], YLabelSize[NT], XTitleOffset[NT], YTitleOffset[NT], XLabelOffset[NT], YLabelOffset[NT];
double Yratio[NY] = { 1., 1., 1., 1. }, Xratio[NX] = { 1., 1., 1., 1. };


TCanvas * canvas;
TPad * canv[NT];
void setupCanvas(const int NX, const int NY);
const double RightMargin  = 0.125;
const double LeftMargin   = 0.125;
const double TopMargin    = 0.125;
const double BottomMargin = 0.095;
const double CanvasSize   = 1400.;

const double Xorigin = 0.085;
const double Yorigin = 0.05;

const int XaxisDivs    = 304;
const double XaxisStop = 0.81;
const double XaxisMin  = -180.;
const double XaxisMax  = 180.;

const int YaxisDivs    = 204;
const double YaxisStop = 0.81;
const double YaxisMin  = 0.04;
const double YaxisMax  = 0.12;

const int myTextAlign = 22;
const int myTextSize  = 26;
const int myTextFont  = 43;
const int myTextColor = 1;

void setupCanvas(const int NX, const int NY)
{
  // --------change as you like following three lines---------
  double rightM = RightMargin, leftM = LeftMargin, topM = TopMargin, bottomM = BottomMargin;
  double gapY = 0., gapX = 0.;
  // ---------------------
  double height[NY], width[NX];
  double deltaW = 1. - leftM - rightM - gapX * (NX - 1);
  double deltaH = 1. - topM - bottomM - gapY * (NY - 1);
  double x1[NX], x2[NX], y1[NY], y2[NY];
  char cid[100];

  // -------------set up for NX x NY plots--------------------
  canvas = new TCanvas("canvas", "plots", 0, 0, CanvasSize * (1. + LeftMargin + RightMargin - TopMargin - BottomMargin),
      CanvasSize); // to keep square shape of individual pads
  canvas->Divide(NX, NY);
  // ---------------------------------
  double xall = 0., yall = 0.;
  for (int i = 0; i < NX; i++) xall += Xratio[i];
  for (int i = 0; i < NY; i++) yall += Yratio[i];
  for (int i = 0; i < NX; i++) width[i] = deltaW * Xratio[i] / xall;
  for (int i = 0; i < NY; i++) height[i] = deltaH * Yratio[i] / yall;
  for (int i = 0; i < NX; i++) {
    if (i == 0) {
      x1[i] = leftM;
      x2[i] = x1[i] + width[i];
    } else {
      x1[i] = x2[i - 1] + gapX;
      x2[i] = x1[i] + width[i];
    }
  }
  for (int i = 0; i < NY; i++) {
    if (i == 0) {
      y1[i] = bottomM;
      y2[i] = y1[i] + height[i];
    } else {
      y1[i] = y2[i - 1] + gapY;
      y2[i] = y1[i] + height[i];
    }
  }
  // ----------------------------------------
  int iplot = 0;
  for (int iy = 0; iy < NY; iy++) {
    for (int ix = 0; ix < NX; ix++) {
      sprintf(cid, "canvas_%d", iplot + 1);
      canv[iplot] = (TPad *) canvas->GetListOfPrimitives()->FindObject(cid);
      canv[iplot]->SetPad(x1[ix], y1[iy], x2[ix], y2[iy]);
      canv[iplot]->SetFillColor(0);
      canv[iplot]->SetFrameFillColor(0);

      canv[iplot]->SetTopMargin(0.);
      canv[iplot]->SetBottomMargin(gapY);
      canv[iplot]->SetLeftMargin(0.);
      canv[iplot]->SetRightMargin(gapX);

      XLabelSize[iplot]   = 0.;
      YLabelSize[iplot]   = 0.;
      XTitleOffset[iplot] = 10.;
      YTitleOffset[iplot] = 10.;
      XLabelOffset[iplot] = 10.;
      YLabelOffset[iplot] = 10.;

      iplot++;
    }
  }
  // ====================end of set up========================
} // setupCanvas

void draw_flow_plots()
{
  gStyle->SetOptTitle(kTRUE);
  gStyle->SetOptStat(0);

  // input and output setup
  TString sInFile = "";
  TString sDummy = "", sFileOutPdf = "";

   sInFile =
   "~/mar19/diphotons/gen5a/systematic_uncertainty/results/parameters_N_MASS_BINS1/histograms_to_pt_y/flow.root";

//  sInFile =
  //  "~/mar19/diphotons/gen5a/systematic_uncertainty/macros/histograms_to_pt_y/flowSim.root";

  sFileOutPdf = "all_flow_plots_exp.pdf";

  TF1 * fDummy = new TF1("fDummy", "1.", -1000., 1000.);

  TF1 * fitFourier = new TF1("fit", funcFourier2, -180, 180, 3);

  TFile * fin = new TFile(sInFile.Data(), "read");

  const int nUsedCentrClass = 3;
  const int firstCentrClass = 0;

  const int nFlowBins = 12;

  const int nUsedRapidityBins   = 4;
  const int firstRapidityBin    = 0;
  const double minRapidity      = 0.1;
  const double binWidthRapidity = 0.2;

  const int nUsedTransMomBins = NPTbins;
  const double minTransMom    = 100.;

  const int firstTransMomBin    = 0;
  const double binWidthTransMom = 200.;
  int selectedPhaseSpaceBins[nUsedCentrClass][NYbins][2] = { // for each centrality and rapidity bin is given lowest and highest used Pt bin
    { { 0, 4 }, {  0, 4 }, {  0, 4 }, {  0, 4 } },           // 0 - 10%
    { { 0, 4 }, {  0, 4 }, {  0, 4 }, {  0, 4 } },           // 10 - 20%
    { { 0, 4 }, {  0, 4 }, {  0, 4 }, {  0, 4 } }// 20 - 30%
  };

  const double nEntriesGoodPS = 50.;

  setupCanvas(NX, NY);
  // --------------------------
  char cid[50];

  sDummy = sFileOutPdf + "[";
  canvas->Print(sDummy.Data());

  TH1D * hwork;

  for (int icc = 0; icc < nUsedCentrClass; icc++) {
    int iplot = 0;
    for (int iycm = 0; iycm < nUsedRapidityBins; iycm++) {
      for (int ipt = 0; ipt < nUsedTransMomBins; ipt++) {
        if (iycm == (nUsedRapidityBins - 1) && ipt == (nUsedTransMomBins - 1)) continue;

        canv[iplot]->cd();
        //  std::cout << icc << " " << ipt << " " << iycm << " " << '\n';


        sDummy = Form("cent%i/standardFlowpt%iy%icent%i", firstCentrClass + icc,
             firstTransMomBin + ipt,firstRapidityBin + iycm,firstCentrClass + icc);

        // sDummy = Form("cent%i/cosDphipt%iy%icent%i", firstCentrClass + icc, firstTransMomBin + ipt,
        //  firstRapidityBin + iycm, firstCentrClass + icc);


        // if ( (ipt < selectedPhaseSpaceBins[icc][iycm][0]) || (ipt > selectedPhaseSpaceBins[icc][iycm][1]) ) {
        //   sDummy = Form("cent%i/flowpt%iy%icent%i", firstCentrClass + icc, firstTransMomBin + ipt,
        //       firstRapidityBin + iycm, firstCentrClass + icc);
        //   hwork = (TH1D *) fin->Get(sDummy.Data() );
        //   hwork->Reset("ICES");
        // } else {
        hwork = (TH1D *) fin->Get(sDummy.Data() );
        // }

        // if (nEntries < nEntriesGoodPS) {
        //  hwork->Reset("ICES");
        // delete hwork->GetListOfFunctions()->FindObject("fit");

        // else {
        //   hwork->Scale(1. / nEntries);
        //   hwork->Fit(fitFourier, "QR");
        // }


        hwork->SetMinimum(YaxisMin);
        hwork->SetMaximum(YaxisMax);
        hwork->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);


        hwork->SetTitle("");
        hwork->GetXaxis()->SetLabelOffset(XLabelOffset[iplot]);
        hwork->GetXaxis()->SetNdivisions(-1 * XaxisDivs);
        hwork->GetYaxis()->SetLabelOffset(YLabelOffset[iplot]);
        hwork->GetYaxis()->SetNdivisions(-1 * YaxisDivs);

        hwork->Draw();

        iplot++;
      }
    }

    // labels:
    TString str = "";
    TLatex label;
    label.SetTextAlign(myTextAlign);
    label.SetTextColor(myTextColor);
    label.SetTextFont(myTextFont);
    label.SetTextSize(myTextSize);
    canvas->cd();
    for (int ix = 0; ix < nUsedRapidityBins; ix++) {
      double x1 = canv[ix]->GetAbsXlowNDC();
      double y1 = canv[ix]->GetAbsYlowNDC();
      double x2 = canv[ix]->GetAbsXlowNDC() + canv[ix]->GetAbsWNDC();
      double y2 = canv[ix]->GetAbsYlowNDC() + canv[ix]->GetAbsHNDC();

      double x = 0.5 * (x1 + x2);
      double y = 1 - 0.5 * TopMargin;




          str = Form("p_{T} #in (%.0f;%.0f) MeV/c", minTransMom + (ix + firstTransMomBin) * binWidthTransMom,
              minTransMom + (ix + firstTransMomBin + 1) * binWidthTransMom);
      label.SetTextAngle(45.);
      label.DrawLatexNDC(x, y, str.Data());
    }


    for (int iy = 0; iy < nUsedTransMomBins; iy++) {
      double x1 = canv[iy * nUsedRapidityBins]->GetAbsXlowNDC();
      double y1 = canv[iy * nUsedRapidityBins]->GetAbsYlowNDC();
      double x2 = canv[iy * nUsedRapidityBins]->GetAbsXlowNDC() + canv[iy * nUsedRapidityBins]->GetAbsWNDC();
      double y2 = canv[iy * nUsedRapidityBins]->GetAbsYlowNDC() + canv[iy * nUsedRapidityBins]->GetAbsHNDC();

      double y = 0.5 * (y1 + y2);
      double x = 1 - 0.5 * RightMargin;

      str =
        Form("y_{cm} #in (%.1f;%.1f)", minRapidity + (iy + firstRapidityBin) * binWidthRapidity,
          minRapidity + (iy + firstRapidityBin + 1) * binWidthRapidity);

      label.SetTextAngle(50.);
      label.DrawLatexNDC(x, y, str.Data());
    }

    canvas->cd();

    canvas->Range(0, 0, 1, 1);
    TGaxis * ycomaxis = new TGaxis(Xorigin, Yorigin, Xorigin, YaxisStop, YaxisMin, YaxisMax, YaxisDivs, "N");
    ycomaxis->SetName("ycomaxis");
    ycomaxis->Draw();
    TGaxis * xcomaxis = new TGaxis(Xorigin, Yorigin, XaxisStop, Yorigin, XaxisMin, XaxisMax, XaxisDivs, "N");
    xcomaxis->SetName("xcomaxis");
    xcomaxis->Draw();

    TLatex title;
    title.SetTextAlign(myTextAlign);
    title.SetTextFont(myTextFont);
    title.SetTextSize(2 * myTextSize);
    title.SetTextColor(myTextColor);

    title.SetTextAngle(90.);
    title.DrawLatexNDC(0.5 * LeftMargin, 1. - 0.5 * (1. - YaxisStop), "#frac{1}{N} #frac{dN}{d#Delta#phi}");

    title.SetTextAngle(0.);
    title.DrawLatexNDC(1. - 0.5 * (1. - XaxisStop), 0.5 * BottomMargin, "#Delta#phi [#circ]");


    canvas->Update();
    canvas->Print(sFileOutPdf.Data());
    canvas->Clear("D");
    canvas->Update();
  }


  // end of printing
  sDummy = sFileOutPdf + "]";
  canvas->Print(sDummy.Data());
} // all_flow_plots
