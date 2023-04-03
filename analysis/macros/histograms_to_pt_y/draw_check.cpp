#include <stdio.h>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooBifurGauss.h"
#include "RooCBShape.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "TFitResult.h"

#include "../GlobalParameters.h"


using namespace std;
using namespace RooFit;

const Float_t minMass = 0;
const Float_t maxMass = 250;


Int_t getSig(TH1D * massAll, TH1D * massCb)
{
  if (massAll->Integral() < 100 || massCb->Integral() < 100) return -1;

  const Float_t minMassBin = massAll->FindBin(minMass);
  const Float_t maxMassBin = massAll->FindBin(maxMass);

  TLatex * latex = new TLatex();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42);
  latex->SetTextAlign(12);


  TString label1 = "#splitline{Signal (#mu #pm 2#sigma) = %s}"
    "{#splitline{Signal / Background = %.2f}{Significance = %.1f}}";

  TString label2 =
    "#splitline{Mean = %.1f}"
    "{Sigma = %.1f}";


  // latex->DrawLatexNDC(labelOffsetX, 0.85, "p+p #sqrt{s} = 3.46 GeV");

  Int_t belowMin                   = 80;
  Int_t belowMax                   = 90;
  Int_t aboveMin                   = 180;
  Int_t aboveMax                   = 200;
  Int_t binBelowMin                = massAll->FindBin(belowMin);
  Int_t binBelowMax                = massAll->FindBin(belowMax);
  Int_t binAboveMin                = massAll->FindBin(aboveMin);
  Int_t binAboveMax                = massAll->FindBin(aboveMax);
  Float_t integralBelowMix         = massCb->Integral(binBelowMin, binBelowMax);
  Float_t integralBelow            = massAll->Integral(binBelowMin, binBelowMax);
  Float_t integralAboveMix         = massCb->Integral(binAboveMin, binAboveMax);
  Float_t integralAbove            = massAll->Integral(binAboveMin, binAboveMax);
  Float_t normFactorIntegralRegion = 1.;

  if (integralBelowMix == 0) normFactorIntegralRegion = integralAbove / integralAboveMix;
  else normFactorIntegralRegion = 0.3 * integralBelow / integralBelowMix + 0.7 * integralAbove / integralAboveMix;

  TH1D * h_scale = (TH1D *) massAll->Clone("h_scale");
  h_scale->Add(massCb);
  h_scale->Divide(massCb);

  TF1 line("line", "pol0", minMass, maxMass);
  h_scale->Fit("line", "EMQ", "", 250, 300);
  Double_t normFactorLineFit = line.GetParameter(0);
  normFactorLineFit -= 1.;
  massAll->GetXaxis()->SetRangeUser(minMass, maxMass);
  massCb->GetXaxis()->SetRangeUser(minMass, maxMass);
  h_scale = (TH1D *) massAll->Clone("h_scale");
  h_scale->GetXaxis()->SetRangeUser(minMass, maxMass);
  massCb->Scale(normFactorLineFit);
  TH1D * massSig = (TH1D *) massAll->Clone("sig");
  massSig->Add(massCb, -1.);

  TF1 * fitSig;
  massSig->GetXaxis()->SetRangeUser(50, 300);
  fitSig = new TF1("fitSig", "gaus", 100, 650);
  fitSig->SetLineColor(kGreen);
  fitSig->SetRange(minMass, maxMass);
  fitSig->SetParameters(massSig->GetMaximum(), 135, 16);
  massSig->Fit(fitSig, "RQEM");

  TFitResultPtr r_mass = massSig->Fit(fitSig, "QEMNS");
  if (!r_mass.Get()) return -1;

  Float_t chi2 = r_mass->Chi2();
  Float_t ndf  = r_mass->Ndf();
  // if (!(chi2 / ndf < 3 && Nsig / sqrt(Nsig + Nbkg) > 5 &&Nbkg >0)) {return 1;}

  Double_t mean = fitSig->GetParameter(1);
  Double_t sigma = fitSig->GetParameter(2);
  Double_t min = mean - 2 * sigma;
  Double_t max = mean + 2 * sigma;
  Int_t minBin = massSig->FindBin(min);
  Int_t maxBin = massSig->FindBin(max);
  Double_t SpBErr, BErr, SpB, S, B, SErr, Sig, SigErr;
  SpB    = massAll->IntegralAndError(minBin, maxBin, SpBErr, "");
  S      = massSig->Integral(minBin, maxBin);
  B      = massCb->IntegralAndError(minBin, maxBin, BErr, "");
  SErr   = sqrt(SpBErr * SpBErr + BErr * BErr);
  Sig    = S / TMath::Sqrt(SpB);
  SigErr = TMath::Sqrt(TMath::Power((SpB + B) / 2. / TMath::Power(SpB, 1.5) * SpBErr,
      2.) + TMath::Power(BErr / TMath::Sqrt(SpB), 2.));

  massSig->DrawClone();
  latex->DrawLatexNDC(0.35, 0.70, Form(label1, Form("%.2e", S), S / B, S / TMath::Sqrt(SpB)));
  latex->DrawLatexNDC(0.35, 0.45, Form(label2, mean, sigma));
  latex->DrawLatexNDC(0.35, 0.35, Form("Normalization = %.2e", normFactorLineFit));
  latex->DrawLatexNDC(0.35, 0.25, Form("chi2 / ndf = %.2f %.2f", chi2, ndf));

  fitSig->DrawClone("same");

  return 0;
} // getSig

// ================================================//================================================

Int_t getSigRooFit(TH1D * massAll, TH1D * massCb)
{
  const Float_t minMassBin = massAll->FindBin(minMass);
  const Float_t maxMassBin = massAll->FindBin(maxMass);
  Int_t BinBelowMin        = massAll->FindBin(80);
  Int_t BinBelowMax        = massAll->FindBin(90);
  Int_t BinAboveMin        = massAll->FindBin(180);
  Int_t BinAboveMax        = massAll->FindBin(200);

  Float_t int_below_mix = massCb->Integral(BinBelowMin, BinBelowMax);
  Float_t int_below     = massAll->Integral(BinBelowMin, BinBelowMax);
  Float_t int_above_mix = massCb->Integral(BinAboveMin, BinAboveMax);
  Float_t int_above     = massAll->Integral(BinAboveMin, BinAboveMax);
  Float_t nbCB = 1;

  if (int_below_mix == 0) nbCB = int_above / int_above_mix;
  else nbCB = 0.3 * int_below / int_below_mix + 0.7 * int_above / int_above_mix;

  //	massCb->Scale(nbCB);
  massAll->GetXaxis()->SetRangeUser(minMass, maxMass);
  massCb->GetXaxis()->SetRangeUser(minMass, maxMass);
  Double_t integralAll = massAll->Integral(minMassBin, maxMassBin);


  RooRealVar x("x", "M_{#gamma#gamma}[MeV/c^{2}]", minMass, maxMass);
  x.setRange("window", 50, 250);

  RooDataHist dhAll("dhAll", "mass spectra of photon pairs in the same event", RooArgList(x), massAll);


  RooDataHist dhBkg("dhBkg", "mass spectra of photon pairs in mixed events", RooArgList(x), massCb);
  Double_t sum1Bkg = dhBkg.sum(0);
  RooHistPdf pdfBkg("pdfBkg", "background PDF", RooArgList(x), dhBkg);
  RooRealVar nBkg("nBkg", "nbackground", integralAll, 1., 2. * integralAll);
  RooExtendPdf epdfBkg("epdfBkg", "epdfBkg", pdfBkg, nBkg, "window");


  RooRealVar mean("mean", "", 135., 120., 180.);
  RooRealVar sigma("sigma", "", 16, 10., 35.);
  RooGaussian pdfSig("pdfSig", "signal p.d.f.", x, mean, sigma);
  RooRealVar nSig("nSig", "nSignal", integralAll, 1., integralAll);
  RooExtendPdf epdfSig("epdfSig", "epdfSig", pdfSig, nSig, "window");

  RooAddPdf apdfAll("apdfAll", "apdfAll", RooArgList(epdfSig, epdfBkg), RooArgList(nSig, nBkg));
  RooFitResult * r = apdfAll.fitTo(dhAll, Range("window"), Save());
  r->Print("v");

  Double_t nSignal = nSig.getVal();
  Double_t nBckg   = nBkg.getVal();

  double nsigVal = epdfSig.expectedEvents(x);
  double nbkgVal = epdfBkg.expectedEvents(x);

  // Double_t fracBkg = nbkgVal / (nsigVal + nbkgVal);
  s Double_t fracBkg = sum2Bkg / sum1Bkg;

  cout << "================\n";
  cout << "Sig / CB: " << nSignal << " / " << nBckg << endl;
  cout << "================\n";


  cout << "================\n";
  cout << "Sig / CB  expectedEvents: " << nsigVal << " / " << nbkgVal << endl;
  cout << "================\n";


  TH1D * massSig = (TH1D *) massAll->Clone();
  massSig->Add(massCb, -fracBkg);
  RooDataHist dhSig("dhSig", "Pion peak", RooArgList(x), massSig);
  RooPlot * xframe = x.frame();
  dhSig.plotOn(xframe, MarkerColor(kRed), LineColor(kRed));
  apdfAll.plotOn(xframe, Components("epdfSig"), LineColor(kRed), LineWidth(1),
    Normalization(dhAll.sumEntries(), RooAbsReal::NumEvent));
  xframe->Draw();


  RooPlot * pdfFrame = x.frame();
  apdfAll.plotOn(pdfFrame, Components(epdfBkg), LineColor(kBlue), LineWidth(2));
  apdfAll.plotOn(pdfFrame, LineColor(kBlack), LineWidth(2));
  apdfAll.plotOn(pdfFrame, Components(epdfSig), LineColor(kRed), LineWidth(2));

  // pdfFrame->Draw();

  Double_t S        = nSig.getVal();
  Double_t B        = nBkg.getVal();
  Double_t SErr     = nSig.getError();
  Double_t BErr     = nBkg.getError();
  Double_t meanVal  = mean.getVal();
  Double_t sigmaVal = sigma.getVal();
  Double_t Sig      = S / TMath::Sqrt(S + B);

  TLatex * latex = new TLatex();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42);
  latex->SetTextAlign(12);


  TString label1 = "#splitline{Signal (#mu #pm 2#sigma) = %s}"
    "{#splitline{Signal / Background = %.2f}{Significance = %.1f}}";

  TString label2 =
    "#splitline{Mean = %.1f}"
    "{Sigma = %.1f}";


  latex->DrawLatexNDC(0.35, 0.70, Form(label1, Form("%.2e", S), S / B, S / TMath::Sqrt(S + B)));
  latex->DrawLatexNDC(0.35, 0.45, Form(label2, meanVal, sigmaVal));
  gPad->Update();
  cout << "================\n";
  cout << Form(label1, Form("%.2e", S), S / B, S / TMath::Sqrt(S + B)) << endl;
  cout << Form(label2, meanVal, sigmaVal) << endl;

  cout << "================\n";


  TString comment;
  comment.ReadLine(cin);
  if (comment == "q") {
    cout << "QUIT\n";
    delete xframe;
    delete pdfFrame;
    return 1;
  }

  delete xframe;
  delete pdfFrame;
  return 0;
} // getSig

Int_t draw_check(TString inputfile = "../../results/parameters/photons_to_histograms/3200A_exp.root",
  TString                outputfile = "draw_check.root")
{
  TFile * in = new TFile(inputfile, "READ");

  inputfile.Replace(0, inputfile.Last('/') + 1, "");
  inputfile.ReplaceAll(".root", "");
  enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
  TH1D * Counter_cent = (TH1D *) in->Get("Counter_cent_cut");
  TFile * out         = new TFile(outputfile, "RECREATE");
  GlobalParameters par("../parameters.conf");
  Int_t MassBins = par.N_MASS_BINS;
  Int_t ptBins   = par.N_PT_BINS;
  Int_t yBins    = par.N_Y_BINS;
  // =======================================================================//
  Float_t MassMin     = 0;
  Float_t MassMax     = 800;
  Float_t massBinSize = (MassMax - MassMin) / Float_t(MassBins);

  Float_t ptMin     = 0;
  Float_t ptMax     = 1200;
  Float_t ptBinSize = (ptMax - ptMin) / Float_t(ptBins);// =50; //MeV

  Float_t yMin     = 0.05;
  Float_t yMax     = 1.05;
  Float_t yBinSize = (yMax - yMin) / Float_t(yBins);// 0.1;

  TCanvas * can = new TCanvas("can", "", 1600, 800);
  can->Divide(2, 1);
  // =======================================================================//
  for (Int_t iCent = 0; iCent < 3; iCent++) {
    TH2D * h_pt_y = new TH2D(Form("h_pt_y_cent%i",
        iCent), ";y_{cm};p_{t}[MeV/c]", yBins, yMin, yMax, ptBins, ptMin, ptMax);
    for (Int_t iY = 0; iY < yBins; iY++) {
      for (Int_t iPt = 0; iPt < ptBins; iPt++) {
        TString hName    = Form("pt%iy%icent%i", iPt, iY, iCent);
        TH1D * massAll   = (TH1D *) in->Get("All/" + hName);
        TH1D * massCb    = (TH1D *) in->Get("Mix/" + hName);
        Float_t good_bin = 1;// (0<h_map->GetBinContent(bin));
        can->cd(1);
        if (0 < getSig(massAll, massCb)) return 2;

        can->Update();
        can->cd(2);
        if (0 < getSigRooFit(massAll, massCb)) return 2;

        can->Update();
      }
    }
  }
  return 0;
} // draw_check
