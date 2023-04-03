#include <stdio.h>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLine.h"
#include "../GlobalParameters.h"
GlobalParameters par("../parameters.conf");


using namespace std;
using namespace RooFit;


Double_t correctionBackground(Double_t * x, Double_t * par)
{
  return par[0] + x[0] * par[1] + par[2] * x[0] * x[0];
}

// Lorentzian Peak function
Double_t gaussian(Double_t * x, Double_t * par)
{
  return par[0] * exp(-0.5 * pow( (x[0] - par[1]) / par[2], 2));
}

Double_t gausPlusPol(Double_t * x, Double_t * par)
{
  return gaussian(x, par) + correctionBackground(x, &par[4]);
}

TH1D * hCb;
Double_t mixedEventBackground(Double_t * x, Double_t * par)
{ // = S/(S+B)
  return hCb->GetBinContent(hCb->GetXaxis()->FindBin(x[0])) * par[0];
}

// Sum of background and peak function
Double_t fitFunction(Double_t * x, Double_t * par)
{
  return gaussian(x, par) + mixedEventBackground(x, &par[3]);
}

Double_t totalFitFunction(Double_t * x, Double_t * par)
{
  return gaussian(x, par) + mixedEventBackground(x, &par[3]) + correctionBackground(x, &par[4]);
}

Float_t getSigNorm(TH1D * massAll, TH1D * massCb)
{
  Int_t belowMin                   = 70;
  Int_t belowMax                   = 85;
  Int_t aboveMin                   = 220;
  Int_t aboveMax                   = 250;
  Int_t binBelowMin                = massAll->FindBin(belowMin);
  Int_t binBelowMax                = massAll->FindBin(belowMax);
  Int_t binAboveMin                = massAll->FindBin(aboveMin);
  Int_t binAboveMax                = massAll->FindBin(aboveMax);
  Float_t integralBelowMix         = massCb->Integral(binBelowMin, binBelowMax);
  Float_t integralBelow            = massAll->Integral(binBelowMin, binBelowMax);
  Float_t integralAboveMix         = massCb->Integral(binAboveMin, binAboveMax);
  Float_t integralAbove            = massAll->Integral(binAboveMin, binAboveMax);
  Float_t normFactorIntegralRegion = 1.;

  normFactorIntegralRegion = integralAbove / integralAboveMix;
  // normFactorIntegralRegion = 0.3 * integralBelow / integralBelowMix + 0.7 * integralAbove / integralAboveMix;
  // return min(integralAbove / integralAboveMix,integralBelow / integralBelowMix);
  TH1D * ratio = (TH1D *) massAll->Clone("ratio");
  ratio->Divide(massCb);
  ratio->GetXaxis()->SetRangeUser(50, 250);
  normFactorIntegralRegion = ratio->GetMinimum();
  return normFactorIntegralRegion;


  // return min(integralAbove / integralAboveMix, integralBelow / integralBelowMix);
} // getSig

Int_t extractYields(TString inputfile,
  TString                   outputfile)
{
  TFile * in = new TFile(inputfile, "READ");

  inputfile.Replace(0, inputfile.Last('/') + 1, "");
  inputfile.ReplaceAll(".root", "");
  enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
  TH1D * Counter_cent = (TH1D *) in->Get("Counter_cent_cut");
  if (Counter_cent == NULL) Counter_cent = (TH1D *) in->Get("Counter_cent");
  TFile * out = new TFile(outputfile, "RECREATE");
  out->cd();

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

  Int_t nBinsScaleFactor = yBins * ptBins;

  TF1 * fit      = new TF1("fit", fitFunction, 50, 250, 4);
  TF1 * fitConst = new TF1("fitConst", "pol0");
  // =======================================================================//
  std::cout << "calculating scales" << '\n';
  TH2D * temp_yield = new TH2D("temp_yield", ";y_{cm};p_{t}[MeV/c]", yBins, yMin, yMax, ptBins, ptMin, ptMax);

  Float_t averageScale[3];
  for (Int_t iCent = 0; iCent < 3; iCent++) {
    TH1D * scaleFactorNorm = new TH1D(Form("scaleFactorNorm%i", iCent), "", nBinsScaleFactor, 0., nBinsScaleFactor);
    Float_t minScale       = 1.0;


    TH1D * integralMassAll = (TH1D *) in->Get("All/pt0y0cent0")->Clone(Form("integralAllcent%i", iCent));
    TH1D * integralMassCb  = (TH1D *) in->Get("Mix/pt0y0cent0")->Clone(Form("integralCbcent%i", iCent));

    for (Int_t iY = 0; iY < yBins; iY++) {
      for (Int_t iPt = 0; iPt < ptBins; iPt++) {
        TString hName  = Form("pt%iy%icent%i", iPt, iY, iCent);
        TH1D * massAll = (TH1D *) in->Get("All/" + hName);
        TH1D * massCb  = (TH1D *) in->Get("Mix/" + hName);
        integralMassAll->Add(massAll);
        integralMassCb->Add(massCb);
        hCb = (TH1D *) massCb->Clone("temp");
        if (massAll->Integral() < 100 || massCb->Integral() < 100) continue;

        Float_t scaleCb = getSigNorm(massAll, massCb);
        hCb->Scale(scaleCb);

        fit->SetParameters(massAll->GetMaximum() * 0.1, 135, 16, 1.0);
        fit->SetParLimits(1, 120., 180.);
        fit->SetParLimits(2, 12., 30.);

        TFitResultPtr r_mass = massAll->Fit(fit, "RQN S");
        Float_t chi2         = 0.;
        Int_t ndf = 0;
        if (r_mass.Get()) {
          chi2 = r_mass->Chi2();
          ndf  = r_mass->Ndf();
        }

        hCb->Scale(fit->GetParameter(3));
        TH1D * sig = (TH1D *) massAll->Clone(massAll->GetName() + (TString) "tempSig");
        sig->Add(hCb, -1);

        Float_t mean = fit->GetParameter(1);
        Float_t sigma = fit->GetParameter(2);
        Double_t min = mean - 2 * sigma;
        Double_t max = mean + 2 * sigma;
        Int_t minBin = sig->FindBin(min);
        Int_t maxBin = sig->FindBin(max);
        Double_t SpBErr, BErr, SpB, S, B, SErr, Sig, SigErr;
        SpB = massAll->IntegralAndError(minBin, maxBin, SpBErr, "");
        S   = sig->Integral(minBin, maxBin);
        Sig = S / TMath::Sqrt(SpB);
        Float_t rapidity = temp_yield->GetXaxis()->GetBinLowEdge(iY + 1);
        Float_t pt       = temp_yield->GetYaxis()->GetBinLowEdge(iPt + 1);

        if (Sig > 20 && mean < 150 && mean > 125 && sigma < 20 &&
          pt < 500 && pt >= 150 && rapidity >= 0.25 && rapidity < 0.55
        )
        {
          if (scaleCb < minScale) minScale = scaleCb;
          scaleFactorNorm->SetBinContent(iY * ptBins + iPt + 1, scaleCb);
          scaleFactorNorm->SetBinError(iY * ptBins + iPt + 1, SpBErr / SpB);
        }
      }
    }
    scaleFactorNorm->Fit(fitConst, "Q");
    scaleFactorNorm->Write();

    averageScale[iCent] = fitConst->GetParameter(0);
    // averageScale[iCent] = minScale;

    // delete scaleFactorNorm;
    Float_t scaleIntegral = getSigNorm(integralMassAll, integralMassCb);
    integralMassCb->Scale(scaleIntegral);
    TH1D * integralSig = (TH1D *) integralMassAll->Clone(Form("integralSigcent%i", iCent));
    integralSig->Add(integralMassCb, -1);
    integralSig->Write();
    averageScale[iCent] = scaleIntegral;
  }


  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TLatex * latex = new TLatex();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42);
  latex->SetTextAlign(12);

  TCanvas * can = new TCanvas("can", "can", 800, 1200);
  can->cd();
  Int_t nXPads = 4;
  Int_t nYPads = 4;
  can->Divide(nXPads, nYPads, 0, 0);
  can->SaveAs(outputfile + "extractionCheck.pdf[");
  Long64_t maxValLeftCanvas = 10000000000;
  TF1 * gaussianFit         = new TF1("gaussianFit", "gaus", 90, 180);
  TF1 * gausAndCorrection   = new TF1("gausAndCorrection", "gaus+pol2(3)", 70, 240);
  TF1 * correction = new TF1("correction", "pol2", 70, 240);
  for (Int_t iCent = 0; iCent < 3; iCent++) {
    TDirectory * cent = out->mkdir(Form("cent%i", iCent));
    cent->cd();
    TH2D * h_pt_y_yield = new TH2D("h_pt_y_yield", ";y_{cm};p_{t}[MeV/c]", yBins, yMin, yMax, ptBins, ptMin, ptMax);
    TH2D * h_pt_y_mean  = new TH2D("h_pt_y_mean", ";y_{cm};p_{t}[MeV/c]", yBins, yMin, yMax, ptBins, ptMin, ptMax);
    TH2D * h_pt_y_sigma = new TH2D("h_pt_y_sigma", ";y_{cm};p_{t}[MeV/c]", yBins, yMin, yMax, ptBins, ptMin, ptMax);


    for (Int_t iY = 0; iY < yBins; iY++) {
      Int_t counterCan = 1;
      can->cd(counterCan);

      for (Int_t iPt = 0; iPt < ptBins; iPt++) {
        TString hName  = Form("pt%iy%icent%i", iPt, iY, iCent);
        TH1D * massAll = (TH1D *) in->Get("All/" + hName);
        TH1D * massCb  = (TH1D *) in->Get("Mix/" + hName);
        if (massAll->Integral() < 100) continue;

        TH1D * massSig = (TH1D *) massAll->Clone(massAll->GetName() + (TString) "sig");
        massCb->Scale(averageScale[iCent]);
        massSig->Add(massCb, -1);
        // if (outputfile.Contains("sim") && iPt < 4) massSig = (TH1D *) in->Get("All/" + hName + "geant");

        gaussianFit->SetParameters(massSig->GetMaximum(), 135, 16.4, 0, 0, 0, 0);
        gaussianFit->SetParLimits(0, massSig->GetMaximum() * 0.2, massSig->GetMaximum() * 2);
        gaussianFit->SetParLimits(1, 125., 155.);
        gaussianFit->SetParLimits(2, 13., 30.);

        TFitResultPtr r_mass = massSig->Fit(gaussianFit, "RQEMSN+");
        gausAndCorrection->SetParameters(gaussianFit->GetParameter(0), gaussianFit->GetParameter(
            1), gaussianFit->GetParameter(2), 0, 0, 0, 0);

        gausAndCorrection->SetParLimits(0, massSig->GetMaximum() * 0.2, massSig->GetMaximum() * 2);
        gausAndCorrection->SetParLimits(1, 125., 155.);
        gausAndCorrection->SetParLimits(2, 13., 30.);
        massSig->Fit(gausAndCorrection, "RQN+");
        massSig->Write();
        correction->SetParameters(gausAndCorrection->GetParameter(3), gausAndCorrection->GetParameter(
            4), gausAndCorrection->GetParameter(5), gausAndCorrection->GetParameter(6));

        TH1D * massSigOriginal = (TH1D *) massSig->Clone(massAll->GetName() + (TString) "sigOriginal");
        if (outputfile.Contains("exp")) massSig->Add(correction, -1);

        Float_t chi2 = 0.;
        Int_t ndf    = 0;
        if (r_mass.Get()) {
          chi2 = r_mass->Chi2();
          ndf  = r_mass->Ndf();
        }
        Float_t gausConst = gaussianFit->GetParameter(0);
        Float_t mean      = gaussianFit->GetParameter(1);
        Float_t sigma     = gaussianFit->GetParameter(2);
        Float_t meanErr   = gaussianFit->GetParError(1);
        Float_t sigmaErr  = gaussianFit->GetParError(2);

        Float_t allSigma = 16.4; // 16.4
        Float_t massPi0  = 135;
        if (outputfile.Contains("sim")) { allSigma = 15.2; massPi0 = mean; }

        Double_t min = massPi0 - 2 * allSigma;
        Double_t max = massPi0 + 2 * allSigma;


        Int_t minBin = massSig->FindBin(min);
        Int_t maxBin = massSig->FindBin(max);


        // for (Int_t binX = 1; binX <= massSig->GetNbinsX(); binX++) {
        //   if (massSig->GetBinContent(binX) < 0)
        //     massSig->SetBinContent(binX, 0);
        // }

        Double_t SpBErr, BErr, SpB, S, S_analytical, B, SErr, Sig, SigErr, S_outside;
        SpB       = massAll->IntegralAndError(minBin, maxBin, SpBErr, "");
        S         = massSig->Integral(minBin, maxBin);
        S_outside = massSig->Integral(maxBin, 2 * maxBin - minBin);

        S_analytical = sigma * gausConst * sqrt(2 * 3.1459) / massBinSize;
        // S =S_analytical;
        B      = massCb->IntegralAndError(minBin, maxBin, BErr, "");
        SErr   = sqrt(SpBErr * SpBErr + BErr * BErr);
        Sig    = S / TMath::Sqrt(SpB);
        SigErr = TMath::Sqrt(TMath::Power((SpB + B) / 2. / TMath::Power(SpB, 1.5) * SpBErr,
            2.) + TMath::Power(BErr / TMath::Sqrt(SpB), 2.));


        if (Sig > 6 && (chi2 / ndf < 8 || Sig > 10) && S / B >= 0.01 && S_outside < S) {
          Int_t bin = h_pt_y_yield->GetBin(iY + 1, iPt + 1);
          //  if (mean > 125. && mean < 150. && sigma > 8. && sigma < 25. && S > 0. && S < SpB->Integral()) {
          h_pt_y_yield->SetBinContent(bin, S);
          h_pt_y_yield->SetBinError(bin, SErr);
          //  }

          h_pt_y_mean->SetBinContent(bin, mean);
          h_pt_y_mean->SetBinError(bin, meanErr);
          h_pt_y_sigma->SetBinContent(bin, sigma);
          h_pt_y_sigma->SetBinError(bin, sigmaErr);

          can->cd(counterCan);
          // if () {
          //  if(counterCan==nXPads*nYPads) {can->SaveAs(outputfile+"extractionCheck.pdf"); counterCan=1;}
          //  else if (counterCan%nXPads!=1) counterCan+=(nXPads-counterCan%nXPads);
          // }

          //	massAll->GetXaxis()->SetRangeUser(40, 240);
          massCb->GetXaxis()->SetRangeUser(40, 240);
          massAll->GetXaxis()->SetRangeUser(120, 150);
          Float_t ymax = massSigOriginal->GetMaximum();
          Float_t ymin = massSigOriginal->GetMinimum();
          massAll->GetXaxis()->SetRangeUser(30, 260);
          massSig->GetXaxis()->SetRangeUser(30, 300);
          //  if (counterCan % nXPads == 1) maxValLeftCanvas = ymax;
          //  else ymax = maxValLeftCanvas;

          // if (counterCan < 5) {
          //  massAll->Scale(1. / 5);
          //  massCb->Scale(1. / 5);
          // }

          massSig->GetYaxis()->SetRangeUser(1.5 * ymin, 1.5 * ymax);
          // if (counterCan < 5) massAll->GetYaxis()->SetRangeUser(-0.02 * ymax, 0.3 * ymax);
          massAll->GetYaxis()->SetTitleOffset(1.2);
          can->cd(counterCan);
          gPad->Clear();
          gPad->SetTickx(2);
          gPad->SetTicky(2);
          massAll->SetTitle(";M_{#gamma#gamma}[MeV/c^{2}];Counts/(5 MeV/c^{2})");
          massAll->SetMarkerStyle(20);
          massSig->SetTitle(";M_{#gamma#gamma}[MeV/c^{2}];Counts/(5 MeV/c^{2})");


          massSig->DrawClone("hist");
          massCb->SetLineColor(kGreen);
          massSig->SetLineColor(kBlue);
          massCb->SetFillStyle(3354);
          massCb->SetFillColorAlpha(massCb->GetLineColor(), 0.35);
          massSigOriginal->SetLineColor(kBlue + 4);
          massSigOriginal->SetLineWidth(2);
          correction->DrawClone("same");
          massSigOriginal->DrawClone("same,hist");
          // massCb->DrawClone("same,bar");
          massSig->DrawClone("same, hist");

          TLine line(30., 0., 300., 0.);
          line.DrawClone();

          TLine line2(167.8, ymin, 167.8, ymax);
          line2.DrawClone();

          TLine line3(102.2, ymin, 102.2, ymax);
          line3.DrawClone();

          gausAndCorrection->SetLineStyle(2);
          gausAndCorrection->SetLineColor(kOrange);
          gausAndCorrection->DrawClone("same");
          TString label1 =
            "#splitline{N (#mu #pm 2#sigma) = %s}"
            "{#splitline{S / B = %.2f}{Sig = %.1f}}";
          TString Signal =
            Form("%.2f x 10^{%d}", S / TMath::Power(10., TMath::Floor(TMath::Log10(S))), TMath::FloorNint(TMath::Log10(
                S)));
          TString label2 = "#splitline{Mean = %.1f MeV/c^{2}}" "{Sigma = %.1f MeV/c^{2}}";


          latex->DrawLatexNDC(0.15, 0.87, Form(label1, Signal.Data(), S / B, S / TMath::Sqrt(SpB)));
          latex->DrawLatexNDC(0.15, 0.62, Form(label2, mean, sigma));
          if (counterCan < 5) latex->DrawLatexNDC(0.15, 0.52, "#times 3");
          TString centValue = Form("%i0-%i0%%", iCent, iCent + 1);
          TString yValue    = Form("%.2f < y_{cm} < %.2f", iY * yBinSize + yMin, (iY + 1) * yBinSize + yMin);
          TString ptValue   = Form("%i < p_{t} < %i MeV/c", (int) (iPt * ptBinSize), (int) ((iPt + 1) * ptBinSize));
          TString label3    = "#splitline{" + centValue + "}{#splitline{" + yValue + "}{" + ptValue + "}}";

          // std::cout << '\n' << '\n' << '\n' << '\n';
          // std::cout << label3 << '\n';
          // std::cout << correction->GetParameter(0) << ' ' << correction->GetParameter(1) << ' '
          //           << correction->GetParameter(2) << ' ' << correction->GetParameter(3) << ' ';
          // std::cout << '\n' << '\n' << '\n' << '\n';


          latex->DrawLatexNDC(0.6, 0.8, label3);
          TLegend * legend = new TLegend(0.58, 0.60, 0.88, 0.73);
          legend->SetBorderSize(0);
          legend->AddEntry(massAll, "same event photon pairs", "l");
          legend->AddEntry(massCb, "mixed event photon pairs", "l");
          legend->AddEntry(massSig, "subtracted #pi^{0} peak", "le");
          legend->AddEntry(gaussianFit, "Gaussian fit of #pi^{0} peak", "l");
          // legend->DrawClone();
          gPad->Update();
          if (counterCan == nXPads * nYPads) { can->SaveAs(outputfile + "extractionCheck.pdf"); counterCan = 0; }
          counterCan++;
        }
      }
    }


    Double_t nEvents = Counter_cent->GetBinContent(iCent + 1);
    h_pt_y_yield->Scale(1. / (yBinSize * ptBinSize));
    h_pt_y_yield->Scale(1. / nEvents);

    // out->cd();
    h_pt_y_yield->Write();
    h_pt_y_mean->Write();
    h_pt_y_sigma->Write();
    delete h_pt_y_yield;
  }
  can->SaveAs(outputfile + "extractionCheck.pdf]");

  delete gaussianFit;
  out->Close();

  return 0;
} // extractYields

Int_t createAccEffMatricesUrqmd(TString inputfile = "", TString urqmdInputFilename = "")
{
  TFile * in         = new TFile(inputfile, "READ");
  TFile * accEffFile = new TFile("urqmdAccEff.root", "RECREATE");
  GlobalParameters par("../parameters.conf");
  TFile * urqmdInputFile = new TFile(urqmdInputFilename, "READ");

  // =======================================================================//
  for (Int_t iCent = 0; iCent < 3; iCent++) {
    TH2D * h_pt_y       = (TH2D *) in->Get(Form("cent%i/h_pt_y_yield", iCent));
    TH2D * h_urqmd_pt_y = (TH2D *) urqmdInputFile->Get(Form("hAll_cent%i", iCent));
    TH2D * hAccEff      = (TH2D *) h_pt_y->Clone(Form("hAccEff_cent%i", iCent));
    hAccEff->Divide(h_urqmd_pt_y);
    accEffFile->cd();
    hAccEff->Write();
  }
  accEffFile->Close();
  return 0;
} // createAccEffMatricesUrqmd

const Float_t massPi0 = 134.98;
Double_t ExpFitFunc_pt(Double_t * xval, Double_t * par)
{
  Double_t c   = par[0];
  Double_t T   = par[1];
  Double_t p_t = xval[0];
  Double_t M   = massPi0;
  Double_t m_t = sqrt(p_t * p_t + M * M);

  return c * p_t * m_t * exp(-m_t / T);
}

Double_t yFit(Double_t * xval, Double_t * par)
{
  Double_t C = par[0];
  Double_t T = par[1];
  Double_t y = xval[0];
  Double_t m = massPi0;

  return C * T * exp(-m * cosh(y) / T) * (m * m + 2 * m * T / cosh(y) + 2 * T * T / pow(cosh(y), 2));
}

Double_t ExpFitFunc_pt2(Double_t * xval, Double_t * par)
{
  Double_t c1 = par[0];
  Double_t T1 = par[1];
  Double_t c2 = par[3];
  Double_t T2 = par[4];

  Double_t p_t = xval[0];
  Double_t M   = massPi0;
  Double_t m_t = sqrt(p_t * p_t + M * M);

  return c1 * p_t * m_t * exp(-m_t / T1) + c2 * p_t * m_t * exp(-m_t / T2);
}

Int_t analyzeYields(TString inputfile = "", TString outputfile = "")
{
  TFile * in  = new TFile(inputfile, "READ");
  TFile * out = new TFile(outputfile, "RECREATE");
  GlobalParameters par("../parameters.conf");
  TFile * accEffFile = new TFile("urqmdAccEff.root", "READ");

  // =======================================================================//
  for (Int_t iCent = 0; iCent < 3; iCent++) {
    TH2D * h_pt_y  = (TH2D *) in->Get(Form("cent%i/h_pt_y_yield", iCent));
    TH2D * hAccEff = (TH2D *) accEffFile->Get(Form("hAccEff_cent%i", iCent));
    h_pt_y->Divide(hAccEff);
    out->cd();
    TH1D * h_y = new TH1D((Form("h_y_cent%i", iCent)), ";y_{cm}", h_pt_y->GetNbinsX(), 0.05, 1.05);
    for (Int_t iY = 1; iY <= h_pt_y->GetNbinsX(); iY++) {
      TH1D * projectionY = h_pt_y->ProjectionY(Form("y%icent%i", iY, iCent), iY, iY);

      if (iCent == 2 && iY == 1) projectionY->SetBinError(projectionY->FindBin(125), 0);
      if (iCent == 0 && iY == 2) projectionY->SetBinError(projectionY->FindBin(125), 0);

      TF1 boltzmannFit("boltzmannFit", ExpFitFunc_pt, 100, 500, 2);
      boltzmannFit.SetParameter(0, projectionY->GetMaximum());
      boltzmannFit.SetParameter(1, 80);
      boltzmannFit.SetLineStyle(2);

      projectionY->Fit("boltzmannFit", "QRNEM+", "", 100, 500);
      boltzmannFit.SetRange(0, 1200);

      // TF1 boltzmannFit("boltzmannFit", ExpFitFunc_pt, 0, 1200, 2);
      // boltzmannFit.SetParameter(0, projectionY->GetMaximum());
      // boltzmannFit.SetParameter(1, 80);
      // // projectionY->GetXaxis()->SetRangeUser(0, 400);
      // // boltzmannFit.SetRange(0, 1200);
      // projectionY->Fit("boltzmannFit", "QREM");
      // boltzmannFit.SetName(Form("boltzFitY%iCent%i", iY, iCent));
      projectionY->GetXaxis()->SetRangeUser(0, 1200);
      TH1D * projectionYWithFitPoints = (TH1D *) projectionY->Clone(Form("y%icent%iWithFitPoints", iY, iCent));
      // boltzmannFit.Write();
      Double_t dN = 0;
      for (Int_t iPt = 1; iPt <= h_pt_y->GetNbinsY(); iPt++) {
        Float_t val = projectionYWithFitPoints->GetBinContent(iPt);
        if (val == 0) {
          projectionYWithFitPoints->SetBinContent(iPt, boltzmannFit.Eval(projectionYWithFitPoints->GetBinCenter(iPt)));
        }
      }
      Double_t C = boltzmannFit.GetParameter(0);
      Double_t T = boltzmannFit.GetParameter(1);
      Double_t dNdYPointErr;
      Double_t dNdYPoint = projectionYWithFitPoints->IntegralAndError(1, h_pt_y->GetNbinsY(), dNdYPointErr, "width");

      // Double_t dN = C * exp(-massPi0/T)*T*(massPi0^2+2*massPi0*T+2*T^2);
      h_y->SetBinContent(iY, dNdYPoint);
      h_y->SetBinError(iY, dNdYPointErr);
    }
    TF1 yieldYFit(Form("y_fit_cent%d0", iCent), yFit, -3, 3, 2);
    yieldYFit.SetParameter(0, h_y->GetMaximum());
    yieldYFit.SetParameter(1, 80);
    yieldYFit.SetRange(0.15, 0.65);
    h_y->Fit(Form("y_fit_cent%d0", iCent), "QREM");
    yieldYFit.SetRange(-3, 3);
    TH1D * h_num = new TH1D((Form("h_num_cent%i", iCent)), ";", 1, 0, 1);
    h_num->SetBinContent(1, yieldYFit.Integral(-3, 3));
  }
  out->cd();
  out->GetList()->Write();

  for (Int_t iCent = 0; iCent < 3; iCent++) {
    TH2D * h_pt_y = (TH2D *) in->Get(Form("cent%i/h_pt_y_yield", iCent));
    h_pt_y->SetName(Form("h_pt_y_cent%i", iCent));
    TH2D * hAccEff = (TH2D *) accEffFile->Get(Form("hAccEff_cent%i", iCent));
    //	h_pt_y->Divide(hAccEff);
    h_pt_y->Write();
  }

  out->Save();
  out->Close();
  accEffFile->Close();

  return 0;
} // analyzeYields

// =====================


// ----------------------------------------------------------------------------------------//


Int_t histograms_to_pt_y(
  TString input_histograms_exp = "/u/prozorov/mar19/diphotons/gen5a/systematic_uncertainty/results/parameters_REF0/photons_to_histograms/3200A_exp.root",
  TString input_histograms_sim = "/u/prozorov/mar19/diphotons/gen5a/systematic_uncertainty/results/parameters_REF0/photons_to_histograms/sim_sim.root",
  TString input_urqmd = "/u/prozorov/mar19/diphotons/gen5a/systematic_uncertainty/results/parameters_REF0/urqmd_to_pt_y/urqmd_pt_y.root")
{
  extractYields(input_histograms_exp, "exp_yields.root");
  extractYields(input_histograms_sim, "sim_yields.root");
  createAccEffMatricesUrqmd("sim_yields.root", input_urqmd);
  analyzeYields("exp_yields.root", "results.root");


  return 0;
}
