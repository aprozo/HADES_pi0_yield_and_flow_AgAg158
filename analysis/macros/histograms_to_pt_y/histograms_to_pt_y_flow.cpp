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


Float_t getNorm(TH1 * All, TH1 * Cb)
{
  Int_t below_min = 80;
  Int_t below_max = 90;
  Int_t above_min = 200;
  Int_t above_max = 230;

  Int_t BinBelowMin     = All->FindBin(below_min);
  Int_t BinBelowMax     = All->FindBin(below_max);
  Int_t BinAboveMin     = All->FindBin(above_min);
  Int_t BinAboveMax     = All->FindBin(above_max);
  Float_t int_below_mix = Cb->Integral(BinBelowMin, BinBelowMax);
  Float_t int_below     = All->Integral(BinBelowMin, BinBelowMax);
  Float_t int_above_mix = Cb->Integral(BinAboveMin, BinAboveMax);
  Float_t int_above     = All->Integral(BinAboveMin, BinAboveMax);
  // Float_t norm = 0.3 * int_below / int_below_mix + 0.7 * int_above / int_above_mix;
  Float_t norm = int_above / int_above_mix;

  return norm;
}

const Float_t massPi0 = 134.98;


// ======================================================================================//
// ==========================FLOW========================================================//
// ======================================================================================//


const Int_t nHarmonics       = 2;
const Double_t thresholdGPSB = 0.1; // threshold of siginificance for considering the phase space bin as good (enough entries to do flow analysis)

enum Method { standard, fitTP, quadrant, nMethods };
inline const TString ToString(Method m)
{
  switch (m) {
      case standard:   return "standard";

      case fitTP:   return "fitTP";

      case quadrant: return "quadrant";

      default:      return "[Unknown]";
  }
}

Double_t flowFourier(Double_t * x_val, Double_t * par)
{
  Double_t x    = x_val[0] * TMath::DegToRad();
  Double_t C    = par[0];
  Double_t v1   = par[1];
  Double_t v2   = par[2];
  Double_t func = C * (1 + 2 * (v1 * cos(x) + v2 * cos(2 * x)));

  return func;
}

TH1D * integrateHists(TH1 ** h, Int_t nBins)
{
  TH1D * Sum = (TH1D *) h[0]->Clone("Integrated" + (TString) h[0]->GetName());

  for (Int_t iBin = 1; iBin < nBins; iBin++) Sum->Add(h[iBin]);
  return Sum;
}

TH1F * hSigBckg;
TH1F * hBckg;
Double_t funcTProfFit(Double_t * x, Double_t * par)
{
  Int_t iBin     = hSigBckg->GetXaxis()->FindBin(x[0]);
  Double_t ratio = hSigBckg->GetBinContent(iBin); // = S/(S+B)
  // Double_t Cb = hBckg->GetBinContent(iBin);

  Double_t value = ratio * par[0] + (1. - ratio)
    * (par[1] + par[2] * x[0] + par[3] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0]); // par[0] = <cos(XX)>, par[1] and par[2] is bckg

  return value;
}

class Flow {
public:

  Flow()
  {
    Geant = NULL;

    Float_t sigma = 16.4; // 15.2 for sim
    Float_t mean  = 135;

    min = mean - 2 * sigma;
    max = mean + 2 * sigma;
  };
  ~Flow()
  {
    // delete Sig;
  };
  // For getting borders of integration//
  void setIntegratedAllCb(TH1F * AllSum, TH1F * CbSum)
  {
    normSum = getNorm(AllSum, CbSum);
    TH1F * SigAll = (TH1F *) AllSum->Clone("sig" + (TString) AllSum->GetName());
    SigAll->Add(CbSum, -normSum);
    TF1 * fit = new TF1("fit", "gaus", 100, 180);
    fit->SetParameters(SigAll->GetMaximum(), 135, 16);
    SigAll->Fit(fit, "RQNEM");
    Float_t mean  = fit->GetParameter(1);
    Float_t sigma = fit->GetParameter(2);
    min = mean - 2 * sigma;
    max = mean + 2 * sigma;

    minBin = SigAll->FindBin(min);
    maxBin = SigAll->FindBin(max);


    delete SigAll;
    delete fit;
  };


  // For getting number of nSignal//
  void setAllCb(TH1F * inAll, TH1F * inCb)
  {
    All = inAll;
    Cb  = inCb;
    Cb->Scale(normSum); // getNorm(All, Cb));

    Sig = (TH1F *) All->Clone("sig" + (TString) All->GetName());
    Sig->Add(Cb, -1);

    TF1 * fit = new TF1("fit", "gaus+pol3(3)", 60, 240);
    fit->SetParameters(Sig->GetMaximum(), 135, 16, 0, 0, 0, 0);
    fit->SetParLimits(1, 125., 155.);
    fit->SetParLimits(2, 8., 25.);

    Sig->Fit(fit, "RQN+");
    TF1 * correction = new TF1("correction", "pol3", 60, 240);
    correction->SetParameters(fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6));
    Sig->Add(correction, -1.);

    //	Int_t minBin = Sig->FindBin(min);
    //	Int_t maxBin = Sig->FindBin(max);

    Double_t err_nAll, err_nCb;


    Double_t nAll = All->IntegralAndError(minBin, maxBin, err_nAll, "");
    Double_t nCb  = Cb->IntegralAndError(minBin, maxBin, err_nCb, "");

    //  if (Geant != NULL) nSig = Geant->IntegralAndError(minBin, maxBin, err_nSig, "");

    // if (Geant != NULL) All->Add(Geant, -1.);
    // nSig     = All->IntegralAndError(minBin, maxBin, err_nAll, ""); // background flow
    // err_nSig = err_nAll;

    //  nSig     = nGeant; // background flow
    //  err_nSig = err_nGeant;

    if (Geant != NULL) {
      nSig     = Geant->Integral(minBin, maxBin);
      err_nSig = TMath::Sqrt(nSig);
    }

    // if (Geant != NULL) nSig = Geant->IntegralAndError(minBin, maxBin, err_nSig, "");


    nSig     = Sig->IntegralAndError(minBin, maxBin, err_nSig, "");
    err_nSig = TMath::Sqrt(err_nAll * err_nAll + err_nCb * err_nCb);


    delete fit;
    delete correction;
  }; // setAllCb


  void fill(TH1D * h[nHarmonics], const Int_t & bin)
  {
    h[0]->SetBinContent(bin + 1, v1);
    h[0]->SetBinError(bin + 1, err_v1);

    h[1]->SetBinContent(bin + 1, v2);
    h[1]->SetBinError(bin + 1, err_v2);
  };


  void fill(TH2D * h[nHarmonics], const Int_t & binX, const Int_t & binY)
  {
    h[0]->SetBinContent(binX + 1, binY + 1, v1);
    h[0]->SetBinError(binX + 1, binY + 1, err_v1);
    h[1]->SetBinContent(binX + 1, binY + 1, v2);
    h[1]->SetBinError(binX + 1, binY + 1, err_v2);
  };


  virtual void evaluate(){ };


  void setGeant(TH1F * in)
  {
    Geant      = in;
    nGeant     = Geant->Integral(min, max);
    err_nGeant = sqrt(nGeant);
  };


  void drawExtraction(TVirtualPad * pad)
  {
    Double_t maxAll;

    pad->SetGrid(1, 1);
    if (SigToCb < 15) maxAll = Sig->GetMaximum();
    else maxAll = Sig->GetMaximum();
    All->GetXaxis()->SetRangeUser(60, 200);
    All->GetYaxis()->SetRangeUser(-0.1 * maxAll, maxAll * 1.2);

    All->SetLineColor(1);
    All->SetLineWidth(3);

    Cb->SetLineColor(4);
    Cb->SetLineWidth(3);
    All->SetTitle("");
    All->DrawClone();
    Cb->DrawClone("same");
    Sig->SetLineWidth(3);
    Sig->SetLineColor(kRed);
    Sig->DrawClone("same");


    TLine * line = new TLine(Sig->GetBinCenter(min), -0.05 * maxAll, Sig->GetBinCenter(min), 1.2 * maxAll);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->DrawClone(" same ");

    TLine * line2 = new TLine(Sig->GetBinCenter(max), -0.05 * maxAll, Sig->GetBinCenter(max), 1.2 * maxAll);
    line2->SetLineColor(kBlack);
    line2->SetLineWidth(2);
    line2->DrawClone(" same ");

    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.07);
    tex.SetLineWidth(3);

    // if (Geant != NULL) {
    //  Geant->SetLineColor(kAzure + 1);
    //  Geant->SetLineWidth(3);
    //  Geant->DrawClone("same");
    // }
    delete line;
    delete line2;
  }; // drawExtraction

  void correctFwResolution(Double_t correction[nHarmonics])
  {
    v1     /= correction[0];
    v2     /= correction[1];
    err_v1 /= correction[0];
    err_v2 /= correction[1];
  };


  inline const Float_t getV1(void) const
  {
    return v1;
  }

  inline const Float_t getV2(void) const
  {
    return v2;
  }

protected:
  Float_t min, max;
  Int_t maxBin, minBin;
  Float_t SigToCb, normSum;
  Double_t v1, v2, err_v1, err_v2, nSig, err_nSig, nGeant, err_nGeant;
  TH1F * All, * Cb, * Sig;
  TH1F * Geant;
  TString name;
};


class StandardFlow : public Flow {
public:
  StandardFlow(TString nameHists, TString title)
  {
    hFlow = new TH1F("standardFlow" + nameHists, title + ";#Delta#phi [#circ];  #frac{1}{N} #frac{dN}{d#Delta#phi}",
        par.N_PHI_BINS, -180, 180);
    flowFit = new TF1("standardFlowFit" + nameHists, flowFourier, -180, 180, 3);
    hFlow->SetMarkerStyle(20);
    hFlow->SetMarkerColor(kBlack);
    hFlow->SetMarkerSize(1);
    hFlow->GetYaxis()->SetRangeUser(0.05, 0.15);
    Geant   = NULL;
    evtFlow = NULL;

    hFlowGeant = new TH1F("standardFlowGeant" + nameHists, "", par.N_PHI_BINS, -180, 180);
    hFlowGeant->SetMarkerStyle(21);
    hFlowGeant->SetMarkerColor(kAzure + 1);
    hFlowGeant->SetMarkerSize(2);
  };
  ~StandardFlow()
  {
    delete hFlow;
    delete flowFit;
    delete hFlowGeant;
  };
  void fillPhiBin(const Int_t & iPhi)
  {
    hFlow->SetBinContent(iPhi + 1, nSig);
    hFlow->SetBinError(iPhi + 1, err_nSig);
    if (Geant != NULL) {
      hFlowGeant->SetBinContent(iPhi + 1, nGeant);
      hFlowGeant->SetBinError(iPhi + 1, err_nGeant);
    }
  };
  void evaluate()
  {
    if (Geant != NULL) {
      hFlowGeant->Scale(1. / hFlowGeant->Integral());
      hFlowGeant->Fit(flowFit, "RQEMS+");
      geant_v1     = flowFit->GetParameter(1);
      geant_err_v1 = flowFit->GetParError(1);
      geant_v2     = flowFit->GetParameter(2);
      geant_err_v2 = flowFit->GetParError(2);
    }

    hFlow->Scale(1. / hFlow->Integral());
    TFitResultPtr r1 = hFlow->Fit(flowFit, "RQEMS+");
    v1      = flowFit->GetParameter(1);
    err_v1  = flowFit->GetParError(1);
    v2      = flowFit->GetParameter(2);
    err_v2  = flowFit->GetParError(2);
    chi2Ndf = r1->Chi2() / r1->Ndf();
  };

  void fillFromGeant(TH1F * h[nHarmonics], const Int_t & bin)
  {
    h[0]->SetBinContent(bin + 1, geant_v1);
    h[0]->SetBinError(bin + 1, geant_err_v1);
    h[1]->SetBinContent(bin + 1, geant_v2);
    h[1]->SetBinError(bin + 1, geant_err_v2);
  };

  void fillFromGeant(TH2D * h[nHarmonics], const Int_t & binX, const Int_t & binY)
  {
    h[0]->SetBinContent(binX + 1, binY + 1, geant_v1);
    h[0]->SetBinError(binX + 1, binY + 1, geant_err_v1);
    h[1]->SetBinContent(binX + 1, binY + 1, geant_v2);
    h[1]->SetBinError(binX + 1, binY + 1, geant_err_v2);
  };

  void setEvt(TH1D * in)
  {
    evtFlow = in;
    evtFlow->Scale(1. / evtFlow->Integral());
    evtFlow->SetLineColor(kBlue);
    evtFlow->SetLineWidth(2);
  };


  void draw(TVirtualPad * pad)
  {
    pad->SetGrid(1, 1);
    hFlow->DrawClone();
    flowFit->DrawClone("same");
    if (Geant != NULL) hFlowGeant->DrawClone("same");
    if (evtFlow != NULL) evtFlow->DrawClone("hist same");

    TLatex * latex = new TLatex();
    latex->SetTextFont(42);
    latex->SetTextAlign(12);
    latex->SetTextSize(0.045);
    latex->DrawLatexNDC(0.43, 0.85, Form("v_{1} = %.2f v_{2} =%.2f", v1, v2));

    auto legend = new TLegend(0.45, 0.69, 0.75, 0.95);
    legend->SetBorderSize(0);

    // if (evtFlow != NULL) legend->AddEntry(evtFlow, "UrQMD", "l");
    if (Geant != NULL) legend->AddEntry(hFlowGeant, "SimulationInfo", "lep");
    legend->AddEntry(hFlow, "Reconstructed", "lep");
    legend->AddEntry(flowFit, "Fit", "l");
  };


  TH1F * getFlowHist()
  {
    return hFlow;
  };


private:
  TH1F * hFlow;
  TF1 * flowFit;
  TH1F * hFlowGeant;
  TH1D * evtFlow;
  Double_t geant_v1, geant_v2, geant_err_v1, geant_err_v2, chi2Ndf;
};


class QuadrantFlow : public Flow {
public:
  void fillQuadPhiBin(const Int_t & iQuadPhi)
  {
    N[iQuadPhi]     = nSig;
    err_N[iQuadPhi] = err_nSig;
  }

  void evaluate()
  {
    Double_t N_0       = N[0];
    Double_t N_90      = N[1];
    Double_t N_180     = N[2];
    Double_t N_270     = N[3];
    Double_t err_N_0   = err_N[0];
    Double_t err_N_90  = err_N[1];
    Double_t err_N_180 = err_N[2];
    Double_t err_N_270 = err_N[3];

    // possible problematic cases:
    if (N_180 == 0) N_180 = 1.;

    Double_t Rin     = N_0 / N_180;
    Double_t err_Rin =
      TMath::Sqrt(TMath::Power( (1. / N_180) * err_N_0,
        2.) + TMath::Power( (-1. * N_0) / (N_180 * N_180) * err_N_180, 2.) );

    // R squeeze out = ( (dN/dDphi)_90 + (dN/dDphi)_270 ) / ( (dN/dDphi)_0 + (dN/dDphi)_180 )
    Double_t Rout     = (N_90 + N_270) / (N_0 + N_180);
    Double_t err_Rout =
      TMath::Sqrt(TMath::Power( (1. / (N_0 + N_180)) * err_N_90,
        2.) + TMath::Power( (1. / (N_0 + N_180)) * err_N_270, 2.)
        + TMath::Power( (-1.) * (N_90 + N_270) / ((N_0 + N_180) * (N_0 + N_180)) * err_N_0, 2.)
        + TMath::Power( (-1.) * (N_90 + N_270) / ((N_0 + N_180) * (N_0 + N_180)) * err_N_180, 2.) );

    // calculate v1 and v2:
    // v2 = (1-Rout)/(2*(1+Rout))
    v2     = (1. - Rout) / (2. * (1. + Rout));
    err_v2 = TMath::Abs(err_Rout / ((1. + Rout) * (1. + Rout)) );

    // v1 = (1+2*v2)*(Rin-1)/(2*(1+Rin))
    v1     = (1. + 2. * v2) * (Rin - 1.) / (2. * (1. + Rin));
    err_v1 = TMath::Sqrt(TMath::Power( (Rin - 1.) / (1. + Rin) * err_v2, 2.)
        + TMath::Power( (1. + 2. * v2) / ((1. + Rin) * (1. + Rin)) * err_Rin, 2.) );
  }; /* evaluate */
private:
  Double_t N[4], err_N[4];
};


class FitTPFlow : public Flow {
public:
  FitTPFlow()
  { };

  ~FitTPFlow()
  {
    delete profileFit;
    delete profileFit2;
  };


  void set(TProfile * cosDphi, TProfile * Cos2Dphi)
  {
    profCosDphi  = cosDphi;
    profCos2Dphi = Cos2Dphi;
    profileFit   = new TF1("profFit" + (TString) profCosDphi->GetName(), funcTProfFit, 80, 250, 5);
    profileFit2  = new TF1("profFit2" + (TString) profCos2Dphi->GetName(), funcTProfFit, 80, 250, 5);
    // profileFit->SetParameters(0.01,0,0,0,0);
    profileFit->SetParLimits(0, -0.2, 0.2);
    profCosDphi->Fit(profileFit, "QERM+");
    v1     = profileFit->GetParameter(0);
    err_v1 = profileFit->GetParError(0);
    profCos2Dphi->Fit(profileFit2, "QERM+");
    v2     = profileFit2->GetParameter(0);
    err_v2 = profileFit2->GetParError(0);
  };

  TProfile * getFlowHist()
  {
    return profCosDphi;
  };
  TProfile * getFlowHist2()
  {
    return profCos2Dphi;
  };


  void draw(TVirtualPad * pad)
  {
    pad->Clear();

    pad->Divide(1, 2);
    pad->cd(1);
    profCosDphi->GetXaxis()->SetRangeUser(60, 300);
    profCosDphi->DrawClone();
    profileFit->DrawClone("same");
    pad->cd(2);
    profCos2Dphi->GetXaxis()->SetRangeUser(60, 300);
    profCos2Dphi->DrawClone();
    profileFit2->DrawClone("same");

    pad->cd(0);

    TLatex * latex = new TLatex();

    latex->SetTextFont(42);
    latex->SetTextAlign(12);
    latex->SetTextSize(0.045);
    latex->DrawLatexNDC(0.43, 0.85, Form("v_{1} = %.2f v_{2} =%.2f", v1, v2));
    pad->Modified();
    pad->Update();
  };


private:
  TProfile * profCosDphi;
  TProfile * profCos2Dphi;
  TF1 * profileFit;
  TF1 * profileFit2;
};


void makeFlowAnalysis(TString inputfile, TString outputfile)
{
  TFile * in   = new TFile(inputfile, "read");
  TFile * out  = new TFile(outputfile + ".root", "recreate");
  Bool_t isSim = kFALSE;

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.2);
  TH1::SetDefaultSumw2();
  if (inputfile.Contains("sim")) { isSim = kTRUE; std::cout << "SIMULATION" << '\n'; }

  // ======================================================================================//
  // ==========================FLOW HISTOGRAMS and SETUP===================================//
  // ======================================================================================//

  Int_t centBins = 3;

  Int_t phiBins      = par.N_PHI_BINS;
  Float_t phiBinSize = 360 / Float_t(phiBins);

  Int_t ptBins      = par.N_PT_FLOW_BINS;
  Float_t ptMin     = 100;
  Float_t ptMax     = 900;
  Float_t ptBinSize = (ptMax - ptMin) / Float_t(ptBins);

  Int_t yBins      = par.N_Y_FLOW_BINS;
  Float_t yMin     = 0.1;
  Float_t yMax     = 0.9;
  Float_t yBinSize = (yMax - yMin) / Float_t(yBins);


  TString methodName[3] = { "standard", "quadrant", "fitTP" };

  TH2D * hFlow[centBins][nMethods][nHarmonics];
  TH1D * hFlowIntegratedPt[centBins][nMethods][nHarmonics];

  TH1D * hFlowVsPt[centBins][yBins][nMethods][nHarmonics], * hFlowVsY[centBins][ptBins][nMethods][nHarmonics];


  for (Int_t iCent = 0; iCent < centBins; iCent++) {
    for (int iMethod = 0; iMethod < nMethods; iMethod++) {
      for (int iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++) {
        hFlow[iCent][iMethod][iHarmonic] =
          new TH2D(Form("cent%imethod%iharmonic%i", iCent, iMethod, iHarmonic), methodName[iMethod]
            + Form(" %i0-%i0%% ;y_{cm};p_{t}, [MeV/c];v%i", iCent, iCent + 1,
            iHarmonic + 1), yBins, yMin, yMax, ptBins, ptMin, ptMax);
        hFlowIntegratedPt[iCent][iMethod][iHarmonic] =
          new TH1D(Form("cent%imethod%iharmonic%iIntegratedPt", iCent, iMethod, iHarmonic), methodName[iMethod]
            + Form(" %i0-%i0%% ;y_{cm};v%i", iCent, iCent + 1, iHarmonic + 1), yBins, yMin, yMax);


        for (Int_t iY = 0; iY < yBins; iY++) {
          hFlowVsPt[iCent][iY][iMethod][iHarmonic] =
            new TH1D(Form("cent%imethod%iharmonic%iy%i", iCent, iMethod, iHarmonic, iY),
              Form(";p_{t}, [MeV/c];v_{%i}", iHarmonic + 1), ptBins, ptMin, ptMax);
        }

        for (Int_t iPt = 0; iPt < ptBins; iPt++) {
          hFlowVsY[iCent][iPt][iMethod][iHarmonic] =
            new TH1D(Form("cent%imethod%iharmonic%ipt%i", iCent, iMethod, iHarmonic, iPt),
              Form(";y ;v_{%i}", iHarmonic + 1), yBins, yMin, yMax);
        }
      }
    }
  }


  // ======================================================================================//
  // ==========================Differential analysis=======================================//
  // ======================================================================================//


  Double_t eventPlaneResolution[centBins][nHarmonics];
  // eventPlaneResolution[0][0] = 0.559; // 0-10 % centrality v1
  // eventPlaneResolution[0][1] = 0.215; // 0-10 % centrality v2
  //
  // eventPlaneResolution[1][0] = 0.694; // 10-20 % centrality v1
  // eventPlaneResolution[1][1] = 0.348; // 10-20 % centrality v2
  //
  // eventPlaneResolution[2][0] = 0.748; // 20-30 % centrality v1
  // eventPlaneResolution[2][1] = 0.417; // 20-30 % centrality v2



    //marvins factors
    eventPlaneResolution[0][0] = 0.54022; // 0-10 % centrality v1
    eventPlaneResolution[0][1] = 0.199263; // 0-10 % centrality v2

    eventPlaneResolution[1][0] = 0.671021; // 10-20 % centrality v1
    eventPlaneResolution[1][1] = 0.322342; // 10-20 % centrality v2

    eventPlaneResolution[2][0] = 0.728786; // 20-30 % centrality v1
    eventPlaneResolution[2][1] = 0.391041; // 20-30 % centrality v2


  TCanvas * standardFlowCanvas       = new TCanvas("standardFlowCanvas", "standardFlowCanvas", 800, 600);
  TCanvas * allPtYBinsCanvas         = new TCanvas("allPtYBinsCanvas", "allPtYBinsCanvas", 800, 600);
  TCanvas * allPtYBinsCanvasProfile  = new TCanvas("allPtYBinsCanvasProfile", "allPtYBinsCanvasProfile", 800, 600);
  TCanvas * allPtYBinsCanvasProfile2 = new TCanvas("allPtYBinsCanvasProfile2", "allPtYBinsCanvasProfile2", 800, 600);
  standardFlowCanvas->SaveAs("check.pdf[");
  TCanvas * checkNdptdy[nMethods];
  for (int iMethod = 0; iMethod < nMethods; iMethod++) {
    checkNdptdy[iMethod] = new TCanvas(Form("method%i", iMethod), "", 1600, 1200);
  }
  checkNdptdy[ (int) standard]->Divide(4, 3, 0, 0);
  allPtYBinsCanvas->Divide(ptBins, yBins, 0, 0);
  allPtYBinsCanvasProfile->Divide(ptBins, yBins, 0, 0);
  allPtYBinsCanvasProfile2->Divide(ptBins, yBins, 0, 0);
  TLatex * latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextAlign(12);
  latex->SetTextSize(0.045);
  allPtYBinsCanvas->SaveAs("checkStandardFlowAllPtY.pdf[");
  allPtYBinsCanvasProfile->SaveAs("checkTProfileFlowAllPtY.pdf[");
  allPtYBinsCanvasProfile2->SaveAs("checkTProfile2FlowAllPtY.pdf[");

  for (Int_t cent = 0; cent < centBins; cent++) {
    out->cd();
    TDirectory * centDir = out->mkdir(Form("cent%i", cent));
    centDir->cd();
    for (Int_t y = 0; y < yBins; y++) {
      TH1F * AllIntegratedPt[phiBins];
      TH1F * CbIntegratedPt[phiBins];

      TH1F * quadAllIntegratedPt[4];
      TH1F * quadCbIntegratedPt[4];

      TProfile * profCosDphiIntegratedPt;
      TProfile * profCos2DphiIntegratedPt;

      TH1F * AllIntegratedPhiIntegratedPt;
      TH1F * CbIntegratedPhiIntegratedPt;

      TString centTitle = Form("  %i0%%-%i0%% ", cent, (cent + 1));
      TString yTitle    = Form(" %.1f < y_{cm} < %.1f ", (y * yBinSize + yMin), ((y + 1) * yBinSize + yMin));

      for (Int_t pt = 0; pt < ptBins; pt++) {
        // Bool_t isGoodPhaseSpaceBin = false;
        // check if phase space have enough pions to do flow analysis
        // if( nSig/TMath::Sqrt(nSig+nCb) > thresholdGPSB ) isGoodPhaseSpaceBin = true;
        // if(!isGoodPhaseSpaceBin) continue;
        //  if (y == yBins - 1) continue;


        TString ptTitle = Form(" %.0f < pt < %.0f MeV/c", (pt * ptBinSize + ptMin), ((pt + 1) * ptBinSize + ptMin));
        TString hName   = Form("pt%iy%icent%i", pt, y, cent);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ==================================== Integration of histograms===========================================//
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////


        TH1F * AllIntegratedPhi = (TH1F *) in->Get((TString) "All/" + hName);
        TH1F * CbIntegratedPhi  = (TH1F *) in->Get((TString) "Mix/" + hName);

        Float_t normSum = getNorm(AllIntegratedPhi, CbIntegratedPhi);
        hSigBckg = (TH1F *) AllIntegratedPhi->Clone((TString) AllIntegratedPhi->GetName() + "hSigToBckg");
        hSigBckg->Add(CbIntegratedPhi, -normSum);
        hSigBckg->Divide(AllIntegratedPhi);

        hBckg = (TH1F *) CbIntegratedPhi->Clone((TString) CbIntegratedPhi->GetName() + "hBckg");

        hBckg->Scale(normSum);

        // integrate histograms in pt
        for (Int_t phi = 0; phi < phiBins; phi++) {
          TH1F * All = (TH1F *) in->Get(Form("All/phi%i", phi) + hName);
          TH1F * Cb  = (TH1F *) in->Get(Form("Mix/phi%i", phi) + hName);
          if (pt == 0) {
            AllIntegratedPt[phi] = (TH1F *) All->Clone(Form("Allphi%iy%icent%i", phi, y, cent));
            CbIntegratedPt[phi]  = (TH1F *) Cb->Clone(Form("Cbphi%iy%icent%i", phi, y, cent));
          } else {
            AllIntegratedPt[phi]->Add(All);
            CbIntegratedPt[phi]->Add(Cb);
          }
        }


        if (pt == 0) {
          AllIntegratedPhiIntegratedPt =
            (TH1F *) AllIntegratedPhi->Clone(Form("AllIntegratedPhiIntegratedPty%icent%i", y, cent));
          CbIntegratedPhiIntegratedPt =
            (TH1F *) CbIntegratedPhi->Clone(Form("CbIntegratedPhiIntegratedPty%icent%i", y, cent));
        } else {
          AllIntegratedPhiIntegratedPt->Add(AllIntegratedPhi);
          CbIntegratedPhiIntegratedPt->Add(CbIntegratedPhi);
        }

        for (Int_t quadPhi = 0; quadPhi < 4; quadPhi++) {
          TH1F * quadAll = (TH1F *) in->Get(Form("All/qphi%i", quadPhi) + hName);
          TH1F * quadCb  = (TH1F *) in->Get(Form("Mix/qphi%i", quadPhi) + hName);
          if (pt ==
            0)
            quadAllIntegratedPt[quadPhi] = (TH1F *) quadAll->Clone(Form("Allqphi%iy%icent%i", quadPhi, y, cent));
          else quadAllIntegratedPt[quadPhi]->Add(quadAll);
          if (pt ==
            0) quadCbIntegratedPt[quadPhi] = (TH1F *) quadCb->Clone(Form("Cbqphi%iy%icent%i", quadPhi, y, cent));
          else quadCbIntegratedPt[quadPhi]->Add(quadCb);
        }


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ==================================== STANDART METHOD=====================================================//
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // if (y == (yBins - 1) && pt == (ptBins - 1)) continue;

        cout << "StandardFlow" << endl;
        StandardFlow standardFlow(hName, centTitle + yTitle + ptTitle);
        standardFlow.setIntegratedAllCb(AllIntegratedPhi, CbIntegratedPhi);
        for (Int_t phi = 0; phi < phiBins; phi++) {
          TH1F * All = (TH1F *) in->Get(Form("All/phi%i", phi) + hName);
          TH1F * Cb  = (TH1F *) in->Get(Form("Mix/phi%i", phi) + hName);
          standardFlow.setAllCb(All, Cb);
          if (isSim) {
            TH1F * Geant = (TH1F *) in->Get(Form("All/phi%i", phi) + hName + "geant");
            standardFlow.setGeant(Geant);
          }
          standardFlow.fillPhiBin(phi);
          checkNdptdy[ (int) standard]->cd(phi + 1);
          standardFlow.drawExtraction(gPad);
        }
        checkNdptdy[ (int) standard]->SaveAs("check.pdf");

        standardFlow.evaluate();
        standardFlowCanvas->cd();
        standardFlow.draw(gPad);
        standardFlowCanvas->SaveAs("check.pdf");

        if (!isSim) standardFlow.correctFwResolution(eventPlaneResolution[cent]);
        standardFlow.fill(hFlow[cent][(int) standard], y, pt);
        standardFlow.fill(hFlowVsY[cent][pt][(int) standard], y);
        standardFlow.fill(hFlowVsPt[cent][y][(int) standard], pt);
        TH1F * flow = (TH1F *) standardFlow.getFlowHist()->Clone("flow" + hName);
        allPtYBinsCanvas->cd(y * yBins + pt + 1);

        gPad->SetTickx(2);
        gPad->SetTicky(2);


        flow->SetMarkerStyle(1);
        flow->GetYaxis()->SetRangeUser(0.06, 0.12);
        flow->GetYaxis()->SetNdivisions(206);
        flow->GetXaxis()->SetNdivisions(207);
        flow->SetTitle("");
        flow->DrawClone();

        TLatex latex2;
        TString label3 = "#splitline{" + centTitle + "}{#splitline{" + yTitle + "}{" + ptTitle + "}}";
        latex2.DrawLatexNDC(0.6, 0.8, label3);

        // flow->Write();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ==================================== QUADRANT METHOD======================================================//
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cout << "quadrantFlow" << endl;
        QuadrantFlow quadrantFlow;
        quadrantFlow.setIntegratedAllCb(AllIntegratedPhi, CbIntegratedPhi);
        checkNdptdy[ (int) quadrant]->Clear();
        checkNdptdy[ (int) quadrant]->Divide(2, 2, 0, 0);
        for (Int_t quadPhi = 0; quadPhi < 4; quadPhi++) {
          TH1F * quadAll = (TH1F *) in->Get(Form("All/qphi%i", quadPhi) + hName);
          TH1F * quadCb  = (TH1F *) in->Get(Form("Mix/qphi%i", quadPhi) + hName);
          quadrantFlow.setAllCb(quadAll, quadCb);
          quadrantFlow.fillQuadPhiBin(quadPhi);
          checkNdptdy[ (int) quadrant]->cd(quadPhi + 1);
          quadrantFlow.drawExtraction(gPad);
        }

        checkNdptdy[ (int) quadrant]->cd(0);
        latex->DrawLatexNDC(0.43, 0.85, Form("v_{1} = %.2f v_{2} =%.2f", quadrantFlow.getV1(), quadrantFlow.getV2()));
        checkNdptdy[ (int) quadrant]->SaveAs("check.pdf");
        quadrantFlow.evaluate();
        if (!isSim) quadrantFlow.correctFwResolution(eventPlaneResolution[cent]);
        quadrantFlow.fill(hFlow[cent][(int) quadrant], y, pt);
        quadrantFlow.fill(hFlowVsY[cent][pt][(int) quadrant], y);
        quadrantFlow.fill(hFlowVsPt[cent][y][(int) quadrant], pt);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ==================================== TPROFILE METHOD======================================================//
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cout << "fitTPFlow" << endl;
        FitTPFlow fitTPFlow;
        TProfile * profCosDphi  = (TProfile *) in->Get("All/cosDphi" + hName);
        TProfile * profCos2Dphi = (TProfile *) in->Get("All/cos2Dphi" + hName);
        checkNdptdy[ (int) fitTP]->cd();
        fitTPFlow.set(profCosDphi, profCos2Dphi);
        fitTPFlow.draw(gPad);
        checkNdptdy[ (int) fitTP]->SaveAs("check.pdf");
        if (!isSim) fitTPFlow.correctFwResolution(eventPlaneResolution[cent]);
        fitTPFlow.fill(hFlow[cent][(int) fitTP], y, pt);
        fitTPFlow.fill(hFlowVsY[cent][pt][(int) fitTP], y);
        fitTPFlow.fill(hFlowVsPt[cent][y][(int) fitTP], pt);

        TProfile * profileFlow2 = (TProfile *) fitTPFlow.getFlowHist2()->Clone("tflow2");
        allPtYBinsCanvasProfile2->cd(y * yBins + pt + 1);
        gPad->SetTickx(2);
        gPad->SetTicky(2);
        profileFlow2->GetYaxis()->SetRangeUser(-0.04, 0.02);
        profileFlow2->SetTitle("");
        profileFlow2->Draw();
        latex2.DrawLatexNDC(0.6, 0.8, label3);

        TProfile * profileFlow = (TProfile *) fitTPFlow.getFlowHist()->Clone("tflow");
        allPtYBinsCanvasProfile->cd(y * yBins + pt + 1);
        gPad->SetTickx(2);
        gPad->SetTicky(2);
        profileFlow->GetYaxis()->SetRangeUser(0, 0.15);
        profileFlow->SetTitle("");
        profileFlow->Draw();
        latex2.DrawLatexNDC(0.6, 0.8, label3);


        if (pt == 0) {
          profCosDphiIntegratedPt =
            (TProfile *) profCosDphi->Clone(Form("profCosDphiIntegratedPty%icent%i", y, cent));
        } else { profCosDphiIntegratedPt->Add(profCosDphi); }
        if (pt == 0) profCos2DphiIntegratedPt =
            (TProfile *) profCos2Dphi->Clone(Form("profCos2DphiIntegratedPty%icent%i", y, cent));
        else profCos2DphiIntegratedPt->Add(profCos2Dphi);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      } // pt loop

      hSigBckg = (TH1F *) AllIntegratedPhiIntegratedPt->Clone(
        (TString) AllIntegratedPhiIntegratedPt->GetName() + "hSigToBckg");
      hSigBckg->Add(CbIntegratedPhiIntegratedPt, -1);
      hSigBckg->Divide(AllIntegratedPhiIntegratedPt);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // ==================================== Integrated PT======================================================//
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      cout << "StandardFlow Integrated" << endl;
      StandardFlow standardFlowIntegratedPt(Form("y%icent%i", y, cent), centTitle + yTitle);
      standardFlowIntegratedPt.setIntegratedAllCb(AllIntegratedPhiIntegratedPt, CbIntegratedPhiIntegratedPt);
      for (Int_t phi = 0; phi < phiBins; phi++) {
        standardFlowIntegratedPt.setAllCb(AllIntegratedPt[phi], CbIntegratedPt[phi]);
        standardFlowIntegratedPt.fillPhiBin(phi);
      }
      standardFlowIntegratedPt.evaluate();
      if (!isSim) standardFlowIntegratedPt.correctFwResolution(eventPlaneResolution[cent]);
      standardFlowIntegratedPt.fill(hFlowIntegratedPt[cent][(int) standard], y);
      TH1F * flowInt = (TH1F *) standardFlowIntegratedPt.getFlowHist()->Clone(Form("flowy%icent%i", y, cent));
      // flowInt->Write();
      // ==================================== QUADRANT METHOD======================================================//
      cout << "QuadrantFlow Integrated" << endl;
      QuadrantFlow quadrantFlowIntegratedPt;
      quadrantFlowIntegratedPt.setIntegratedAllCb(AllIntegratedPhiIntegratedPt, CbIntegratedPhiIntegratedPt);
      for (Int_t quadPhi = 0; quadPhi < 4; quadPhi++) {
        quadrantFlowIntegratedPt.setAllCb(quadAllIntegratedPt[quadPhi], quadCbIntegratedPt[quadPhi]);
        quadrantFlowIntegratedPt.fillQuadPhiBin(quadPhi);
      }
      quadrantFlowIntegratedPt.evaluate();
      quadrantFlowIntegratedPt.fill(hFlowIntegratedPt[cent][(int) quadrant], y);
      // ==================================== TPROFILE METHOD======================================================//
      cout << "fitTPFlow Integrated" << endl;
      FitTPFlow fitTPFlowIntegratedPt;
      fitTPFlowIntegratedPt.set(profCosDphiIntegratedPt, profCos2DphiIntegratedPt);
      fitTPFlowIntegratedPt.fill(hFlowIntegratedPt[cent][(int) fitTP], y);
    } // y loop
    allPtYBinsCanvas->SaveAs("checkStandardFlowAllPtY.pdf");
    allPtYBinsCanvasProfile->SaveAs("checkTProfileFlowAllPtY.pdf");
    allPtYBinsCanvasProfile2->SaveAs("checkTProfile2FlowAllPtY.pdf");
  } // cent loop

  allPtYBinsCanvas->SaveAs("checkStandardFlowAllPtY.pdf]");
  allPtYBinsCanvasProfile->SaveAs("checkTProfileFlowAllPtY.pdf]");
  allPtYBinsCanvasProfile2->SaveAs("checkTProfile2FlowAllPtY.pdf]");
  out->cd();
  out->GetList()->Write();


  standardFlowCanvas->SaveAs("check.pdf]");

  for (int iMethod = 0; iMethod < nMethods; iMethod++) {
    TFile * outMethod = new TFile(outputfile + Form("%i.root", iMethod), "recreate");
    outMethod->cd();
    for (Int_t iCent = 0; iCent < centBins; iCent++) {
      for (int iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++) {
        hFlow[iCent][iMethod][iHarmonic]->SetName(Form("cent%iharmonic%i", iCent, iHarmonic));
        hFlowIntegratedPt[iCent][iMethod][iHarmonic]->SetName(Form("cent%iharmonic%iIntegratedPt", iCent, iHarmonic));
        hFlow[iCent][iMethod][iHarmonic]->Write();
        hFlowIntegratedPt[iCent][iMethod][iHarmonic]->Write();

        for (Int_t iY = 0; iY < yBins; iY++) {
          hFlowVsPt[iCent][iY][iMethod][iHarmonic]->Write();
        }

        for (Int_t iPt = 0; iPt < ptBins; iPt++) {
          hFlowVsY[iCent][iPt][iMethod][iHarmonic]->Write();
        }
      }
    }
    outMethod->Close();
  }


  out->Save();
} // pi0_flow macro

Int_t histograms_to_pt_y_flow(
  TString input_histograms = "/u/prozorov/mar19/diphotons/gen5a/systematic_uncertainty/results/parameters_REF0/photons_to_histograms_flow/3200A_exp.root"
  // TString input_histograms = "/u/prozorov/mar19/diphotons/gen5a/systematic_uncertainty/results/parameters_REF0/photons_to_histograms_flow/sim_sim.root"
)
{
  makeFlowAnalysis(input_histograms, "flow");
  return 0;
}
