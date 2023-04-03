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
const Float_t maxMass = 350;

// ================================================//================================================

Double_t getSigRooFit(TH1D * massAll, TH1D * massCb, Double_t & nSigExtracted, Double_t & nCbExtracted)
{
	massAll->GetXaxis()->SetRangeUser(minMass, maxMass);
	massCb->GetXaxis()->SetRangeUser(minMass, maxMass);
	Double_t integralAll = massAll->Integral();


	RooRealVar x("x", "M_{#gamma#gamma}[MeV/c^{2}]", minMass, maxMass);
	x.setRange("window", 50, 250);

	RooDataHist dhAll("dhAll", "mass spectra of photon pairs in the same event", RooArgList(x), massAll);


	RooDataHist dhBkg("dhBkg", "mass spectra of photon pairs in mixed events", RooArgList(x), massCb);
	Double_t sum1Bkg = dhBkg.sum(0);
	RooHistPdf pdfBkg("pdfBkg", "background PDF", RooArgList(x), dhBkg);
	RooRealVar nBkg("nBkg", "nbackground", integralAll, 1., 2. * integralAll);
	RooExtendPdf epdfBkg("epdfBkg", "epdfBkg", pdfBkg, nBkg, "window");


	RooRealVar mean("mean", "", 135., 120., 180.);
	RooRealVar sigma("sigma", "", 16, 10., 30.);
	RooGaussian pdfSig("pdfSig", "signal p.d.f.", x, mean, sigma);
	RooRealVar nSig("nSig", "nSignal", integralAll, 1., integralAll);
	RooExtendPdf epdfSig("epdfSig", "epdfSig", pdfSig, nSig, "window");

	RooAddPdf apdfAll("apdfAll", "apdfAll", RooArgList(epdfSig, epdfBkg), RooArgList(nSig, nBkg));
	RooFitResult * r = apdfAll.fitTo(dhAll, Range("window"), Save());

	Double_t nSignal = nSig.getVal();
	Double_t nBckg   = nBkg.getVal();

	Double_t nsigVal = epdfSig.expectedEvents(x);
	Double_t nbkgVal = epdfBkg.expectedEvents(x);

	//Double_t fracBkg = nbkgVal / (nsigVal + nbkgVal);
	Double_t fracBkg = nBckg / (nSignal + nBckg);

  nSigExtracted=nSignal;
  nCbExtracted=nBckg;

	return fracBkg;
} // getSig

Int_t draw_check2(TString inputfile = "../../results/parameters_BEST/photons_to_histograms/3200A_exp.root",
                  TString outputfile = "draw_check.root")
{
	TFile * in = new TFile(inputfile, "READ");

	inputfile.Replace(0, inputfile.Last('/') + 1, "");
	inputfile.ReplaceAll(".root", "");
	enum Counters_cent { cent0_10, cent10_20, cent20_30, cent30_40, all_cents, cNumCounters_cent };
	TH1D * Counter_cent = (TH1D *) in->Get("Counter_cent");
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
	//
	// TCanvas * can = new TCanvas("can", "", 1600, 800);
	// can->Divide(2, 1);
	// =======================================================================//

	Double_t fracInt[3];
	Double_t nSigExtractedInt[3];
	Double_t nCbExtractedInt[3];

	Double_t frac[3][yBins][ptBins];
	Double_t nSigExtracted[3][yBins][ptBins];
	Double_t nCbExtracted[3][yBins][ptBins];
	for (Int_t iCent = 0; iCent < 3; iCent++) {
		TH2D * h_pt_y = new TH2D(Form("h_pt_y_cent%i",
		                              iCent), ";y_{cm};p_{t}[MeV/c]", yBins, yMin, yMax, ptBins, ptMin, ptMax);

		TH1D * massAllIntegrated   = (TH1D *) in->Get( Form("All/cent%i", iCent));
		TH1D * massCbIntegrated    = (TH1D *) in->Get(Form("Mix/cent%i", iCent));
		fracInt[iCent]=getSigRooFit(massAllIntegrated, massCbIntegrated,nSigExtractedInt[iCent],nCbExtractedInt[iCent]);
		cout<<" ratio = " <<fracInt[iCent]<<endl;

		for (Int_t iY = 0; iY < yBins; iY++) {
			for (Int_t iPt = 0; iPt < ptBins; iPt++) {
				TString hName    = Form("pt%iy%icent%i", iPt, iY, iCent);
				TH1D * massAll   = (TH1D *) in->Get("All/" + hName);
				TH1D * massCb    = (TH1D *) in->Get("Mix/" + hName);
				Float_t good_bin = 1;// (0<h_map->GetBinContent(bin));

				frac[iCent][iY][iPt]=getSigRooFit(massAll, massCb,nSigExtracted[iCent][iY][iPt],nCbExtracted[iCent][iY][iPt]);
				cout<<" y"<<iY<<" pt"<<iPt<<" ratio = " <<frac[iCent][iY][iPt]<<endl;


				// can->Update();
			}
		}
	}





	for (Int_t iCent = 0; iCent < 3; iCent++) {

		cout<<" ratio = " <<fracInt[iCent]<<endl;
		cout<<" signal = "<<nSigExtractedInt[iCent]<<"  cb = "<<nCbExtractedInt[iCent]<<endl;
		for (Int_t iY = 0; iY < yBins; iY++) {
			for (Int_t iPt = 0; iPt < ptBins; iPt++) {
				cout<<" y"<<iY<<" pt"<<iPt<<" ratio = " <<frac[iCent][iY][iPt]<<endl;
				cout<<" signal = "<<nSigExtracted[iCent][iY][iPt]<<"  cb = "<<nCbExtracted[iCent][iY][iPt]<<endl;


			}
		}
	}







	return 0;
} // draw_check
