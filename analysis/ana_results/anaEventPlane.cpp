#include <stdio.h>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#define FACTOR 0.797884561

Double_t SphericalBesselI(Int_t order, Double_t arg) // compute half-integer modified Bessel functions
{                                                    // order >0, for order>5/2, interpolation is used (improve this by using recurrence!!!)
	if (order < 0) return 0;

	if (arg < 1e-7) return 0;

	switch (order) {
	case 0:   // 1/2
		return FACTOR * sqrt(arg) * TMath::SinH(arg) / arg;

		break;

	case 1:   // 3/2
		return FACTOR * sqrt(arg) * ( -TMath::SinH(arg) / (arg * arg) + TMath::CosH(arg) / arg );

		break;

	case 2:   // 5/2
		return FACTOR * sqrt(arg)
		       * ( (3. / (arg * arg) + 1.) * TMath::SinH(arg) / arg - 3. * TMath::CosH(arg) / (arg * arg) );

		break;

	default:
		break;
	}

	return 0.5 * (TMath::BesselI(order, arg) + TMath::BesselI(order + 1, arg)); // use average of integer-order Bessel
}

Double_t calcEventPlaneResolutionFactor(TH1F * h, Int_t n) // compute Ollitrault correction factor
{
	// Int_t i1 = h->FindBin(-90.);
	Int_t i2 = h->FindBin(90.);
	Int_t i3 = h->FindBin(180.);

	Double_t rat = ( h->Integral(i2, i3) ) / h->Integral(1, i3);

	Double_t chisq = -2. * TMath::Log(2. * rat);
	Double_t chi   = sqrt(chisq);

	Int_t n1 = (n - 1) / 2;
	Int_t n2 = (n + 1) / 2;
	Double_t sumBessel;

	if (n % 2 == 1) sumBessel = TMath::BesselI(n1, 0.5 * chisq) + TMath::BesselI(n2, 0.5 * chisq); // integer order
	// half-integer order approximated via interpolation of nearest integer-order modified Bessel functions
	// (Exact solution is not implemented in Root!), checked ok vs. Ollitrault's plot in arXiv:nucl-ex/9711003v2
	//  else sumBessel = 0.5*( TMath::BesselI(n1,0.5*chisq) + 2*TMath::BesselI(n2,0.5*chisq) + TMath::BesselI(n2+1,0.5*chisq) ) ;
	else sumBessel = SphericalBesselI(n1, 0.5 * chisq) + SphericalBesselI(n2, 0.5 * chisq); // exact solution for order 1/2, 3/2, 5/2

	Float_t fac = sqrt(TMath::Pi()) / 2. * chi * exp(-0.5 * chisq) * sumBessel;

	return fac;
}

void anaEventPlane(TString input = "eventPlane.root")
{
	TFile * in  = new TFile(input, "READ");
	TFile * out = new TFile("eventPlaneResolution.root", "RECREATE");

	TH1F * hEventPlaneCharaAll = (TH1F *) in->Get("hEventPlaneCharaAll");
	TH1F * hEventPlaneInfo     = (TH1F *) in->Get("hEventPlaneInfo");

	const Int_t nCentralityClasses = 4;
	TH1F * hDeltaABcent10[nCentralityClasses];
	TH1F * hDeltaABcent5[nCentralityClasses * 2];
	TH1F * hEventPlaneChara[nCentralityClasses];

	out->cd();
	gStyle->SetOptStat(0);
	TH1F * hResFactor[4];
	TH1F * hResFactorK5[4];
	for (Int_t n = 0; n < 3; n++) {
		hResFactor[n] =new TH1F(Form("hResFactor%i", n + 1),";Centrality [%];R_{n}", 4, 0, 40);
		hResFactorK5[n] = new TH1F(Form("hResFactor%ik5", n + 1),Form("Event Plane Resolution Factor (n=%i);Centrality [%];R_{n}", n + 1), 8, 0, 40);


		hResFactor[n]->SetMarkerStyle(20);
		hResFactor[n]->SetLineColor(n + 2);
		hResFactorK5[n]->SetMarkerStyle(24);
		hResFactorK5[n]->SetMarkerColor(n + 2);
		hResFactorK5[n]->SetLineColor(n + 2);

		for (Int_t centBin = 0; centBin < nCentralityClasses; centBin++) {
			hDeltaABcent10[centBin]   = (TH1F *) in->Get(Form("hDeltaABcent10_%d", centBin));
			hEventPlaneChara[centBin] = (TH1F *) in->Get(Form("hEventPlaneChara%d", centBin));

			hResFactor[n]->SetBinContent(centBin + 1,calcEventPlaneResolutionFactor(hDeltaABcent10[centBin], n + 1));
		}


		for (Int_t centBin = 0; centBin < nCentralityClasses * 2; centBin++) {
			hDeltaABcent5[centBin] = (TH1F *) in->Get(Form("hDeltaABcent5_%d", centBin));
			hResFactorK5[n]->SetBinContent(centBin + 1,calcEventPlaneResolutionFactor(hDeltaABcent5[centBin], n + 1));
		}
	}


	TCanvas * can = new TCanvas("can", "can", 800, 600);
	can->cd();
	TLegend * leg = new TLegend(0.13, 0.66, 0.28, 0.88);
	for (Int_t n = 0; n < 3; n++) {
		leg->AddEntry(hResFactor[n], Form("n=%i", n + 1), "l");
		hResFactor[n]->GetYaxis()->SetRangeUser(0, 1);
		hResFactor[n]->Draw("hist text90 same");
		hResFactorK5[n]->Draw("P same");
	}
	leg->Draw("same");
	can->SaveAs("eventPlaneResolution.pdf");

	TCanvas * can2 = new TCanvas("can2", "can2", 800, 600);
	can2->cd();

	Double_t ymax  = hEventPlaneChara[0]->GetMaximum();
	TLegend * leg2 = new TLegend(0.67, 0.67, 0.89, 0.89);
	for (Int_t centBin = 0; centBin < nCentralityClasses; centBin++) {
		hEventPlaneChara[centBin]->GetYaxis()->SetTitleOffset(1.5);

		leg2->AddEntry(hEventPlaneChara[centBin], Form("%i0-%i0 %%", centBin, centBin + 1), "l");
		Double_t scale = ymax / hEventPlaneChara[centBin]->GetMaximum();
		hEventPlaneChara[centBin]->Scale(scale);
		hEventPlaneChara[centBin]->GetYaxis()->SetRangeUser(45E5, 49E5);
		hEventPlaneChara[centBin]->SetLineColor(centBin + 6);
		hEventPlaneChara[centBin]->Draw("hist same");
	}
	leg2->Draw("same");
	can2->SaveAs("eventPlanes.pdf");

	out->cd();
	out->GetList()->Write();
	out->Save();
} // anaEventPlane
