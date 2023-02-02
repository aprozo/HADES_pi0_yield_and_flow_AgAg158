#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include <iostream>
#include <stdio.h>





Double_t fitEnergyCorr(Double_t* x, Double_t* par)
{
	//if(x[0]<=par[5])
	return par[0]+par[1]/sqrt(x[0]) + par[2]/x[0] +par[3]/pow(x[0],1.5)+par[4]/pow(x[0],2);
//	else return (par[0]-par[2]/par[5]-2*par[3]/pow(par[5],1.5)-3*par[4]/pow(par[5],2)) + (par[1]+2*par[2]/sqrt(par[5])+3*par[3]/par[5] +4*par[4]/pow(par[5],1.5))/sqrt(x[0]);
}




void fit()
{
	TFile * in = new TFile("checkHistogramsIteration0parN1.0.root", "read");
	gStyle->SetOptStat(0);

	TH1D * h = (TH1D*)in->Get("hEnergyCorrection");
	h->GetXaxis()->SetRangeUser(50,800);
	h->GetYaxis()->SetRangeUser(0.9,2);
	TF1 * fit = new TF1("fit", fitEnergyCorr, 80, 800,6);
	fit->SetParameters(9.93671e-01,-1.42255e+00,2.38587e+01,-1.10000e+02,1.57192e+03,500);

	// EXT PARAMETER                APPROXIMATE
	// NO.   NAME      VALUE            ERROR
	//  1  p0           9.93671e-01   2.34832e-01
	//  2  p1          -1.42255e+00   4.21771e+00
	//  3  p2           2.38587e+01   6.19008e+01
	//  4  p3          -1.10000e+02   7.70508e+02
	//  5  p4           1.57192e+03   8.70551e+03


	h->Fit(fit,"REM");
	TCanvas* can = new TCanvas("can", "", 800, 600);
	can->cd(); h->SetMarkerStyle(20);
	h->Draw("P");
  h->SetTitle("");
	h->GetXaxis()->SetTitle("E_{#gamma}, MeV");
	h->GetYaxis()->SetTitle("m_{PDG}/m_{meas}");
	fit->Draw("same");
  can->SaveAs("nonlinearity.pdf");
}
