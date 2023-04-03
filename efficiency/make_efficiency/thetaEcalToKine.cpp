#include <iostream>
#include <map>
#include <stdio.h>
#include "TH3F.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPad.h"
#include "TTree.h"
#include "TChain.h"
#include "smoothMatrix.C"


void thetaEcalToKine(void)
{
	gStyle->SetOptStat(0);       // switch off stat box
	TFile *out= new TFile("thetaKineToEcal.root","recreate");
	TFile *fileIn= new TFile("../sim_photon_ttree_gen5a.root","Read");
	TChain *treeIn=(TChain*)fileIn->Get("nt");
	// TH2D * h2=new TH2D("h2","h2",400,-2,3,100,10,50);
	// treeIn->Draw("theta_ecal:(0.954*theta_ecal+1-180*theta_kine/3.14)>>h2","emc_cell>=10","colz");
	// h2->FitSlicesX(0,0,-1,0);
	// TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");
	// TF1 *fit = new TF1("fit","pol1",15,50);
	// h2_1->Fit("fit","REM");


	TH2D * phi=new TH2D("phi","phi",400,-2,3,17,0,17);
	treeIn->Draw("columnEMC:(phi_ecal-phi_kine)>>phi","","colz");
	phi->FitSlicesX(0,0,-1,0);
	TH1D *phi_1 = (TH1D*)gDirectory->Get("phi_1");


	out->cd();
	// h2->Write();
	// h2_1->Write();

	phi->Write();
	phi_1->Write();


}
