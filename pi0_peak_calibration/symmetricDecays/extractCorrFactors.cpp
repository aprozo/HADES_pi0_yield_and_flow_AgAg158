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


#include "CorrFactors.h"
#include "EcalCell.h"
#include "exp_setup_cells.h"


using namespace std;
//


Double_t fitEnergyCorr(Double_t* x, Double_t* par)
{
	if(x[0]<=par[5]) return par[0]*x[0]+par[1]*sqrt(x[0])+par[2] + par[3]/sqrt(x[0]) + par[4]/x[0];
	else return (par[0]-par[2]/par[5]-2*par[3]/pow(par[5],1.5) -3*par[4]/pow(par[5],2))*x[0] + (par[1]+2*par[2]/sqrt(par[5])+3*par[3]/par[5] +4*par[4]/pow(par[5],1.5))*sqrt(x[0]);
}




void extractCorrFactors()
{
	Float_t parN=1.0;
	TFile * out = new TFile(Form("checkHistogramsIteration%iparN%.1f.root", 0, 1.0), "RECREATE");

	gROOT->SetBatch(1);


	CorrFactors corrFactors[6][163];
	//readCorrFactors("AllCorrFactors.root", corrFactors);

	TFile * inFile = new TFile("pi0PeakHists.root", "read");

	vector<EcalCell> ecalCells;
	EcalCell allCells(0, 0);
	allCells.Initialize();


	// for (Int_t iSector = 0; iSector < 6; iSector++) {
	// 	for (Int_t iCell = 0; iCell < 163; iCell++) {
	// 		EcalCell a(iSector, iCell);
	// 		if (isGoodCell(iSector, iCell)) a.Initialize();
	// 		ecalCells.push_back(a);
	// 	}
	// }

	// TDirectory * saved = gDirectory;
	// TDirectory * sub   = out->mkdir("Cells");
	//
	// for (Int_t iSector = 0; iSector < 6; iSector++) {
	// 	for (Int_t iCell = 0; iCell < 163; iCell++) {
	// 		if (!isGoodCell(iSector, iCell)) continue;
	// 		cout << "iSector " << iSector << " iCell " << iCell << endl;
	// 		ecalCells[iSector * 163 + iCell].ReadHists(inFile);
	// 		ecalCells[iSector * 163 + iCell].Evaluate(iteration, corrFactors[iSector][iCell], parN);
	// 		sub->cd();
	// 		ecalCells[iSector * 163 + iCell].Write();
	// 		saved->cd();
	// 	}
	// }

	allCells.ReadHists(inFile);
	allCells.Evaluate(0, corrFactors[0][0], parN);
	out->cd();
	allCells.Write();

	// TH1D * h= allCells.getEnergyCorrection();
	// TF1 * fit = new TF1("fit", fitEnergyCorr,6,80, 800);
	// fit->SetParameters(1,-0.25,0.65,-0.47,0.12,500);
	// h->Fit(fit,"RQEM");
	//
	// fit->Write();



	out->Save();
//	out->Close();
} /* fill_hists_diphoton */
