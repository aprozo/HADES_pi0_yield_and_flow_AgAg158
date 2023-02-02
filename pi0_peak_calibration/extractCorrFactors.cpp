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


void extractCorrFactors(Int_t iteration = 0, Float_t parN = 1.0)
{
	TFile * out = new TFile(Form("checkHistogramsIteration%iparN%.1f.root", iteration, parN), "RECREATE");

	gROOT->SetBatch(1);


	CorrFactors corrFactors[6][163];
	if (iteration!=0) readCorrFactors("AllCorrFactors.root", corrFactors);

	TFile * inFile;

//	if (iteration==0) inFile= new TFile("pi0PeakHistsIteration0.root", "read");
//	else
   inFile= new TFile("pi0PeakHists.root", "read");

	vector<EcalCell> ecalCells;
	EcalCell allCells(0, 0);
	allCells.Initialize();


	for (Int_t iSector = 0; iSector < 6; iSector++) {
		for (Int_t iCell = 0; iCell < 163; iCell++) {
			EcalCell a(iSector, iCell);
			if (isGoodCell(iSector, iCell)) a.Initialize();
			ecalCells.push_back(a);
		}
	}

	TDirectory * saved = gDirectory;
	TDirectory * sub   = out->mkdir("Cells");

	for (Int_t iSector = 0; iSector < 6; iSector++) {
		for (Int_t iCell = 0; iCell < 163; iCell++) {
			if (!isGoodCell(iSector, iCell)) continue;
			cout << "iSector " << iSector << " iCell " << iCell << endl;
			ecalCells[iSector * 163 + iCell].ReadHists(inFile);
			ecalCells[iSector * 163 + iCell].Evaluate(iteration, corrFactors[iSector][iCell], parN);
			sub->cd();
			ecalCells[iSector * 163 + iCell].Write();
			saved->cd();
		}
	}

	allCells.ReadHists(inFile);
	allCells.Evaluate(iteration, corrFactors[0][0], parN);
	out->cd();
	allCells.Write();


	writeCorrFactors("AllCorrFactors.root", corrFactors);
	writeCorrFactors(Form("AllCorrFactorsIteration%iparN%.1f.root", iteration, parN), corrFactors);
	out->Save();
	out->Close();
} /* fill_hists_diphoton */
