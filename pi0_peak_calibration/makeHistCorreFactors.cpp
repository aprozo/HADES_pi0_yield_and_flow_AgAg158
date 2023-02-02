#include "TFile.h"
#include "CorrFactors.h"
#include "TH1.h"
#include "TString.h"

void makeHistCorreFactors()
{
	CorrFactors corrFactors[6][163];

	TFile * out       = new TFile("AllCorrFactors.root", "RECREATE");
	TH2D * AllCorrFactors =new TH2D("AllCorrFactors","Correction factor ; sector*163+cell ; Old energy, MeV",  (163*6), 0, (163*6), nEnergyBins, energyBinAxis);


	TFile * in       = new TFile("allEnergiesCorrFactors.root", "read");
	TH1D * allEnergiesCorrFactors =(TH1D*)in->Get("allEnergiesCorrFactors");

	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			Double_t val=allEnergiesCorrFactors->GetBinContent(sector*163+cell+1);
			if(val<0.9&&val>0.89999) val=1.;
			for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++)
				AllCorrFactors->SetBinContent(sector*163+cell+1,iEnergyBin+1,val);
		}
	}
	out->cd();
	AllCorrFactors->Write();
	out->Save();
	out->Close();
}
