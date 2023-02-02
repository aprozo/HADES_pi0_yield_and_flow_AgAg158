#ifndef __CORRFACTORS__
#define  __CORRFACTORS__


#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"


const Int_t nEnergyBins   = 1;
const Float_t minEnergy = 0.;
const Float_t maxEnergy = 1500.;

const Double_t energyBinAxis[nEnergyBins+1] ={minEnergy,maxEnergy}; //100*(1.3)^n


const Double_t averageCorrFactors[nEnergyBins]={1.};



//const Double_t energyBinAxis[nEnergyBins+1] = {minEnergy,maxEnergy};

Int_t getEnergyBin(Float_t energy){
	for (size_t i = 0; i < nEnergyBins; i++) {
		if (energy >=energyBinAxis[i] && energy<energyBinAxis[i+1]) return i;
	}
	if (energy>=energyBinAxis[nEnergyBins]) return (nEnergyBins-1);
	return 0;
}

Double_t getAverageCorrection(Int_t bin){
	return averageCorrFactors[bin];
}

struct CorrFactors {
	CorrFactors(){
		for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
			val[iEnergyBin] =1.;
		}
	};
	Float_t val[nEnergyBins];
};


void readCorrFactors(const TString &inputFile, CorrFactors corr[6][163])
{

	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
				corr[sector][cell].val[iEnergyBin] =1.0;
			}
		}
	}

	TFile * in       = new TFile(inputFile, "read");
	if (in->IsZombie()) return;
	TH2D * AllCorrFactors =(TH2D*)in->Get("AllCorrFactors");
	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			//if(AllCorrFactors->GetBinContent(sector*163+cell+1,1)<0.5) AllCorrFactors->SetBinContent(sector*163+cell+1,1, AllCorrFactors->GetBinContent(sector*163+cell+1,1));
			for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++) {
				corr[sector][cell].val[iEnergyBin] =AllCorrFactors->GetBinContent(sector*163+cell+1,0);
			}
		}
	}
}

void writeCorrFactors(const TString &outputFile, CorrFactors corr[6][163])
{
	TFile * out       = new TFile(outputFile, "recreate");
	TH2D * AllCorrFactors =new TH2D("AllCorrFactors","Correction factor ; sector*163+cell ; Old energy, MeV",  (163*6), 0, (163*6), nEnergyBins, energyBinAxis);
	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			for (Int_t iEnergyBin = 0; iEnergyBin < nEnergyBins; iEnergyBin++)
				AllCorrFactors->SetBinContent(sector*163+cell+1,iEnergyBin+1,corr[sector][cell].val[iEnergyBin]);
		}
	}
	out->cd();

	AllCorrFactors->Write();
	out->Save();
}



#endif
