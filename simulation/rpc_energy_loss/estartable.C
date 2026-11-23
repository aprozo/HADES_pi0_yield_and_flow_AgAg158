#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCutG.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "TColor.h"
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "/u/sudol/unigen/v2.1/drawStyle.C"
using namespace std;

static Float_t emomentum[28];
static Float_t eloss[28];

void read_eloss_table()
{
	TString FileName("aluminiumenergyloss.txt");
	ifstream inFile1;
	inFile1.open(FileName);
	if (inFile1.fail()) {
		cout << "-E- main : File " << FileName
		     << " does not exist!" << endl;
		return;
	}

	Char_t line1[1000];
	TString strLine1;
	do {
		inFile1.getline(line1,1000);
		strLine1 = line1;
	} while(!strLine1.Contains("MeV       MeV cm2/g"));

	for (Int_t i=0; i<29; i++) {
		inFile1.getline(line1,1000);
		sscanf(line1," %e %e",&emomentum[i],&eloss[i]);
		printf("%e \t %e \n",emomentum[i],eloss[i]);
		if (inFile1.eof()) break;
	}
}

void estartable(void)
{
	drawStyle();
	gStyle->SetOptTitle(0);
	Float_t aldensity=2.6989;// g/cm^3
	read_eloss_table();
	TCanvas *can  = new TCanvas("can","can",800,600);
	TH1F  *energy_loss[10];
	TLegend *leg =plotLegend("left_top","");
	TString inputFile    = "hist_lepton_energy_loss.root";
	TFile    *inFile      = new TFile(inputFile,"read");
	TProfile *true_enloss = (TProfile*) inFile->Get("energy_loss_pfx");
	setLegendEntry(leg,true_enloss,"sim data ",0.5,"pL",0.05);
	for (Int_t i=0; i<10; i++) {
		energy_loss[i]= new TH1F(Form("energy_loss_vs_emom_al_%d_mm",(i+1)),";p_{e} [MeV/c];E_{loss} [MeV]",200,0,1000);
		energy_loss[i]->SetMarkerColor(kOrange-9+2*i);
		energy_loss[i]->SetMarkerStyle(8);
		energy_loss[i]->GetXaxis()->SetTitleOffset(0.75);
		energy_loss[i]->GetXaxis()->SetTitleSize(0.055);
		energy_loss[i]->GetYaxis()->SetTitleOffset(0.8);
		energy_loss[i]->GetYaxis()->SetTitleSize(0.055);
		//energy_loss[i]->SetMarkerSize(2);
		if(i==9 || i==4) setLegendEntry(leg,energy_loss[i],Form("estar %d mm",i+1),1,"pL",0.05);
	}

	for (Int_t i=0; i<10; i++) {
		for (Int_t k=2; k<28; k++)
		{
			Int_t bin=energy_loss[i]->GetXaxis()->FindBin(emomentum[k]);
			energy_loss[i]->SetBinContent(bin,eloss[k]*aldensity*i/10);
		}
	}

	can->cd();
	energy_loss[0]->GetYaxis()->SetRangeUser(0.001,120);


	energy_loss[9]->Draw("same P");
	energy_loss[4]->Draw("same P");



	true_enloss->Rebin(2);
	true_enloss->SetMarkerColor(kBlack);
	true_enloss->SetMarkerStyle(8);
//true_enloss->SetMarkerSize(2);
	true_enloss->Draw("same");
	leg->Draw("same");
	can->SaveAs("EnergyLossRpcEstar.pdf");
}
