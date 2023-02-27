#include "TObject.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TString.h"
#include "TSystem.h"
#include "TF1.h"
#include <iostream>
#include "TStyle.h"
#include "TCanvas.h"
//#include "exp_setup_cells.h" // some useful functions for ECAL cell operations


using namespace std;


void fitHists()
{
	//TFile * in = new TFile("photonsBeta.root", "read");
	//TFile * in = new TFile("photonsBeta200BeamBins.root", "read");
	TFile * in = new TFile("photonsBetaAfterCut.root", "read");
	TFile * out = new TFile("betaBeamMeanSigma.root", "recreate");
	out->cd();

	TH1D *  nEventsBeam= (TH1D *)in->Get("nEventsBeam");
	Int_t nBeamBins=nEventsBeam->GetNbinsX();
	TProfile * averageMean = new TProfile("averageMean",";cell",1200,0,1200,0.9,1.2);
	TProfile * averageSigma = new TProfile("averageSigma",";cell",1200,0,1200);
	TH2F * betaBeamMean = new TH2F("mean",";beam time, days;cell",nBeamBins,68,88,1200,0,1200);
	TH2F * betaBeamSigma = new TH2F("sigma",";beam time, days; cell",nBeamBins,68,88,1200,0,1200);
	TH2D * beamCells = new TH2D("beamCells","Integral of photons within 1 sigma [1/event];beam time, days;sector*200+cell; Integral/events",nBeamBins,68,88,1200,0,1200);
	gStyle->SetOptStat(0);
	beamCells->GetZaxis()->SetRangeUser(10,1E5);


	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			TString histName=Form("sector%i_cell%i_beam_vs_beta",sector,cell);
			TH2F* betaBeam = (TH2F *) in->Get(histName);
			if(betaBeam->GetEntries()==0) continue;
			for (Int_t binX = 1; binX <=betaBeam->GetNbinsX(); binX++) {
				Int_t nEvents = nEventsBeam->GetBinContent(binX);
				for (Int_t binY = 1; binY <=betaBeam->GetNbinsY(); binY++) {
					Double_t val = betaBeam->GetBinContent(binX,binY);
					if(nEvents>0) val/=nEvents;
					betaBeam->SetBinContent(binX,binY,val);
				}
			}

			std::cout << "sector " <<sector<<" cell "<<cell<< '\n';

			for (Int_t bin = 1; bin <= betaBeam->GetNbinsX(); bin++) {
				TH1D * projY = betaBeam->ProjectionY("px", bin, bin);
				if (projY->GetEntries() ==0) continue;
				TF1 * ff = new TF1("ff", "gaus", 0.9, 1.3);
				//ff->SetParameters(projY->GetMaximum(),1., 0.05);

				projY->Fit("ff","QR");
				Float_t mean  = ff->GetParameter(1);
				Float_t sigma = ff->GetParameter(2);
				delete ff;

				if ((mean > 0.9 ||  mean < 1.2) && sigma < 0.2) {
					betaBeamMean->SetBinContent(bin,sector*200+cell,mean);
					betaBeamSigma->SetBinContent(bin,sector*200+cell,sigma);
					averageMean->Fill(sector*200+cell,mean);
					averageSigma->Fill(sector*200+cell,sigma);
				}
				else {
					mean=0;
					sigma=0;
				}

				Double_t integral =  projY->Integral(projY->FindBin(mean-sigma),projY->FindBin(mean+sigma));
				beamCells->SetBinContent(bin,sector*200+cell,integral);


			} //loop over  beamtime bin


			betaBeam->Write();


		}//loop over  cells
	}//loop over  bsectors
	betaBeamMean->Write();
	betaBeamSigma->Write();
	averageMean->Write();
	averageSigma->Write();
	beamCells->Write();


TCanvas * can =new TCanvas("can","",800,600);
can->cd();
can->SetLogz();
beamCells->GetYaxis()->SetRangeUser(320,1170);
beamCells->Draw("colz");
can->SaveAs("IntegralCellsAfter.pdf");

} // fitHists
