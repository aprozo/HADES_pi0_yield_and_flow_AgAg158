#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TColor.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TLine.h"
#include "TLatex.h"
#include "exp_setup_cells.h"
#include "EcalCellPar.h"


const Float_t energyAtOffset = 20;


// if(pid.beamtime<88.46) AgAg158
// if(pid.beamtime>88.46&&pid.beamtime<90.8) AgAg123@2500A
// if(pid.beamtime>90.8)	AgAg123@200A

Double_t twcFunc(Double_t * x, Double_t * par)
{
	if(x[0]<par[5]) return par[0]+par[1]*exp(sqrt(par[2]*log(par[3]+par[4]*x[0])));
	else return par[6] + par[7]*x[0];
}


// ===========================================================
// Data structures
// ============================================================


typedef struct {
	Float_t sector, cell, cell2, tof, theta, phi, richQa,
	        emc_mq, mom, beta, charge, mm_qa, energy, maxenergy, time, time2, size, emc_theta,
	        emc_phi, emc_path, distRpcEmc, timeCorr, beamtime, TOT, TOT2,
	        dthetadphi_qa, dtime, chi2, raw_time, raw_time2, goodrich, sys;
} PID;
PID pid;


// ===========================================================
// EcalCell Parent Class
// ============================================================
class EcalCell {
protected:
	TH2D * hist;
	TProfile * prof;
	TF1 * fitOld;
	TF1 * fitNew;
	TString cellName;
	Int_t sector, cell;
	Float_t totOffset;


	Params oldParams, newParams;
public:

	EcalCell(const Int_t &newSector, const Int_t &newCell) : sector(newSector), cell(newCell)
	{
		cellName = Form("Sector%i_Cell%i", sector, cell);
	}

	virtual void setTotThrehold(const Float_t & _totOffset){
		totOffset = _totOffset;
	}

	void setParameters(const EcalCellPar & _pars){
		newParams = _pars.newParams; oldParams = _pars.oldParams;
	}

	EcalCellPar getParameters(){
		return { oldParams, newParams };
	}

	virtual void fill(){
	}

	virtual void calibrate(){
	}

	virtual void draw(){
	}
}; // end of a class

// ================================================================================//
// ================================================================================//

class EcalCellEnergy: public EcalCell {
public:
	EcalCellEnergy(const Int_t &newSector, const Int_t &newCell) : EcalCell(newSector, newCell){
		Int_t nbinx           = 350;
		Float_t x_min         = 250;
		Float_t x_max         = 250;
		Int_t nbiny           = 200;
		Float_t y_min         = -100;
		Float_t y_max         = 1500;
		TString fitFuncEnergy = "[0] + exp([1]+[2] * x) ";

		hist   = new TH2D(cellName + "Energy", cellName + "Energy; ToT; mom", nbinx, x_min, x_max, nbiny, y_min, y_max);// tot vs mom
		fitOld = new TF1(cellName + "fitOld", fitFuncEnergy, x_min, x_max);
	}

	void setTotThrehold(const Float_t & _totOffset) override
	{
		totOffset = _totOffset;
		TString fitFuncEnergy2Par = (TString) Form(" %.2f + [0] + exp([1]) *( exp([2] * x) - exp([2]* %.2f) )",
		                                           energyAtOffset, totOffset);
		fitNew = new TF1(cellName + "fitNew", fitFuncEnergy2Par, 250, 1500);
		fitNew->SetParLimits(0, -10, 5);
		fitNew->SetParLimits(1, 0, 10);
		fitNew->SetParLimits(2, 0, 0.02);
	}

	void fill() override
	{
		TChain * nt = new TChain(Form("Sector%iCell%i", sector, cell));

		nt->Add(Form("/lustre/hebe/hades/user/prozorov/mar19/leptonttree/splitted/tree_Sector%i_Cell%i.root", sector,
		             cell));
		PID pid;
		nt->SetBranchAddress("pid", &pid);
		Long_t nEntry = nt->GetEntries();

		for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry) {
			nt->GetEntry(iEntry);

			Float_t energy1, energy2;
			energy1 = pid.maxenergy;
			if (pid.size == 1) energy2 = 0;
			else energy2 = pid.energy - pid.maxenergy;


			if (pid.size == 1 && (pid.emc_mq < 0.6 || (pid.mom > 600 && pid.emc_mq < 1.3)) && pid.emc_path > 90 &&
			    pid.richQa < 1.2)
			{
				if (!(pid.TOT > (totOffset + 150) || pid.beamtime > 90.8)) continue;

				Float_t energydiff        = (100 * abs(pid.mom - (energy1 + energy2)) / pid.mom);
				Float_t within_resolution = energydiff * sqrt(pid.mom / 1000.) / 3;
				if (within_resolution > 10) continue;

				hist->Fill(pid.TOT, pid.mom);
			}
		}


		for (int i = 1; i <= hist->GetYaxis()->GetNbins(); i++) {
			Double_t mom = ((TAxis *) hist->GetYaxis())->GetBinCenter(i);
			TH1D * projectioAtMom = (TH1D *) hist->ProjectionX("projectioAtMom", i - 5, i + 5);
			Double_t totMax       = ((TAxis *) projectioAtMom->GetXaxis())->GetBinCenter(projectioAtMom->GetMaximumBin());

			for (int j = 1; j <= hist->GetXaxis()->GetNbins(); j++) {
				if (hist->GetBinContent(j, i) > 0 && mom > 600 &&
				    abs(totMax - ((TAxis *) hist->GetXaxis())->GetBinCenter(j)) < 100)
					hist->SetBinContent(j, i, hist->GetBinContent(j, i) * 10);
			}
		}
		delete nt;
	} // fillHist

	void calibrate() override
	{
		fitOld->SetParameter(0, oldParams.adc0);
		fitOld->SetParameter(1, oldParams.adc1);
		fitOld->SetParameter(2, oldParams.adc2);

		fitNew->SetParameter(1, oldParams.adc1);
		fitNew->SetParameter(2, oldParams.adc2);
		fitNew->SetRange(totOffset - 10, 1500);

		hist->Fit(fitNew, "REQNM");
		cout << "Sector " << sector << " Cell " << cell << " done!" << endl;

		newParams.adc0 = energyAtOffset + fitNew->GetParameter(0) - exp(fitNew->GetParameter(1) ) * exp(fitNew->GetParameter(
																											2) * totOffset);
		newParams.adc1 = fitNew->GetParameter(1);
		newParams.adc2 = fitNew->GetParameter(2);
	}

	void drawCrossPoint()
	{
		TLine * lh = new TLine(250, energyAtOffset, 1500, energyAtOffset);
		TLine * lv = new TLine(totOffset, -100, totOffset, 1500);

		lh->SetLineColor(kBlack);
		lv->SetLineColor(kBlack);
		lh->SetLineWidth(2);
		lv->SetLineWidth(2);
		lh->Draw("same");
		lv->Draw("same");
	}

	void draw() override
	{
		gStyle->SetOptStat(0);
		gPad->SetGrid(1);
		hist->SetTitle("");
		hist->Draw("colz");
		drawCrossPoint();
		fitOld->SetLineColor(kGreen + 2);
		fitOld->Draw("same");
		fitNew->Draw("same");

		TLegend * legend = new TLegend(0.1, 0.7, 0.42, 0.85);
		//legend->AddEntry(fitOld, "old Fit", "L");
		//legend->AddEntry(fitNew, "new Fit 2 par(ToT treshold)", "L");
		legend->SetLineWidth(0);
	//	legend->Draw("same");
	}
};


// ================================================================================//
// ================================================================================//


class EcalCellTime: public EcalCell {
public:
	EcalCellTime(const Int_t &newSector, const Int_t &newCell) : EcalCell(newSector, newCell){
		Int_t nbinx   = 300;
		Float_t x_min = 0;
		Float_t x_max = 1200;
		Int_t nbiny   = 150;
		Float_t y_min = -10;
		Float_t y_max = 20;

		hist = new TH2D(cellName + "Time", cellName + "Energy; Energy, MeV; t_{ECAL}-t_{RPC} - #Delta t, ns", nbinx, x_min, x_max,
		                nbiny, y_min, y_max);


		//	TString fitFuncTimeNew = "[0]+[1]*exp(sqrt([2]*log([3]+[4]*x)))  ";

		TString fitFuncTimeNew = "[0]+[1]/sqrt(x-[2])";


		fitNew = new TF1(cellName + "fitNew", fitFuncTimeNew, x_min, x_max);
		//fitNew = new TF1(cellName + "fitNew", twcFunc, x_min, x_max,8);


		fitNew->SetParameter(0, 0);
		fitNew->SetParameter(1, 35);
		fitNew->SetParameter(2, energyAtOffset);
		fitNew->SetParameter(3, 3);
		fitNew->SetParLimits(3, 0,10);
		fitNew->SetParameter(4, energyAtOffset);

		TString fitFuncTimeOld = "[0]  +[1]/sqrt(x-[2])";
		fitOld = new TF1(cellName + "fitOld", fitFuncTimeOld, x_min, x_max);
		fitOld->SetParameter(0, oldParams.twc0);
		fitOld->SetParameter(1, oldParams.twc1);
		fitOld->SetParameter(2, oldParams.twc2);
	}


	void fill() override
	{
		TChain * nt = new TChain(Form("Sector%iCell%i", sector, cell));

		nt->Add(Form("/lustre/hebe/hades/user/prozorov/mar19/leptonttree/splitted/tree_Sector%i_Cell%i.root", sector,
		             cell));
		PID pid;
		nt->SetBranchAddress("pid", &pid);
		Long_t nEntry = nt->GetEntries();
		for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry) {
			nt->GetEntry(iEntry);

			Float_t time      = pid.raw_time + oldParams.tdc1;
			Float_t newEnergy = newParams.adc0 + exp(newParams.adc1 + newParams.adc2 * pid.TOT);
			// Float_t oldTwc  = oldParams.twc0 + oldParams.twc1 / TMath::Sqrt(pid.TOT - oldParams.twc2);
			// Float_t newTwc = newParams.twc0 + newParams.twc1 / TMath::Sqrt(pid.TOT - newParams.twc2);
			// // time = time - newTwc;
			Float_t energydiff        = (100 * abs(pid.mom - newEnergy) / pid.mom);
			Float_t within_resolution = energydiff * sqrt(pid.mom / 1000.) / 3;
			if (within_resolution > 10) continue;

			if (pid.size==1 && (pid.emc_mq < 1.5 || (newEnergy> 400 && pid.emc_mq < 0.8) )
			    &&pid.emc_path >0 && pid.distRpcEmc >0 && pid.sys==0) // &&pid.distRpcEmc < 350)
				hist->Fill(newEnergy, time - pid.tof - pid.timeCorr);
		}
		delete nt;
	} // fillHist


	void calibrate() override
	{
		prof = hist->ProfileX("prof" + cellName, 1, -1, "g"); // gives the errors of corresponding x value


		for (Int_t binX = 0.5*(hist->GetXaxis()->GetNbins()); binX < hist->GetXaxis()->GetNbins(); binX++) {
			for (Int_t binY = 0; binY < hist->GetYaxis()->GetNbins(); binY++) {
				hist->SetBinContent(binX,binY,  hist->GetBinContent(binX,binY));
			}
		}


		fitNew->SetRange(energyAtOffset+5, 1200);
		hist->Fit(fitNew, "RQNEM");

		newParams.twc0 = fitNew->GetParameter(0);
		newParams.twc1 = fitNew->GetParameter(1);
		newParams.twc2 = fitNew->GetParameter(2);
//		newParams.twc3 = fitNew->GetParameter(3);
//		newParams.twc4 = fitNew->GetParameter(4);


		if (!isGoodCell(sector, cell)) {
			newParams.twc0 = 0;
			newParams.twc1 = 0;
			newParams.twc2 = 0;
//			newParams.twc3 = 0;
//			newParams.twc4 = 0;
		}
	}










	//
	// void calibrate() override
	// {
	// 	TString fitFuncTimeNew = "[0]+[1]*exp(sqrt([2]*log([3]+[4]*x)))";
	//
	// 	TF1 * firstFit = new TF1(cellName + "firstFit", fitFuncTimeNew, 0, 1200);
	// 	firstFit->SetParameter(0, 15);
	// 	firstFit->SetParameter(1, -2);
	// 	firstFit->SetParameter(2, 0.5);
	// 	firstFit->SetParameter(3, -150);
	// 	firstFit->SetParameter(4, 50);
	// 	firstFit->SetParameter(5, 600);
	//
	// 	firstFit->SetParLimits(1, -50, 0);
	// 	firstFit->SetParLimits(2, 0, 3);
	// 	firstFit->SetParLimits(4, 0.01, 100);
	//
	// 	firstFit->SetRange(energyAtOffset+5, 600);
	// 	hist->Fit(firstFit, "RQNEM");
	// 	firstFit->SetRange(energyAtOffset+5, 1200);
	//
	// 	prof = hist->ProfileX("prof" + cellName, 1, -1, "g"); // gives the errors of corresponding x value
	// 	TH1D* proj=   hist->ProjectionY("py", hist->GetXaxis()->FindBin(700),hist->GetXaxis()->FindBin(1200));
	// 	TF1 *gausFit = new TF1("gausFit","gaus",-10,20);
	// 	proj->Fit(gausFit,"RQN");
	// 	Float_t asymp =  gausFit->GetParameter(1);
	// 	fitNew->SetParameter(6, asymp);
	// 	Int_t bin;
	//
	// 	for (bin = 0.3*(prof->GetXaxis()->GetNbins()); bin < prof->GetXaxis()->GetNbins(); bin++) {
	// 		Float_t fitVal = firstFit->Eval(prof->GetBinCenter(bin));
	// 		if((fitVal - asymp)<=0) break;
	// 	}
	//
	// 	fitNew->SetParameter(7, -0.01);
	// 	fitNew->SetParameter(5, prof->GetBinCenter(bin));
	// 	fitNew->SetParLimits(5, prof->GetBinCenter(bin) -10,  prof->GetBinCenter(bin) + 10);
	//
	//
	// 	fitNew->SetParameter(0, firstFit->GetParameter(0));
	// 	fitNew->SetParameter(1, firstFit->GetParameter(1));
	// 	fitNew->SetParameter(2, firstFit->GetParameter(2));
	// 	fitNew->SetParameter(3, firstFit->GetParameter(3));
	// 	fitNew->SetParameter(4, firstFit->GetParameter(4));
	//
	//
	// 	fitNew->SetRange(energyAtOffset+5, 1200);
	// 	hist->Fit(fitNew, "RQNEM");
	//
	// 	// TProfile * check = (TProfile *)prof->Clone("check");
	// 	// // check->Add(fitNew,-1.);
	// 	//
	// 	// prof->Rebin(4);
	// 	// for (bin = 0.8*(prof->GetXaxis()->GetNbins()-1); bin > prof->FindBin(600); bin--) {
	// 	// 	Float_t fitVal = fitNew->Eval(prof->GetBinCenter(bin));
	// 	// 	if(prof->GetBinContent(bin-1)< fitVal && prof->GetBinContent(bin)> fitVal ) break;
	// 	// }
	// 	// Float_t par5 = prof->GetBinCenter(bin);
	// 	// Float_t par6 = prof->GetBinContent(bin);
	// 	// //	delete check;
	// 	//
	// 	// // fitNew->FixParameter(0,fitNew->GetParameter(0));
	// 	// // fitNew->FixParameter(1,fitNew->GetParameter(1));
	// 	// // fitNew->FixParameter(2,fitNew->GetParameter(2));
	// 	// // fitNew->FixParameter(3,fitNew->GetParameter(3));
	// 	// // fitNew->FixParameter(4,fitNew->GetParameter(4));
	// 	// //
	// 	//
	// 	// fitNew->SetParameter(5, par5);
	// 	// //fitNew->SetParLimits(5, par5-50,par5+50);
	// 	// fitNew->SetParameter(6, par6);
	// 	// //fitNew->SetParLimits(6, par6-1,par6+1);
	// 	// //fitNew->SetParameter(7, -0.001);
	// 	// //fitNew->SetParLimits(7,0,-0.01);
	//
	// 	newParams.twc0 = fitNew->GetParameter(0);
	// 	newParams.twc1 = fitNew->GetParameter(1);
	// 	newParams.twc2 = fitNew->GetParameter(2);
	// 	newParams.twc3 = fitNew->GetParameter(3);
	// 	newParams.twc4 = fitNew->GetParameter(4);
	//
	//
	// 	if (!isGoodCell(sector, cell)) {
	// 		newParams.twc0 = 0;
	// 		newParams.twc1 = 0;
	// 		newParams.twc2 = 0;
	// 		newParams.twc3 = 0;
	// 		newParams.twc4 = 0;
	// 	}
	// }






	void draw() override
	{
		hist->GetXaxis()->SetTitleOffset(0.9);
		hist->GetYaxis()->SetTitleOffset(0.9);
		hist->GetXaxis()->SetTitleSize(0.05);
		hist->GetYaxis()->SetTitleSize(0.05);
		gStyle->SetOptStat(0);
		gPad->SetGrid(1);
		hist->SetTitle("");
		hist->Draw("colz");
		fitNew->Draw("same"); fitNew->SetLineWidth(1);
		Float_t chi2Ndf = fitNew->GetChisquare()/fitNew->GetNDF();
		//fitOld->Draw("same");
		fitOld->SetLineColor(kGreen);       fitOld->SetLineWidth(1);
		prof->Draw("same");

		TLegend * legend = new TLegend(0.1, 0.7, 0.42, 0.85);
		//	legend->AddEntry(fitOld, "old Fit", "L");
		legend->AddEntry(fitNew, "new Fit", "L");
		legend->SetLineWidth(0);
		//legend->Draw("same");

		TLatex * latex = new TLatex();
		latex->SetTextSize(0.045);
		latex->SetTextFont(42);
		latex->SetTextAlign(12);
		Float_t labelOffsetX = 0.43; // 0.43 Lambda, K0s 0.53 Hypertriton, XiM
		TString label1 ="a_{0} =%.2f a_{1} =%.2f a_{2} =%.2f";
		latex->DrawLatexNDC(labelOffsetX, 0.70, Form(label1,newParams.twc0,newParams.twc1,newParams.twc2));
	//	latex->DrawLatexNDC(labelOffsetX, 0.80, Form("chi2/Ndf = %.2f",chi2Ndf));
	}
};


typedef struct {
	Float_t sector, cell,  tof, emc_mq, mom, multEmc,
	        beta, charge, energy, time, emc_path,distRpcEmc,timeCorr,
	        beamtime, TOT, qualityDThDPh, qualityDTime, raw_time, sys;
}
PIDCharged;



class EcalCellTimeCharged: public EcalCellTime {
public:
	EcalCellTimeCharged(const Int_t &newSector, const Int_t &newCell) : EcalCellTime(newSector, newCell){
	};

	void fill() override
	{
		TChain * nt = new TChain(Form("Sector%iCell%i", sector, cell));
		nt->Add(Form("/lustre/hebe/hades/user/prozorov/mar19/ChargedParticlesTime/splitted/tree_Sector%i_Cell%i.root", sector, cell));

		PIDCharged pid;
		nt->SetBranchAddress("pid", &pid);
		Long_t nEntry = nt->GetEntries();
		for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry) {
			nt->GetEntry(iEntry);

			Float_t time      = pid.raw_time + oldParams.tdc1;
			Float_t newEnergy = newParams.adc0 + exp(newParams.adc1 + newParams.adc2 * pid.TOT);
			Float_t weight =1.;
			if (pid.emc_path > 90 &&  pid.emc_mq < 0.5 && pid.beta>0.5  && pid.distRpcEmc>0)
				hist->Fill(newEnergy, time - pid.tof-pid.timeCorr, weight);

		}
		delete nt;
	} // fillHist

};
