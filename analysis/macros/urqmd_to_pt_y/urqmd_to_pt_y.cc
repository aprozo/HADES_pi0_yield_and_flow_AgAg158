#include "iostream"
#include "TString.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TChain.h"
#include "../GlobalParameters.h"

struct Pid { Float_t pt, y, phi, centBin; };


Int_t urqmd_to_pt_y(Int_t maxEntries = 100000000)
{
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	GlobalParameters par("../parameters.conf");

	Double_t minPt = 0.;
	Double_t maxPt = 1200.; // MeV
	Double_t minY  = 0.05;
	Double_t maxY  = 1.05;
	TFile * out = new TFile("urqmd_pt_y.root", "recreate");
	TH2D * yieldPtY[3];

	for (Int_t iCent= 0; iCent< 3; iCent++) {
		yieldPtY[iCent]= new TH2D(Form("hAll_cent%i", iCent), ";y_{cm} ;p_{t}, [MeV/c]", par.N_Y_BINS, minY, maxY, par.N_PT_BINS, minPt, maxPt);
	}
	Float_t yBinSize = yieldPtY[0]->GetXaxis()->GetBinWidth(1);
	Float_t ptBinSize = yieldPtY[0]->GetYaxis()->GetBinWidth(1);

	TFile in(par.PATH_TO_URQMD_DATA, "read");
	TH1D * Counter_cent = (TH1D *) in.Get("EventNum");

	TChain * nt = new TChain("nt");
	nt->Add(par.PATH_TO_URQMD_DATA);

	Pid pid;
	nt->SetBranchAddress("pid", &pid.pt);
	Int_t nentries = (Int_t) nt->GetEntries();
	fprintf(stderr, "Entries: %d\n", nentries);
	for (Int_t entryLoop = 0; entryLoop < nentries; entryLoop++) {
		if (entryLoop % 50000000 == 0) printf(" Entry number: %d\n", entryLoop);
		nt->GetEntry(entryLoop);
		Int_t cent=(int)(pid.centBin-1);
		yieldPtY[cent]->Fill(pid.y-par.MIDRAPIDITY, pid.pt);
	}
	out->cd();

	for (Int_t iCent= 0; iCent< 3; iCent++) {
		Double_t nEvents= Counter_cent->GetBinContent(iCent+2);
		yieldPtY[iCent]->Scale(1./(yBinSize*ptBinSize));
		yieldPtY[iCent]->Scale(1. / nEvents);
		yieldPtY[iCent]->Write();
	}
	out->Save();
	out->Close();
	return 0;
} // urqmd_to_pt_y

int main(int argc, char ** argv)
{
	std::cout << "urqmd_to_pt_y()\n";

	if (argc == 1) return urqmd_to_pt_y(1000000);

	if (argc == 2) return urqmd_to_pt_y(atoi(argv[1]));

	std::cout << "ERROR: Wrong number of arguments!\n urqmd_to_pt_y nMaxEvents=1000000\n";
	return 2;
}
