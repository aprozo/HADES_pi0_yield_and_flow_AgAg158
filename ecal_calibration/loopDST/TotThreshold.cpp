#include "Includes.h"
#include "exp_setup_cells.h"
#include "EcalCellPar.h"



void analysis(TString inFile, TString outputFile, Int_t nev = 1000000.)
{
	Bool_t createHades = kTRUE; // kFALSE;  // kTRUE = create HADES new
	HLoop loop(createHades);  // kTRUE : create Hades  (needed to work with standard eventstructure)
	gHades->setBeamTimeID(HADES::kMar19);
	Bool_t ret = kFALSE;

	if (inFile.Contains(",")) {
		ret = loop.addMultFiles(inFile); // file1,file2,file3
	} else {
		ret = loop.addFiles(inFile); // myroot*.root
	}
	if (ret == 0) {
		cout << "READBACK: ERROR : cannot find inputfiles : " << inFile.Data() << endl;
		return;
	}
	// if(!loop.setInput("")) { // all input categories
	if (!loop.setInput("-*,+HStart2Hit,+HEmcCal")) {
		cout << "READBACK: ERROR : cannot read input !" << endl;
		exit(1);
	} // read all categories
	loop.printCategories();
	loop.printChain();

	TFile * out = new TFile(outputFile, "RECREATE");
	out->cd();
	Hades * gHades = new Hades;
	gHades->setQuietMode(2);

	HCategory * catEmcCal          = loop.getCategory("HEmcCal");
	HCategory * catStart2Hit       = loop.getCategory("HStart2Hit");
	HIterator *iterStart = (HIterator*)catStart2Hit->MakeIterator("native");

	Array<EcalCellPar> ecalCells = initialize();
	readParams("mar19EmcCalPar.txt", ecalCells);

	// ==========================================HISTOGRAMS======================================================================//
	TH2F * hTotCell    = new TH2F("hTotCell", "hTotCell", 1200, 0, 1200, 800, 100, 900);

	TH2F * hTotCellCleaned    = new TH2F("hTotCellCleaned", "hTotCellCleaned", 1200, 0, 1200, 800, 100, 900);
	TDirectory * subD = out->mkdir("TotVsTime");
	subD->cd();

	TH2F * hTimeTot[6][163];
	TH1D * 	hTot[6][163];


	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			if(isGoodCell(sector,cell))
				hTimeTot[sector][cell] = new TH2F(Form("sector%i_cell%i", sector, cell), Form("sector%i_cell%i;TOT,ns;Time,ns", sector, cell), 300, 200, 700, 500, -20, 600);
		}
	}
	// =======================================================================================================================================//

	Int_t nEntries = loop.getEntries();
	gStyle->SetOptStat(0);
	printf("nEntries = %i\n", nEntries);
	if (nEntries > nev && nev > 0) nEntries = nev;


	for (Int_t i = 0; i < nEntries; i++) {
		if (loop.nextEvent(i) <= 0) { cout << " end recieved " << endl; break; } // last event reached
		HTool::printProgress(i, nEntries, 1, "Analyze pairs :");

		Float_t start_time = 0.0;
		HStart2Hit *start = 0;
		iterStart->Reset();
		while ((start=(HStart2Hit *)iterStart->Next()) != 0) {
			start_time =start->getTime();
		}


		for (Int_t i = 0; i < catEmcCal->getEntries(); i++) { // loop over cal level
			HEmcCal * emc_cal = (HEmcCal *) catEmcCal->getObject(i);
			Int_t pos         = emc_cal->getCell();
			Int_t cell        = systemMap[pos];
			Int_t sector      = emc_cal->getSector();
			Float_t time      = emc_cal->getTime()-start_time ;
			Float_t tot       = emc_cal->getWidth();

			if(!isGoodCell(sector,cell)) continue;

			hTimeTot[sector][cell]->Fill(time, tot);
			if (time < 40 && time > -20) {
				hTotCell->Fill(sector * 200 + cell, tot);
			}
		}//loop over EmcCal
	}//loop over Events



	for (Int_t sector = 0; sector < 6; sector++) {
		for (Int_t cell = 0; cell < 163; cell++) {
			if(!isGoodCell(sector,cell)) continue;
			TH2F * h        = hTimeTot[sector][cell];
			TH1D * proj     = hTimeTot[sector][cell]->ProjectionY("proj",150,300);
			Float_t average = proj->GetBinLowEdge(proj->GetMaximumBin());
			Float_t min     = average - 20;
			Float_t max     = average + 20;

			for (Int_t biny = 1; biny < h->GetYaxis()->GetNbins(); biny++) {
				if (h->GetYaxis()->GetBinLowEdge(biny) < min || h->GetYaxis()->GetBinLowEdge(biny) > max) {
					for (Int_t binx = 1; binx < h->GetXaxis()->GetNbins(); binx++) {
						h->SetBinContent(binx, biny, 0);
					}
				}
			}
			hTot[sector][cell] = h->ProjectionX("Tot_"+(TString)h->GetName());

		}
	}    //filter histograms


	out->cd();

	out->GetList()->Write();




	out->Save();
	out->Close();
	return;
} /* fill_ttree_diphoton */
