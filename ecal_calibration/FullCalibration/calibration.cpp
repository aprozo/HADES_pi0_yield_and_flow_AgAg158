#include "exp_setup_cells.h"
#include "EcalCell.h"
#include "EcalCellPar.h"

using namespace std;

void calibration()
{
  EcalCellPar ecalCellPars[6][163];
  readParams ("mar19EmcCalParNewEnergy.txt", ecalCellPars, "old");
  //readParams ("mar19EmcCalParNewTime.txt",   ecalCellPars, "new");

  TFile fTotThresholds("ToTtreshold_hist.root", "read");
  TH1 * hTotThresholds = (TH1 *) fTotThresholds.Get("hToTtreshold");

  int N = 2;           // number of canvases (Npages needed +1); 4x4->52; 3x3->92
  int row = 1, col = 1; // number of plots in roows and cols on canvas
  TCanvas * canvas[N];

  for (Int_t i = 0; i < N; i++) {
    canvas[i] = new TCanvas(Form("canvas%i", i), Form("canvas%i", i), 800, 600);
    canvas[i]->Divide(row, col);
  }
  Int_t counter = 0;
  // Int_t sector = 2;
  // Int_t cell = 109;

  for (Int_t sector = 0; sector < 2; sector++) {
    for (Int_t cell = 142; cell < 143; cell++) {
      if (!isGoodCell(sector, cell)) continue;
      cout << "sector " << sector << " cell " << cell << endl;
      EcalCellTime module(sector, cell);
      //EcalCellEnergy module(sector, cell);
      module.setParameters(ecalCellPars[sector][cell]);
      Double_t totTreshold = hTotThresholds->GetBinContent(sector * 200 + cell + 1);
      module.setTotThrehold(totTreshold);
      module.fill();
      module.calibrate();
      ecalCellPars[sector][cell] = module.getParameters();

      canvas[(int) (counter / (row * col))]->cd(counter % (row * col) + 1);
      module.draw();
      counter++;
    }
  }

  writeParams("mar19EmcCalParNewTime.txt", ecalCellPars);

  TString outputfile = "check";

  canvas[0]->SaveAs(outputfile + ".pdf[");
  for (Int_t i = 0; i < N - 1; i++) canvas[i]->SaveAs(outputfile + ".pdf");
  canvas[N - 1]->SaveAs(outputfile + ".pdf]");
} // calibration
