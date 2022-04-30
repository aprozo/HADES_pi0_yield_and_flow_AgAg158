#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TFitResultPtr.h"


// #include "exp_setup_cells.h"

using namespace std;

void getTotThresholds()
{
  TFile * file    = new TFile("TotThreshold.root");
  TH2D * hTotCell = (TH2D *) file->Get("hTotCell");
  TH1D * projs[6][163];                           // to save them to file and have a look
  TF1 * fit = new TF1("fit", "landau", 200, 600); // loose range to make sure ge will have fill peak inside
  TFitResultPtr fitStatus(0);                     // to make sure fit was done properly
  TH1D * hToTtreshold =
    new TH1D("hToTtreshold", "ToT corresponding to 20 MeV energy;200*sector+cell;tot", 1200, 0, 1200);

  Int_t binsize = hTotCell->GetXaxis()->GetNbins();
  Int_t sector, cell;
  Float_t maxbin, maxbinval = 0, minbinval = 0, tresholdval;

  for (Int_t i = 0; i < binsize; i++) {
    sector = (int) (i / 200.);
    cell   = i % 200;
    if (cell > 162) continue;  // do not store/fit those empty ones

    projs[sector][cell] = (TH1D *) hTotCell->ProjectionY(TString::Format("proj_%d_%d", sector, cell), i + 1, i + 1);
    projs[sector][cell]->GetXaxis()->SetRangeUser(200, 600);
    maxbin    = projs[sector][cell]->GetMaximumBin();
    maxbinval = projs[sector][cell]->GetBinCenter(maxbin);
    minbinval = projs[sector][cell]->GetBinCenter(maxbin - 50); // more or less background value

    if (projs[sector][cell]->GetEntries() > 0 && projs[sector][cell]->GetBinContent(maxbin) > 100) {
      fit->SetRange(maxbinval - 30, maxbinval + 5);
      fit->SetParameter(1, maxbinval);
      fitStatus = projs[sector][cell]->Fit(fit, "RQ", "");
      fit->SetRange(maxbinval - 50, maxbinval + 10); // extend the range a bit so we can find this "half height" value

      if (fitStatus >= 0 && fit->GetParameter(2) < 30) tresholdval = fit->GetX((fit->GetMaximum() - minbinval) / 4);  // get ToT value at half of gauss's height
      else tresholdval = 0;
    } else { tresholdval = 0; }
    hToTtreshold->SetBinContent(i + 1, tresholdval);
  }
  TFile * fileTot = new TFile("ToTtreshold_hist.root", "recreate");
  fileTot->cd();
  hToTtreshold->Write();
  fileTot->Save();
  fileTot->Close();

  // debug purposes
  TFile * fileProj = new TFile("ProjectionsDebug.root", "recreate");
  fileProj->cd();
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 163; j++) {
      projs[i][j]->Write();
    }
  }
  fileProj->Save();
  fileProj->Close();
} // getToTtreshold
