#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "THStack.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEntryList.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include "TFitResult.h"

#include "exp_setup_cells.h"

using namespace std;


typedef struct {
  Float_t sector, cell, cell2, tof, theta, phi, richQa, emc_mq, mom,
    beta, charge, mm_qa, energy, maxenergy, time, time2, size, emc_theta, emc_phi,
    emc_path, candpart_length, emc_length, beamtime, TOT, TOT2, dthetadphi_qa, dtime, chi2, raw_time, raw_time2,
    goodrich, sys;
}
PID;


template <typename Cell>
using Array = vector<vector<Cell> >;

struct Params {
  Double_t tdc0 = 1, tdc1 = 0;
  Double_t adc0 = 0, adc1 = 1, adc2 = 0;
  Double_t twc0 = 0, twc1 = 0, twc2 = 0;
};


class EcalModule {
public:
  EcalModule(const Int_t &newSector, const Int_t &newCell)
    :
    sector(newSector),
    cell(newCell)
  {
    Int_t nbinx = 150;
    Int_t nbiny = 150;

    name  = Form("Sector%i_Cell%i", sector, cell);
    hTime = new TH2D(name + "Time", name + "Time", nbinx, 0, 1000, nbiny, -10, 20); // tot vs time-tof
    hTot  = new TH1F(name + "tot", name + "tot", nbinx, 0, 600);                    // tot vs time-tof
    //  TString fitFuncTime = "[0]  +  [1]/sqrt(x-[2])";

    TString fitFuncTime = "[0]+[1]*exp(sqrt([2]*log([3] +[4]*x)))";
    fitTime = new TF1(name + "FitTime", fitFuncTime, 0, 1000);
  }

  void setTotThrehold(const Float_t & TotThreshold){ offsetTot = TotThreshold; }

  void fillHist()
  {
    TChain * nt = new TChain(Form("Sector%iCell%i", sector, cell));

    nt->Add(Form("/lustre/hebe/hades/user/prozorov/mar19/leptonttree/splitted/tree_Sector%i_Cell%i.root", sector,
      cell));

    PID pid;
    nt->SetBranchAddress("pid", &pid);
    Long_t nEntry = nt->GetEntries();
    for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry) {
      nt->GetEntry(iEntry);

      // Float_t twc  = oldParams.twc0 + oldParams.twc1 / TMath::Sqrt(pid.TOT - oldParams.twc2);
      Float_t time = pid.time;
      time = pid.raw_time + oldParams.tdc1;
      // Float_t newTwc = newParams.twc0 + newParams.twc1 / TMath::Sqrt(pid.TOT - newParams.twc2);
      // // time = time - newTwc;

      if (pid.emc_path > 0 && (pid.emc_mq < 2.5 || (pid.TOT < 600 && pid.emc_mq < 0.8) ) &&
        pid.sys == 0) hTime->Fill(pid.energy, time - pid.tof - 0.7);
      if (time - pid.tof > 5 && pid.sys == 0 && pid.emc_mq < 3.) hTot->Fill(pid.TOT);
    }


    Float_t max = hTot->GetMaximum();
    for (Int_t binx = 1; binx <= hTot->GetNbinsX(); binx++) {
      if (hTot->GetBinContent(binx) > 0.9 * max) {
        offsetTot = hTot->GetBinLowEdge(binx);
        break;
      }
    }


    delete nt;
  } // fillHist

  ;

  void calibrate()
  {
    cout << name << endl;
    pfx = hTime->ProfileX("pfx" + name, 1, -1, "g"); // gives the errors of corresponding tot
    fitTime->SetParameter(0, 0);
    fitTime->SetParameter(1, -0.2);
    fitTime->SetParameter(2, 2);
    fitTime->SetParameter(3, -100);
    fitTime->SetParameter(4, 0.8);

    fitTime->SetParLimits(0, -10, 20);
    fitTime->SetParLimits(1, 0, -0.8);
    fitTime->SetParLimits(2, 0, 2);
    fitTime->SetParLimits(3, -200, -20);
    fitTime->SetParLimits(4, 0.01, 1);

    hTime->Fit(fitTime, "RQN");

    newParams.twc0 = fitTime->GetParameter(0);
    newParams.twc1 = fitTime->GetParameter(1);
    newParams.twc2 = fitTime->GetParameter(2);

    if (!isGoodCell(sector, cell)) {
      newParams.twc0 = 0;
      newParams.twc1 = 0;
      newParams.twc2 = 0;
    }
  } // calibrate

  void drawTime()
  {
    hTime->GetYaxis()->SetTitle("time-tof, ns");
    hTime->GetXaxis()->SetTitle("TOT, ns");
    hTime->GetXaxis()->SetTitleOffset(0.9);
    hTime->GetYaxis()->SetTitleOffset(0.9);
    hTime->GetXaxis()->SetTitleSize(0.05);
    hTime->GetYaxis()->SetTitleSize(0.05);
    gStyle->SetOptStat("n");
    gPad->SetGrid(1);
    hTime->Draw("colz");
    fitTime->Draw("same");
    pfx->Draw("same");
    // hMostProbTime->Draw("samePE");
  }

  Params oldParams, newParams;
private:


  TH2D * hTime;
  TProfile * pfx;
  TH1F * hTot;
  TF1 * fitTime;
  TF1 * oldEnergy;
  TString name;
  Float_t offsetTot;
  Int_t sector, cell;
}; // end of a class

istream& operator >> (istream &is, Params &par)
{
  if (is) {
    is >> par.tdc0 >> par.tdc1 >> par.adc0 >> par.adc1 >> par.adc2 >> par.twc0 >> par.twc1 >> par.twc2;
  }
  return is;
}

ostream& operator << (ostream &is, Params &par)
{
  if (is) {
    is << setw(15) << par.tdc0 << "\t"
       << setw(15) << par.tdc1 << "\t"
       << setw(15) << par.adc0 << "\t"
       << setw(15) << par.adc1 << "\t"
       << setw(15) << par.adc2 << "\t"
       << setw(15) << par.twc0 << "\t"
       << setw(15) << par.twc1 << "\t"
       << setw(15) << par.twc2;
  }
  return is;
}

void readParams(const TString &inputFile, Array<EcalModule> &ecalCells)
{
  ifstream inFile(inputFile);

  if (inFile.fail()) {
    cout << "-E- main : File " << inputFile << " does not exist!" << endl;
    return;
  }
  string line;
  do {
    getline(inFile, line);
  } while (line != "[EmcCalPar]");
  Params params;
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      inFile >> sector >> cell >> params;
      ecalCells[sector][cell].oldParams = params;
      ecalCells[sector][cell].newParams = params;
      if (inFile.eof())
        break;
    }
  }
}

void readNewParams(const TString &inputFile, Array<EcalModule> &ecalCells)
{
  ifstream inFile(inputFile);

  if (inFile.fail()) {
    cout << "-E- main : File " << inputFile << " does not exist!" << endl;
    return;
  }
  string line;
  do {
    getline(inFile, line);
  } while (line != "[EmcCalPar]");

  Params params;
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      inFile >> sector >> cell >> params;
      ecalCells[sector][cell].newParams = params;
      if (inFile.eof())
        break;
    }
  }
}

void writeParams(const TString &outputFile, Array<EcalModule> &ecalCells)
{
  ofstream outFile(outputFile);

  if (outFile.fail()) {
    cout << "-E- main : File " << outputFile << " cannot be written!" << endl;
    return;
  }
  outFile << "##############################################################################\n"
          << "##############continue;################################################################\n"
          << "# Calibration parameters for the Emc\n"
          << "# Format:\n"
          << "# sector(0...5)   position(0...162)   TdcSlope  TdcOffset AdcOffset1\n"
          << "# AdcOffset2 AdcSlope   Twc0  Twc1  Twc2 Twc3 Twc4\n"
          << "##############################################################################\n"
          << "[EmcCalPar]\n";
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      outFile << sector << "\t" << cell << "\t" << ecalCells[sector][cell].newParams << "\n";
    }
  }
}

Array<EcalModule> initialize()
{
  Array<EcalModule> ecalCells;

  ecalCells.reserve(6);
  for (Int_t sector = 0; sector < 6; sector++) {
    vector<EcalModule> oneSector;
    oneSector.reserve(163);
    for (Int_t cell = 0; cell < 163; cell++) {
      EcalModule a(sector, cell);
      oneSector.push_back(a);
    }
    ecalCells.push_back(oneSector);
  }
  return ecalCells;
  // gStyle->SetOptStat("n");
}

void calibrateCells(void)
{
  Array<EcalModule> ecalCells = initialize();

  readParams("../EnergyCalibration/mar19EmcCalPar_new.txt", ecalCells);
  TFile fTotThresholds("ToTtreshold_hist.root", "read");
  TH1 * hTotThresholds = (TH1 *) fTotThresholds.Get("hToTtreshold");
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      if (!isGoodCell(sector, cell)) continue;
      Double_t totTreshold = hTotThresholds->GetBinContent(sector * 200 + cell + 1);
      ecalCells[sector][cell].setTotThrehold(totTreshold);
      ecalCells[sector][cell].fillHist(ecalCells);
      ecalCells[sector][cell].calibrate();
    }
  }
  // ===============================================================================================//
  // ==========================================DRAWING==============================================//
  // ===============================================================================================//
  int N = 32;           // number of canvases (Npages needed +1); 4x4->52; 3x3->92
  int row = 3, col = 3; // number of plots in roows and cols on canvas
  TCanvas * canvas[N];


  for (Int_t i = 0; i < N; i++) {
    canvas[i] =
      new TCanvas(Form("canvas%i", i), Form("canvas%i", i), 10, 10, 3000, 3000);
    canvas[i]->Divide(row, col);
  }
  Int_t counter = 0;
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      if (!isGoodCell(sector, cell)) continue;
      canvas[(int) (counter / (row * col))]->cd(counter % (row * col) + 1);
      ecalCells[sector][cell].drawEnergy();
      counter++;
    }
  }

  canvas[0]->SaveAs("timeCalibration.pdf[");
  for (Int_t i = 0; i < N - 1; i++) canvas[i]->SaveAs("timeCalibration.pdf");
  canvas[0]->SaveAs("timeCalibration.pdf]");

  // ===============================================================================================//
  // ===============================================================================================//
  // ===============================================================================================//

  writeParams("mar19EmcCalPar_time.txt", ecalCells);
} // calibrateCells

int timeCalibration()
{
  calibrateCells();
  return 0;
}
