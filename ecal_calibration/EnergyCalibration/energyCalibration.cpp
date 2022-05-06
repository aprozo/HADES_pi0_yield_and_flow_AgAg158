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
#include "TFitResultPtr.h"
#include "TLine.h"

#include "exp_setup_cells.h"

using namespace std;


// if(pid.beamtime<88.46) AgAg158
// if(pid.beamtime>88.46&&pid.beamtime<90.8) AgAg123@2500A
// if(pid.beamtime>90.8)	AgAg123@200A


// ===========================================================
// Data structures
// ============================================================
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

// ===========================================================
// ECalModule Class
// ============================================================
class EcalModule {
private:
  TH2D * hEnergy;  // size 1 with artificial dot (forcing fit to pass it)
  TF1 * oldEnergy; // old parameters fit
  TF1 * fitEnergy; // new fit, artificial dot used


  TString name;
  Int_t sector, cell;
  Double_t totDot; // artificial dot

public:
  Params oldParams, newParams;
  TH2D * hEnergySize2;

  EcalModule(const Int_t &newSector, const Int_t &newCell) : sector(newSector), cell(newCell)
  {
    Int_t nbinx = 350;
    Int_t nbiny = 200;
    TString fitFuncEnergy = "[0] + exp([1]+[2] * x) ";

    name      = Form("Sector%i_Cell%i", sector, cell);
    hEnergy   = new TH2D(name + "Energy", name + "Energy; ToT; mom", nbinx, 250, 1500, nbiny, -100, 1500); // tot vs mom
    oldEnergy = new TF1(name + "OldEnergy", fitFuncEnergy, 200, 1500);
  }

  void setTotThrehold(const Float_t & TotThreshold)
  {
    totDot = TotThreshold;
    TString fitFuncEnergy2Par =
      (TString) Form("20 + [0] + exp([1]) *( exp([2] * x) - exp([2]* %.d) )", TotThreshold);
    fitEnergy = new TF1(name + "fitEnergy", fitFuncEnergy2Par, 200, 1500);
    fitEnergy->SetParLimits(0, -10, 5);
    fitEnergy->SetParLimits(1, 0, 10);
    fitEnergy->SetParLimits(2, 0, 0.02);
  }





  void fillHist(Array<EcalModule> &ecalCells)
  {
    TChain * nt = new TChain(Form("Sector%iCell%i", sector, cell));

    nt->Add(Form("/lustre/hebe/hades/user/prozorov/mar19/leptonttree/splitted/tree_Sector%i_Cell%i.root", sector,
      cell));
    PID pid;
    nt->SetBranchAddress("pid", &pid);
    Long_t nEntry         = nt->GetEntries();
    Double_t cutOffBorder = 150;

    for (Int_t iEntry = 0; iEntry < nEntry; ++iEntry) {
      nt->GetEntry(iEntry);
      Float_t energy1, energy2;
      energy1 = pid.maxenergy;
      if (pid.size == 1) energy2 = 0;
      else energy2 = pid.energy - pid.maxenergy;


      if ((pid.emc_mq < 0.6 || (pid.mom > 600 && pid.emc_mq < 1.3)) && pid.emc_path > 90 && pid.richQa < 1.2) {
        if (!(pid.TOT > (totDot + 150) || pid.beamtime > 90.8)) continue;
        Float_t energydiff        = (100 * abs(pid.mom - (energy1 + energy2)) / pid.mom);
        Float_t within_resolution = energydiff * sqrt(pid.mom / 1000.) / 3;
        if (within_resolution > 10) continue;

        if (pid.size == 1) hEnergy->Fill(pid.TOT, pid.mom);
        //  else if (energy1 > 400 && pid.size == 2 && abs(pid.time - pid.time2) < 1.5) {
        //   ecalCells[sector][pid.cell2].hEnergySize2->Fill(pid.TOT2, pid.mom - energy1);
        // }
      }
    }


    for (int i = 1; i <= hEnergy->GetYaxis()->GetNbins(); i++) {
      Double_t centerY  = ((TAxis *) hEnergy->GetYaxis())->GetBinCenter(i);
      TH1D * maxValProj = (TH1D *) hEnergy->ProjectionX("maxValProj", i - 5, i + 5);
      Int_t maxXbin     = maxValProj->GetMaximumBin();
      Double_t centerYX = ((TAxis *) maxValProj->GetXaxis())->GetBinCenter(maxXbin);

      for (int j = 1; j <= hEnergy->GetXaxis()->GetNbins(); j++) {
        if (hEnergy->GetBinContent(j, i) > 0 && centerY > 600 &&
          abs(centerYX - ((TAxis *) hEnergy->GetXaxis())->GetBinCenter(j)) < 100)
          hEnergy->SetBinContent(j, i, hEnergy->GetBinContent(j, i) * 15);
      }
    }


    delete nt;
  } // fillHist

  void calibrate()
  {
    oldEnergy->SetParameter(0, oldParams.adc0);
    oldEnergy->SetParameter(1, oldParams.adc1);
    oldEnergy->SetParameter(2, oldParams.adc2);

    // fitEnergy->SetParameter(0, oldParams.adc0);
    fitEnergy->SetParameter(1, oldParams.adc1);
    fitEnergy->SetParameter(2, oldParams.adc2);
    fitEnergy->SetRange(totDot - 10, 1500);

    hEnergy->Fit(fitEnergy, "REQNM");
    cout << "Sector " << sector << " Cell " << cell << " done!" << endl;

    newParams.adc0 = 20 + fitEnergy->GetParameter(0) - exp(fitEnergy->GetParameter(1) )
      * exp(fitEnergy->GetParameter(2) * totDot);
    newParams.adc1 = fitEnergy->GetParameter(1);
    newParams.adc2 = fitEnergy->GetParameter(2);
  }

  void drawCrossPoint()
  {
    TLine * lh = new TLine(250, 20, 1500, 20);
    TLine * lv = new TLine(totDot, -100, totDot, 1500);

    lh->SetLineColor(kBlack);
    lv->SetLineColor(kBlack);
    lh->SetLineWidth(2);
    lv->SetLineWidth(2);
    lh->Draw("same");
    lv->Draw("same");
  }

  void drawEnergy()
  {
    gStyle->SetOptStat("n");
    gPad->SetGrid(1);

    hEnergy->Draw("colz same");
    drawCrossPoint();
    oldEnergy->SetLineColor(kGreen + 2);
    oldEnergy->Draw("same");
    fitEnergy->Draw("same");

    TLegend * legend = new TLegend(0.1, 0.7, 0.42, 0.85);
    legend->AddEntry(oldEnergy, "old Fit", "L");
    legend->AddEntry(fitEnergy, "new Fit 2 par(ToT treshold)", "L");

    legend->SetLineWidth(0);
    legend->Draw("same");
  }

  void saveCrosscheck()
  {
    hEnergy->Write();
  }
}; // end of a class

// ======================================================================
// Other functions
// ======================================================================

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
}

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

void writeParams(const TString &outputFile, Array<EcalModule> &ecalCells)
{
  ofstream outFile(outputFile);

  if (outFile.fail()) {
    cout << "-E- main : File " << outputFile << " cannot be written!" << endl;
    return;
  }

  outFile << "##############################################################################\n"
          << "##############################################################################\n"
          << "# Calibration parameters for the Emc\n"
          << "# Format:\n"
          << "# sector(0...5)   position(0...162)   TdcSlope  TdcOffset AdcOffset1\n"
          << "# AdcOffset2 AdcSlope   Twc0  Twc1  Twc2"
          << "##############################################################################\n"
          << "[EmcCalPar]\n";
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      outFile << sector << "\t" << cell << "\t" << ecalCells[sector][cell].newParams << "\n";
    }
  }
}

void calibrateCells(void)
{
  Array<EcalModule> ecalCells = initialize();

  readParams("mar19EmcCalPar.txt", ecalCells);
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
  int N = 52;           // number of canvases (Npages needed +1); 4x4->52; 3x3->92
  int row = 3, col = 3; // number of plots in roows and cols on canvas
  TCanvas * canvasEnergy[N];


  for (Int_t i = 0; i < N; i++) {
    canvasEnergy[i] =
      new TCanvas(Form("canvasEnergy%i", i), Form("canvasEnergy%i", i), 10, 10, 3000, 3000);
    canvasEnergy[i]->Divide(row, col);
  }
  Int_t counter = 0;
  for (Int_t sector = 0; sector < 6; sector++) {
    for (Int_t cell = 0; cell < 163; cell++) {
      if (!isGoodCell(sector, cell)) continue;
      canvasEnergy[(int) (counter / (row * col))]->cd(counter % (row * col) + 1);
      ecalCells[sector][cell].drawEnergy();
      counter++;
    }
  }

  canvasEnergy[0]->SaveAs("energyCalibration.pdf[");
  for (Int_t i = 0; i < N - 1; i++) canvasEnergy[i]->SaveAs("energyCalibration.pdf");
  canvasEnergy[0]->SaveAs("energyCalibration.pdf]");

  // ===============================================================================================//
  // ===============================================================================================//
  // ===============================================================================================//

  writeParams("mar19EmcCalPar_new.txt", ecalCells);
} // calibrateCells

// ======================================================================
// Main
// ======================================================================

void calibrateOneCell(void)
{
  Array<EcalModule> ecalCells = initialize();

  readParams("mar19EmcCalPar.txt", ecalCells);
  TFile fTotThresholds("ToTtreshold_hist.root", "read");
  TH1 * hTotThresholds = (TH1 *) fTotThresholds.Get("hToTtreshold");
  Int_t sector         = 4;
  Int_t cell = 44;

  Double_t totTreshold = hTotThresholds->GetBinContent(sector * 200 + cell + 1);
  ecalCells[sector][cell].setTotThrehold(totTreshold);
  ecalCells[sector][cell].fillHist(ecalCells);
  ecalCells[sector][cell].calibrate();

  // ===============================================================================================//
  // ==========================================DRAWING==============================================//
  // ===============================================================================================//
  int N = 1;            // number of canvases (Npages needed +1); 4x4->52; 3x3->92
  int row = 1, col = 1; // number of plots in roows and cols on canvas
  TCanvas * canvasEnergy[N];


  for (Int_t i = 0; i < N; i++) {
    canvasEnergy[i] =
      new TCanvas(Form("canvasEnergy%i", i), Form("canvasEnergy%i", i), 10, 10, 3000, 3000);
    canvasEnergy[i]->Divide(row, col);
  }
  Int_t counter = 0;


  canvasEnergy[(int) (counter / (row * col))]->cd(counter % (row * col) + 1);
  ecalCells[sector][cell].drawEnergy();
  counter++;

  canvasEnergy[0]->SaveAs("energyCalibration.pdf[");
  canvasEnergy[0]->SaveAs("energyCalibration.pdf");
  canvasEnergy[0]->SaveAs("energyCalibration.pdf]");

  // ===============================================================================================//
  // ===============================================================================================//
  // ===============================================================================================//
} // calibrateCells

// ======================================================================
// Main
// ======================================================================


int energyCalibration()
{
  // calibrateOneCell();
  calibrateCells();
  return 0;
}
