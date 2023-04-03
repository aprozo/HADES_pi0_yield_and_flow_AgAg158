#include <iostream>
#include <map>
#include <stdio.h>
#include "TH3F.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPad.h"
#include "TTree.h"
#include "TChain.h"
#include "smoothMatrix.C"

using namespace std;

/*
 * typedef struct {
 * Float_t mult_tofrpc,cent_bin,charge,emc_match,sector,emc_cell,theta_rpc,phi_rpc,mom,beta,emc_size,
 * emc_energy,emc_time,emc_mq,emc_path,mm_qa,within_resolution,beamtime;}
 * PID_exp_lepton;
 *
 * typedef struct {
 * Float_t imppar,cent_bin,charge,emc_match,sector,emc_cell,theta_rpc,phi_rpc,mom,beta,mom_kine,cflag,emc_size,
 * emc_energy,emc_time,emc_energy_corr,emc_mq,emc_path,mm_qa,geantid,within_resolution:emc_energy_centr;}
 * PID_sim_lepton;
 *
 *
 *
 * static PID_exp_lepton pid_exp_lepton;
 * static PID_sim_lepton pid_sim_lepton;
 *
 */
enum  { sim_lepton, exp_lepton };
enum  { electron, positron, all_leptons, nParticle };


Int_t particle_color[3] = { kAzure + 1, kGreen + 1, kBlack };

void setup_hist_draw(TH1 * hist, const Int_t & particle)
{
  hist->SetMarkerStyle(20);
  hist->SetLineColor(particle_color[particle]);
  hist->SetMarkerColor(particle_color[particle]);
  hist->SetLineWidth(3);
  hist->GetYaxis()->SetTitle("Eff");

  hist->GetYaxis()->SetRangeUser(0, 1);
  if (hist->GetNbinsX() == 40) {
    hist->GetXaxis()->SetTitle("theta, deg");
    hist->GetXaxis()->SetRangeUser(10, 50);
  } else if (hist->GetNbinsX() == 12) {
    hist->GetXaxis()->SetTitle("mom, MeV");
    hist->GetXaxis()->SetRangeUser(0, 700);
  }
}

void analizeLeptonTtreesDraw(void)
{
  gStyle->SetOptStat(0); // switch off stat box
  TFile * fileIn[2];
  const int nQuantity = 1;
  TH1::SetDefaultSumw2;
  const int nEvents = 1000000;
  fileIn[sim_lepton] = new TFile("../sim_photon_ttree_gen5a.root", "Read");
  fileIn[exp_lepton] = new TFile("../exp_lepton_ttree_gen5a.root", "Read");

  TChain * treeIn[2];
  treeIn[0] = (TChain *) fileIn[0]->Get("nt");
  treeIn[1] = (TChain *) fileIn[1]->Get("nt");
  TString quantity[nQuantity]  = { "mom" };
  Int_t binNumber[nQuantity]   = { 10 };// , 40 };
  Int_t binMinLimit[nQuantity] = { 0 };// , 10 };
  Int_t binMaxLimit[nQuantity] = { 1000 };// , 50 };

  TString fileNames[2] = { "sim_lepton", "exp_lepton" };
  TFile * out = new TFile("effLeptons.root", "RECREATE");

  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////


  TString gen_cuts[2];

  gen_cuts[exp_lepton] = " mom>100&&theta_rpc>=12&&theta_rpc<=45 ";
  gen_cuts[sim_lepton] = "  mom>100&&theta_rpc>=12&&theta_rpc<=45 "; // same as in the experiment
  // gen_cuts[sim_lepton]="pid.mom>20&&pid.theta_rpc>12&&pid.theta_rpc<45&&((pid.cflag&64)==64)";//maximum possible efficiency


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////-----------------/////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////
  TString emc_cuts[2];
  // emc_cuts[exp_lepton] = " within_resolution<6.9 && emc_path>0 && (emc_mq <2. || abs(qualityDTime+0.35)<0.5)";
  // emc_cuts[sim_lepton] = " within_resolution<5.5 && emc_path>0 && (emc_mq <2. || abs(qualityDTime)<0.5)"; // same as in the experiment
  // emc_cuts[exp_lepton] = "within_resolution<6.9&& emc_path>0 ";
  // emc_cuts[sim_lepton] = "within_resolution<5.5&& emc_path>0 "; // same as in the experiment

   emc_cuts[exp_lepton] = "within_resolution<6.9 && emc_mq <10.0 && emc_mq >=0 ";
   emc_cuts[sim_lepton] = "within_resolution<5.5 && emc_mq <10.0 && emc_mq >=0 "; // same as in the experiment
  // emc_cuts[sim_lepton]="pid.true_emc_hit_track";                                 //maximum possible efficiency


  //////////////////////////////////////////////////////////////////////////////////////////////

  const Int_t centBins = 3;
  TCanvas * can[nQuantity][centBins];
  TH1D * denominator[nQuantity][2][centBins];
  TH1D * numerator[nQuantity][2][centBins];

  for (Int_t cent = 2; cent < 3; cent++) {
    TString general_cuts = Form("cent_bin<3", cent);
    // TString general_cuts = Form("cent_bin>=0", cent);
    for (Int_t SimExp = 0; SimExp < 2; SimExp++) {
      for (Int_t iQuantity = 0; iQuantity < nQuantity; iQuantity++) {
        TString hist_name, drawString, cuts;

        hist_name = fileNames[SimExp] + quantity[iQuantity] + "_denominator" + Form("_cent%i", cent);
        denominator[iQuantity][SimExp][cent] =
          new TH1D(hist_name, hist_name, binNumber[iQuantity], binMinLimit[iQuantity], binMaxLimit[iQuantity]);


        drawString = quantity[iQuantity] + ">>" + hist_name;
        cuts       = general_cuts + "&&" + gen_cuts[SimExp];


        treeIn[SimExp]->Draw(drawString, cuts, "");

        hist_name = fileNames[SimExp] + quantity[iQuantity] + "_efficiency" + Form("_cent%i", cent);
        numerator[iQuantity][SimExp][cent] = new TH1D(hist_name, hist_name, binNumber[iQuantity],
            binMinLimit[iQuantity], binMaxLimit[iQuantity]);
        drawString = quantity[iQuantity] + ">>" + hist_name;
        cuts       = general_cuts + "&&" + gen_cuts[SimExp] + "&&" + emc_cuts[SimExp];

        cout << drawString << "  " << cuts << endl;
        treeIn[SimExp]->Draw(drawString, cuts, "");

         numerator[iQuantity][SimExp][cent]->Divide(denominator[iQuantity][SimExp][cent]);
         //setup_hist_draw(numerator[iQuantity][SimExp][cent], SimExp);
        numerator[iQuantity][SimExp][cent]->GetYaxis()->SetRangeUser(0,1);
        numerator[iQuantity][SimExp][cent]->SetLineColor(SimExp==0?kRed:kBlue);
        numerator[iQuantity][SimExp][cent]->SetLineWidth(2);
        numerator[iQuantity][SimExp][cent]->SetMarkerStyle(20);
        numerator[iQuantity][SimExp][cent]->SetTitle("0-30%;E [MeV]; #epsilon");


      }
    }
    Float_t effExp = (Float_t) numerator[0][exp_lepton][cent]->GetEntries()
      / (Float_t) denominator[0][exp_lepton][cent]->GetEntries();
    Float_t effSim = (Float_t) numerator[0][sim_lepton][cent]->GetEntries()
      / (Float_t) denominator[0][sim_lepton][cent]->GetEntries();
    cout << "Cent " << cent << "   efficiency Exp :  " << effExp << endl;
    cout << "Cent " << cent << "   efficiency Sim :  " << effSim << endl;
    cout << "Cent " << cent << "   efficiency correction :  " << effExp / effSim << endl;
  }


  out->cd();
  for (Int_t cent = 2; cent < 3; cent++) {
    for (Int_t iQuantity = 0; iQuantity < nQuantity; iQuantity++) {
      for (Int_t SimExp = 0; SimExp < 2; SimExp++) {
        numerator[iQuantity][SimExp][cent]->Write();
        denominator[iQuantity][SimExp][cent]->Write();
      }
    }
  }
  out->Save();
  out->Close();

} // analizeLeptonTtreesDraw
