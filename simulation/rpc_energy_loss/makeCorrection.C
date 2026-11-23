#include <iostream>
#include <map>
#include <stdio.h>
#include <TTree.h>  
#include <TSystem.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TString.h>
#include <TTree.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TProfile.h>


int makeCorrection(void)
{

    TFile *in = new TFile("photon_correction.root","read");

    TH2F *corr_2D_thphi;
    TH2F *dev_2D_thphi;
    TH2F *corr_2D_cellEn;
    TH2F *dev_2D_cellEn;
    TH1F *corr_1D_cell;
    TH1F *dev_1D_cell;

    corr_2D_thphi = (TH2F*)in->Get("devMomThetaVsPhi");
    dev_2D_thphi = (TH2F*)in->Get("devMomThetaVsPhiN");
    corr_2D_thphi->Divide(dev_2D_thphi);
    //corr_2D_thphi->Scale(1./100);
    corr_2D_thphi->SetNameTitle("energy_correction2D_ThetaPhi","energy_correction2D_ThetaPhi");

    corr_2D_cellEn = (TH2F*)in->Get("devMomCellVsEne_noABS");
    dev_2D_cellEn = (TH2F*)in->Get("devMomCellVsEneN");
    corr_2D_cellEn->Divide(dev_2D_cellEn);
    corr_2D_cellEn->Scale(1./100);
    corr_2D_cellEn->SetNameTitle("energy_correction2D_CellEnergy","energy_correction2D_CellEnergy");

    corr_1D_cell = (TH1F*)in->Get("devMomCell_noABS");
    dev_1D_cell = (TH1F*)in->Get("devMomCellN");
    corr_1D_cell->Divide(dev_1D_cell);
    //corr_1D_cell->Scale(1./100);
    corr_1D_cell->SetNameTitle("energy_correction1D_CellEnergy","energy_correction2D_CellEnergy");

    TFile *out = new TFile("my_correction_energy_ag1650ag.root","recreate");
    out->cd();
    corr_2D_thphi->Write();
    corr_2D_cellEn->Write();
    corr_1D_cell->Write();
    return 0;
}
