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

typedef struct { Float_t eventnumber,downflag,trigdecision,trigbit,mult_track, mult_tofrpc,mult_tof,mult_start,start_time,start_strip,mult_emc, system,sector,metahit,tof,tof_cor,tof_tof,tof_rpc,xpos,dx_tofmdc,   eloss,leftamp,rightamp,theta,phi,seeninrich,emc_mq,mom,beta, beta_cor,chi2,charge,mass,isused,chi2_mdc0,metamatch_quality,  ishadron,emc_cell,emc_energy,emc_time,dedx_mdc_comb,zprime,id,  geantid,cflag,geantparentid,mom_geant,imppar,mult_wall,RP,RPAB,RP_geant;} PID;

//typedef struct { Float_teventnumber,imppar,id,mass,opangle,mom,theta,phi,y,pt,sec1,mom1,theta1,phi1,emc_ene1,sec2,mom2,theta2,phi2,emc_ene2;}PID; 


void fill_hist (void)
{
  TFile *corrInput                      = new TFile("corretion_energy_ag1650ag.root","read");
  //TH1F  *corr1D_energy                  = (TH1F*) corrInput->Get("energy_correction1D_Cell");
  TH2F  *corr2D_energy                   = (TH2F*) corrInput->Get("energy_correction2D_CellEnergy");
  
  TH2F  *energy_loss                    = new TH2F("energy_loss_vs_emom",";p_{e} [MeV]",100,0,1000,100,-500,500);
  TH2F  *energy_correction              = new TH2F("energy_correction_vs_emom",";p_{e} [MeV]",100,0,1000,60,0,300);
  TH2F  *energy_loss_theta              = new TH2F("energy_loss_vs_theta", "Electron mom - Energy in emc, MeV",200,-800,800,50,10,80);
  TH2F  *energy_loss_phi                = new TH2F("energy_loss_vs_phi",   "Electron mom - Energy in emc, MeV",180,0,360,200,-800,800);
  TH2F  *energy_vs_mom                  = new TH2F("energy_vs_mom",";p_{e} [MeV]",200,0,1000,200,0,1000);
  TH2F  *correnergy_vs_mom              = new TH2F("correnergy_vs_mom",";p_{e} [MeV]",200,0,1000,200,0,1000);
  TProfile  *energylossprx; 
  TCanvas *c1                           = new TCanvas("c1","c1",10,10,1200,1200);
  TCanvas *c2                           = new TCanvas("c2","c2",10,10,1200,1200);
  TCanvas *c3                           = new TCanvas("c3","c3",10,10,1200,1200);
  TCanvas *c4                           = new TCanvas("c4","c4",10,10,1200,1200);
  TCanvas *c5                           = new TCanvas("c5","c5",10,10,1200,1200);
  TCanvas *c6                           = new TCanvas("c6","c6",10,10,1200,1200);  
  TString filename = "Ag_Ag_1650MeV_10Mevents_lepton_ttree.root";
  //TString filename = "/lustre/hades/user/hadestof/prop17/200/root/Ag_Ag_1650MeV_1Mevents_diphoton_ttree_new.root.root";

  
  TFile  *inFile = new TFile(filename,"read");
  TString outputFile = "hist_lepton_energy_loss.root";
  TChain *nt = (TChain*)inFile->Get("nt");

    PID pid;
    nt->SetBranchAddress("pid" , &pid.eventnumber);
    Long64_t nentries = nt->GetEntries();
    TFile *out    = new TFile(outputFile,"RECREATE");
    corr2D_energy->RebinY(4);
     corr2D_energy->Scale(1/4.);
     for (Long64_t i = 1; i <= nentries ; i++)
    {
            nt->GetEntry(i);
	    //if(pid.emc_cell!=7) continue;
      
	    if(!(pid.isused>0&&pid.seeninrich>0&&(pid.geantid==2||pid.geantid==3)&&pid.beta>0.95&&pid.emc_mq>0&&pid.emc_mq<2.5))continue;

      //if(pid.mom_geant<170) continue;
      Int_t binX = corr2D_energy->GetXaxis()->FindBin((Int_t)pid.emc_cell);
      Int_t binY = corr2D_energy->GetYaxis()->FindBin(pid.emc_energy);
      
      Float_t corr1D = corr2D_energy->GetBinContent(binX,binY);

      if(!TMath::Finite(corr1D)) corr1D=0.;      
       if(corr1D!=0)
      {
      Float_t enCorr = pid.emc_energy*(1+corr1D);
      //if ()continue;
      
      energy_vs_mom              ->Fill(pid.mom_geant,      pid.emc_energy);      
      energy_loss                ->Fill(pid.mom_geant,      pid.mom_geant-enCorr);
      correnergy_vs_mom          ->Fill(pid.mom_geant,      enCorr);
      energy_correction          ->Fill(pid.mom_geant,      pid.emc_energy*corr1D);      }
      
   }
     


 
     
     c1->cd();
     gPad->SetGrid();
     gPad->SetLogz();
     energy_vs_mom->GetYaxis()->SetTitle("Energy in emc, MeV");
     energy_vs_mom->Draw("colz");
     c1->SaveAs("energy_vs_mom.png","png");
     
     c2->cd();
     gPad->SetGrid();
     gPad->SetLogz();
     energy_loss->GetYaxis()->SetTitle("Energy loss in RPC, MeV");
     energy_loss->Draw("colz");
     energy_loss->ProfileX("energy_loss_pfx",40,100)->Draw("same");
     
     c2->SaveAs("energy_loss.png","png");

     c3->cd();
     gPad->SetGrid();
     gPad->SetLogz();
     correnergy_vs_mom->GetYaxis()->SetTitle("Corrected Energy in emc, MeV");
     correnergy_vs_mom->Draw("colz");
     c3->SaveAs("correnergy_vs_mom.png","png");
     

          c6->cd();
     gPad->SetGrid();
     gPad->SetLogz();
     energy_correction->GetYaxis()->SetTitle("Correction, MeV");
     energy_correction->Draw("colz");
     c6->SaveAs("energy_correction.png","png");
     
out->cd();
     energy_vs_mom->Write();
     energy_loss->Write();
     correnergy_vs_mom->Write();
     energy_correction->Write();
     energy_loss->ProfileX("energy_loss_pfx",40,100)->Write();
out->Save();
 out->Close();
}
