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
   typedef struct {
   Float_t mult_tofrpc,cent_bin,charge,emc_match,sector,emc_cell,theta_rpc,phi_rpc,mom,beta,emc_size,
   emc_energy,emc_time,emc_mq,emc_path,mm_qa,within_resolution,beamtime;}
   PID_exp_lepton;

   typedef struct {
   Float_t imppar,cent_bin,charge,emc_match,sector,emc_cell,theta_rpc,phi_rpc,mom,beta,mom_kine,cflag,emc_size,
   emc_energy,emc_time,emc_energy_corr,emc_mq,emc_path,mm_qa,geantid,within_resolution:emc_energy_centr;}
   PID_sim_lepton;

   typedef struct {
   Float_t imppar,cent_bin,emc_hit,sector,emc_cell,theta_kine,phi_kine,mom_kine,emc_size,
   emc_energy,emc_time,emc_beta,emc_energy_corr,emc_to_vert_dist,within_resolution,geantid,matchedCells,matchedPart;}
   PID_sim_photon;

   static PID_exp_lepton pid_exp_lepton;
   static PID_sim_lepton pid_sim_lepton;
   static PID_sim_photon pid_sim_photon;
 */
enum  {sim_photon, sim_lepton, exp_lepton,        nFiles};
enum  {electron, positron, all_leptons, photon,   nParticle};


Int_t particle_color[3]={kAzure+1,kGreen+1,kBlack};

void setup_hist_draw(TH1 * hist,const Int_t & particle){
	hist->SetMarkerStyle(20);
	hist->SetLineColor  (particle_color[particle]);
	hist->SetMarkerColor(particle_color[particle]);
	hist->SetLineWidth(3);
	hist->GetYaxis()->SetTitle("Eff");
	//hist->Smooth();
	hist->GetYaxis()->SetRangeUser(0,1);
	if(hist->GetNbinsX()==40) {
		hist->GetXaxis()->SetTitle("theta, deg");
		hist->GetXaxis()->SetRangeUser(10,50);
	}
	else if(hist->GetNbinsX()==12) {
		hist->GetXaxis()->SetTitle("mom, MeV");
		hist->GetXaxis()->SetRangeUser(0,700);
	}
}

void setup_hist_draw(TH2 * hist){
	hist->GetZaxis()->SetRangeUser(0,1);
	hist->GetZaxis()->SetTitle("Eff");
	hist->GetYaxis()->SetTitle("theta, deg");
	hist->GetXaxis()->SetTitle("mom, MeV");


}

void ttree_make_plots_efficiency(void)
{
	gStyle->SetOptStat(0);       // switch off stat box
	Int_t nEvents =10000000;
	Int_t nSector =6;
	TFile *fileIn[3];
	const int nParams=2;


	fileIn[sim_photon] = new TFile("../sim_photon_ttree_gen5a.root","Read");
	fileIn[sim_lepton] = new TFile("../lepton_ttree_sim_gen3_x6stat_06Aug.root");
	fileIn[exp_lepton] = new TFile("../exp_leptons_ttree_all_days.root","Read");

	TChain *treeIn[nFiles];
	for (Int_t f=0; f<nFiles; f++)  {   treeIn[f]=(TChain*)fileIn[f]->Get("nt");}

	TString params    [nFiles][nParams] ={{"mom_kine","theta_kine"},{"mom","theta_rpc"},{"mom","theta_rpc"}};
	TString params2d  [nFiles] ={"theta_kine:mom_kine","theta_rpc:mom","theta_rpc:mom"};
	TString param_name[nParams] ={"_mom","_theta"};
	Int_t param_bins[nParams] ={ 12, 40};
	Int_t param_low [nParams] ={ 0, 10};
	Int_t param_high[nParams] ={ 700, 50};

	TString fileNames[nFiles]={"sim_photon","sim_lepton","exp_lepton"};
	TCanvas *can[nParams];
	TCanvas *can2d[nFiles];

	TFile *out= new TFile("efficiencyEcalAllSizes.root","RECREATE");


	TH1D *denominator[nParams][nFiles][nSector];
	TH1D * numerator[nParams][nFiles][nSector];

	TH2D *denominator2d[nFiles][nSector];
	TH2D * numerator2d[nFiles][nSector];


	for (Int_t sector=0; sector<nSector; sector++)  {
		for (Int_t par=0; par<nParams; par++)  {
			can[par] = new TCanvas(params[1][par],params[1][par],1200,1000);
			for (Int_t f=0; f<nFiles; f++)  {
				TString hist_name=Form("sector%d",sector)+fileNames[f]+param_name[par]+"_denominator";
				denominator[par][f][sector]=new TH1D(hist_name,hist_name,param_bins[par],param_low[par],param_high[par]);
				hist_name=Form("sector%d",sector)+fileNames[f]+param_name[par]+"_efficiency";
				numerator[par][f][sector]=new TH1D(hist_name,hist_name,param_bins[par],param_low[par],param_high[par]);
			}
		}
	}

	for (Int_t sector=0; sector<nSector; sector++)  {
		for (Int_t f=0; f<nFiles; f++)  {
			can2d[f] = new TCanvas(fileNames[f],fileNames[f],1200,1000);
			TString hist_name=Form("sector%d",sector)+fileNames[f]+"_mom_theta_denominator";
			denominator2d[f][sector]=new TH2D(hist_name,hist_name,param_bins[0],param_low[0],param_high[0],param_bins[1],param_low[1],param_high[1]);
			hist_name=Form("sector%d",sector)+fileNames[f]+"_mom_theta_efficiency";

			numerator2d[f][sector]=new TH2D(hist_name,hist_name,param_bins[0],param_low[0],param_high[0],param_bins[1],param_low[1],param_high[1]);
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////

	TString general_cuts="cent_bin==0";
	TString gen_cuts[nFiles];
	gen_cuts[sim_photon]="mom_kine>=100";
	gen_cuts[exp_lepton]="mom>60&&theta_rpc>=12&&theta_rpc<=45";
	gen_cuts[sim_lepton]="mom>60&&theta_rpc>=12&&theta_rpc<=45";        //same as in the experiment
	//gen_cuts[sim_lepton]="pid.mom>20&&pid.theta_rpc>12&&pid.theta_rpc<45&&((pid.cflag&64)==64)";//maximum possible efficiency


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////-----------------/////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	TString general_emc_cuts="within_resolution>0";
	//cuts for leptons and photons separately///////////////////////
	/////
	TString emc_cuts[nFiles];
	emc_cuts[sim_photon]="(matchedCells==0)&&(matchedPart==0)&&emc_beta>0.95&&emc_beta<1.05  && within_resolution<5.5&&emc_cell>=10";
	emc_cuts[exp_lepton]="emc_match && emc_mq <2 && within_resolution<6.9 && emc_path>0";
	emc_cuts[sim_lepton]="emc_match && emc_mq <2 && within_resolution<5.5 && emc_path>0";     //same as in the experiment
	//emc_cuts[sim_lepton]="pid.true_emc_hit_track";                                 //maximum possible efficiency


	//////////////////////////////////////////////////////////////////////////////////////////////



	TH1D* ratio[2];
	for (Int_t par=0; par<nParams; par++)  {
		can[par]->cd();
		for (Int_t f=0; f<nFiles; f++)  {
			for (Int_t sector=0; sector<nSector; sector++)  {
				if (sector==0||sector==3) continue;
				numerator[par][f][sector]->Sumw2();
				denominator[par][f][sector]->Sumw2();
				TString hist_name=Form("sector%d",sector)+fileNames[f]+param_name[par]+"_denominator";
				TString drawString=params[f][par]+">>"+ hist_name;
				TString cuts=general_cuts+"&&"+gen_cuts[f]+"&&"+Form("sector==%d",sector);
				treeIn[f]->Draw(drawString,cuts,"",nEvents);


				hist_name=Form("sector%d",sector)+fileNames[f]+param_name[par]+"_efficiency";
				drawString=params[f][par]+">>"+ hist_name;
				cuts=general_cuts+"&&"+gen_cuts[f]+"&&"+general_emc_cuts+"&&"+emc_cuts[f]+"&&"+Form("sector==%d",sector);
				treeIn[f]->Draw(drawString,cuts,"",nEvents);

				numerator[par][f][sector]->Divide(denominator[par][f][sector]);
				setup_hist_draw(numerator[par][f][sector],f);
			}
			//ratio[par]=(TH1D*)nominator[par][2]->Clone("mom_lepton_ratio");
			//ratio[par]->Divide(nominator[par][1]);

			//ratio[par]->Draw("same");
		}


	}

	for (Int_t f=0; f<nFiles; f++)  {
		can2d[f]->cd();
		for (Int_t sector=0; sector<nSector; sector++)  {
			if (sector==0||sector==3) continue;
			TString hist_name=Form("sector%d",sector)+fileNames[f]+"_mom_theta_denominator";
			TString drawString=params2d[f]+">>"+ hist_name;
			TString cuts=general_cuts+"&&"+gen_cuts[f]+"&&"+Form("sector==%d",sector);
			treeIn[f]->Draw(drawString,cuts,"",nEvents);

			hist_name=Form("sector%d",sector)+fileNames[f]+"_mom_theta_efficiency";
			drawString=params2d[f]+">>"+ hist_name;
			cuts=general_cuts+"&&"+gen_cuts[f]+"&&"+general_emc_cuts+"&&"+emc_cuts[f]+"&&"+Form("sector==%d",sector);
			treeIn[f]->Draw(drawString,cuts,"",nEvents);
			numerator2d[f][sector]->Divide(denominator2d[f][sector]);

			setup_hist_draw(numerator2d[f][sector]);
		}
	}


	TLegend *leg=new TLegend(0.27,0.10,0.52,0.26);
	leg->SetEntrySeparation(0.1);
	for (Int_t f=0; f<nFiles; f++)  {
		leg->AddEntry(numerator[0][f][2],fileNames[f],"lep");
	}

	for (Int_t par=0; par<nParams; par++)  {
		can[par]->cd();
		for (Int_t f=0; f<nFiles; f++)  {
			numerator[par][f][2]->Draw("same");
		}
		leg->Draw("same");
	}


	for (Int_t f=0; f<nFiles; f++)  {
		can2d[f]->cd();
		numerator2d[f][2]->Draw("colz");
	}

	out->cd();


	for (Int_t sector=0; sector<nSector; sector++)  {
		if (sector==0||sector==3) continue;
		TString hist_name;
		for (Int_t f=0; f<nFiles; f++)  {
			for (Int_t par=0; par<nParams; par++)  {
				hist_name=Form("sector%d",sector)+fileNames[f]+param_name[par]+"_efficiency";
				numerator[par][f][sector]->Write();
			}
			hist_name=Form("sector%d",sector)+fileNames[f]+"_mom_theta_efficiency";
			numerator2d[f][sector]->Write();

		}
	}




}
