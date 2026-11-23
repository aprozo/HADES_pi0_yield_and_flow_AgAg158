
//----------------------------------------------------//
//           tutaj znajduje sie input i output        //
// /lustre/hades/user/sudol/HADES_iQMD_1M/root    //
//----------------------------------------------------//

#define iqmd_cxx
#include "iqmd.h"
#include <TH2.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include "htool.h"

#include <stdio.h>
#include <iostream>
#include <vector>

#include "my_dalitz.h"

using namespace std;

void iqmd::Loop(Int_t nevt )
{
//   In a ROOT session, you can do:
//      Root > .L iqmd.C
//      Root > iqmd t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

        setup_my_dalitz();

        FILE * pFile;

        Float_t impact_bins[]  = {9.18, 7.95, 6.5, 4.6, 0.};// latest version from WIKI
        TFile *fout = new TFile("hists_iqmd_pi0_test.root","recreate");
        fout->cd();

        TH1F *cosheli1 = new TH1F("heli1","heli1",100,-1,1);
        TH1F *cosheli2 = new TH1F("heli2","heli2",100,-1,1);
        TH2F *iqmd_pcosTheCM[4];
        TH2F *mom_pi0dil = new TH2F("mom_pi0dil","mom_pi0dil",100,0,1,100,0,1);
        TH2F *phi_pi0dil = new TH2F("phi_pi0dil","phi_pi0dil",100,-3,3,100,-3,3);
        TH2F *the_pi0dil = new TH2F("the_pi0dil","the_pi0dil",100,0,6,100,0,6);
        TH2F *mom12 = new TH2F("mom12","mom12",100,0,1,100,0,1);
        TH2F *phi12 = new TH2F("phi12","phi12",100,-3,3,100,-3,3);
        TH2F *the12 = new TH2F("the12","the12",100,0,6,100,0,6);
        TH2F *costhe12 = new TH2F("costhe12","costhe12",100,-1,1,100,-1,1);
        TH1F *p_lep = new TH1F("p_lep","p_lep",100,0,1);
        TH1F *m_diele = new TH1F("m_diele","m_diele",100,0,1);
        TH1F *phiang_pi0 = new TH1F("phi_pi0","phi_pi0",100,0,6);
        TH1F *the_pi0 = new TH1F("the_pi0","the_pi0",100,0,6);
        TH2F *dTheta_dPhi_pi0_diele = new TH2F("dTheta_dPhi_pi0_diele","dTheta_dPhi_pi0_diele",100,-TMath::Pi(),TMath::Pi(),100,TMath::Pi(),TMath::Pi());
        TH2F *hPtY_deuter = new TH2F("hPtY_deuter","hPtY_deuter",110,-2.2,2.2,100,0,2000);
        TH2F *hPtY_carbon = new TH2F("hPtY_carbon","hPtY_carbon",110,-2.2,2.2,100,0,2000);
        TProfile2D *pv1PtY_pi0 = new TProfile2D("pv1PtY_pi0","pv1PtY_pi0",11,-1.1,1.1,10,0,1000);
        TProfile2D *pv2PtY_pi0 = new TProfile2D("pv2PtY_pi0","pv2PtY_pi0",11,-1.1,1.1,10,0,1000);
        TProfile2D *pv1PtY_diele = new TProfile2D("pv1PtY_diele","pv1PtY_diele",11,-1.1,1.1,10,0,1000);
        TProfile2D *pv2PtY_diele = new TProfile2D("pv2PtY_diele","pv2PtY_diele",11,-1.1,1.1,10,0,1000);
        for(Int_t i=0; i<4; i++) {
                iqmd_pcosTheCM[i] = new TH2F(Form("iqmd_pcosTheCM_mult%i",i),Form("iqmd_pcosTheCM_mult%i",i),60,0,1200,40,-1,1);
        }

        if (fChain == 0) return;

        Long64_t nentries = fChain->GetEntriesFast();
        if(nevt>0 && nevt < nentries) nentries = nevt;

        Int_t fileCount =0;
        Int_t evtCount;
        Long64_t nbytes = 0, nb = 0;

//   for (Long64_t jentry=0; jentry<5000;jentry++) {
        for (Long64_t jentry=0; jentry<nentries; jentry++) {

                HTool::printProgress(jentry,nentries,1,"Analyze :");

                Long64_t ientry = LoadTree(jentry);
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);   nbytes += nb;

                Int_t ip=4;       //Hades binning
                if (b <= impact_bins[0] && b > impact_bins[1]) {ip=0;}
                if (b <= impact_bins[1] && b > impact_bins[2]) {ip=1;}
                if (b <= impact_bins[2] && b > impact_bins[3]) {ip=2;}
                if (b <= impact_bins[3] && b > impact_bins[4]) {ip=3;}

                //  if(ip >= 4) continue;

                evtCount++;

                if(jentry%1000==0) {
                        evtCount = 1;
                        if(jentry!=0) fclose (pFile);
                        pFile = fopen (Form("/lustre/hades/user/harabasz/flowtest/evt/iqmd_test_evt_%i.evt",fileCount),"w");
                        fileCount++;
                }

                Int_t npi0 = 0;
                for(Int_t p =0; p<npart; p++) {
                        if (imod[p] == 7) { ++npi0; }
                }

                fprintf (pFile, "%i  %i  1.25  %f  1 \n",evtCount, npart + 2*npi0, b );

                Float_t w = 1.0;
                for(Int_t p =0; p<npart; p++) {
                        //	  if(imod[p] !=7 && imod[p] !=8) continue;

//	   if(chrg[p] < 3) cout << "imod : "<<imod[p] <<endl;
                        if(chrg[p]== 2 && mass[p] <= 3) imod[p] = 49.;
                        if(chrg[p]== 2 && mass[p] >  3) imod[p] = 47.;
                        if(chrg[p]== 3 && mass[p] <= 6) imod[p] = 61.;
                        if(chrg[p]== 3 && mass[p] >  6) imod[p] = 62.;
                        if(chrg[p]== 4 && mass[p] <= 7) imod[p] = 63.;
                        if(chrg[p]== 4 && mass[p] >  8) imod[p] = 64.;
                        if(chrg[p]== 5 && mass[p] <=10) imod[p] = 65.;
                        if(chrg[p]== 5 && mass[p] > 11) imod[p] = 66.;
                        if(chrg[p]== 6) imod[p] = 67.;
                        if(chrg[p]== 7) imod[p] = 68.;
                        if(chrg[p]== 8) imod[p] = 69.;
                        if(chrg[p]== 9) imod[p] = 70.;
                        if(chrg[p]==10) imod[p] = 71.;
                        if(chrg[p]==11) imod[p] = 72.;
                        if(chrg[p]==12) imod[p] = 73.;
                        if(chrg[p] > 12) imod[p] = chrg[p] +  61.;


                        TLorentzVector p_cm;
//	   p_cm.Boost(0.,0.,-0.629); // go to NN CM system for 1230MeV

                        // note: 108 protons + ; T = 1.58 GeV; restEnergy = 0.938 GeV





                        p_cm.SetXYZT(px[p],py[p],pz[p],etot[p]);
                        if (imod[p] == 45) {
                                hPtY_deuter->Fill(p_cm.Rapidity()-0.74,p_cm.Perp()*1000);
                        }
                        if (imod[p] == 67) {
                                hPtY_carbon->Fill(p_cm.Rapidity()-0.74,p_cm.Perp()*1000);
                                cout << "Rest mass of carbon: " << p_cm.M() << endl;
                        }
                        if (imod[p] == 7) {

                                TLorentzVector pi0 = p_cm;

                                pv1PtY_pi0->Fill(p_cm.Rapidity()-0.74,p_cm.Perp()*1000,TMath::Cos(p_cm.Phi())); // Psi_RP is always 0 in this data sample
                                pv2PtY_pi0->Fill(p_cm.Rapidity()-0.74,p_cm.Perp()*1000,TMath::Cos(2*p_cm.Phi()));

                                Float_t px;
                                Float_t theta_pi0 = p_cm.Theta();
                                Float_t phi_pi0 = p_cm.Phi();
                                std::vector<TLorentzVector> leptons;
                                int s = 0;
                                do {
                                        leptons = my_dalitz(p_cm);
                                        px = leptons[0].Px();
                                        // cout << s << " " << px << endl;
                                        ++s;
                                } while (!TMath::Finite(px));
                                p_cm = leptons[0];
                                fprintf (pFile, "%f  %f  %f  %f  %i  %2.1f \n", p_cm.E(),p_cm.Px(),p_cm.Py(),p_cm.Pz(), 2, w);
                                p_cm = leptons[1];
                                fprintf (pFile, "%f  %f  %f  %f  %i  %2.1f \n", p_cm.E(),p_cm.Px(),p_cm.Py(),p_cm.Pz(), 3, w);
                                p_lep->Fill(p_cm.Vect().Mag());
                                p_cm = leptons[2];
                                fprintf (pFile, "%f  %f  %f  %f  %i  %2.1f \n", p_cm.E(),p_cm.Px(),p_cm.Py(),p_cm.Pz(), 1, w);
                                w += 1.0;

                                TLorentzVector diele = leptons[0] + leptons[1];
                                Float_t theta_diele = diele.Theta();
                                Float_t phi_diele = diele.Phi();

                                pv1PtY_diele->Fill(diele.Rapidity()-0.74,diele.Perp()*1000,TMath::Cos(diele.Phi())); // Psi_RP is always 0 in this data sample
                                pv2PtY_diele->Fill(diele.Rapidity()-0.74,diele.Perp()*1000,TMath::Cos(2*diele.Phi()));

                                m_diele->Fill(diele.M());
                                dTheta_dPhi_pi0_diele->Fill((phi_diele-phi_pi0)*TMath::Sin(theta_pi0),theta_diele-theta_pi0);
                                mom12->Fill(leptons[0].Vect().Mag(),leptons[1].Vect().Mag());
                                phi12->Fill(leptons[0].Phi(),leptons[1].Phi());
                                the12->Fill(leptons[0].Theta(),leptons[1].Theta());
                                costhe12->Fill(TMath::Cos(leptons[0].Theta()),TMath::Cos(leptons[1].Theta()));
                                mom_pi0dil->Fill(diele.Vect().Mag(),pi0.Vect().Mag());
                                phi_pi0dil->Fill(diele.Phi(),pi0.Phi());
                                the_pi0dil->Fill(diele.Theta(),pi0.Theta());
                                phiang_pi0->Fill(phi_pi0);
                                the_pi0->Fill(theta_pi0);

                                TLorentzVector dilep_cm(diele);
//               dilep_cm.Boost(0.,0.,-0.629); // go to NN CM system for 1230MeV
                                dilep_cm.Boost(-pi0.BoostVector()); // go to NN CM system for 1230MeV
                                Float_t thetaCM = dilep_cm.Theta()*TMath::RadToDeg();

                                bool swapHeli = (rnd.Uniform() > 0.5);

                                TLorentzVector ev4_cm(leptons[0]);
                                TLorentzVector pv4_cm(leptons[1]);
//               ev4_cm.Boost(0.,0.,-0.629);
//               pv4_cm.Boost(0.,0.,-0.629);
                                ev4_cm.Boost(-pi0.BoostVector());
                                pv4_cm.Boost(-pi0.BoostVector());
                                TLorentzVector ev4_dilep(ev4_cm);
                                TLorentzVector pv4_dilep(pv4_cm);
                                ev4_dilep.Boost(-dilep_cm.BoostVector());
                                pv4_dilep.Boost(-dilep_cm.BoostVector());
                                Float_t heli1 = ev4_dilep.Angle(dilep_cm.Vect());
                                Float_t heli2 = pv4_dilep.Angle(dilep_cm.Vect());
                                if (swapHeli) {
                                        Float_t temp = heli1;
                                        heli1 = heli2;
                                        heli2 = temp;
                                }
                                Float_t cosHeli1 = TMath::Cos(heli1);
                                Float_t cosHeli2 = TMath::Cos(heli2);

                                // cosheli1->Fill(TMath::Cos(pv4_dilep.Theta()));
                                // cosheli2->Fill(TMath::Cos(ev4_dilep.Theta()));
                                cosheli1->Fill(cosHeli1);
                                cosheli2->Fill(cosHeli2);

                                // cout << "helicity c: " << TMath::Cos(leptons[0].Theta()) << endl;;
                                // cout << "helicity d: " << TMath::Cos(leptons[1].Theta()) << endl;

                        } else {

                                fprintf (pFile, "%f  %f  %f  %f  %i  1.0 \n", p_cm.E(),p_cm.Px(),p_cm.Py(),p_cm.Pz(), (int)imod[p]);
                        }
                } // loopa po czastkach
        } // loopa po evencie
        fclose (pFile);

        fout->cd();
        fout->GetList()->Write();
/*
   iqmd_pcosTheCM[0]->Write();
   iqmd_pcosTheCM[1]->Write();
   iqmd_pcosTheCM[2]->Write();
   iqmd_pcosTheCM[3]->Write();
   p_lep->Write();
   dTheta_dPhi_pi0_diele->Write();
   pv1PtY->Write();
   pv2PtY->Write();
 */
        fout->Close();


}
