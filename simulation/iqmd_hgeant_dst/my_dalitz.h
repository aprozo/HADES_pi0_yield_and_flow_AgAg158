#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"
#include "hphysicsconstants.h"

TRandom3 rnd;
TF1 f_heli("f_heli","1+x*x",-1,1);
TH1F *invmasdil_mult0;


void setup_my_dalitz() {
//========= Macro generated from object: invmasdil_mult0/
//========= by ROOT version5.34/34
        Double_t xAxis1[39] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.11, 0.13, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.48, 0.51, 0.54, 0.57, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.399};

        invmasdil_mult0 = new TH1F("invmasdil_mult0","",38, xAxis1);
        invmasdil_mult0->SetBinContent(1,5.45068e+09);
        invmasdil_mult0->SetBinContent(2,1.74017e+09);
        invmasdil_mult0->SetBinContent(3,9.55862e+08);
        invmasdil_mult0->SetBinContent(4,6.15152e+08);
        invmasdil_mult0->SetBinContent(5,4.1709e+08);
        invmasdil_mult0->SetBinContent(6,2.84479e+08);
        invmasdil_mult0->SetBinContent(7,1.90492e+08);
        invmasdil_mult0->SetBinContent(8,1.22066e+08);
        invmasdil_mult0->SetBinContent(9,7.32671e+07);
        invmasdil_mult0->SetBinContent(10,2.8452e+07);
        invmasdil_mult0->SetBinContent(11,3.26355e+06);
        invmasdil_mult0->SetBinContent(12,6450);
        invmasdil_mult0->SetBinError(1,744176);
        invmasdil_mult0->SetBinError(2,417153);
        invmasdil_mult0->SetBinError(3,309170);
        invmasdil_mult0->SetBinError(4,248022);
        invmasdil_mult0->SetBinError(5,204228);
        invmasdil_mult0->SetBinError(6,168665);
        invmasdil_mult0->SetBinError(7,138019);
        invmasdil_mult0->SetBinError(8,110483);
        invmasdil_mult0->SetBinError(9,85596.2);
        invmasdil_mult0->SetBinError(10,37717.4);
        invmasdil_mult0->SetBinError(11,12774.1);
        invmasdil_mult0->SetBinError(12,567.891);
        invmasdil_mult0->SetEntries(1e+08);
}

std::vector<TLorentzVector> my_dalitz(TLorentzVector pi0) {
        Double_t phi_X_g = rnd.Uniform(-TMath::Pi(),TMath::Pi());
        Double_t theta_X_g = TMath::ACos(rnd.Uniform(-1,1));

        Double_t m_g = 0;

        Double_t M = HPhysicsConstants::mass(7)/1000; // pi0 mass in GeV
        Double_t E_X_g = 0.5*M; // PDG (18.16)
        Double_t p_X_g = E_X_g;

        TVector3 v3_gamma_X;
        v3_gamma_X.SetMagThetaPhi(p_X_g, theta_X_g, phi_X_g);

        TLorentzVector photon1_X(v3_gamma_X, p_X_g);
        TLorentzVector photon2_X(-v3_gamma_X, p_X_g); // opposite direction, E = cp for real photon

        // cout << "theta 1a : " << e1_ee.Theta() << endl;
        // cout << "theta 2a : " << e2_ee.Theta() << endl;
        // cout << "t1-t2-pi : " << e1_ee.Theta() + e2_ee.Theta() - TMath::Pi() <<  endl;

        TVector3 newUz;
        newUz.SetMagThetaPhi(1,gamma_X.Theta(),gamma_X.Phi());

        e1_ee.RotateUz(newUz);
        e2_ee.RotateUz(newUz);
        e1_ee.Boost(gamma_X.BoostVector());
        e2_ee.Boost(gamma_X.BoostVector());


        e1_ee.Boost(pi0.BoostVector());
        e2_ee.Boost(pi0.BoostVector());
        photon_X.Boost(pi0.BoostVector());

        std::vector<TLorentzVector> leptons;
        leptons.push_back(e1_ee);
        leptons.push_back(e2_ee);
        leptons.push_back(photon_X);

        if (!TMath::Finite(e1_ee.X())) {
                cout << "- " << E_X_g << " " << p_X_g << " " << m_g << " " << px << " " << e1_ee.Px() << " " << p_e1_ee << " " << theta_e_ee << " " << phi_e_ee << endl;
        }
        return leptons;

}
