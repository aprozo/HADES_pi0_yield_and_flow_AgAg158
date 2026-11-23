#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TH2F.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// ===========================Acceptance==================================================//
TH2F * acceptance_mar19;
Bool_t isInAcc_mar19(const Double_t phi, const Double_t theta)
{
        UInt_t bin = acceptance_mar19->FindBin(180 + phi, theta);
        Float_t w  = acceptance_mar19->GetBinContent(bin);
        if (w > 0.8&&theta>12&&theta<47) return kTRUE;

        return kFALSE;
}

// ========================================================================================//


struct Header {
        Int_t eventN;
        Int_t nParticles;
        Double_t colEnergy;
        Double_t impPar;
};

struct Particle {

        Particle(){
        }
        Double_t Px, Py, Pz, E;
        Int_t id;
};

struct Pid {
        Pid(){
        }

        Pid(Particle p1,Particle p2)
        {
          TLorentzVector  photon1(p1.Px, p1.Py, p1.Pz,p1.E);
          TLorentzVector  photon2(p2.Px, p2.Py, p2.Pz,p2.E);
          TLorentzVector  pi0=photon1+photon2;
          //photon1.Boost(-pi0.BoostVector()); // go to CMS
          //photon2.Boost(-pi0.BoostVector()); // go to CMS
          pt      = pi0.Pt();// in Mev
          y       = pi0.Rapidity();

          theta1  = 180 * photon1.Theta() / TMath::Pi();
          phi1    = 180 * photon1.Phi() / TMath::Pi();
          theta2  = 180 * photon2.Theta() / TMath::Pi();
          phi2    = 180 * photon2.Phi() / TMath::Pi();
          if (isInAcc_mar19(phi1, theta1)&&isInAcc_mar19(phi2, theta2)) isInAcc = 1;
          else isInAcc = 0;

        }

        Float_t pt, y, isInAcc, theta1, phi1, theta2, phi2;
};


struct Event {
        Header head;
        Pid pi0;
};

// =======================================INPUT============================================//
// ========================================================================================//
istream& operator >> (istream& stream, Header& head)
{
        Float_t w;
        if (stream) {
                stream >> head.eventN >> head.nParticles >> head.colEnergy >> head.impPar >> w;
        }
        return stream;
}

istream& operator >> (istream& stream, Particle& particle)
{ if (stream) {
          Float_t temp1,temp2,temp3,temp4;
          stream >> particle.E >> particle.Px >> particle.Py >> particle.Pz >> particle.id
                 >>temp1>>temp2>>temp3>>temp4;
          particle.E*=1000;
          particle.Px*=1000;
          particle.Py*=1000;
          particle.Pz*=1000;
        }
        return stream;
}

istream& operator >> (istream& stream, Event& event)
{
        if (stream) {
                stream >> event.head;
                Particle photon1,photon2;
                stream >> photon1;
                stream >> photon2;
                Pid pi0(photon1,photon2);
                event.pi0=pi0;
        }
        return stream;
}

// ========================================================================================//



void analyze_evt_file()
{

   TString outputFile="output.root";
   Int_t nev=1000000;
   TFile * acc = new TFile("acceptance_gen4.root", "Read");
   acceptance_mar19 = (TH2F *) acc->Get("Acceptance2D_expsetup");

    TFile * output = new TFile(outputFile, "recreate");
    output->cd();
    Pid pid;
    TTree nt("nt", "pid");
    nt.Branch("pid", &pid, "pt:y:isInAcc:theta1:phi1:theta2:phi2");

for (Int_t i=0;i<7;i++){
    ifstream in(Form(Form("pi0_10mil_evts_%i.evt",i)));
    if (!in) cout << "An error occurred opening the file" << endl;

    while (!in.eof()&&in) {

      Event event;
      in >> event;
      pid = event.pi0;
      nt.Fill();
    }
    in.close();
  }
    nt.Write();

} // makeTtree
