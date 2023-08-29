#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"

#include "TString.h"
#include "TSystem.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "addFiles.h"
using namespace std;




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

	Pid(Particle p)
	{
		TLorentzVector pion(p.Px, p.Py, p.Pz, p.E);

		pt      = pion.Pt();        // in Mev
		y       = pion.Rapidity();
		phi     = pion.Phi()*TMath::RadToDeg();
	}

	Float_t pt, y, phi, centBin;
};


struct Event {
	Event()
	{
		particles.clear();
	};

	Header head;
	vector<Pid> particles;
	Int_t nPi0;
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
{
	Float_t w;

	if (stream) {
		stream >> particle.E >> particle.Px >> particle.Py >> particle.Pz >> particle.id >> w;
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
		event.nPi0 = 0;
		for (Int_t i = 0; i < event.head.nParticles; i++) {
			Particle particle;
			stream >> particle;
			if (particle.id == 7) {
				event.nPi0++;
				Pid pi0(particle);
				event.particles.push_back(pi0);
			}
		}
	}
	return stream;
}

// ========================================================================================//


Float_t impact_bins[] = { 0., 3.8, 5.4, 6.6, 7.6 }; // latest version from WIKI

Int_t getCentBin(const Double_t & b)
{
	Int_t impBin = 5;     // Hades binning

	if (b >= impact_bins[0] && b < impact_bins[1]) impBin = 1;
	if (b >= impact_bins[1] && b < impact_bins[2]) impBin = 2;
	if (b >= impact_bins[2] && b < impact_bins[3]) impBin = 3;
	if (b >= impact_bins[3] && b < impact_bins[4]) impBin = 4;
	return impBin;
}

void makeTtreePi0Acceptance(TString inFile, TString outputFile, Int_t nev=100000000.)
{
	vector<TString> fileArray=addMultFiles(inFile);
	TFile * output = new TFile(outputFile, "recreate");
	output->cd();
	Pid pid;
	TTree nt("nt", "pid");
	nt.Branch("pid", &pid, "pt:y:phi:centBin");

	enum Counters_cent {    cent00_10,    cent10_20,    cent20_30, cent30_40,   all_cents,  cNumCounters_cent};
	TH1D* Counter_cent = new TH1D( "nEventsCentrality", "nEventsCentrality", cNumCounters_cent, 0, cNumCounters_cent);


	Counter_cent->GetXaxis()->SetBinLabel( 1, "cent00_10" );
	Counter_cent->GetXaxis()->SetBinLabel( 2, "cent10_20" );
	Counter_cent->GetXaxis()->SetBinLabel( 3, "cent20_30" );
	Counter_cent->GetXaxis()->SetBinLabel( 4, "cent30_40" );
	Counter_cent->GetXaxis()->SetBinLabel( 5, "all_cent" );


	for (Int_t i = 0; i < fileArray.size(); i++) {

		ifstream in(fileArray.at(i).Data());

		if (!in) cout << "An error occurred opening the file" << endl;
		if (i % 1 == 0) cout << i << endl;
		while (!in.eof()&&in) {
			Event event;
			in >> event;
			pid.centBin = getCentBin(event.head.impPar);
			if (pid.centBin>4) continue;
			Counter_cent->Fill(pid.centBin);
			Counter_cent->Fill(all_cents);
			for (Int_t j = 0; j < event.nPi0; j++) {
				Pid pion = event.particles[j];
				pion.centBin = pid.centBin;
				pid = pion;
				nt.Fill();
			}
		}
		in.close();
	}

	nt.Write();
	Counter_cent->Write();
} // makeTtree
