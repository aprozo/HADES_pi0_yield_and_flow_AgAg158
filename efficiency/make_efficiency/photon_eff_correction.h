#include "TH3.h"
#include "TFile.h"
#include "TH1.h"

TFile *FileEff= new TFile("efficiencyEcalAllSizes.root","Read");
TH3* photon_eff=(TH3*)FileEff->Get("efficiencyPhiThetaEnergy");

Double_t getEfficiency( const Float_t & phi, const Float_t &  theta, Float_t mom)
{
	if (mom >=650) mom=649;
	if (mom<=150) mom=152;

	Double_t eff= photon_eff->Interpolate(phi,theta,mom);    // factor arising from correcting by exp/sim leptons ratio = 0.86
	return (eff>0&&eff<1) ? eff : 1;
}
