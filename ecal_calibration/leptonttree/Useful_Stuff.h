#include "hrichhit.h"
#include "hrichcal.h"

// -------------------------------------------------------------------
// -------------------------------------------------------------------
// This file contains all used functions for analysis macro
// - Lepton selection Cut functions
// - efficiency correction functions
// -------------------------------------------------------------------
// For systematic studies use (one of) the following variations:
// - Change of the Comversion supression Cut from 25 Cals to 35
// - Change of RICH matching cut to momentum independent at 2 degree
// - Change of M2 cut to momentum independent at 10k MeV2
// -------------------------------------------------------------------
// -------------------------------------------------------------------


Bool_t isConversionCutDiff(HParticleCand* cand, HCategory* RingCat, HCategory* CalCat)
{
    // Function to remove candidates matched with ring
    // with many other Cals around
    // a) Double Ring - Conversion
    // b) Blob - unphysical
    // Cut is based on difference between Ring_Cals and Cals in 3Rad area around Ring center

    Int_t richIndex = cand->getRichInd();
    HRichHit* ring = static_cast<HRichHit*>(RingCat->getObject(richIndex));
    Double_t Posx = ring->getRich700CircleCenterX();
    Double_t Posy = ring->getRich700CircleCenterY();
    Double_t Radius = ring->getRich700CircleRadius();

    Int_t N_ofCals = ring->getRich700NofRichCals();
    Int_t Cal_Counter = 0;
    Int_t size = CalCat->getEntries();
    HRichCal* Cal =0;
    for (Int_t i=0; i < size; i++) {
	Cal = HCategoryManager::getObject(Cal,CalCat,i);
        Double_t x_coord = (Cal->getCol()-96)*6.625;       // translate from col/row in mm
	Double_t y_coord = (Cal->getRow()-96)*6.625;

        Double_t dist_center = sqrt((Posx-x_coord)*(Posx-x_coord)+(Posy-y_coord)*(Posy-y_coord));
        if (dist_center<(2*Radius)) Cal_Counter++;
    }
    //cout << Radius << "  " << Cal_Counter << "    " << N_ofCals << endl;
    if ((Cal_Counter-N_ofCals)<15) return false;
    return true;
}


Bool_t isConversion(HParticleCand* cand, HCategory* RingCat, HCategory* CalCat)
{
    // Function to remove candidates matched with ring
    // with many other Cals around
    // a) Double Ring - Conversion
    // b) Blob - unphysical
    // Cut is only based on Cals in 3Rad area around ring center - double peak structure in spectrum for tweo rings

    Int_t richIndex = cand->getRichInd();
    HRichHit* ring = static_cast<HRichHit*>(RingCat->getObject(richIndex));
    Double_t Posx = ring->getRich700CircleCenterX();
    Double_t Posy = ring->getRich700CircleCenterY();
    Double_t Radius = ring->getRich700CircleRadius();

    Int_t Cal_Counter = 0;
    Int_t size = CalCat->getEntries();
    HRichCal* Cal =0;
    for (Int_t i=0; i < size; i++) {
	//Cal = HCategoryManager::getObject(Cal,CalCat,i);
        HRichCal* Cal = static_cast<HRichCal*>(CalCat->getObject(i));
        Double_t x_coord = (Cal->getCol()-96)*6.625;       // translate from col/row in mm
	Double_t y_coord = (Cal->getRow()-96)*6.625;

        Double_t dist_center = sqrt((Posx-x_coord)*(Posx-x_coord)+(Posy-y_coord)*(Posy-y_coord));
        if (dist_center<(2*Radius)) Cal_Counter++;
    }
    if (Cal_Counter<25) return false;    // 25 original, 30 for consisteny check
    return true;
}





Bool_t isOARing(HParticleCand* cand, HCategory* RingCat){
    Bool_t pairfound = false;
    Double_t Radius_Cut = 18.7473;
    Int_t richIndex = cand -> getRichInd();
    HRichHit* ring = static_cast<HRichHit*>(RingCat->getObject(richIndex));

    Int_t size = RingCat->getEntries();
    for(Int_t j = 0; j < size && pairfound == false; j++)
    {
	if (j == richIndex) continue;
	HRichHit* ringref = static_cast<HRichHit*>(RingCat->getObject(j));
        if (ringref->getRich700CircleRadius()<Radius_Cut) continue;
	Double_t Theta1 = (TMath::Pi() / 180.0) * ring->getTheta();
	Double_t Theta2 = (TMath::Pi() / 180.0) * ringref->getTheta();
        Double_t Phi1 = (TMath::Pi() / 180.0) * ring->getPhi();
        Double_t Phi2 = (TMath::Pi() / 180.0) * ringref->getPhi();


	Double_t opening_angle = (180.0 / TMath::Pi()) * TMath::ACos(TMath::Sin(Theta1)*TMath::Sin(Theta2)*TMath::Cos(Phi1 - Phi2)+TMath::Cos(Theta1)*TMath::Cos(Theta2));

        //cout << "oa:  " <<  opening_angle << endl;
	if (opening_angle > 9) continue;
	else pairfound = true;
    }
    //cout << pairfound << endl;
    return pairfound;
}







