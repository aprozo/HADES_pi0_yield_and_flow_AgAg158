
#include "EcalCellPar.h"
#include "CorrFactors.h"

CorrFactors corrFactors[6][163];
EcalCellPar ecalCellPars[6][163];

void initializeEcalEnergyCorrections(){
	readCorrFactors("AllCorrFactors.root",corrFactors);
	readParams("mar19EmcCalParNewEnergy.txt", ecalCellPars, "old");
}

void correctNonlinearity(Float_t & energy){
	energy*= 0.99367-1.42255/sqrt(energy) + 23.8587/energy -110/pow(energy,1.5)+1571.92/pow(energy,2);
}
void correctBeamDependence(Float_t & energy, const Float_t & beamTime){
	energy*= (155.932-0.261425*beamTime)/134.977;
}

void recalibrate(HEmcCluster * cls, HCategory * emcCalCat,Float_t beamTime=80.1568327436)
{
	Float_t energySum = 0.;
	Int_t sector      = cls->getSector();
	Int_t clusterSize = cls->getNCells();
	for (int i = 0; i < clusterSize; i++) {
		UChar_t c1 = cls->getCell(i);
		Int_t cell = getSystemMap(c1);
		HLocation loc;
		loc.set(2, 0, 0); // set location of cluster
		loc[0] = sector;
		loc[1] = c1;
		HEmcCal * cal = (HEmcCal *) emcCalCat->getObject(loc);
		Float_t TOT   = cal->getWidth();
		Params newParams = ecalCellPars[sector][cell].newParams;
		Float_t newEnergy = newParams.adc0 + exp(newParams.adc1 + newParams.adc2 * TOT); // recalibration of a single module energy within cluster
		energySum += newEnergy;
	}
	energySum*= corrFactors[sector][getSystemMap(cls->getCell())].val[0]; // correct the whole cluster for each cell coefficient
	correctNonlinearity(energySum);
	correctBeamDependence(energySum,beamTime);
	cls->setEnergy(energySum);
}
