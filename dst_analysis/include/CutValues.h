#pragma once

class CutValues
{

private:
	// Double_t E_beam = 1.58; //nucleon energy in AGeV
	// Double_t m_0 = 0.93827; //nucleon mass in GeV
	// Double_t  Energy = E_beam+m_0;
	// Double_t  p_z =sqrt(Energy*Energy-m_0*m_0);
	// Double_t y1 = 0.5*TMath::Log((Energy+p_z)/(Energy-p_z)); //Beam Rapiidity
	// Double_t y0 = 0.5*y1;                                    //Mid  Rapiidity

  const Double_t midrapidity = 0.822;

	Int_t nPtBins = 24;
  Int_t ptMin  = 0;
  Int_t ptMax  = 1200;

	Int_t nYBins = 8;
  Float_t yMin = 0.15;
  Float_t yMax = 0.95;

	Int_t centMin  = 0;
	Int_t centMax  = 30;
	Int_t nCentralityBins = 3;

	Int_t MassBins = 400;

	Int_t phiBins = 12;

	Int_t ptBinSize   = (ptMax - ptMin) / nPtBins;
	Float_t yBinSize  = (yMax - yMin) / nYBins;



	void setDefaultValues() {

		fOnlySector2and4= kFALSE;
		fOnlySize1= kFALSE;

		// Ecal Gamma Cuts
		fPhotonLowThetaCut= 16.5;
		fUpperEcalBetaCut= 1.3;
		fLowerEcalBetaCut= 0.8;
		fUpperEcalEnergyCut= 1500.;
		fLowerEcalEnergyCut= 100.;
		fEcalCandMatchQuality= 2.5;
	}


public:
	CutValues() {
		setDefaultValues();
	}

	~CutValues(){

	};
	Bool_t fOnlySector2and4;
	Bool_t fOnlySize1;

	//Ecal Photon Cuts
	Double_t fPhotonLowThetaCut;
	Double_t fUpperEcalBetaCut;
	Double_t fLowerEcalBetaCut;
	Double_t fUpperEcalEnergyCut;
	Double_t fLowerEcalEnergyCut;
	Double_t fEcalCandMatchQuality;

	Int_t getNPtBins(){
		return nPtBins;
	}
	Int_t getNYBins(){
		return nYBins;
	}

	Double_t getMinPt(){
		return ptMin;
	}
	Double_t getMaxPt(){
		return ptMax;
	}
	Double_t getMinY(){
		return yMin;
	}
	Double_t getMaxY(){
		return yMax;
	}

  Double_t getMidrapidity(){
    return midrapidity;
  }

	Bool_t getOnlySector2and4(){
		return fOnlySector2and4;
	}

	Bool_t getOnlySize1(){
		return fOnlySize1;
	}

	Double_t getEcalCandMatchQuality(){
		return fEcalCandMatchQuality;
	}

	void setFeb22Cuts()
	{
	}

	// std::vector<HParticleCandSim*> getGoodCandidates(HCategory* candCat){
	// 	std::vector<HParticleCandSim*> GoodCandidatesVec;
	// 	Int_t nofPartCands = candCat->getEntries();
	// 	GoodCandidatesVec.reserve(nofPartCands);
	// 	for(Int_t iP = 0; iP < nofPartCands; iP++) {
	// 		HParticleCandSim* cand = static_cast<HParticleCandSim*>(candCat->getObject(iP));
	// 		if(!cand)
	// 			continue;
	// 		if(!cand->isFlagBit(kIsUsed) )
	// 			//If nonexistent or duplicate continue
	// 			continue;
  //
	// 		if( !cand->isFlagAND(4,
	// 		                     Particle::kIsAcceptedHitInnerMDC,
	// 		                     Particle::kIsAcceptedHitOuterMDC,
	// 		                     Particle::kIsAcceptedHitMETA,
	// 		                     Particle::kIsAcceptedRK )
	// 		    || cand->getInnerSegmentChi2() < 0
	// 		    || cand->getChi2()             > 1000 )
	// 			//If cand is from bad fit or bad track/meta
	// 			continue;
  //
	// 		if(cand->getCharge()==0)
	// 			continue;
  //
	// 		if(!HParticleTool::isGoodMetaCell(cand,4,true) )
	// 			//If cand is not from good meta cell
	// 			continue;
  //
	// 		if((cand->getMetaMatchQuality() > 4.0
	// 		    || cand->getMetaMatchQuality() < 0.0))
	// 			//If cand meta matching is bad
	// 			continue;
  //
	// 		if(cand->getMomentumOrg() > 4000)
	// 			//If momentum too large, skip particle... Since Gen4
	// 			continue;
  //
	// 		GoodCandidatesVec.push_back(cand);
	// 	}
	// 	return GoodCandidatesVec;
	// }
  //
	// std::vector<HParticleCandSim*> getLeptonCandidates(vector<HParticleCandSim*>  const & GoodCands, HCategory* RingCat, HVertex& vertexClust, Double_t RMCut, Double_t RadiusCut, Double_t MassCut){
	// 	std::vector<HParticleCandSim*> LeptonCandidatesVec;
  //
  //
	// 	for (UInt_t iGoodCands=0; iGoodCands<GoodCands.size(); iGoodCands++) {
	// 		// if(GoodCands.at(iGoodCands)->getCharge()!=-1)continue;
  //
	// 		if(!GoodCands.at(iGoodCands)->isFlagBit(kIsLepton)) { continue; }
  //
	// 		if(GoodCands.at(iGoodCands)->getMass()>MassCut) continue;
  //
	// 		Int_t richIndex = GoodCands.at(iGoodCands)->getRichInd();
	// 		HRichHitSim* ring = (richIndex >= 0) ? static_cast<HRichHitSim*>(RingCat->getObject(richIndex)) : NULL;
	// 		Double_t Radius = ring->getRich700CircleRadius();
	// 		if(Radius<RadiusCut) continue;
  //
	// 		if(GoodCands.at(iGoodCands)->getMass()>MassCut) continue;
  //
	// 		Double_t dPhiCorr,dThetaCorr;
	// 		Double_t corRM = getCorrectedRM(GoodCands.at(iGoodCands), ring, vertexClust, dPhiCorr, dThetaCorr);
	// 		//Double_t corRM =1;
	// 		if(corRM>RMCut) continue;
  //
	// 		LeptonCandidatesVec.push_back(GoodCands.at(iGoodCands));
	// 	}
  //
	// 	return LeptonCandidatesVec;
	// }
  //
  //
	// TH1F* GetCutValueHistogram()
	// {
	// 	TH1F* hCutValues = new TH1F("hCutValues","hCutValues; ",9,-0.5,8.5);
	// 	string labelsX[9] = {"fPhotonLowThetaCut","fUpperEcalBetaCut","fLowerEcalBetaCut","fUpperEcalEnergyCut","fLowerEcalEnergyCut","fEcalCandMatchQuality"};
	// 	for (Int_t i = 1; i <= 9; i++) {
	// 		hCutValues->GetXaxis()->SetBinLabel(i, labelsX[i-1].c_str());
	// 	}
  //
	// 	hCutValues->SetBinContent(3,fPhotonLowThetaCut);
	// 	hCutValues->SetBinContent(4,fUpperEcalBetaCut);
	// 	hCutValues->SetBinContent(5,fLowerEcalBetaCut);
	// 	hCutValues->SetBinContent(6,fUpperEcalEnergyCut);
	// 	hCutValues->SetBinContent(7,fLowerEcalEnergyCut);
	// 	hCutValues->SetBinContent(8,fEcalCandMatchQuality);
  //
	// 	return hCutValues;
	// }
  //
	// TH1F* GetBoolHistogram()
	// {
	// 	TH1F*  hBoolValues= new TH1F("hBoolValues","hBoolValues; ",8,-0.5,7.5);
	// 	string labelsX[8] = {"fOnlyFillMixerWithRecLambda","fOnlySector2and4","fOnlySector2","fOnlySize1","setMass","acceptPipP","fUseRelPimMassCut"};
	// 	for (Int_t i = 1; i <= 8; i++) {
	// 		hBoolValues->GetXaxis()->SetBinLabel(i, labelsX[i-1].c_str());
	// 	}
  //
	// 	hBoolValues->SetBinContent(3,(int)fOnlySector2and4);
	// 	hBoolValues->SetBinContent(5,(int)fOnlySize1);
  //
  //
	// 	return hBoolValues;
	// }

};
