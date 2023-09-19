#include <stdio.h>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include <iomanip> // std::setw
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TLatex.h"

#include "../GlobalParameters.h"

#include "TPolyLine.h"
#include "TPolyMarker.h"
#define TCOND(a, b, c) ((a) ? (b) : (c))
template <typename T, typename U>
U TernaryComp(T v, T c, U rs, U re, U rl)
{
	return TCOND(v < c, rs, TCOND(v > c, rl, re));
} ///< return the appropriate value for v < c (rs), v == c (re) and v > c (rl)

void set_nice_color5()
{
	new TColor(2001, (66. / 255.), (98. / 255.), (255. / 255.)); // blue
	new TColor(2002, (255. / 255.), (89. / 255.), (74. / 255.)); // red
	new TColor(2003, (25. / 255.), (170. / 255.), (25. / 255.)); // grenn
	new TColor(2004, (255. / 255.), (166. / 255.), (33. / 255.));
	new TColor(2005, (0. / 255.), (170. / 255.), (255. / 255.));
	new TColor(2006, (153. / 255.), (0. / 255.), (153. / 255.));
}

Int_t hexToColorIndex(Int_t hex)
{
	Int_t ri = (hex & 0xff0000) >> 16;
	Int_t gi = (hex & 0xff00) >> 8;
	Int_t bi = (hex & 0xff);
	Int_t ci = TColor::GetColor(hex);

	new TColor(ci, ri / 255., gi / 255., bi / 255.);
	return ci;
}

Double_t DrawHADESLogoText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t o)
{
	auto x = [&](Double_t x)
	{
		return (x / 1208.);
	};
	auto y = [&](Double_t y)
	{
		return (y / -312. + 1.);
	}; // Inverted because pixel coordinates are also inverted!

	TVirtualPad *vpCur = gPad;
	TPad *pHADESLogo = new TPad(Form("p%sHADESLogo", vpCur->GetName()), "", x1, y1, x2, y2);

	pHADESLogo->Draw();
	pHADESLogo->cd();

	TPolyLine *plHADESLogo = new TPolyLine();
	plHADESLogo->SetFillStyle(1001);
	plHADESLogo->SetFillColor(kBlack);
	plHADESLogo->SetLineWidth(0.);

	auto DrawPolyLine = [&](Int_t n, Double_t *x, Double_t *y, Option_t *option)
	{
		plHADESLogo->DrawPolyLine(n, x, y, option);
		delete[] x;
		delete[] y;
	};

	DrawPolyLine(4, new Double_t[4]{x(0.), x(1208.), x(1208.), x(0.)},
				 new Double_t[4]{y(0.), y(0.), y(58. + o), y(58. + o)}, "f"); // Upper Bar
	DrawPolyLine(4, new Double_t[4]{x(0.), x(1208.), x(1208.), x(0.)},
				 new Double_t[4]{y(254. - o), y(254. - o), y(312.), y(312.)}, "f"); // Lower Bar

	DrawPolyLine(4, new Double_t[4]{x(10. - o), x(78. + o), x(78. + o), x(10. - o)},
				 new Double_t[4]{y(68. - o), y(68. - o), y(244. + o), y(244. + o)}, "f"); // H Left Part
	DrawPolyLine(8,
				 new Double_t[8]{x(142. - o), x(214. + o), x(214. + o), x(142. - o), x(142. - o), x(94. - o), x(94. - o),
								 x(142. - o)},
				 new Double_t[8]{y(68. - o), y(68. - o), y(244. + o), y(244. + o), y(179. + o), y(179. + o), y(128. - o),
								 y(128. - o)},
				 "f"); // H Right Part

	DrawPolyLine(11,
				 new Double_t[11]{x(240. - o), x(330. - o), x(394. + o), x(490. + o), x(412. - o), x(404. - o), x(340. - o),
								  x(360. - o), x(382. - o), x(364.), x(314. + o)},
				 new Double_t[11]{y(244. + o), y(68. - o), y(68. - o), y(244. + o), y(244. + o), y(225. + o), y(225. + o),
								  y(183. - o), y(183. - o), y(144. + o), y(244. + o)},
				 "f"); // A

	DrawPolyLine(4, new Double_t[4]{x(516. - o), x(588. + o), x(588. + o), x(516. - o)},
				 new Double_t[4]{y(68. - o), y(68. - o), y(244. + o), y(244. + o)}, "f"); // D Left Part
	DrawPolyLine(29,
				 new Double_t[29]{x(604. - o), x(641.), x(670. + o), x(691. + o), x(709. + o), x(723. + o), x(731. + o), x(734. + o), x(734. + o), x(730. + o), x(722. + o), x(709. + o), x(693. + o), x(677. + o), x(653.), x(604. - o),
								  x(604. - o), x(626.), x(645. - o), x(654. - o), x(662. - o), x(668. - o), x(671. - o), x(671. - o), x(667. - o), x(659. - o), x(643. - o), x(627.), x(604. - o)},
				 new Double_t[29]{y(68. - o), y(71. - o), y(78. - o), y(87. - o), y(99. - o), y(115. - o), y(132.), y(147.), y(166.), y(186.), y(204. + o), y(220. + o), y(231. + o), y(238. + o), y(244. + o), y(244. + o), y(203. - o),
								  y(202. - o), y(197. - o), y(192. - o), y(184. - o), y(172.), y(160.), y(153.), y(143.),
								  y(134. + o), y(125. + o), y(120. + o), y(116. + o)},
				 "f"); // D Right Part

	DrawPolyLine(4, new Double_t[4]{x(764. - o), x(834. + o), x(834. + o), x(764. - o)},
				 new Double_t[4]{y(68. - o), y(68. - o), y(244. + o), y(244. + o)}, "f"); // E Left Part
	DrawPolyLine(4, new Double_t[4]{x(850. - o), x(962. + o), x(962. + o), x(850. - o)},
				 new Double_t[4]{y(68. - o), y(68. - o), y(116. + o), y(116. + o)}, "f"); // E Upper Right Part
	DrawPolyLine(4, new Double_t[4]{x(850. - o), x(945. + o), x(945. + o), x(850. - o)},
				 new Double_t[4]{y(132. - o), y(132. - o), y(180. + o), y(180. + o)}, "f"); // E Middle Right Part
	DrawPolyLine(4, new Double_t[4]{x(850. - o), x(962. + o), x(962. + o), x(850. - o)},
				 new Double_t[4]{y(196. - o), y(196. - o), y(244. + o), y(244. + o)}, "f"); // E Lower Right Part

	DrawPolyLine(13,
				 new Double_t[13]{x(1110. - o), x(1144.), x(1161. + o), x(1173. + o), x(1182. + o), x(1188. + o), x(1192. + o),
								  x(1192. + o), x(1142. - o), x(1135. - o), x(1128. - o), x(1119.), x(1110. - o)},
				 new Double_t[13]{y(68. - o), y(68. - o), y(72. - o), y(79. - o), y(87. - o), y(96. - o), y(109.), y(119. + o),
								  y(119. + o), y(113. + o), y(109. + o), y(107. + o), y(106. + o)},
				 "f"); // S Upper Right Part
	DrawPolyLine(12,
				 new Double_t[12]{x(992. - o), x(1060. + o), x(1070. + o), x(1079.), x(1094. + o), x(1094. + o), x(1036. - o),
								  x(1023. - o), x(1015. - o), x(1007. - o), x(999. - o), x(992. - o)},
				 new Double_t[12]{y(190. - o), y(190. - o), y(201. - o), y(205. - o), y(207. - o), y(244. + o), y(244. + o),
								  y(236. + o), y(229. + o), y(221. + o), y(211. + o), y(198.)},
				 "f"); // S Lower Left Part
	DrawPolyLine(48,
				 new Double_t[48]{x(1094. + o), x(1094. + o), x(1084.), x(1075. + o), x(1069. + o), x(1069. + o), x(1077. + o),
								  x(1093. + o), x(1117. + o), x(1145. + o), x(1169. + o), x(1182. + o), x(1193. + o), x(1197. + o), x(1198. + o),
								  x(1198. + o), x(1197. + o), x(1192. + o), x(1182. + o), x(1171.), x(1157.), x(1110. - o),
								  x(1110. - o), x(1122.), x(1131. - o), x(1135. - o), x(1135. - o), x(1132. - o), x(1126. - o), x(1116. - o),
								  x(1101. - o),
								  x(1078. - o), x(1061. - o), x(1046. - o), x(1032. - o), x(1022. - o), x(1015. - o), x(1005. - o), x(1000. - o),
								  x(998. - o), x(998. - o), x(1000. - o), x(1005. - o), x(1014. - o), x(1026. - o), x(1037. - o), x(1049.), x(1064.)},
				 new Double_t[48]{y(68. - o), y(106. + o), y(107. + o), y(110. + o), y(114. + o), y(120. - o), y(125. - o),
								  y(130. - o), y(134. - o), y(138. - o), y(145. - o), y(152. - o), y(162. - o), y(172.), y(178.), y(205.), y(211.), y(221. + o), y(232. + o), y(239. + o), y(244. + o), y(244. + o), y(207. - o), y(207. - o), y(203. - o),
								  y(197.), y(194.), y(187. + o), y(182. + o), y(178. + o), y(176. + o), y(174. + o), y(171. + o), y(167. + o),
								  y(162. + o), y(157. + o), y(152. + o), y(141. + o), y(132.), y(124.), y(109.), y(102.), y(93. - o), y(84. - o),
								  y(77. - o), y(73. - o), y(70. - o), y(68. - o)},
				 "f"); // S Center Part

	delete plHADESLogo;
	vpCur->cd();
	return x2 - x1;
} // DrawHADESLogoText

////////////////////////////////////////////////////////////////////////////////
/// Draws the HADES logo text in the current pad at the NDC position [x,y].
/// The width and height are calculated based on the parameter s and the parameter
/// a defines the alignment of the logo with respect to [x,y]. Both are used
/// identical to the corresponding parameters for text objects (See TAttText).
/// Furthermore, a different aspect ratio ar can be specified.
/// For more details see also the actual implementation of Tools::DrawHADESLogoText().

Double_t DrawHADESLogoText(Double_t x, Double_t y, Double_t s, Int_t a, Double_t o, Double_t ar)
{
	Double_t padw = gPad->XtoPixel(gPad->GetX2()), padh = gPad->YtoPixel(gPad->GetY1());
	Float_t hpx = s * TCOND(s < 3., TMath::Min(padw, padh), 1.), wpx = hpx * ar;

	return DrawHADESLogoText(TernaryComp(a / 10, 2, x, x - wpx / padw / 2., x - wpx / padw),
							 TernaryComp(a % 10, 2, y, y - hpx / padh / 2., y - hpx / padh),
							 TernaryComp(a / 10, 2, x + wpx / padw, x + wpx / padw / 2., x),
							 TernaryComp(a % 10, 2, y + hpx / padh, y + hpx / padh / 2., y), o);
}

using namespace std;

vector<TString> getFileNames(TString separatedList)
{
	vector<TString> out;

	separatedList.ReplaceAll("../parameters_files/parameters_", "");
	separatedList.ReplaceAll(".conf", "");
	TObjArray *arr = separatedList.Tokenize(" ");
	Int_t n = arr->GetEntries();
	if (n == 0)
	{
		cout << Form("No file in input string '%s' found !", separatedList.Data()) << endl;
	}
	else
	{
		for (Int_t i = 0; i < n; i++)
		{
			TString name = ((TObjString *)arr->At(i))->GetString();
			name.ReplaceAll(" ", "");
			out.push_back(name);
		}
	}
	delete arr;
	return out;
}

template <typename T>
T *computeSystematics(T ***h_diff, std::vector<vector<Int_t>> sysGroupVariationIndex,
					  std::vector<TString> sysGroupName)
{
	T *output = (T *)h_diff[0][0]->Clone((TString) "total_systemtics_" + (TString)h_diff[0][0]->GetName());
	TDirectory *current = gDirectory;

	// TDirectory * dir=current->mkdir((TString)h_diff[0][0]->GetName());
	// dir->cd();

	for (Int_t iY = 1; iY <= output->GetNbinsY(); iY++)
	{
		for (Int_t iX = 1; iX <= output->GetNbinsX(); iX++)
		{
			Int_t bin = output->GetBin(iX, iY);
			Double_t total_sys_err = 0.;
			Int_t nGroups = sysGroupVariationIndex.size();
			TH1D *allGroupsSystematics = new TH1D((TString)h_diff[0][0]->GetName() + Form("y%ipt%i_sys", iX,
																						  iY),
												  "", nGroups, 0, nGroups); // put all variations into a single 1d histogram to average the error with const fit

			for (Int_t iGroup = 0; iGroup < sysGroupVariationIndex.size(); iGroup++)
			{
				allGroupsSystematics->GetXaxis()->SetBinLabel(iGroup + 1, sysGroupName[iGroup]);

				// ================one systematics group error (e.g. BETA_CUT1 and BETA_CUT2)====================//
				Int_t nVariations = sysGroupVariationIndex[iGroup].size();
				TDirectory *dirGroup = (TDirectory *)current->Get(sysGroupName[iGroup]);
				if (!dirGroup)
					dirGroup = current->mkdir(sysGroupName[iGroup]);
				dirGroup->cd();
				TH1D *singleGroupVariations = new TH1D((TString)h_diff[0][0]->GetName() + Form("y%ipt%i", iX,
																							   iY),
													   "", nVariations, 0,
													   nVariations); // put all variations into a single 1d histogram to average the error with const fit
				for (Int_t iVariation = 0; iVariation < nVariations; iVariation++)
				{
					singleGroupVariations->GetXaxis()->SetBinLabel(iVariation + 1,
																   Form("%i", sysGroupVariationIndex[iGroup][iVariation]));
					Double_t diffVal = abs(h_diff[iGroup][iVariation]->GetBinContent(bin)); // difference between reference and current variation, 0 - for ref

					if (h_diff[iGroup][iVariation] != NULL && diffVal != 1 && diffVal < 0.9)
						singleGroupVariations->SetBinContent(iVariation + 1, diffVal);
					singleGroupVariations->SetBinError(iVariation + 1, h_diff[iGroup][iVariation]->GetBinError(bin));
				}

				TF1 *fitConstant = new TF1("fitConstant", "pol0", 0, 20); // fit for 1d histogram for computeulating systemtical error within one group of changes, maximum of 20 changes
				fitConstant->SetParameter(0, singleGroupVariations->GetMaximum());
				singleGroupVariations->Fit("fitConstant", "QRMNE");
				singleGroupVariations->Write();
				Double_t sys_err = fitConstant->GetParameter(0);
				delete singleGroupVariations;
				// ===========================================================================================//
				total_sys_err += sys_err * sys_err; // add single group error to total
				allGroupsSystematics->SetBinContent(iGroup + 1, sys_err);
			} // iGroup
			total_sys_err = sqrt(total_sys_err);
			current->cd();
			TDirectory *allGroupsDir = (TDirectory *)current->Get(h_diff[0][0]->GetName() + (TString) "_sys_bin");
			if (!allGroupsDir)
				allGroupsDir = current->mkdir(h_diff[0][0]->GetName() + (TString) "_sys_bin");
			allGroupsDir->cd();
			allGroupsSystematics->Write();
			current->cd();
			output->SetBinContent(bin, total_sys_err);
		}
	}

	return output;
} // computeSystematics

// ===============================================================================================//

template <typename T>
void setSystematicErrors(T *valueHist, T *errorHist)
{
	for (Int_t iY = 1; iY <= valueHist->GetNbinsY(); iY++)
	{
		for (Int_t iX = 1; iX <= valueHist->GetNbinsX(); iX++)
		{
			Int_t bin = valueHist->GetBin(iX, iY);
			valueHist->SetBinError(bin, valueHist->GetBinContent(bin) * errorHist->GetBinContent(bin));
		}
	}
}

template <typename T>
T calculate(T a, T b)
{
	return a + b;
}

template <typename T>
T *compute(TFile **in, std::vector<vector<Int_t>> sysGroupVariationIndex, std::vector<TString> sysGroupName,
		   TString hName, const Int_t &indexRef)
{
	// ========================define array of histograms difference===========//
	T *h_ref = (T *)in[indexRef]->Get(hName);

	if (h_ref == NULL)
	{
		std::cout << (TString) "no histogram" + hName + in[indexRef]->GetName() << '\n';
		return NULL;
	}
	Int_t nGroups = sysGroupVariationIndex.size();
	T **h_diff[nGroups];

	for (Int_t iGroup = 0; iGroup < nGroups; iGroup++)
	{
		Int_t nVariations = sysGroupVariationIndex[iGroup].size();
		h_diff[iGroup] = new T *[nVariations];
		for (Int_t iVariation = 0; iVariation < nVariations; iVariation++)
		{
			Int_t index = sysGroupVariationIndex[iGroup][iVariation];
			h_diff[iGroup][iVariation] = (T *)in[index]->Get(hName)->Clone(hName);
			if (h_diff[iGroup][iVariation] == NULL)
			{
				std::cout << (TString) "no histogram" + hName + in[index]->GetName() << '\n';
				return NULL;
			}
			h_diff[iGroup][iVariation]->Add(h_ref, -1.); // subtract fo
			h_diff[iGroup][iVariation]->Divide(h_ref);
		} // iVariation
	}	  // iGroup
		  // ========================================================================//

	T *h_sys_errors = computeSystematics(h_diff, sysGroupVariationIndex, sysGroupName);
	h_sys_errors->SetName(h_sys_errors->GetName() + (TString) "_only_errors");
	h_sys_errors->Write();
	T *h_out = (T *)h_ref->Clone(h_ref->GetName() + (TString) "_sys");
	setSystematicErrors(h_out, h_sys_errors);
	delete h_sys_errors;
	return h_out;
}

// ===============================================================================================//
// ===============================================================================================//

void ana_results(TString inputList, TString outputfile = "systematics.root")
{
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	TString nameRef = "REF0";
	//	nameRef      = "ENERGY_CUT2";

	Int_t flowMethod = 0; // 0 -standard, 1 - quadrants , 2 -tprofile
	TString nameMethod[] = {"standard", "quarters", "tprofile"};

	Int_t centBins = 3;

	Int_t ptBins = 4;
	Float_t ptMin = 100;
	Float_t ptMax = 900;
	Float_t ptBinSize = (ptMax - ptMin) / Float_t(ptBins);

	Int_t yBins = 4;
	Float_t yMin = 0.1;
	Float_t yMax = 0.9;
	Float_t yBinSize = (yMax - yMin) / Float_t(yBins);

	std::vector<TString> parSetsNames = getFileNames(inputList); // read names of different systematic variations
	Int_t nParametersSets = parSetsNames.size();				 // total number of files
	if (nParametersSets == 1)
		nameRef = parSetsNames[0];
	TFile *in[nParametersSets];
	TFile *inFlow[nParametersSets + 3];
	std::vector<GlobalParameters> par;

	std::map<TString, vector<Int_t>> parSetMap;		   // from name to file index
	std::vector<vector<Int_t>> sysGroupVariationIndex; // dimensions of variation by groups with file index
	std::vector<TString> sysGroupName;				   //
	Int_t indexRef = 0;								   // reference file number

	for (Int_t i = 0; i < nParametersSets; i++)
	{
		TString filename = "../results/parameters_" + parSetsNames[i] + "/histograms_to_pt_y/results.root";
		in[i] = new TFile(filename, "READ");

		TString filenameFlow = "../results/parameters_" + parSetsNames[i] + Form("/histograms_to_pt_y/flow%i.root", flowMethod);
		inFlow[i] = new TFile(filenameFlow, "READ");

		TString name = parSetsNames[i];
		if (name == nameRef)
			indexRef = i;
		name.Chop();																	// chop last index for puting it into the map, e.g. BETA_CUT1 -> BETA_CUT
		parSetMap[name].push_back(i);													// save the reference index
		par.push_back({"../parameters_files/parameters_" + parSetsNames[i] + ".conf"}); // remeber the paramtere file , currently not needed
	}

	// ==================make a vector out of map=====================//
	for (auto const &elem : parSetMap)
	{
		sysGroupName.push_back(elem.first);
		sysGroupVariationIndex.push_back(elem.second);
		// =====check Map=====/
		std::cout << "Name of group: " << elem.first << endl;
		for (auto const &index : elem.second)
			std::cout << in[index]->GetName() << std::endl;
		std::cout << std::endl;
		// =====check Map=====/
	}

	//	==================add 3 flow methods to groups=====================//
	std::vector<vector<Int_t>> sysGroupVariationIndexFlow = sysGroupVariationIndex; // dimensions of variation by groups with file index
	std::vector<TString> sysGroupNameFlow = sysGroupName;

	for (Int_t iMethod = 0; iMethod < 3; iMethod++)
	{
		TString filename = "../results/parameters_" + parSetsNames[indexRef] + Form("/histograms_to_pt_y/flow%i.root", iMethod);
		inFlow[nParametersSets + iMethod] = new TFile(filename, "READ");
		TString name = nameMethod[iMethod]; // parSetsNames[indexRef] + Form("method%i", iMethod);
	}

	sysGroupVariationIndexFlow.push_back({nParametersSets});
	sysGroupNameFlow.push_back(nameMethod[0]);

	sysGroupVariationIndexFlow.push_back({nParametersSets + 1});
	sysGroupNameFlow.push_back(nameMethod[1]);

	sysGroupVariationIndexFlow.push_back({nParametersSets + 2});
	sysGroupNameFlow.push_back(nameMethod[2]);

	//	==================flow map check=====================//
	std::cout << "========================================= " << endl;

	for (Int_t i = 0; i < sysGroupNameFlow.size(); i++)
	{
		std::cout << "Name of group: " << sysGroupNameFlow[i] << endl;
		for (Int_t j = 0; j < sysGroupVariationIndexFlow[i].size(); j++)
		{
			std::cout << "Index: " << sysGroupVariationIndexFlow[i][j] << endl;
		}
	}

	Int_t colors[] =
		{kOrange, kPink + 2, kGreen, kAzure + 5, kGray, kTeal + 2, kCyan, kYellow, kRed, kBlue, kMagenta};

	// Draw yield hists //
	std::cout << "==========Draw yield vs Y============== " << endl;
	TCanvas *can = new TCanvas("can", "can", 800, 600);
	can->SaveAs("check.pdf[");
	TCanvas *canY[centBins];
	// =======================================================================//
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		canY[iCent] = new TCanvas(Form("canY%i", iCent), Form("canY%i", iCent), 158, 28, 900, 918);
		canY[iCent]->cd();

		auto legend = new TLegend(0.2, 0.15, 0.45, 0.45);
		Int_t groupCounter = 0;
		for (auto const &elem : parSetMap)
		{
			// std::map<TString, vector<Int_t> > parSetMap;
			Int_t variationCounter = 0;
			for (auto const &index : elem.second)
			{
				TH1D *h_y = (TH1D *)in[index]->Get(Form("h_y_cent%i", iCent));
				if (h_y == NULL)
					cout << "no histogram h_y in file " << in[index]->GetName() << endl;
				h_y->GetYaxis()->SetTitle("dN/dy [1/event]");
				if (iCent == 0)
					h_y->GetYaxis()->SetRangeUser(0., 10);
				else if (iCent == 1)
					h_y->GetYaxis()->SetRangeUser(0., 5);
				else
					h_y->GetYaxis()->SetRangeUser(0., 3);
				h_y->SetLineColor(colors[groupCounter] + variationCounter);
				legend->AddEntry(h_y, parSetsNames[index], "l");

				if (index == indexRef)
					h_y->SetLineColor(kBlack);
				h_y->SetLineWidth(2);
				h_y->DrawClone("same hist");
				variationCounter++;
			}
			groupCounter++;
		}
		legend->Draw("same");

		canY[iCent]->SaveAs("check.pdf");
		canY[iCent]->Clear();
	}

	std::cout << "==========Draw yield vs Pt============== " << endl;
	TLatex *latex = new TLatex();

	latex->SetTextFont(42);
	latex->SetTextAlign(12);
	latex->SetTextSize(0.045);

	TCanvas *canPtY[centBins];
	// =======================================================================//
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		canPtY[iCent] = new TCanvas(Form("canPtY%i", iCent), Form("canPtY%i", iCent), 158, 28, 900, 918);
		canPtY[iCent]->cd();

		for (Int_t iY = 0; iY <= 8; iY++)
		{
			auto legend2 = new TLegend(0.7, 0.55, 0.85, 0.95);
			TString yTitle = Form("%.2f < y_{cm} < %.2f", iY * 0.1 + 0.15, (iY + 1) * 0.1 + 0.15);
			Int_t groupCounter = 0;
			for (auto const &elem : parSetMap)
			{
				// std::map<TString, vector<Int_t> > parSetMap;
				Int_t variationCounter = 0;
				for (auto const &index : elem.second)
				{
					TH2D *h_pt_y = (TH2D *)in[index]->Get(Form("h_pt_y_cent%i", iCent));
					if (h_pt_y ==
						NULL)
						std::cout << "not found! index = " << index << "name file = " << in[index]->GetName() << '\n';
					TH1D *h_y = h_pt_y->ProjectionY(Form("py%i%i%i", groupCounter, variationCounter, iCent), iY + 1, iY + 1);
					h_y->GetYaxis()->SetNdivisions(504);
					h_y->GetYaxis()->SetTitleOffset(1.4);
					h_y->SetLineColor(colors[groupCounter] + variationCounter);
					h_y->GetYaxis()->SetTitle("d^{2}N/p_{t}dy [1/event]");
					if (index == indexRef)
						h_y->SetLineColor(kBlack);
					legend2->AddEntry(h_y, parSetsNames[index], "l");
					h_y->SetLineWidth(2);
					h_y->GetXaxis()->SetRangeUser(0.00001, 800);
					h_y->GetYaxis()->SetRangeUser(0, 1.05 * h_y->GetMaximum());
					if (groupCounter == 0)
						h_y->DrawClone("hist");
					else
						h_y->DrawClone("same hist");
					variationCounter++;
				}
				groupCounter++;
			}

			for (auto const &elem : parSetMap)
			{
				for (auto const &index : elem.second)
				{
					if (index == indexRef)
					{
						TH2D *h_pt_y = (TH2D *)in[index]->Get(Form("h_pt_y_cent%i", iCent));
						TH1D *h_y = h_pt_y->ProjectionY(Form("pyRef%i", iCent), iY + 1, iY + 1);
						h_y->SetLineColor(kBlack);
						h_y->DrawClone("same hist");
					}
				}
			}

			latex->DrawLatexNDC(0.2, 0.77, yTitle);
			latex->DrawLatexNDC(0.2, 0.85, Form("%i0-%i0%%", iCent, (iCent + 1)));
			legend2->Draw("same");
			canPtY[iCent]->SaveAs("check.pdf");
		}
	}

	can->SaveAs("check.pdf]");
	//	Draw flow hists //

	Int_t model_colors[] = {kOrange, kPink, kAzure + 1, kOrange + 2, kTeal + 2};

	Int_t method_colors[] = {kOrange, kPink + 2, kGreen, kAzure + 5, kGray, kTeal + 2, kCyan, kOrange + 2, kPink, kAzure + 1, kMagenta};

	can->SaveAs("checkFlow_" + nameMethod[flowMethod] + ".pdf[");
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		for (Int_t iHarmonic = 0; iHarmonic < 2; iHarmonic++)
		{
			for (Int_t iPt = 0; iPt < 4; iPt++)
			{
				auto leg = new TLegend(0.65, 0.55, 0.85, 0.86);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				can->cd();

				for (Int_t iGroup = sysGroupNameFlow.size() - 3; iGroup < sysGroupNameFlow.size(); iGroup++)
				{
					for (Int_t iVariation = 0; iVariation < sysGroupVariationIndexFlow[iGroup].size(); iVariation++)
					{
						Int_t index = sysGroupVariationIndexFlow[iGroup][iVariation];
						TH2D *hist = (TH2D *)inFlow[index]->Get(Form("cent%iharmonic%i", iCent, iHarmonic));

						TH1D *projection = hist->ProjectionX(Form("iPt%i%i%i%i%i", iPt, iGroup, iVariation, iCent, iHarmonic), iPt + 1, iPt + 1);

						TString centTitle = Form("%i0-%i0%% ", iCent, (iCent + 1));
						TString ptTitle = Form("%.0f<p_{t}<%.0f MeV/c", (iPt * ptBinSize + ptMin), ((iPt + 1) * ptBinSize + ptMin));

						latex->DrawLatexNDC(0.25, 0.24, ptTitle);
						latex->DrawLatexNDC(0.25, 0.3, centTitle);
						projection->SetTitle("");
						projection->GetYaxis()->SetTitle(Form("v%i", iHarmonic + 1));
						projection->GetYaxis()->SetTitleOffset(1.2);
						projection->SetMarkerColor(model_colors[iGroup] + iVariation);
						projection->SetLineColor(model_colors[iGroup] + iVariation);
						projection->SetMarkerStyle(20 + iGroup - sysGroupNameFlow.size() + 3);
						projection->SetMarkerSize(2);

						//	if (index == (nParametersSets + flowMethod)) projection->SetLineColor(kBlack);
						if (index != indexRef)
							leg->AddEntry(projection, sysGroupNameFlow[iGroup], "lp");
						projection->SetLineWidth(2);
						projection->GetYaxis()->SetRangeUser(-0.25, 0.25);
						if ((iGroup + iVariation) == sysGroupNameFlow.size() - 3)
							projection->DrawClone("PE");
						else
							projection->DrawClone("same PE");
						//	if (index == (nParametersSets + flowMethod)) projection->DrawClone("same PE");
					}
				}
				leg->Draw("same");

				Double_t logoW = DrawHADESLogoText(0.16, 0.84, 40., 12, -1., 3.);
				latex->DrawLatexNDC(0.16 + logoW, 0.84, "Work in progress");
				latex->DrawLatexNDC(0.15, 0.77, "Ag+Ag #sqrt{s_{NN}} = 2.55 GeV");
				// latex->DrawLatexNDC(0.25, 0.3, Form("%i0 - %i0 %%", iCent, (iCent + 1)));

				// leg->Clear();

				can->Modified();
				can->Update();
				can->SaveAs("checkFlow_" + nameMethod[flowMethod] + ".pdf");
			}
		}
	}

	can->SaveAs("checkFlow_" + nameMethod[flowMethod] + ".pdf]");

	// ===============================================================//
	TFile *out = new TFile(outputfile, "RECREATE");
	const Int_t nHarmonics = 2;
	TH2D *h_pt_y_sys[centBins]; // number of pi0 vs rapidity and pt
	TH1D *h_y_sys[centBins];	// number of pi0 vs rapidity
	TH1D *h_num_sys[centBins];	// number of pi0 integrated

	TH2D *h_flow_sys[centBins][nHarmonics];				 // number of pi0 integrated
	TH1D *h_flow_ptIntegrated_sys[centBins][nHarmonics]; // number of pi0 integrated

	out->cd();
	// =======================================================================//
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		TString centName = Form("cent%i0_%i0", iCent, (iCent + 1));
		std::cout << '\n'
				  << "centrality" << centName << '\n'
				  << '\n';

		TDirectory *centDir = out->mkdir(centName);
		centDir->cd();
		h_pt_y_sys[iCent] = compute<TH2D>(in, sysGroupVariationIndex, sysGroupName, Form("h_pt_y_cent%i", iCent), indexRef); // not yet in files
		h_pt_y_sys[iCent]->Write();
		TH2D *h_pt_y = (TH2D *)in[indexRef]->Get(Form("h_pt_y_cent%i", iCent));
		h_pt_y->Write();

		h_y_sys[iCent] = compute<TH1D>(in, sysGroupVariationIndex, sysGroupName, Form("h_y_cent%i", iCent), indexRef);
		h_y_sys[iCent]->Write();
		TH1D *h_y = (TH1D *)in[indexRef]->Get(Form("h_y_cent%i", iCent));
		h_y->Write();

		h_num_sys[iCent] = compute<TH1D>(in, sysGroupVariationIndex, sysGroupName, Form("h_num_cent%i", iCent), indexRef);
		h_num_sys[iCent]->Write();
		TH1D *h_num = (TH1D *)in[indexRef]->Get(Form("h_num_cent%i", iCent));
		h_num->Write();

		for (Int_t iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++)
		{
			h_flow_sys[iCent][iHarmonic] = compute<TH2D>(inFlow, sysGroupVariationIndexFlow, sysGroupNameFlow, Form("cent%iharmonic%i", iCent, iHarmonic), indexRef);
			Int_t bin1 = h_flow_sys[iCent][iHarmonic]->FindBin(0.8, 600);
			h_flow_sys[iCent][iHarmonic]->SetBinContent(bin1, 10000);
			h_flow_sys[iCent][iHarmonic]->SetBinError(bin1, 0);
			Int_t bin2 = h_flow_sys[iCent][iHarmonic]->FindBin(0.6, 800);
			h_flow_sys[iCent][iHarmonic]->SetBinContent(bin2, 10000);
			h_flow_sys[iCent][iHarmonic]->SetBinError(bin2, 0);
			Int_t bin3 = h_flow_sys[iCent][iHarmonic]->FindBin(0.8, 800);
			h_flow_sys[iCent][iHarmonic]->SetBinContent(bin3, 10000);
			h_flow_sys[iCent][iHarmonic]->SetBinError(bin3, 0);

			h_flow_sys[iCent][iHarmonic]->Write();

			TH2D *h_flow = (TH2D *)inFlow[indexRef]->Get(Form("cent%iharmonic%i", iCent, iHarmonic));

			h_flow->SetBinContent(bin1, 10000);
			h_flow->SetBinError(bin1, 0);

			h_flow->SetBinContent(bin2, 10000);
			h_flow->SetBinError(bin2, 0);

			h_flow->SetBinContent(bin3, 10000);
			h_flow->SetBinError(bin3, 0);

			h_flow->Write();

			h_flow_ptIntegrated_sys[iCent][iHarmonic] = compute<TH1D>(inFlow, sysGroupVariationIndexFlow, sysGroupNameFlow, Form("cent%iharmonic%iIntegratedPt", iCent, iHarmonic), indexRef);
			h_flow_ptIntegrated_sys[iCent][iHarmonic]->Write();
			TH1D *h_flow_ptIntegrated = (TH1D *)inFlow[indexRef]->Get(Form("cent%iharmonic%iIntegratedPt", iCent, iHarmonic));
			h_flow_ptIntegrated->Write();
		}
	} // iCent
	exit;
} // ana_results
