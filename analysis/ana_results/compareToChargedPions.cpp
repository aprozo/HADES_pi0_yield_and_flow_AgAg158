// #include "drawStyle.C"
#include "Includes.h"
#include "TProfile2D.h"

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

void setHist(TH1D *hist)
{
	hist->SetMarkerStyle(20);
	hist->SetMarkerSize(1.5);
	hist->SetLineColor(1);
	hist->SetMarkerColor(1);
	hist->SetMinimum(-0.25);
	hist->SetMaximum(0.25);
	hist->GetYaxis()->SetTitleOffset(0.95);
	hist->GetYaxis()->SetTitleSize(0.05);
	hist->SetFillStyle(0);
}

void drawStyle()
{
	// root -l
	// .L drawStyle.C++

	// gSystem->Load("/misc/galatyuk/tools/drawStyle_C.so");
	// drawStyle();

	cout << " -- SET DRAW OPTIONS LOADED --" << endl;

	// defining canvas options ------------------------
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);

	// defining size of pad ---------------------------
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.1);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadLeftMargin(0.2);
	gStyle->SetPadBorderMode(0);

	// switch options of the histogram ----------------
	//    gStyle->SetOptStat(0);         // switch off stat box
	//    gStyle->SetErrorX(0);          // switch off error in x
	gStyle->SetTitleStyle(0);	   // make title box transparent
	gStyle->SetTitleBorderSize(0); // make border of box disapear
}

void setCanvas(TCanvas *canvas)
{
	gStyle->SetLineStyleString(22, "80 18 12 18 12 12");
	canvas->SetFillColor(10);
	canvas->SetBorderMode(0);
	canvas->SetBorderSize(0);
	canvas->SetTickx();
	canvas->SetTicky();
	canvas->SetFrameLineWidth(2);
	canvas->SetFrameBorderMode(0);
	canvas->SetFrameBorderSize(0);
	canvas->SetLeftMargin(0.1);
	canvas->SetRightMargin(0.2);
	canvas->SetTopMargin(0.1);
	canvas->SetBottomMargin(0.1);
	canvas->ToggleEventStatus();
	//    canvas->Range(-194.483,-10.3682,1041.38,-2.08469);
}

void setTH1(TH1 *hist, TString xAxisTitle = " ", TString yAxisTitle = " ", Int_t Ndivision = 505,
			Int_t marker_style = 20, Float_t marker_size = 0.8, Int_t color = 1)
{
	hist->GetXaxis()->SetTitle(xAxisTitle);
	hist->GetYaxis()->SetTitle(yAxisTitle);
	hist->GetXaxis()->SetTitleSize(0.06);
	hist->GetYaxis()->SetTitleSize(0.06);
	hist->GetXaxis()->SetTitleFont(42);
	hist->GetYaxis()->SetTitleFont(42);
	hist->GetXaxis()->SetNdivisions(Ndivision);
	hist->GetYaxis()->SetNdivisions(Ndivision);
	hist->GetYaxis()->SetTitleOffset(1.3);
	hist->GetXaxis()->SetTitleOffset(1.1);
	//    hist->GetXaxis()->CenterTitle();
	//    hist->GetYaxis()->CenterTitle();
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);
	// hist->GetXaxis()->SetRangeUser(0., 650.);
	hist->SetMarkerStyle(marker_style);
	hist->SetMarkerSize(marker_size);
	hist->SetMarkerColor(color);
	hist->SetLineColor(color);
	hist->SetLineWidth(1);
	hist->SetFillColorAlpha(color, 0.25);
	//  hist->Sumw2();
	//    hist->SetMinimum(5e-5);
	//    hist->SetMaximum(1e2);
}

TH1D *makeMirror(TH1D *input)
{
	Int_t nBins = input->GetNbinsX();
	TH1D *mirror = new TH1D(input->GetName() + (TString) "Mirror", "",
							nBins, -(input->GetXaxis()->GetXmax()), -(input->GetXaxis()->GetXmin()));

	for (Int_t iY = 1; iY <= nBins; iY++)
	{
		Double_t val = input->GetBinContent(iY);
		Double_t error = input->GetBinError(iY);
		mirror->SetBinContent(nBins - iY + 1, val);
		mirror->SetBinError(nBins - iY + 1, error);
	}
	return mirror;
}

Double_t yFit(Double_t *xval, Double_t *par)
{
	Double_t C = par[0];
	Double_t T = par[1];
	Double_t y = xval[0];
	Double_t m = 134.98;

	return C * T * exp(-m * cosh(y) / T) * (m * m + 2 * m * T / cosh(y) + 2 * T * T / pow(cosh(y), 2));
}

Double_t ExpFitFunc_pt(Double_t *xval, Double_t *par)
{
	Double_t c = par[0];
	Double_t T = par[1];
	Double_t p_t = xval[0];
	Double_t M = 134.98;
	Double_t m_t = sqrt(p_t * p_t + M * M);

	return c * p_t * m_t * exp(-m_t / T);
}

Double_t ExpFitFunc_pt2(Double_t *xval, Double_t *par)
{
	Double_t c1 = par[0];
	Double_t T1 = par[1];
	Double_t c2 = par[3];
	Double_t T2 = par[4];

	Double_t p_t = xval[0];
	Double_t M = 134.98;
	Double_t m_t = sqrt(p_t * p_t + M * M);

	return c1 * p_t * m_t * exp(-m_t / T1) + c2 * p_t * m_t * exp(-m_t / T2);
}

Double_t temperatureFitFunc(Double_t *xval, Double_t *par)
{
	Double_t C = par[0];
	Double_t y = xval[0];

	return C / cosh(y);
}

void compareToChargedPions(TString inputfile = "systematics.root", TString outputfile = "compareToChargedPions.root")
{
	drawStyle();
	gStyle->SetErrorX(0.25);

	set_nice_color5();
	Int_t color_BK[] = {2001, 2002, 2003, 2004, 2005, 2006};
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(55);

	TFile *out = new TFile(outputfile, "RECREATE");
	TFile *in = new TFile(inputfile, "READ");

	TFile *inChargedPionsY = new TFile("chargedPionsHists/pi0_from_pip_pim.root", "READ");
	TCanvas *can = new TCanvas("can", "can", 800, 600);
	setCanvas(can);
	can->SaveAs("pi0ChargedPions.pdf[");
	TLatex *latex = new TLatex();

	latex->SetTextFont(42);
	latex->SetTextAlign(12);

	TLegend *leg = new TLegend(0.25, 0.12, 0.65, 0.33, NULL, "brNDC");
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetNColumns(3);
	leg->SetHeader("       ");

	gStyle->SetOptStat(0);
	for (Int_t iCent = 0; iCent < 3; iCent++)
	{
		TString centName = Form("cent%i0_%i0/", iCent, (iCent + 1));
		std::cout << '\n'
				  << "centrality" << centName << '\n'
				  << '\n';

		TH1D *sys = (TH1D *)in->Get(centName + Form("h_y_cent%i_sys", iCent));
		TH1D *results = (TH1D *)in->Get(centName + Form("h_y_cent%i", iCent));
		delete sys->GetListOfFunctions()->FindObject(Form("y_fit_cent%i0", iCent));
		delete results->GetListOfFunctions()->FindObject(Form("y_fit_cent%i0", iCent));
		sys->SetBinContent(sys->FindBin(0.15), 0);
		results->SetBinContent(sys->FindBin(0.15), 0);
		sys->SetBinError(sys->FindBin(0.15), 0);
		results->SetBinError(sys->FindBin(0.15), 0);

		TH1D *resultsMirror = makeMirror(results);
		TH1D *sysMirror = makeMirror(sys);

		TH1D *chargedPions = (TH1D *)inChargedPionsY->Get(Form("dNdy_%i0_%i0", iCent, iCent + 1));

		// setTH1(TH1 *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndivision=505, Int_t marker_style =20, Float_t marker_size =0.8, Int_t color=1)

		setTH1(chargedPions, "#it{y}_{c.m.}", "d#it{N}/d#it{y}_{c.m.} [1/event]", 505, 20 + iCent, 1.2, kGray + 1);
		setTH1(results, "#it{y}_{c.m.}", "d#it{N}/d#it{y}_{c.m.} [1/event]", 505, 20 + iCent, 1.2, 2001 + iCent);
		setTH1(sys, " ", " ", 505, 20 + iCent, 1.2, 2001 + iCent);
		setTH1(resultsMirror, " ", " ", 505, 24 + iCent, 1.2, 2001 + iCent);
		setTH1(sysMirror, " ", " ", 505, 24 + iCent, 1.2, 2001 + iCent);
		//	results->GetXaxis()->SetCanExtend(kTRUE);
		// results->ExtendAxis(-1.5,results->GetXaxis());

		results->GetYaxis()->SetRangeUser(0, 7.5);
		results->GetXaxis()->SetRangeUser(-1.5, 1.5);
		results->GetXaxis()->SetTitleOffset(0.6);
		results->GetYaxis()->SetTitleOffset(0.6);
		TF1 *yieldYFit = new TF1(Form("y_fit_cent%d", iCent), yFit, -1.5, 1.5, 2);
		yieldYFit->SetLineWidth(1);
		yieldYFit->SetLineStyle(9);
		yieldYFit->SetLineColor(2001 + iCent);
		yieldYFit->SetParameter(0, results->GetMaximum());
		yieldYFit->SetParameter(1, 80);
		yieldYFit->SetRange(0.25, 0.75);
		results->Fit(yieldYFit, "N");
		yieldYFit->SetRange(-1.5, 1.5);

		// leg->AddEntry(results, Form(" %i0-%i0 %%", iCent, iCent + 1), "p");
		leg->AddEntry((TObject *)0, Form(" %i0-%i0 %%", iCent, iCent + 1), "");
		leg->AddEntry(results, "    ", "p");
		leg->AddEntry(chargedPions, "    ", "p");

		yieldYFit->GetYaxis()->SetRangeUser(0., 7.5);
		yieldYFit->SetTitle(";#it{y}_{c.m.};d#it{N}/d#it{y}_{c.m.} [1/event]");
		yieldYFit->GetXaxis()->SetTitleOffset(0.6);
		yieldYFit->GetYaxis()->SetTitleOffset(0.6);
		if (iCent == 0)
			yieldYFit->Draw();
		else
			yieldYFit->Draw("same");
		chargedPions->DrawClone("same E P");
		results->DrawClone("same P E");
		yieldYFit->DrawClone("same");
		sys->DrawClone("same E2");
		sysMirror->DrawClone("same E2");

		resultsMirror->DrawClone("same P E");
	} // iCent
	leg->Draw("same");
	latex->SetTextSize(0.045);
	latex->DrawLatexNDC(0.43, 0.85, "Ag+Ag #sqrt{s_{NN}} = 2.55 GeV");
	latex->SetTextSize(0.035);
	latex->DrawLatexNDC(0.42, 0.28, "#pi^{0} #rightarrow #gamma #gamma ");
	latex->DrawLatexNDC(0.52, 0.28, " (#pi^{+} + #pi^{-})/2");

	can->SaveAs("pi0ChargedPions.pdf");

	gStyle->SetErrorX(0.25);
	for (Int_t iCent = 0; iCent < 3; iCent++)
	{
		TString centName = Form("cent%i0_%i0/", iCent, (iCent + 1));
		std::cout << '\n'
				  << "centrality" << centName << '\n'
				  << '\n';

		TCanvas *yieldPt = new TCanvas(Form("yieldPtCent%i", iCent), "yieldPt", 800, 600);
		setCanvas(yieldPt);
		yieldPt->cd();

		auto legend = new TLegend(0.62, 0.52, 0.82, 0.77);
		legend->SetBorderSize(0);
		legend->SetFillColor(0);
		legend->SetFillStyle(0);
		legend->SetHeader("       ");
		TH2D *h_pt_y = (TH2D *)in->Get(centName + Form("h_pt_y_cent%i_sys", iCent));
		Float_t yMin = h_pt_y->GetXaxis()->GetXmin();
		Float_t yBinSize = h_pt_y->GetXaxis()->GetBinWidth(1);

		for (Int_t iY = 1; iY <= h_pt_y->GetNbinsX(); iY++)
		{
			TString yTitle = Form("%.2f < y_{cm} < %.2f", iY * yBinSize + yMin, (iY + 1) * yBinSize + yMin);

			TH1D *neutralPi0 = (TH1D *)h_pt_y->ProjectionY(Form("py%i", iY), iY, iY);

			// if (iCent == 2 && iY == 1) neutralPi0->SetBinError(neutralPi0->FindBin(125), 0);
			// if (iCent == 0 && iY == 2) neutralPi0->SetBinError(neutralPi0->FindBin(125), 0);

			neutralPi0->GetYaxis()->SetTitle("dN/dp_{t}[1/event]");
			neutralPi0->GetYaxis()->SetTitleOffset(1.25);

			// TF1 bFit1("bFit1", ExpFitFunc_pt, 50, 850, 2);
			// bFit1.SetParNames("C1(y)", "T_B1(y)");
			// TF1 bFit2("bFit2", ExpFitFunc_pt, 50, 850, 2);
			// bFit1.SetParNames("C1(y)", "T_B1(y)");
			// TF1 boltzmannFit("boltzmannFit", ExpFitFunc_pt2, 100, 800, 4);
			// boltzmannFit.SetParNames("C1(y)", "T_B1(y)", "C2(y)", "T_B2(y)");
			//
			// bFit1.SetParameter(0, neutralPi0->GetMaximum());
			// bFit1.SetParameter(1, 60);
			// bFit1.SetParError(1, 60 * 0.1);
			// bFit2.SetParameter(0, neutralPi0->GetMaximum() * 0.5);
			// bFit2.SetParameter(1, 90);
			// bFit2.SetParError(1, 90 * 0.1);
			//
			// neutralPi0->Fit("bFit1", "QRNE", "", 150, 300);
			// neutralPi0->Fit("bFit2", "QRNE", "", 450, 850);
			//
			// Float_t C1 = bFit1.GetParameter(0);
			// Float_t C2 = bFit2.GetParameter(0);
			// Float_t T1 = bFit1.GetParameter(1);
			// Float_t T2 = bFit2.GetParameter(1);
			//
			// boltzmannFit.SetParError(0, C1 * 0.1);
			// boltzmannFit.SetParError(1, T1 * 0.05);
			// boltzmannFit.SetParError(2, C2 * 0.1);
			// boltzmannFit.SetParError(3, T2 * 0.05);
			//
			// boltzmannFit.SetParameters(bFit1.GetParameter(0), bFit1.GetParameter(1), bFit2.GetParameter(0), bFit2.GetParameter(
			//     1));
			// //  neutralPi0->Fit("boltzmannFit", "QRN", "", 200, 600);
			//
			// // Float_t t2_range = 5;
			// boltzmannFit.SetParLimits(1, T1 - 10, T1 + 10);
			// boltzmannFit.SetParLimits(3, T2 - 5, T2 + 5);
			//
			//
			// neutralPi0->Fit("boltzmannFit", "QRNEM", "", 150, 850);
			// //  boltzmannFit.ReleaseParameter(0);
			// // boltzmannFit.ReleaseParameter(1);
			// //  boltzmannFit.ReleaseParameter(2);
			// // boltzmannFit.ReleaseParameter(3);
			//
			//
			// neutralPi0->Fit("boltzmannFit", "RNEM", "", 150, 850);
			//
			//
			// boltzmannFit.SetRange(0, 1200);

			TF1 boltzmannFit("boltzmannFit", ExpFitFunc_pt, 100, 600, 2);
			boltzmannFit.SetParameter(0, neutralPi0->GetMaximum());
			boltzmannFit.SetParameter(1, 80);
			boltzmannFit.SetLineStyle(2);

			neutralPi0->Fit("boltzmannFit", "QRNEM+", "", 100, 400);
			boltzmannFit.SetRange(0, 1200);

			Double_t C = boltzmannFit.GetParameter(0);
			Double_t T = boltzmannFit.GetParameter(1);

			neutralPi0->SetMarkerStyle(20);
			neutralPi0->SetMarkerColor(kRed);
			neutralPi0->GetYaxis()->SetRangeUser(0.000001, 1.5 * neutralPi0->GetMaximum());

			neutralPi0->SetTitle(centName);
			neutralPi0->GetXaxis()->SetRangeUser(0, 800);
			neutralPi0->Draw("E0,P");
			// boltzmannFit.DrawClone("same");
			TFile *inChargedPionsPtY = new TFile(Form("chargedPionsHists/pions_pt_%i0_%i0.root", iCent, iCent + 1), "READ");

			if (iY <= 7)
			{
				TH1D *chargedPip = (TH1D *)inChargedPionsPtY->Get(Form("hPtSpectra_CorrectedExpPiP_YBin%i", 8 + iY));
				TH1D *chargedPim = (TH1D *)inChargedPionsPtY->Get(Form("hPtSpectra_CorrectedExpPiM_YBin%i", 8 + iY));
				TH1D *chargedPi0 = (TH1D *)inChargedPionsPtY->Get(Form("hPtSpectra_CorrectedExpPi0_YBin%i", 8 + iY));
				chargedPi0->SetMarkerColor(kAzure + 1);
				chargedPi0->SetMarkerStyle(20);
				chargedPim->SetMarkerStyle(21);
				chargedPip->SetMarkerStyle(25);
				// chargedPi0->SetMarkerSize(2);
				// chargedPim->SetMarkerSize(2);
				// chargedPip->SetMarkerSize(2);
				//	gStyle->SetErrorX(0);
				chargedPim->SetTitle(centName);
				chargedPim->Draw("same,E0,P,X0");
				chargedPip->Draw("same,E0,P,X0");
				chargedPi0->Draw("same,E0,P,X0");

				//
				// bFit1.SetParameter(0, chargedPi0->GetMaximum());
				// bFit1.SetParameter(1, 60);
				//
				// bFit2.SetParameter(0, chargedPi0->GetMaximum());
				// bFit2.SetParameter(1, 95);
				//
				// chargedPi0->Fit("bFit1", "QRN","", 0,250);
				// chargedPi0->Fit("bFit2", "QRN","", 250,800);
				//
				// Float_t C1 = bFit1.GetParameter(0);
				// Float_t C2 = bFit2.GetParameter(0);
				// Float_t T1 = bFit1.GetParameter(1);
				// Float_t T2 = bFit2.GetParameter(1);
				//
				// boltzmannFit.SetParError(3, T1*0.05);
				// boltzmannFit.SetParError(1, T1*0.05);
				// boltzmannFit.SetParError(0, C1*0.1);
				// boltzmannFit.SetParError(2, C2*0.1);
				//
				// boltzmannFit.SetParameters(bFit1.GetParameter(0), bFit1.GetParameter(
				//                 1), bFit2.GetParameter(0), bFit2.GetParameter(1));
				//
				//
				// chargedPi0->Fit("boltzmannFit", "QRN","", 200,600);
				//
				// Float_t t2_range = 3;
				// boltzmannFit.SetParLimits(3, T2-t2_range, T2+t2_range);
				// chargedPi0->Fit("boltzmannFit", "QRNEM","", 200,600);
				// boltzmannFit.ReleaseParameter(3);
				// chargedPi0->Fit("boltzmannFit", "QRNEM","", 200,600);
				// boltzmannFit.SetLineColor(kAzure+8);
				// boltzmannFit.DrawClone("same");

				if (iY == 1)
				{
					legend->AddEntry(chargedPim, "#pi^{-}", "ep");
					legend->AddEntry(chargedPip, "#pi^{+}", "ep");
					legend->AddEntry(chargedPi0, "(#pi^{-} + #pi^{+})/2 ", "ep");
					legend->AddEntry(neutralPi0, "#pi^{0}", "ep");
				}
			}
			TLatex *latex = new TLatex();
			latex->SetTextSize(0.035);
			latex->SetTextFont(42);
			latex->SetTextAlign(12);
			latex->DrawLatexNDC(0.6, 0.77, yTitle);
			latex->DrawLatexNDC(0.6, 0.85, Form("%i0-%i0%%", iCent, (iCent + 1)));
			legend->Draw("same");
			yieldPt->SaveAs("pi0ChargedPions.pdf");
		}
	}

	////log scale//////////

	for (Int_t iCent = 0; iCent < 3; iCent++)
	{
		TString centName = Form("cent%i0_%i0/", iCent, (iCent + 1));
		TCanvas *allYields = new TCanvas(Form("allYieldsCent%i", iCent), "allYields", 800, 600);
		setCanvas(allYields);
		allYields->cd();

		auto legend = new TLegend(0.52, 0.62, 0.76, 0.87);
		legend->SetBorderSize(0);
		legend->SetFillColor(0);
		legend->SetFillStyle(0);
		legend->SetHeader("       ");

		auto legend5 = new TLegend(0.62, 0.42, 0.78, 0.65);
		legend5->SetBorderSize(0);
		legend5->SetFillColor(0);
		legend5->SetFillStyle(0);
		legend5->SetHeader("       ");

		TH2D *h_pt_y = (TH2D *)in->Get(centName + Form("h_pt_y_cent%i_sys", iCent));
		Float_t yMin = h_pt_y->GetXaxis()->GetXmin();
		Float_t yMax = h_pt_y->GetXaxis()->GetXmax();
		Float_t yBinSize = h_pt_y->GetXaxis()->GetBinWidth(1);
		Float_t yBins = h_pt_y->GetNbinsX();

		for (Int_t iY = 1; iY <= 7; iY++)
		{
			TH1D *neutralPi0 = (TH1D *)h_pt_y->ProjectionY(Form("py%i%i", iY, iCent), iY, iY);
			setHist(neutralPi0);
			neutralPi0->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy_{cm}[1/event]");
			neutralPi0->GetYaxis()->SetTitleOffset(0.95);
			neutralPi0->Scale(pow(10, iY - 1));
			if (iY == 4 && iCent == 1)
				neutralPi0->SetBinError(1, neutralPi0->GetBinError(1) * 0.8);
			TF1 *boltzmannFit = new TF1(Form("boltzmannFit%i%i", iY, iCent), ExpFitFunc_pt, 0, 1200, 2);
			boltzmannFit->SetParameter(0, neutralPi0->GetMaximum());
			boltzmannFit->SetParameter(1, 80);
			boltzmannFit->SetLineStyle(9);
			boltzmannFit->SetLineWidth(4);
			boltzmannFit->SetLineColor(kOrange - 7 + iY);

			neutralPi0->Fit(Form("boltzmannFit%i%i", iY, iCent), "QRNEM+");
			Double_t C = boltzmannFit->GetParameter(0);
			Double_t T = boltzmannFit->GetParameter(1);
			Double_t TErr = boltzmannFit->GetParError(1);

			Float_t xmin = neutralPi0->GetBinCenter(neutralPi0->FindFirstBinAbove(0) - 1);
			Float_t xmax = neutralPi0->GetBinCenter(neutralPi0->FindLastBinAbove(0) + 1);
			boltzmannFit->SetRange(xmin, xmax);
			neutralPi0->SetMarkerStyle(20);
			neutralPi0->SetMarkerSize(1);
			neutralPi0->SetLineWidth(2);
			neutralPi0->SetFillColorAlpha(kOrange - 7 + iY, 0.25);
			neutralPi0->SetLineColor(kOrange - 7 + iY);
			neutralPi0->SetMarkerColor(kOrange - 7 + iY);
			neutralPi0->GetYaxis()->SetRangeUser(10E-6, 10E5);

			gPad->SetLogy();
			neutralPi0->Draw("E2,P same");
			boltzmannFit->Draw("same");

			TLatex *latex = new TLatex();
			latex->SetTextSize(0.035);
			latex->SetTextFont(42);
			latex->SetTextAlign(12);
			latex->DrawLatexNDC(0.2, 0.2, Form("%i0-%i0%%", iCent, (iCent + 1)));

			TString yTitle = Form("%.2f < y_{cm} < %.2f (#upoint10^{%i})", iY * yBinSize + yMin, (iY + 1) * yBinSize + yMin,
								  iY - 1);
			legend->AddEntry(neutralPi0, yTitle, "ep");

			TFile *inChargedPionsPtY = new TFile(Form("chargedPionsHists/pions_pt_%i0_%i0.root", iCent, iCent + 1), "READ");

			if (iY <= 7)
			{
				TH1D *chargedPip = (TH1D *)inChargedPionsPtY->Get(Form("hPtSpectra_CorrectedExpPiP_YBin%i", 8 + iY));
				TH1D *chargedPim = (TH1D *)inChargedPionsPtY->Get(Form("hPtSpectra_CorrectedExpPiM_YBin%i", 8 + iY));
				TH1D *chargedPi0 = (TH1D *)inChargedPionsPtY->Get(Form("hPtSpectra_CorrectedExpPi0_YBin%i", 8 + iY));
				chargedPi0->SetMarkerColor(kAzure + 1);
				chargedPi0->SetMarkerStyle(20);
				chargedPim->SetMarkerStyle(21);
				chargedPip->SetMarkerStyle(25);
				chargedPim->Scale(pow(10, iY - 1));
				chargedPip->Scale(pow(10, iY - 1));
				chargedPi0->Scale(pow(10, iY - 1));

				chargedPim->SetTitle(centName);
				chargedPim->Draw("same,E2,P,X0");
				chargedPip->Draw("same,E2,P,X0");
				chargedPi0->Draw("same,E2,P,X0");

				if (iY == 1)
				{
					legend5->AddEntry(chargedPim, "#pi^{-}", "ep");
					legend5->AddEntry(chargedPip, "#pi^{+}", "ep");
					legend5->AddEntry(chargedPi0, "(#pi^{-} + #pi^{+})/2 ", "ep");
					legend5->AddEntry(neutralPi0, "#pi^{0}", "ep");
				}
			}
			neutralPi0->Fit(Form("boltzmannFit%i%i", iY, iCent), "QRNEM+");
			boltzmannFit->SetLineWidth(4);

			neutralPi0->Draw("E2,P same");
			boltzmannFit->Draw("same");
			legend5->Draw("same");
		}
		Double_t logoW = DrawHADESLogoText(0.21, 0.84, 40., 12, -1., 3.);
		latex->DrawLatexNDC(0.21 + logoW, 0.84, "Work in progress");
		latex->DrawLatexNDC(0.25, 0.77, "Ag+Ag #sqrt{s_{NN}} = 2.55 GeV");
		// latex->DrawLatexNDC(0.2, 0.7, centName);
		// latex->DrawLatexNDC(0.25, 0.3, Form("%i0 - %i0 %%", iCent, (iCent + 1)));

		legend->Draw("same");
		allYields->SaveAs("pi0ChargedPions.pdf");
	}

	can->SaveAs("pi0ChargedPions.pdf]");

	TH1D *temperatureExtraction[3];

	for (Int_t iCent = 0; iCent < 3; iCent++)
	{
		TString centName = Form("cent%i0_%i0/", iCent, (iCent + 1));
		TCanvas *allYields = new TCanvas(Form("allYieldsCent%i", iCent), "allYields", 800, 600);
		setCanvas(allYields);
		allYields->cd();

		auto legend = new TLegend(0.57, 0.62, 0.82, 0.87);
		legend->SetBorderSize(0);
		legend->SetFillColor(0);
		legend->SetFillStyle(0);
		legend->SetHeader("       ");
		TH2D *h_pt_y = (TH2D *)in->Get(centName + Form("h_pt_y_cent%i_sys", iCent));
		Float_t yMin = h_pt_y->GetXaxis()->GetXmin();
		Float_t yMax = h_pt_y->GetXaxis()->GetXmax();
		Float_t yBinSize = h_pt_y->GetXaxis()->GetBinWidth(1);
		Float_t yBins = h_pt_y->GetNbinsX();

		temperatureExtraction[iCent] =
			new TH1D(Form("temperatureExtraction%i", iCent), ";y_{cm};T[MeV]", yBins, yMin, yMax);

		for (Int_t iY = 1; iY <= yBins; iY++)
		{
			TH1D *neutralPi0 = (TH1D *)h_pt_y->ProjectionY(Form("py%i%i", iY, iCent), iY, iY);
			setHist(neutralPi0);
			neutralPi0->GetYaxis()->SetTitle("d^{2}N/dp_{t}dy_{cm}[1/event]");
			neutralPi0->GetYaxis()->SetTitleOffset(0.9);
			neutralPi0->Scale(pow(10, iY - 1));
			if (iY == 4 && iCent == 1)
				neutralPi0->SetBinError(1, neutralPi0->GetBinError(1) * 0.8);
			if (iY == yBins && iCent == 2)
				neutralPi0->SetBinError(8, neutralPi0->GetBinError(8) * 0.5);
			TF1 *boltzmannFit = new TF1(Form("boltzmannFit%i%i", iY, iCent), ExpFitFunc_pt, 0, 1200, 2);
			boltzmannFit->SetParameter(0, neutralPi0->GetMaximum());
			boltzmannFit->SetParameter(1, 80);
			boltzmannFit->SetLineStyle(2);
			boltzmannFit->SetLineColor(kOrange - 7 + iY);

			neutralPi0->Fit(Form("boltzmannFit%i%i", iY, iCent), "QRNEM+");
			Double_t C = boltzmannFit->GetParameter(0);
			Double_t T = boltzmannFit->GetParameter(1);
			Double_t TErr = boltzmannFit->GetParError(1);
			if (iY < 8)
			{
				temperatureExtraction[iCent]->SetBinContent(iY, T);
				temperatureExtraction[iCent]->SetBinError(iY, TErr);
			}

			Float_t xmin = neutralPi0->GetBinCenter(neutralPi0->FindFirstBinAbove(0));
			Float_t xmax = neutralPi0->GetBinCenter(neutralPi0->FindLastBinAbove(0));
			boltzmannFit->SetRange(xmin, xmax);
			neutralPi0->SetMarkerStyle(20);
			neutralPi0->SetLineWidth(3);
			neutralPi0->SetFillColorAlpha(kOrange - 7 + iY, 0.25);
			neutralPi0->SetLineColor(kOrange - 7 + iY);
			neutralPi0->SetMarkerColor(kOrange - 7 + iY);
			neutralPi0->GetYaxis()->SetRangeUser(10E-6, 10E8);
			gPad->SetLogy();
			neutralPi0->Draw("E2,P same");
			boltzmannFit->Draw("same");

			TLatex *latex = new TLatex();
			latex->SetTextSize(0.035);
			latex->SetTextFont(42);
			latex->SetTextAlign(12);
			latex->DrawLatexNDC(0.2, 0.2, Form("%i0-%i0%%", iCent, (iCent + 1)));

			TString yTitle = Form("%.2f < y_{cm} < %.2f (#upoint10^{%i})", iY * yBinSize + yMin, (iY + 1) * yBinSize + yMin,
								  iY - 1);
			legend->AddEntry(neutralPi0, yTitle, "ep");
		}
		legend->Draw("same");
		Double_t logoW = DrawHADESLogoText(0.21, 0.84, 40., 12, -1., 3.);
		latex->DrawLatexNDC(0.21 + logoW, 0.84, "Work in progress");
		latex->DrawLatexNDC(0.2, 0.77, "Ag+Ag #sqrt{s_{NN}} = 2.55 GeV");
		latex->DrawLatexNDC(0.2, 0.7, centName);
		latex->DrawLatexNDC(0.25, 0.3, Form("%i0 - %i0 %%", iCent, (iCent + 1)));

		allYields->SaveAs("boltzmannFits.pdf");
	}

	{
		TCanvas *allYieldsTemperature = new TCanvas("allYieldsTemperature", "allYieldsTemperature", 800, 600);
		setCanvas(allYieldsTemperature);
		allYieldsTemperature->cd();

		TLatex *latex = new TLatex();
		latex->SetTextSize(0.035);
		latex->SetTextFont(42);
		latex->SetTextAlign(12);
		auto legend = new TLegend(0.64, 0.62, 0.8, 0.87);
		legend->SetBorderSize(0);
		legend->SetFillColor(0);
		legend->SetFillStyle(0);
		legend->SetHeader("       ");

		for (Int_t iCent = 0; iCent < 3; iCent++)
		{
			TH1D *temperatureExtractionMirror = makeMirror(temperatureExtraction[iCent]);
			TF1 *temperatureFit = new TF1(Form("temperatureFit%i", iCent), "[0]/cosh(x)", -1., 1.);
			setTH1(temperatureExtraction[iCent], "#it{y}_{c.m.}", "Slope [MeV]", 505, 20, 1.2, 2001 + iCent);
			setTH1(temperatureExtractionMirror, "#it{y}_{c.m.}", "Slope [MeV]", 505, 24, 1.2, 2001 + iCent);
			temperatureExtraction[iCent]->GetXaxis()->SetRangeUser(-1.5, 1.5);
			temperatureFit->SetRange(0.15, 0.65);
			temperatureExtraction[iCent]->Fit(Form("temperatureFit%i", iCent), "QRNEM+");
			temperatureFit->SetRange(-1., 1.);
			temperatureFit->SetLineColor(2001 + iCent);
			temperatureFit->SetTitle(";y_{cm};Slope [MeV]");
			temperatureFit->DrawClone(iCent == 0 ? "" : "same");
			temperatureExtraction[iCent]->DrawClone("same");
			temperatureExtractionMirror->DrawClone("same");
			legend->AddEntry(temperatureExtraction[iCent], Form("%i0-%i0%%", iCent, iCent + 1), "ep");

			latex->DrawLatexNDC(0.3, 0.4 - iCent * 0.1,
								Form("T_{eff %i0-%i0%%} = %.1f #pm %.1f MeV", iCent, (iCent + 1), temperatureFit->GetParameter(0),
									 temperatureFit->GetParError(0)));
		}
		legend->Draw("same");
		allYieldsTemperature->SaveAs("temperature.pdf");
	}

	// ======================================================FLOW plots===============================================================//

	Int_t flowMethod = 0;
	Int_t nHarmonics = 2;
	Int_t centBins = 3;

	Int_t ptBins = 4;
	Float_t ptMin = 100;
	Float_t ptMax = 900;
	Float_t ptBinSize = (ptMax - ptMin) / Float_t(ptBins);

	Int_t yBins = 4;
	Float_t yMin = 0.1;
	Float_t yMax = 0.9;
	Float_t yBinSize = (yMax - yMin) / Float_t(yBins);

	TFile *flowInput = new TFile("systematics.root", "read");

	// READING comparison
	const int nCompares = 3;

	TString nameChargedPions[nCompares] = {"#pi^{-}", "#pi^{+}", "(#pi^{-}+#pi^{+})/2"};
	TFile *fileChargedPions = new TFile("chargedPionsHists/ChargedPionFlow_Sasha_Binning.root", "read");

	// READING MODELS
	const int nModels = 3;
	TString filenameModels[nModels] = {"1580_urqmd_eos", "1580_urqmd_cascade", "1580_gibuu"};
	TString titleModel[nModels] = {"UrQMD EoS", "UrQMD cascade", "GiBUU"};

	TFile *fin[nModels];
	for (int j = 0; j < nModels; j++)
		fin[j] = new TFile("models/agag_" + filenameModels[j] + ".root", "read");

	can->cd();
	can->SaveAs("flowChargedPionComparison.pdf[");
	TLegend *leg2 = new TLegend(0.67, 0.65, 0.92, 0.88);
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
	Int_t model_colors[] = {kAzure, kPink, kGreen - 3};

	latex->SetTextSize(0.045);
	latex->SetTextFont(42);
	latex->SetTextAlign(12);

	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		for (int iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++)
		{
			TString centName = Form("centrality %i0 - %i0 %%", iCent, (iCent + 1));
			TH2D *hFlow = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i_sys", iCent, iCent + 1, iCent,
													  iHarmonic));
			TH1D *hFlowIntegratedPt =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt_sys", iCent, iCent + 1, iCent,
											iHarmonic));

			TH2D *hFlowStat = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i", iCent, iCent + 1, iCent,
														  iHarmonic));
			TH1D *hFlowIntegratedPtStat =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt", iCent, iCent + 1, iCent,
											iHarmonic));
			TProfile2D *chargedPionsFlow2D[nCompares];

			chargedPionsFlow2D[0] =
				(TProfile2D *)fileChargedPions->Get(Form("pim2D_v%i_corr_centr%i", iHarmonic + 1, iCent)); // PIM
			chargedPionsFlow2D[0]->SetTitle(Form(";y_{cm};p_{t},[MeV/c];v%i", iHarmonic + 1));
			chargedPionsFlow2D[1] =
				(TProfile2D *)fileChargedPions->Get(Form("pip2D_v%i_corr_centr%i", iHarmonic + 1, iCent)); // PIP
			chargedPionsFlow2D[1]->SetTitle(Form(";y_{cm};p_{t},[MeV/c];v%i", iHarmonic + 1));
			chargedPionsFlow2D[2] = (TProfile2D *)chargedPionsFlow2D[0]->Clone("Average"); // PIM+PIP
			chargedPionsFlow2D[2]->Add(chargedPionsFlow2D[1]);

			for (size_t yOrPt = 0; yOrPt < 2; yOrPt++)
			{
				for (int iBin = 0; iBin < 4; iBin++)
				{
					TH1D *projection;
					TH1D *projectionStat;
					if (yOrPt == 0)
					{
						projection = hFlow->ProjectionX("pX", iBin + 1, iBin + 1);
						projectionStat = hFlowStat->ProjectionX("pXstat", iBin + 1, iBin + 1);
					}
					else
					{
						projection = hFlow->ProjectionY("pY", iBin + 1, iBin + 1);
						projectionStat = hFlowStat->ProjectionY("pYstat", iBin + 1, iBin + 1);
					}

					setHist(projection);
					setHist(projectionStat);
					// projectionStat->SetErrorX(0.5);
					projection->SetLineColor(kBlack);

					leg2->Clear();
					projection->GetYaxis()->SetTitle(Form("v_{%i}", iHarmonic + 1));
					projection->Draw("E2");

					leg2->AddEntry(projection, "#pi^{0}", "lep");

					TH1D *chargedPionsFlowY[nCompares];
					for (int iCompare = 0; iCompare < nCompares - 1; iCompare++)
					{
						if (yOrPt == 0)
						{
							chargedPionsFlowY[iCompare] =
								chargedPionsFlow2D[iCompare]->ProjectionX(Form("comp%i", iCompare), iBin + 1, (iBin + 1));
						}
						else
						{
							chargedPionsFlowY[iCompare] =
								chargedPionsFlow2D[iCompare]->ProjectionY(Form("comp%i", iCompare), iBin + 1, (iBin + 1));
						}
						TH1D *hist = chargedPionsFlowY[iCompare];
						setHist(hist);
						hist->SetLineColor(2001 + iCompare);
						hist->SetMarkerColor(2001 + iCompare);
						hist->DrawClone("P E2 same");
						leg2->AddEntry(hist, nameChargedPions[iCompare], "lep");
					}
					projection->Draw("same E2");
					projectionStat->Draw("same");
					leg2->Draw();
					TString name;
					if (yOrPt == 0)
					{
						name = Form("%.0f<p_{t}<%.0f MeV/c", (iBin * ptBinSize + ptMin), ((iBin + 1) * ptBinSize + ptMin));
					}
					else
					{
						name = Form("%.1f< y_{cm} <%.1f", (iBin * yBinSize + yMin), ((iBin + 1) * yBinSize + yMin));
					}
					latex->DrawLatexNDC(0.25, 0.24, name);

					Double_t logoW = DrawHADESLogoText(0.16, 0.84, 40., 12, -1., 3.);
					latex->DrawLatexNDC(0.16 + logoW, 0.84, "Work in progress");
					latex->DrawLatexNDC(0.17, 0.77, "Ag+Ag #sqrt{s_{NN}} = 2.55 GeV");
					latex->DrawLatexNDC(0.17, 0.7, centName);

					//
					// //////////////////add models //////////////////////
					//
					//
					//
					TProfile2D *hModelPi02D = (TProfile2D *)fin[0]->Get(Form("pi0/centrclass_%i/hV%i_ip2_icc%i", iCent, iHarmonic + 1, iCent));
					TString nameProjection = Form("%i%i%i%i", iCent, iHarmonic, yOrPt, iBin);
					TH1D *hModelFlow;
					if (yOrPt == 0)
					{
						hModelFlow = hModelPi02D->ProjectionX(nameProjection + "pi0", iBin + 1, iBin + 1);
					}
					else
					{
						TH1D *hModelPi0Projection = hModelPi02D->ProjectionY(nameProjection + "pi0_oldpt", iBin + 6, iBin + 6);
						hModelFlow = new TH1D(nameProjection + "pi0", "", 4, 100, 900);
						for (size_t binX = 1; binX <= 4; binX++)
						{
							hModelFlow->SetBinContent(binX, hModelPi0Projection->GetBinContent(binX));
							hModelFlow->SetBinError(binX, hModelPi0Projection->GetBinError(binX));
						}
					}
					// TH1D * hModelFlow      = hModelPi02D->ProjectionX(nameProjection + "pi0", iBin + 1, iBin + 1);
					hModelFlow->SetLineColor(kBlack);
					hModelFlow->DrawClone("HIST L same");
					hModelFlow->SetFillColorAlpha(hModelFlow->GetLineColor(), 0.25);
					hModelFlow->DrawClone("E3 same");

					TProfile2D *hModelPiP2D =
						(TProfile2D *)fin[0]->Get(Form("pip/centrclass_%i/hV%i_ip4_icc%i", iCent, iHarmonic + 1, iCent));
					//  TH1D * hModelPiP = hModelPiP2D->ProjectionX(nameProjection + "pip", iBin + 1, iBin + 1);
					TH1D *hModelPiP;
					if (yOrPt == 0)
					{
						hModelPiP = hModelPiP2D->ProjectionX(nameProjection + "pip", iBin + 1, iBin + 1);
					}
					else
					{
						TH1D *hModelPiPProjection = hModelPiP2D->ProjectionY(nameProjection + "pip_oldpt", iBin + 6, iBin + 6);
						hModelPiP = new TH1D(nameProjection + "pip", "", 4, 100, 900);
						for (size_t binX = 1; binX <= 4; binX++)
						{
							hModelPiP->SetBinContent(binX, hModelPiPProjection->GetBinContent(binX));
							hModelPiP->SetBinError(binX, hModelPiPProjection->GetBinError(binX));
						}
					}

					hModelPiP->SetLineColor(2002);
					hModelPiP->DrawClone("HIST L same");
					hModelPiP->SetFillColorAlpha(hModelPiP->GetLineColor(), 0.25);
					hModelPiP->DrawClone("E3 same");

					TProfile2D *hModelPiM2D =
						(TProfile2D *)fin[0]->Get(Form("pim/centrclass_%i/hV%i_ip3_icc%i", iCent, iHarmonic + 1, iCent));
					// TH1D * hModelPiM = hModelPiM2D->ProjectionX(nameProjection + "pim", iBin + 1, iBin + 1);

					TH1D *hModelPiM;
					if (yOrPt == 0)
					{
						hModelPiM = hModelPiM2D->ProjectionX(nameProjection + "pim", iBin + 1, iBin + 1);
					}
					else
					{
						TH1D *hModelPiMProjection = hModelPiM2D->ProjectionY(nameProjection + "pim_oldpt", iBin + 6, iBin + 6);
						hModelPiM = new TH1D(nameProjection + "pim", "", 4, 100, 900);
						for (size_t binX = 1; binX <= 4; binX++)
						{
							hModelPiM->SetBinContent(binX, hModelPiMProjection->GetBinContent(binX));
							hModelPiM->SetBinError(binX, hModelPiMProjection->GetBinError(binX));
						}
					}
					hModelPiM->SetLineColor(2001);
					hModelPiM->DrawClone("HIST L same");
					hModelPiM->SetFillColorAlpha(hModelPiM->GetLineColor(), 0.25);
					hModelPiM->DrawClone("E3 same");

					can->SaveAs("flowChargedPionComparison.pdf");
				}
			}
		}
	}

	can->SaveAs("flowChargedPionComparison.pdf]");

	// ===================================================================================================================//
	// ===================================================compare to models============================================== //
	// ===================================================================================================================//

	can->SaveAs("modelComparison.pdf[");
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		for (int iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++)
		{
			TString centName = Form("%i0 - %i0 %%", iCent, (iCent + 1));
			std::cout << centName << "   v" << (iHarmonic + 1) << '\n';
			TH2D *hFlow = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i_sys", iCent, iCent + 1, iCent,
													  iHarmonic));
			TH1D *hFlowIntegratedPt =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt_sys", iCent, iCent + 1, iCent,
											iHarmonic));

			TH2D *hFlowStat = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i", iCent, iCent + 1, iCent,
														  iHarmonic));
			TH1D *hFlowIntegratedPtStat =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt", iCent, iCent + 1, iCent, iHarmonic));
			for (int iBin = 0; iBin < 4; iBin++)
			{
				can->cd();
				TH1D *projection;
				TH1D *projectionStat;
				projection = hFlow->ProjectionX("pX", iBin + 1, iBin + 1);
				projectionStat = hFlowStat->ProjectionX("pXstat", iBin + 1, iBin + 1);
				setHist(projection);
				setHist(projectionStat);
				projection->SetLineColor(kBlack);

				leg2->Clear();
				projection->GetYaxis()->SetTitle(Form("v_{%i}", iHarmonic + 1));
				//	projection->Draw("E2");

				//	leg2->AddEntry(projection, "data", "lep");
				//
				// projection->Draw("same E2");
				// projectionStat->Draw("same");

				for (int iModel = 0; iModel < nModels; iModel++)
				{
					// Float_t ptValMin = (iBin* ptBinSize + ptMin)*0.001;
					// Float_t ptValMax = ((iBin+1) * ptBinSize + ptMin)*0.001;
					leg2->Clear();
					TProfile2D *hModelPi02D =
						(TProfile2D *)fin[iModel]->Get(Form("pi0/centrclass_%i/hV%i_ip2_icc%i", iCent, iHarmonic + 1, iCent));
					TString nameProjection = Form("%i%i%i%i", iCent, iHarmonic, iModel, iBin);
					TH1D *hModelFlow = hModelPi02D->ProjectionX(nameProjection + "pi0", iBin + 1, iBin + 1);
					hModelFlow->SetLineColor(2003);
					hModelFlow->GetYaxis()->SetRangeUser(-0.25, 0.25);
					hModelFlow->GetYaxis()->SetTitle(Form("v_{%i}", iHarmonic + 1));
					hModelFlow->GetXaxis()->SetTitle("y_{cm}");
					hModelFlow->GetYaxis()->SetTitleOffset(1.2);
					// if (iModel == 0)
					hModelFlow->DrawClone("HIST L");
					// else hModelFlow->DrawClone("HIST L same");
					hModelFlow->SetFillColorAlpha(hModelFlow->GetLineColor(), 0.25);
					hModelFlow->DrawClone("E3 same");
					leg2->AddEntry(hModelFlow, titleModel[iModel] + " #pi^{0}", "f");

					TProfile2D *hModelPiP2D =
						(TProfile2D *)fin[iModel]->Get(Form("pip/centrclass_%i/hV%i_ip4_icc%i", iCent, iHarmonic + 1, iCent));
					TH1D *hModelPiP = hModelPiP2D->ProjectionX(nameProjection + "pip", iBin + 1, iBin + 1);
					hModelPiP->SetLineColor(2002);
					hModelPiP->DrawClone("HIST L same");
					hModelPiP->SetFillColorAlpha(hModelPiP->GetLineColor(), 0.25);
					hModelPiP->DrawClone("E3 same");
					leg2->AddEntry(hModelPiP, titleModel[iModel] + " #pi^{+}", "f");

					TProfile2D *hModelPiM2D =
						(TProfile2D *)fin[iModel]->Get(Form("pim/centrclass_%i/hV%i_ip3_icc%i", iCent, iHarmonic + 1, iCent));
					TH1D *hModelPiM = hModelPiM2D->ProjectionX(nameProjection + "pim", iBin + 1, iBin + 1);
					hModelPiM->SetLineColor(2001);
					hModelPiM->DrawClone("HIST L same");
					hModelPiM->SetFillColorAlpha(hModelPiM->GetLineColor(), 0.25);
					hModelPiM->DrawClone("E3 same");
					leg2->AddEntry(hModelPiM, titleModel[iModel] + " #pi^{-}", "f");

					leg2->Draw();
					TString name;
					name = Form("%.0f<p_{t}<%.0f MeV/c", (iBin * ptBinSize + ptMin), ((iBin + 1) * ptBinSize + ptMin));

					latex->DrawLatexNDC(0.25, 0.24, name);
					latex->DrawLatexNDC(0.25, 0.3, centName);

					can->SaveAs("modelComparison.pdf");
				}
			} // iBin
		}	  // iHarmonic
	}		  // iCent
	can->SaveAs("modelComparison.pdf]");

	can->SaveAs("modelComparisonPi0.pdf[");
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		for (int iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++)
		{
			TString centName = Form("%i0 - %i0 %%", iCent, (iCent + 1));
			std::cout << centName << "   v" << (iHarmonic + 1) << '\n';
			TH2D *hFlow = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i_sys", iCent, iCent + 1, iCent,
													  iHarmonic));
			TH1D *hFlowIntegratedPt =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt_sys", iCent, iCent + 1, iCent,
											iHarmonic));

			TH2D *hFlowStat = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i", iCent, iCent + 1, iCent,
														  iHarmonic));
			TH1D *hFlowIntegratedPtStat =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt", iCent, iCent + 1, iCent, iHarmonic));
			for (int iBin = 0; iBin < 4; iBin++)
			{
				can->cd();
				TH1D *projection;
				TH1D *projectionStat;
				projection = hFlow->ProjectionX("pX", iBin + 1, iBin + 1);
				projectionStat = hFlowStat->ProjectionX("pXstat", iBin + 1, iBin + 1);
				setHist(projection);
				setHist(projectionStat);
				projection->SetLineColor(kBlack);

				leg2->Clear();
				projection->GetYaxis()->SetTitle(Form("v_{%i}", iHarmonic + 1));
				projection->Draw("E2");
				leg2->AddEntry(projection, "data", "lep");
				projection->Draw("same E2");
				projectionStat->Draw("same");

				for (int iModel = 0; iModel < nModels; iModel++)
				{
					// Float_t ptValMin = (iBin* ptBinSize + ptMin)*0.001;
					// Float_t ptValMax = ((iBin+1) * ptBinSize + ptMin)*0.001;

					TProfile2D *hModelPi02D =
						(TProfile2D *)fin[iModel]->Get(Form("pi0/centrclass_%i/hV%i_ip2_icc%i", iCent, iHarmonic + 1, iCent));
					TString nameProjection = Form("%i%i%i%i", iCent, iHarmonic, iModel, iBin);
					TH1D *hModelFlow = hModelPi02D->ProjectionX(nameProjection + "pi0", iBin + 1, iBin + 1);
					hModelFlow->SetLineColor(2001 + iModel);
					hModelFlow->GetYaxis()->SetRangeUser(-0.25, 0.25);
					hModelFlow->GetYaxis()->SetTitle(Form("v_{%i}", iHarmonic + 1));
					hModelFlow->GetXaxis()->SetTitle("y_{cm}");
					hModelFlow->GetYaxis()->SetTitleOffset(1.2);

					hModelFlow->DrawClone("HIST L same");
					hModelFlow->SetFillColorAlpha(hModelFlow->GetLineColor(), 0.25);
					hModelFlow->DrawClone("E3 same");
					leg2->AddEntry(hModelFlow, titleModel[iModel], "f");
				}
				leg2->Draw();
				TString name;
				name = Form("%.0f<p_{t}<%.0f MeV/c", (iBin * ptBinSize + ptMin), ((iBin + 1) * ptBinSize + ptMin));

				latex->DrawLatexNDC(0.25, 0.24, name);
				latex->DrawLatexNDC(0.25, 0.3, centName);

				can->SaveAs("modelComparisonPi0.pdf");
			} // iBin
		}	  // iHarmonic
	}		  // iCent
	can->SaveAs("modelComparisonPi0.pdf]");

	leg = new TLegend(0.52, 0.65, 0.75, 0.8);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);

	can->SaveAs("flowPi0.pdf[");
	for (Int_t iCent = 0; iCent < centBins; iCent++)
	{
		for (int iHarmonic = 0; iHarmonic < nHarmonics; iHarmonic++)
		{
			TString centName = Form("%i0 - %i0 %%", iCent, (iCent + 1));
			TH2D *hFlow =
				(TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i_sys", iCent, iCent + 1, iCent, iHarmonic));
			TH1D *hFlowIntegratedPt =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt_sys", iCent, iCent + 1, iCent,
											iHarmonic));

			TH2D *hFlowStat = (TH2D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%i", iCent, iCent + 1, iCent,
														  iHarmonic));
			TH1D *hFlowIntegratedPtStat =
				(TH1D *)flowInput->Get(Form("cent%i0_%i0/cent%iharmonic%iIntegratedPt", iCent, iCent + 1, iCent, iHarmonic));

			for (size_t yOrPt = 0; yOrPt < 2; yOrPt++)
			{
				leg->Clear();
				for (int iBin = 0; iBin < 3; iBin++)
				{
					can->cd();
					TH1D *projection;
					TH1D *projectionStat;
					if (yOrPt == 0)
					{
						projection = hFlow->ProjectionX(Form("pX%i%i", iBin, yOrPt), iBin + 1, iBin + 1);
						projectionStat = hFlowStat->ProjectionX(Form("pX%i%istat", iBin, yOrPt), iBin + 1, iBin + 1);
					}
					else
					{
						projection = hFlow->ProjectionY(Form("pY%i%i", iBin, yOrPt), iBin + 1, iBin + 1);
						projectionStat = hFlowStat->ProjectionY(Form("pY%i%istat", iBin, yOrPt), iBin + 1, iBin + 1);
					}

					setHist(projection);
					setHist(projectionStat);
					projection->SetLineColor(2001 + iBin);
					projection->SetMarkerColor(2001 + iBin);

					projectionStat->SetLineColor(2001 + iBin);
					projectionStat->SetMarkerColor(2001 + iBin);

					projection->GetYaxis()->SetTitle(Form("v_{%i}", iHarmonic + 1));
					if (iBin == 0)
						projection->Draw("E2");
					else
						projection->Draw("E2 same");
					projectionStat->Draw("same");

					// latex->DrawLatexNDC(0.25, 0.3, Form("%i0 - %i0 %%", iCent, (iCent + 1)));

					TString name;
					if (yOrPt == 0)
					{
						name = Form("%.0f<p_{t}<%.0f MeV/c", (iBin * ptBinSize + ptMin), ((iBin + 1) * ptBinSize + ptMin));
					}
					else
					{
						name = Form("%.1f< y_{cm} <%.1f", (iBin * yBinSize + yMin), ((iBin + 1) * yBinSize + yMin));
					}

					leg->AddEntry(projection, name, "lep");
					// latex->DrawLatexNDC(0.25, 0.3, centName);
				} // iBin
				leg->Draw("same");
				Double_t logoW = DrawHADESLogoText(0.16, 0.84, 40., 12, -1., 3.);
				latex->DrawLatexNDC(0.16 + logoW, 0.84, "Work in progress");
				latex->DrawLatexNDC(0.17, 0.77, "Ag+Ag #sqrt{s_{NN}} = 2.55 GeV");
				latex->DrawLatexNDC(0.17, 0.7, centName);
				can->SaveAs("flowPi0.pdf");
			} // yOrPt
		}	  // iHarmonic
	}		  // iCent
	can->SaveAs("flowPi0.pdf]");
} // ana_results
