scaleFactorNorm0->SetTitle(";a.u. phase space bin p_{t}-y; normalization factor");
gStyle->SetOptStat(0)
TCanvas * can = new TCanvas("can", "can", 800, 600)


  scaleFactorNorm2->SetMarkerStyle(20)
  scaleFactorNorm2->SetMarkerColor(kOrange + 1)
  scaleFactorNorm2->SetMarkerSize(2)

  scaleFactorNorm1->SetMarkerColor(kBlue + 1)
  scaleFactorNorm1->SetMarkerStyle(20)
  scaleFactorNorm1->SetMarkerSize(2)

  scaleFactorNorm0->SetMarkerStyle(20)
  scaleFactorNorm0->SetMarkerColor(kGreen + 1)
  scaleFactorNorm0->SetMarkerSize(2)

  TF1 * fit0 = scaleFactorNorm0->GetFunction("fitConst");
fit0->SetLineColor(scaleFactorNorm0->GetMarkerColor());
fit0->SetLineStyle(9);

TF1 * fit1 = scaleFactorNorm1->GetFunction("fitConst");
fit1->SetLineColor(scaleFactorNorm1->GetMarkerColor());
fit1->SetLineStyle(9);

TF1 * fit2 = scaleFactorNorm2->GetFunction("fitConst");
fit2->SetLineColor(scaleFactorNorm2->GetMarkerColor());
fit2->SetLineStyle(9);


TLegend * legend = new TLegend(0.57, 0.15, 0.77, 0.36);
legend->SetBorderSize(0);
legend->AddEntry(scaleFactorNorm0, "00-10%", "lep");
legend->AddEntry(scaleFactorNorm1, "10-20%", "lep");
legend->AddEntry(scaleFactorNorm2, "20-30%", "lep");
scaleFactorNorm0->GetXaxis()->SetRangeUser(70, 140);
scaleFactorNorm0->GetYaxis()->SetRangeUser(0.007, 0.014);
scaleFactorNorm0->Draw("P hist")
scaleFactorNorm1->Draw("same P hist")
scaleFactorNorm2->Draw("same P hist")

fit0->Draw("same");
fit1->Draw("same");
fit2->Draw("same");

legend->Draw("same")
can->SaveAs("NormalizationFactorCheck.pdf")
