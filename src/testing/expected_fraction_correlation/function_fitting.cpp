{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();

	TF1 *ff1 = new TF1("ff1", "1/(1 + TMath::Exp(x - TMath::Log(56))) - 0.5", -0.1, 4.5);
	ff1->GetXaxis()->SetRangeUser(-0.1,4.5);
	ff1->GetYaxis()->SetRangeUser(-0.1,1.5);
	ff1->GetXaxis()->SetLimits(-0.1,4.5);
	ff1->GetYaxis()->SetLimits(-0.1,1.5);
	ff1->Draw();

	TF1 *ff2 = new TF1("ff2", "2*(1/(1 + TMath::Exp(x - TMath::Log(56))) - 0.5)", -0.1, 4.5);
	ff2->Draw("same");

	TF1 *ff3 = new TF1("ff3", "3*(1/(1 + TMath::Exp(x - TMath::Log(56))) - 0.5)", -0.1, 4.5);
	ff3->Draw("same");

	TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
	c2->SetGrid();

	TF1 *ff4 = new TF1("ff4", "(1/(1 + TMath::Exp(x - TMath::Log(56)))) - 0.5", -0.1, 4.5);
	ff4->GetXaxis()->SetRangeUser(-0.1,4.5);
	ff4->GetYaxis()->SetRangeUser(-0.1,1.5);
	ff4->GetXaxis()->SetLimits(-0.1,4.5);
	ff4->GetYaxis()->SetLimits(-0.1,1.5);
	ff4->Draw();

	TF1 *ff5 = new TF1("ff5", "(TMath::Erf(-(x - TMath::Log(56))/2.))", -0.1, 4.5);
	ff5->SetLineColor(4);
	ff5->Draw("same");

/*	TF1 *ff5b = new TF1("ff5b", "TMath::Log(TMath::Erf(-(x - TMath::Log(56))/2.))", -0.1, 4.5);
	ff5b->SetLineColor(4);
	ff5b->SetLineStyle(9);
	ff5b->Draw("same");*/

	TF1 *ff6 = new TF1("ff6", "(TMath::TanH(-(x - TMath::Log(56))/2.))", -0.1, 4.5);
	ff6->SetLineColor(2);
	ff6->Draw("same");

	TF1 *ff7 = new TF1("ff7", "-(x - TMath::Log(56))/(TMath::Sqrt(1 + (x - TMath::Log(56))*(x - TMath::Log(56))/4.))/2.", -0.1, 4.5);
	ff7->SetLineColor(3);
	ff7->Draw("same");

	TF1 *ff8 = new TF1("ff8", "TMath::ATan(-(x - TMath::Log(56)))/(TMath::Pi()/2.)", -0.1, 4.5);
	ff8->SetLineColor(6);
	ff8->Draw("same");

/*	TF1 *ff5 = new TF1("ff5", "1*(1/(1 + TMath::Exp((x - TMath::Log(56))/2)) - 0.5)", -0.1, 4.5);
	ff5->Draw("same");

	TF1 *ff6 = new TF1("ff6", "1*(1/(1 + TMath::Exp((x - TMath::Log(56))/3)) - 0.5)", -0.1, 4.5);
	ff6->Draw("same");*/

	TCanvas *c3 = new TCanvas("c3", "", 1200, 800);
	c3->SetGrid();

	double x1[] = { 1.*TMath::Log(1) + 0.*TMath::Log(56), 1.*TMath::Log(4) + 0.*TMath::Log(56), 1.*TMath::Log(16) + 0.*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)}; 
	double x2[] = { 0.75*TMath::Log(1) + 0.25*TMath::Log(56), 0.75*TMath::Log(4) + 0.25*TMath::Log(56), 0.75*TMath::Log(16) + 0.25*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double x3[] = { 0.5*TMath::Log(1) + 0.5*TMath::Log(56), 0.5*TMath::Log(4) + 0.5*TMath::Log(56), 0.5*TMath::Log(16) + 0.5*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double x4[] = { 0.25*TMath::Log(1) + 0.75*TMath::Log(56), 0.25*TMath::Log(4) + 0.75*TMath::Log(56), 0.25*TMath::Log(16) + 0.75*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};

	double y1[] = { 1.02785, 0.91522, 0.49253, -0.00639};
	double y2[] = { 0.75360, 0.67983, 0.36582, -0.00639};
	double y3[] = { 0.50140, 0.45472, 0.24379, -0.00639};
	double y4[] = { 0.25431, 0.23379, 0.12598, -0.00639};

	TGraph *gr1 = new TGraph(4, x1, y1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(2);
	gr1->SetName("pure_600");
	gr1->GetXaxis()->SetRangeUser(-0.1,4.5);
	gr1->GetYaxis()->SetRangeUser(-0.1,1.5);
	gr1->GetXaxis()->SetLimits(-0.1,4.5);
	gr1->GetYaxis()->SetLimits(-0.1,1.5);
	gr1->GetXaxis()->SetTitle("Expected proton fraction");
	gr1->GetYaxis()->SetTitle("Proton fraction");
	gr1->Draw("AP");

	TGraph *gr2 = new TGraph(4, x2, y2);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(3);
	gr2->SetName("75-25_Fe-bgd");
	gr2->Draw("SAME;P");

	TGraph *gr3 = new TGraph(4, x3, y3);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(4);
	gr3->SetName("50-50_Fe-bgd");
	gr3->Draw("SAME;P");

	TGraph *gr4 = new TGraph(4, x4, y4);
	gr4->SetMarkerStyle(20);
	gr4->SetMarkerColor(6);
	gr4->SetName("25-75_Fe-bgd");
	gr4->Draw("SAME;P");

/*	TF1 *fitfunc = new TF1("fitfunc", "[0]*(1/(1 + TMath::Exp((x - TMath::Log(56))/[1])) - 0.5)", -0.1, 4.5);
	fitfunc->SetParameter(0,2.);
	fitfunc->SetParameter(1,1.);
	gr1->Fit("fitfunc");
	gr1->GetFunction("fitfunc")->SetLineWidth(1);
	gr1->GetFunction("fitfunc")->SetLineColor(1);
	gr2->Fit("fitfunc");
	gr2->GetFunction("fitfunc")->SetLineWidth(1);
	gr2->GetFunction("fitfunc")->SetLineColor(1);
	gr3->Fit("fitfunc");
	gr3->GetFunction("fitfunc")->SetLineWidth(1);
	gr3->GetFunction("fitfunc")->SetLineColor(1);
	gr4->Fit("fitfunc");
	gr4->GetFunction("fitfunc")->SetLineWidth(1);
	gr4->GetFunction("fitfunc")->SetLineColor(1);

	c3->SaveAs("function_fit_sigmoidal.pdf");*/

	TF1 *fitfunc = new TF1("fitfunc", "[0]*(TMath::Erf((TMath::Log(56) - x)/[1]))", -0.1, 4.5);
	fitfunc->SetParameter(0,1./2.);
	fitfunc->SetParameter(1,2.);
	gr1->Fit("fitfunc");
	gr1->GetFunction("fitfunc")->SetLineWidth(1);
	gr1->GetFunction("fitfunc")->SetLineColor(1);
	gr2->Fit("fitfunc");
	gr2->GetFunction("fitfunc")->SetLineWidth(1);
	gr2->GetFunction("fitfunc")->SetLineColor(1);
	gr3->Fit("fitfunc");
	gr3->GetFunction("fitfunc")->SetLineWidth(1);
	gr3->GetFunction("fitfunc")->SetLineColor(1);
	gr4->Fit("fitfunc");
	gr4->GetFunction("fitfunc")->SetLineWidth(1);
	gr4->GetFunction("fitfunc")->SetLineColor(1);

	c3->SaveAs("function_fit_erf.pdf");

/*	TF1 *fitfunc = new TF1("fitfunc", "[0]*(-(x - TMath::Log(56))/(TMath::Sqrt(1 + (x - TMath::Log(56))*(x - TMath::Log(56))/(2*[1])))/[1])", -0.1, 4.5);
	fitfunc->SetParameter(0,1.);
	fitfunc->SetParameter(1,2.);
	gr1->Fit("fitfunc");
	gr1->GetFunction("fitfunc")->SetLineWidth(1);
	gr1->GetFunction("fitfunc")->SetLineColor(1);
	gr2->Fit("fitfunc");
	gr2->GetFunction("fitfunc")->SetLineWidth(1);
	gr2->GetFunction("fitfunc")->SetLineColor(1);
	gr3->Fit("fitfunc");
	gr3->GetFunction("fitfunc")->SetLineWidth(1);
	gr3->GetFunction("fitfunc")->SetLineColor(1);
	gr4->Fit("fitfunc");
	gr4->GetFunction("fitfunc")->SetLineWidth(1);
	gr4->GetFunction("fitfunc")->SetLineColor(1);

	c3->SaveAs("function_fit_sqrt.pdf");*/

/*	TF1 *fitfunc = new TF1("fitfunc", "[0]*TMath::ATan(-(x - TMath::Log(56))/[1])/(TMath::Pi()/2.)", -0.1, 4.5);
	fitfunc->SetParameter(0,1./2.);
	fitfunc->SetParameter(1,2.);
	gr1->Fit("fitfunc");
	gr1->GetFunction("fitfunc")->SetLineWidth(1);
	gr1->GetFunction("fitfunc")->SetLineColor(1);
	gr2->Fit("fitfunc");
	gr2->GetFunction("fitfunc")->SetLineWidth(1);
	gr2->GetFunction("fitfunc")->SetLineColor(1);
	gr3->Fit("fitfunc");
	gr3->GetFunction("fitfunc")->SetLineWidth(1);
	gr3->GetFunction("fitfunc")->SetLineColor(1);
	gr4->Fit("fitfunc");
	gr4->GetFunction("fitfunc")->SetLineWidth(1);
	gr4->GetFunction("fitfunc")->SetLineColor(1);

	c3->SaveAs("function_fit_arctan.pdf");*/
}
