{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	// 11 = p-O mixture, 12 = He-O mixture, 13 = p-He mixture
	double x11[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(16), 0.75*TMath::Log(1)  + 0.25*TMath::Log(16), 0.5*TMath::Log(1)  + 0.5*TMath::Log(16), 0.25*TMath::Log(1)  + 0.75*TMath::Log(16), 0.*TMath::Log(1)  + 1.*TMath::Log(16) };
	double x12[] = { 1.*TMath::Log(4)  + 0.*TMath::Log(16), 0.75*TMath::Log(4)  + 0.25*TMath::Log(16), 0.5*TMath::Log(4)  + 0.5*TMath::Log(16), 0.25*TMath::Log(4)  + 0.75*TMath::Log(16), 0.*TMath::Log(4)  + 1.*TMath::Log(16) };
	double x13[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(4) , 0.75*TMath::Log(1)  + 0.25*TMath::Log(4) , 0.5*TMath::Log(1)  + 0.5*TMath::Log(4) , 0.25*TMath::Log(1)  + 0.75*TMath::Log(4) , 0.*TMath::Log(1)  + 1.*TMath::Log(4)  };
	double xerr[] = { 0., 0., 0., 0. , 0., };

	double y11[] =   { 1.02785, 0.87109, 0.74186, 0.61165, 0.49253 };
	double y11le[] = { 0.01372, 0.01552, 0.01529, 0.02006, 0.02404 };
	double y11he[] = { 0.01267, 0.01550, 0.01667, 0.02110, 0.02343 };
	double y12[] =   { 0.91522, 0.79732, 0.69518, 0.59112, 0.49253 };
	double y12le[] = { 0.01642, 0.01698, 0.01494, 0.02073, 0.02404 };
	double y12he[] = { 0.01416, 0.01591, 0.01729, 0.02134, 0.02343 };
	double y13[] =   { 1.02785, 0.98224, 0.95279, 0.92899, 0.91522 };
	double y13le[] = { 0.01372, 0.01396, 0.01167, 0.01481, 0.01642 };
	double y13he[] = { 0.01267, 0.01231, 0.00934, 0.01286, 0.01416 };

	nrp = 5;

	TGraphAsymmErrors *gr11 = new TGraphAsymmErrors(nrp, x11, y11, xerr, xerr, y11le, y11he);
	gr11->SetMarkerStyle(20);
	gr11->SetMarkerColor(2);
	gr11->SetName("p-O_data");
	gr11->GetXaxis()->SetRangeUser(-0.1,4.5);
	gr11->GetYaxis()->SetRangeUser(-0.1,1.6);
	gr11->GetXaxis()->SetLimits(-0.1,4.5);
	gr11->GetYaxis()->SetLimits(-0.1,1.6);
	gr11->GetXaxis()->SetTitle("lnA");
	gr11->GetYaxis()->SetTitle("MVA fraction");
	gr11->Draw("AP");

	TGraphAsymmErrors *gr12 = new TGraphAsymmErrors(nrp, x12, y12, xerr, xerr, y12le, y12he);
	gr12->SetMarkerStyle(20);
	gr12->SetMarkerColor(3);
	gr12->SetName("He-O_data");
	gr12->Draw("SAME;P");

	TGraphAsymmErrors *gr13 = new TGraphAsymmErrors(nrp, x13, y13, xerr, xerr, y13le, y13he);
	gr13->SetMarkerStyle(20);
	gr13->SetMarkerColor(4);
	gr13->SetName("p-He_data");
	gr13->Draw("SAME;P");

	cout << endl << "Fitting lnA_p-Fe_others.cpp ------------------------------------------" << endl;
	TF1 *fit1 = new TF1("fit1", "([0]-[1])*((x - TMath::Log(16))/(TMath::Log(1) - TMath::Log(16))) + [1]", -0.1, 4.5);
	fit1->SetParameter(0, 0.1);
	fit1->SetParameter(1, 0.);
	gr11->Fit("fit1");
	gr11->GetFunction("fit1")->SetLineWidth(1);
	gr11->GetFunction("fit1")->SetLineColor(1);
	TF1 *fit1 = new TF1("fit1", "([0]-[1])*((x - TMath::Log(16))/(TMath::Log(4) - TMath::Log(16))) + [1]", -0.1, 4.5);
	fit1->SetParameter(0, 0.1);
	fit1->SetParameter(1, 0.);
	gr12->Fit("fit1");
	gr12->GetFunction("fit1")->SetLineWidth(1);
	gr12->GetFunction("fit1")->SetLineColor(1);
	TF1 *fit1 = new TF1("fit1", "([0]-[1])*((x - TMath::Log(4))/(TMath::Log(1)  - TMath::Log(4)))  + [1]", -0.1, 4.5);
	fit1->SetParameter(0, 0.1);
	fit1->SetParameter(1, 0.);
	gr13->Fit("fit1");
	gr13->GetFunction("fit1")->SetLineWidth(1);
	gr13->GetFunction("fit1")->SetLineColor(1);
/*	gr11->Fit("pol1");
	gr11->GetFunction("pol1")->SetLineWidth(1);
	gr11->GetFunction("pol1")->SetLineColor(1);
	gr12->Fit("pol1");
	gr12->GetFunction("pol1")->SetLineWidth(1);
	gr12->GetFunction("pol1")->SetLineColor(1);
	gr13->Fit("pol1");
	gr13->GetFunction("pol1")->SetLineWidth(1);
	gr13->GetFunction("pol1")->SetLineColor(1);*/

	TLine *l1 = new TLine(TMath::Log(1), -0.1, TMath::Log(1), 1.6);
	l1->SetLineWidth(1);
	l1->SetLineColor(6);
	l1->SetLineStyle(1);
	l1->Draw("same");
	TLine *l2 = new TLine(TMath::Log(4), -0.1, TMath::Log(4), 1.6);
	l2->SetLineWidth(1);
	l2->SetLineColor(6);
	l2->SetLineStyle(1);
	l2->Draw("same");
	TLine *l3 = new TLine(TMath::Log(16), -0.1, TMath::Log(16), 1.6);
	l3->SetLineWidth(1);
	l3->SetLineColor(6);
	l3->SetLineStyle(1);
	l3->Draw("same");
	TLine *l4 = new TLine(TMath::Log(56), -0.1, TMath::Log(56), 1.6);
	l4->SetLineWidth(1);
	l4->SetLineColor(6);
	l4->SetLineStyle(1);
	l4->Draw("same");

	c1->SaveAs("fit-lnA_p-Fe_others.pdf");
}
