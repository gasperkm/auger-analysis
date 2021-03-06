{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	// 11 = p-O mixture, 12 = He-O mixture, 13 = p-He mixture
	double x11[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(16), 0.75*TMath::Log(1)  + 0.25*TMath::Log(16), 0.5*TMath::Log(1)  + 0.5*TMath::Log(16), 0.25*TMath::Log(1)  + 0.75*TMath::Log(16), 0.*TMath::Log(1)  + 1.*TMath::Log(16) };
	double x12[] = { 1.*TMath::Log(4)  + 0.*TMath::Log(16), 0.75*TMath::Log(4)  + 0.25*TMath::Log(16), 0.5*TMath::Log(4)  + 0.5*TMath::Log(16), 0.25*TMath::Log(4)  + 0.75*TMath::Log(16), 0.*TMath::Log(4)  + 1.*TMath::Log(16) };
	double x13[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(4) , 0.75*TMath::Log(1)  + 0.25*TMath::Log(4) , 0.5*TMath::Log(1)  + 0.5*TMath::Log(4) , 0.25*TMath::Log(1)  + 0.75*TMath::Log(4) , 0.*TMath::Log(1)  + 1.*TMath::Log(4)  };
	double xerr[] = { 0., 0., 0., 0. , 0., };

	double y11[] =   { 1.46308, 1.31434, 1.21514, 1.09778, 1.06299 };
	double y11le[] = { 0.13022, 0.09059, 0.09485, 0.08487, 0.11278 };
	double y11he[] = { 0.12181, 0.08422, 0.08946, 0.07720, 0.11121 };
	double y12[] =   { 1.39969, 1.26855, 1.18273, 1.08298, 1.06299 };
	double y12le[] = { 0.12279, 0.08669, 0.09429, 0.08298, 0.11278 };
	double y12he[] = { 0.11333, 0.08143, 0.08962, 0.07576, 0.11121 };
	double y13[] =   { 1.46308, 1.40887, 1.38972, 1.37787, 1.39969 };
	double y13le[] = { 0.13022, 0.09540, 0.09952, 0.09302, 0.12279 };
	double y13he[] = { 0.12181, 0.08910, 0.09531, 0.08628, 0.11333 };

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

	cout << endl << "Fitting lnA_O-Fe_others.cpp ------------------------------------------" << endl;
	gr11->Fit("pol1");
	gr11->GetFunction("pol1")->SetLineWidth(1);
	gr11->GetFunction("pol1")->SetLineColor(1);
	gr12->Fit("pol1");
	gr12->GetFunction("pol1")->SetLineWidth(1);
	gr12->GetFunction("pol1")->SetLineColor(1);
	gr13->Fit("pol1");
	gr13->GetFunction("pol1")->SetLineWidth(1);
	gr13->GetFunction("pol1")->SetLineColor(1);

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

	c1->SaveAs("fit-lnA_O-Fe_others.pdf");
}
