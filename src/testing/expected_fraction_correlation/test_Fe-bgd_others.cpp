{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp = 3;

	// Points are: 0 = 50% sig + 50% bgd, 1 = 75% sig + 25% bgd, 2 = 25% sig + 75% bgd
	// x1, x11, x21: sig = p, bgd = O
	// x2, x12, x22: sig = He, bgd = O
	// x3, x13, x23: sig = p, bgd = He

	// p/Fe treatment
	double x1[] = { 0.65561, 0.82780, 0.48341 };
	double x2[] = { 0.48341, 0.56951, 0.39732 };
	double x3[] = { 0.82780, 0.91390, 0.74171 };

	double y1[] = { 0.74186, 0.87109, 0.61165 };
	double y2[] = { 0.69518, 0.79732, 0.59112 };
	double y3[] = { 0.95279, 0.98224, 0.92899 };

	TGraph *gr1 = new TGraph(nrp, x1, y1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(2);
	gr1->SetName("p-sig_O-bgd_p-Fe");
	gr1->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr1->GetYaxis()->SetRangeUser(-0.1,1.5);
	gr1->GetXaxis()->SetLimits(-0.1,1.1);
	gr1->GetYaxis()->SetLimits(-0.1,1.5);
	gr1->GetXaxis()->SetTitle("Expected proton fraction");
	gr1->GetYaxis()->SetTitle("Proton fraction");
	gr1->Draw("AP");

	TGraph *gr2 = new TGraph(nrp, x2, y2);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(3);
	gr2->SetName("He-sig_O-bgd_p-Fe");
	gr2->Draw("SAME;P");

	TGraph *gr3 = new TGraph(nrp, x3, y3);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(4);
	gr3->SetName("p-sig_He-bgd_p-Fe");
	gr3->Draw("SAME;P");

	// He/Fe treatment
	double x11[] = { 0.65561, 0.82780, 0.48341 };
	double x12[] = { 0.48341, 0.56951, 0.39732 };
	double x13[] = { 0.82780, 0.91390, 0.74171 };

	double y11[] = { 0.82663, 0.94663, 0.71108 };
	double y12[] = { 0.78339, 0.88815, 0.69677 };
	double y13[] = { 1.02850, 1.05370, 1.00954 };

	TGraph *gr11 = new TGraph(nrp, x11, y11);
	gr11->SetMarkerStyle(21);
	gr11->SetMarkerColor(2);
	gr11->SetMarkerSize(0.9);
	gr11->SetName("p-sig_O-bgd_He-Fe");
	gr11->Draw("SAME;P");

	TGraph *gr12 = new TGraph(nrp, x12, y12);
	gr12->SetMarkerStyle(21);
	gr12->SetMarkerColor(3);
	gr12->SetMarkerSize(0.9);
	gr12->SetName("He-sig_O-bgd_He-Fe");
	gr12->Draw("SAME;P");

	TGraph *gr13 = new TGraph(nrp, x13, y13);
	gr13->SetMarkerStyle(21);
	gr13->SetMarkerColor(4);
	gr13->SetMarkerSize(0.9);
	gr13->SetName("p-sig_He-bgd_He-Fe");
	gr13->Draw("SAME;P");

	// O/Fe treatment
	double x21[] = { 0.65561, 0.82780, 0.48341 };
	double x22[] = { 0.48341, 0.56951, 0.39732 };
	double x23[] = { 0.82780, 0.91390, 0.74171 };

	double y21[] = { 1.21514, 1.31434, 1.09778 };
	double y22[] = { 1.18273, 1.26855, 1.08298 };
	double y23[] = { 1.38972, 1.40887, 1.37787 };

	TGraph *gr21 = new TGraph(nrp, x21, y21);
	gr21->SetMarkerStyle(22);
	gr21->SetMarkerColor(2);
	gr21->SetName("p-sig_O-bgd_O-Fe");
	gr21->Draw("SAME;P");

	TGraph *gr22 = new TGraph(nrp, x22, y22);
	gr22->SetMarkerStyle(22);
	gr22->SetMarkerColor(3);
	gr22->SetName("He-sig_O-bgd_O-Fe");
	gr22->Draw("SAME;P");

	TGraph *gr23 = new TGraph(nrp, x23, y23);
	gr23->SetMarkerStyle(22);
	gr23->SetMarkerColor(4);
	gr23->SetName("p-sig_He-bgd_O-Fe");
	gr23->Draw("SAME;P");

	cout << endl << "Fitting test_Fe-bgd.cpp ---------------------------------------------" << endl;
	gr1->Fit("pol1");
	gr1->GetFunction("pol1")->SetLineWidth(1);
	gr1->GetFunction("pol1")->SetLineColor(1);
	gr2->Fit("pol1");
	gr2->GetFunction("pol1")->SetLineWidth(1);
	gr2->GetFunction("pol1")->SetLineColor(1);
	gr3->Fit("pol1");
	gr3->GetFunction("pol1")->SetLineWidth(1);
	gr3->GetFunction("pol1")->SetLineColor(1);

	gr11->Fit("pol1");
	gr11->GetFunction("pol1")->SetLineWidth(1);
	gr11->GetFunction("pol1")->SetLineColor(1);
	gr12->Fit("pol1");
	gr12->GetFunction("pol1")->SetLineWidth(1);
	gr12->GetFunction("pol1")->SetLineColor(1);
	gr13->Fit("pol1");
	gr13->GetFunction("pol1")->SetLineWidth(1);
	gr13->GetFunction("pol1")->SetLineColor(1);

	gr21->Fit("pol1");
	gr21->GetFunction("pol1")->SetLineWidth(1);
	gr21->GetFunction("pol1")->SetLineColor(1);
	gr22->Fit("pol1");
	gr22->GetFunction("pol1")->SetLineWidth(1);
	gr22->GetFunction("pol1")->SetLineColor(1);
	gr23->Fit("pol1");
	gr23->GetFunction("pol1")->SetLineWidth(1);
	gr23->GetFunction("pol1")->SetLineColor(1);

/*	TLine *l1 = new TLine(-0.1, 0.55941, 1.1, 0.55941);
	l1->SetLineWidth(1);
	l1->SetLineColor(6);
	l1->SetLineStyle(9);
	l1->Draw("same");
	TLine *l2 = new TLine(-0.1, 0.67013, 1.1, 0.67013);
	l2->SetLineWidth(1);
	l2->SetLineColor(6);
	l2->SetLineStyle(9);
	l2->Draw("same");
	TLine *l3 = new TLine(-0.1, 1.01530, 1.1, 1.01530);
	l3->SetLineWidth(1);
	l3->SetLineColor(6);
	l3->SetLineStyle(9);
	l3->Draw("same");*/

	c1->SaveAs("fit_Fe-bgd_others.pdf");
}
