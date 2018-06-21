{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	double x1[] = { 0.65561, 0.48341, 0.82780 }; 
	double x2[] = { 0.82780, 0.56951, 0.91390 };
	double x3[] = { 0.48341, 0.39732, 0.74171 };

	double y1[] = { 1.21514, 1.18273, 1.38972 };
        double y2[] = { 1.31434, 1.26855, 1.40887 };
	double y3[] = { 1.09778, 1.08298, 1.37787 };

	nrp = 3;

	TGraph *gr1 = new TGraph(nrp, x1, y1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(2);
	gr1->SetName("50sig-50bgd");
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
	gr2->SetName("75sig-25bgd");
	gr2->Draw("SAME;P");

	TGraph *gr3 = new TGraph(nrp, x3, y3);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(4);
	gr3->SetName("25sig-75bgd");
	gr3->Draw("SAME;P");

	TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
	c2->SetGrid();

	// 11 = p-O mixture, 12 = He-O mixture, 13 = p-He mixture
	// 0 = 100% sig, 1 = 50-50, 2 = 75-25, 3 = 25-75, 4 = 100% bgd
	double x11[] = { 1., 0.65561, 0.82780, 0.48341, 0.31122 };
	double x12[] = { 0.65561, 0.48341, 0.56951, 0.39732, 0.31122 };
	double x13[] = { 1., 0.82780, 0.91390, 0.74171, 0.65561 };

	double y11[] = { 1.46308, 1.21514, 1.31434, 1.09778, 1.06299 };
	double y12[] = { 1.39969, 1.18273, 1.26855, 1.08298, 1.06299 };
	double y13[] = { 1.46308, 1.38972, 1.40887, 1.37787, 1.39969 };

	nrp = 5;

	TGraph *gr11 = new TGraph(nrp, x11, y11);
	gr11->SetMarkerStyle(22);
	gr11->SetMarkerColor(2);
	gr11->SetName("p-O_data");
	gr11->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr11->GetYaxis()->SetRangeUser(-0.1,1.5);
	gr11->GetXaxis()->SetLimits(-0.1,1.1);
	gr11->GetYaxis()->SetLimits(-0.1,1.5);
	gr11->GetXaxis()->SetTitle("Expected proton fraction");
	gr11->GetYaxis()->SetTitle("Proton fraction");
	gr11->Draw("AP");

	TGraph *gr12 = new TGraph(nrp, x12, y12);
	gr12->SetMarkerStyle(22);
	gr12->SetMarkerColor(3);
	gr12->SetName("He-O_data");
	gr12->Draw("SAME;P");

	TGraph *gr13 = new TGraph(nrp, x13, y13);
	gr13->SetMarkerStyle(22);
	gr13->SetMarkerColor(4);
	gr13->SetName("p-He_data");
	gr13->Draw("SAME;P");

	cout << endl << "Fitting test_O-Fe_others.cpp ----------------------------------------" << endl;
	gr11->Fit("pol1");
	gr11->GetFunction("pol1")->SetLineWidth(1);
	gr11->GetFunction("pol1")->SetLineColor(1);
	gr12->Fit("pol1");
	gr12->GetFunction("pol1")->SetLineWidth(1);
	gr12->GetFunction("pol1")->SetLineColor(1);
	gr13->Fit("pol1");
	gr13->GetFunction("pol1")->SetLineWidth(1);
	gr13->GetFunction("pol1")->SetLineColor(1);

	c2->SaveAs("fit_O-Fe_others.pdf");

/*	TCanvas *c3 = new TCanvas("c3", "", 1200, 800);
	c3->SetGrid();

	double x21[3], x22[3], x23[3];
	double y21[3], y22[3], y23[3];

	for(int i = 0; i < nrp; i++)
	{
		x21[i] = x11[i] - (1 - TMath::Log(16)/TMath::Log(56));
		y21[i] = y11[i];
		x22[i] = x12[i] - (1 - TMath::Log(16)/TMath::Log(56));
		y22[i] = y12[i];
		x23[i] = x13[i] - (1 - TMath::Log(4)/TMath::Log(56));
		y23[i] = y13[i];
	}

	TGraph *gr21 = new TGraph(nrp, x21, y21);
	gr21->SetMarkerStyle(22);
	gr21->SetMarkerColor(2);
	gr21->SetName("p-O_data");
	gr21->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr21->GetYaxis()->SetRangeUser(-0.1,1.5);
	gr21->GetXaxis()->SetLimits(-0.1,1.1);
	gr21->GetYaxis()->SetLimits(-0.1,1.5);
	gr21->GetXaxis()->SetTitle("Expected proton fraction (shifted)");
	gr21->GetYaxis()->SetTitle("Proton fraction (shifted)");
	gr21->Draw("AP");

	TGraph *gr22 = new TGraph(nrp, x22, y22);
	gr22->SetMarkerStyle(22);
	gr22->SetMarkerColor(3);
	gr22->SetName("He-O_data");
	gr22->Draw("SAME;P");

	TGraph *gr23 = new TGraph(nrp, x23, y23);
	gr23->SetMarkerStyle(22);
	gr23->SetMarkerColor(4);
	gr23->SetName("p-He_data");
	gr23->Draw("SAME;P");

	cout << endl << "Fitting test_O-Fe_others.cpp (shifted) ------------------------------" << endl;
	gr21->Fit("pol1");
	gr21->GetFunction("pol1")->SetLineWidth(1);
	gr21->GetFunction("pol1")->SetLineColor(1);
	gr22->Fit("pol1");
	gr22->GetFunction("pol1")->SetLineWidth(1);
	gr22->GetFunction("pol1")->SetLineColor(1);
	gr23->Fit("pol1");
	gr23->GetFunction("pol1")->SetLineWidth(1);
	gr23->GetFunction("pol1")->SetLineColor(1);

	c3->SaveAs("fit_O-Fe_others-shift.pdf");*/
}
