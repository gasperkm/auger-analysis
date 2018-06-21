{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	double x1[] = { 1., 0.65561, 0.31122, 0.}; 
//	double x2[] = { 1., 0.65561, 0.31122, 0.};
	double x3[] = { 0.5, 0.32780, 0.15561};
	double x4[] = { 0.75, 0.49171, 0.23341};
	double x5[] = { 0.25, 0.16390, 0.07780};

	double y1[] = { 1.07134, 0.741378, 0.00337782, -0.341724 };
//        double y2[] = { 1.03657, 0.90446, 0.50348, -0.00420};
	double y3[] = { 0.341808, 0.21745, -0.141539 };
	double y4[] = { 0.658748, 0.469511, -0.0690796 };
	double y5[] = { 0.00712702, -0.0518741, -0.23107 };

	nrp = 4;

	TGraph *gr1 = new TGraph(nrp, x1, y1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(1);
	gr1->SetName("pure_600");
	gr1->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr1->GetYaxis()->SetRangeUser(-0.3,1.5);
	gr1->GetXaxis()->SetLimits(-0.1,1.1);
	gr1->GetYaxis()->SetLimits(-0.3,1.5);
	gr1->GetXaxis()->SetTitle("Expected proton fraction");
	gr1->GetYaxis()->SetTitle("Proton fraction");
	gr1->Draw("AP");

/*	TGraph *gr2 = new TGraph(nrp, x2, y2);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(1);
	gr2->SetName("pure_300");
	gr2->Draw("SAME;P");*/

	nrp = 3;

	TGraph *gr3 = new TGraph(nrp, x3, y3);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(2);
	gr3->SetName("50-50_O-bgd");
	gr3->Draw("SAME;P");

	TGraph *gr4 = new TGraph(nrp, x4, y4);
	gr4->SetMarkerStyle(20);
	gr4->SetMarkerColor(4);
	gr4->SetName("75-25_O-bgd");
	gr4->Draw("SAME;P");

	TGraph *gr5 = new TGraph(nrp, x5, y5);
	gr5->SetMarkerStyle(20);
	gr5->SetMarkerColor(3);
	gr5->SetName("25-75_O-bgd");
	gr5->Draw("SAME;P");

	TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
	c2->SetGrid();

	double x11[] = { 1.,  0.5, 0.75, 0.25 };
	double x12[] = { 0.65561, 0.32780, 0.49171, 0.16390 };
	double x13[] = { 0.31122, 0.15561, 0.23341, 0.07780 };

	double y11[] = { 1.07134, 0.341808, 0.658748, 0.00712702 };
	double y12[] = { 0.741378, 0.21745, 0.469511, -0.0518741 };
	double y13[] = { 0.00337782, -0.141539, -0.0690796, -0.23107 };

	nrp = 4;

	TGraph *gr11 = new TGraph(nrp, x11, y11);
	gr11->SetMarkerStyle(20);
	gr11->SetMarkerColor(2);
	gr11->SetName("p-sig");
	gr11->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr11->GetYaxis()->SetRangeUser(-0.3,1.5);
	gr11->GetXaxis()->SetLimits(-0.1,1.1);
	gr11->GetYaxis()->SetLimits(-0.3,1.5);
	gr11->GetXaxis()->SetTitle("Expected proton fraction");
	gr11->GetYaxis()->SetTitle("Proton fraction");
	gr11->Draw("AP");

	TGraph *gr12 = new TGraph(nrp, x12, y12);
	gr12->SetMarkerStyle(20);
	gr12->SetMarkerColor(3);
	gr12->SetName("He-sig");
	gr12->Draw("SAME;P");

	TGraph *gr13 = new TGraph(nrp, x13, y13);
	gr13->SetMarkerStyle(20);
	gr13->SetMarkerColor(4);
	gr13->SetName("O-sig");
	gr13->Draw("SAME;P");

	cout << endl << "Fitting test_p-O.cpp ------------------------------------------------" << endl;
	gr11->Fit("pol1");
	gr11->GetFunction("pol1")->SetLineWidth(1);
	gr11->GetFunction("pol1")->SetLineColor(1);
	gr12->Fit("pol1");
	gr12->GetFunction("pol1")->SetLineWidth(1);
	gr12->GetFunction("pol1")->SetLineColor(1);
	gr13->Fit("pol1");
	gr13->GetFunction("pol1")->SetLineWidth(1);
	gr13->GetFunction("pol1")->SetLineColor(1);

	c2->SaveAs("fit_p-O.pdf");
}
