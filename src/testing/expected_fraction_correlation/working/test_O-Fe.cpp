{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	double x1[] = { 1., 0.65561, 0.31122, 0. }; 
//	double x2[] = { 1., 0.65561, 0.31122, 0. };
	double x3[] = { 0.5, 0.32780, 0.15561 };
	double x4[] = { 0.75, 0.49171, 0.23341 };
	double x5[] = { 0.25, 0.16390, 0.07780 };

	double y1[] = { 1.46308, 1.39969, 1.06299, -0.01271 };
//        double y2[] = { 1.49835, 1.42090, 1.15386, -0.01329 };
	double y3[] = { 0.71115, 0.67875, 0.50417 };
	double y4[] = { 1.07554, 1.02974, 0.74695 };
	double y5[] = { 0.37499, 0.36019, 0.26296 };

	nrp = 4;

	TGraph *gr1 = new TGraph(nrp, x1, y1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(1);
	gr1->SetName("pure_600");
	gr1->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr1->GetYaxis()->SetRangeUser(-0.1,1.5);
	gr1->GetXaxis()->SetLimits(-0.1,1.1);
	gr1->GetYaxis()->SetLimits(-0.1,1.5);
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
	gr3->SetName("50-50_Fe-bgd");
	gr3->Draw("SAME;P");

	TGraph *gr4 = new TGraph(nrp, x4, y4);
	gr4->SetMarkerStyle(20);
	gr4->SetMarkerColor(4);
	gr4->SetName("75-25_Fe-bgd");
	gr4->Draw("SAME;P");

	TGraph *gr5 = new TGraph(nrp, x5, y5);
	gr5->SetMarkerStyle(20);
	gr5->SetMarkerColor(3);
	gr5->SetName("25-75_Fe-bgd");
	gr5->Draw("SAME;P");

	TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
	c2->SetGrid();

/*	double x11[] = { 1., 1., 0.5, 0.75, 0.25 };
	double x12[] = { 0.65561, 0.65561, 0.32780, 0.49171, 0.16390 };
	double x13[] = { 0.31122, 0.31122, 0.15561, 0.23341, 0.07780 };

	double y11[] = { 1.46308, 1.49835, 0.71115, 1.07554, 0.37499 };
	double y12[] = { 1.39969, 1.42090, 0.67875, 1.02974, 0.36019 };
	double y13[] = { 1.06299, 1.15386, 0.50417, 0.74695, 0.26296 };*/

	double x11[] = { 1., 0.5, 0.75, 0.25 };
	double x12[] = { 0.65561, 0.32780, 0.49171, 0.16390 };
	double x13[] = { 0.31122, 0.15561, 0.23341, 0.07780 };

	double y11[] = { 1.46308, 0.71115, 1.07554, 0.37499 };
	double y12[] = { 1.39969, 0.67875, 1.02974, 0.36019 };
	double y13[] = { 1.06299, 0.50417, 0.74695, 0.26296 };

	nrp = 4;

	TGraph *gr11 = new TGraph(nrp, x11, y11);
	gr11->SetMarkerStyle(22);
	gr11->SetMarkerColor(2);
	gr11->SetName("p-sig");
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
	gr12->SetName("He-sig");
	gr12->Draw("SAME;P");

	TGraph *gr13 = new TGraph(nrp, x13, y13);
	gr13->SetMarkerStyle(22);
	gr13->SetMarkerColor(4);
	gr13->SetName("O-sig");
	gr13->Draw("SAME;P");

	cout << endl << "Fitting test_O-Fe.cpp -----------------------------------------------" << endl;
	gr11->Fit("pol1");
	gr11->GetFunction("pol1")->SetLineWidth(1);
	gr11->GetFunction("pol1")->SetLineColor(1);
	gr12->Fit("pol1");
	gr12->GetFunction("pol1")->SetLineWidth(1);
	gr12->GetFunction("pol1")->SetLineColor(1);
	gr13->Fit("pol1");
	gr13->GetFunction("pol1")->SetLineWidth(1);
	gr13->GetFunction("pol1")->SetLineColor(1);

	c2->SaveAs("fit_O-Fe.pdf");
}
