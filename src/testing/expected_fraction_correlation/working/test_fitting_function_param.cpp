{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	// 0 = ln1, 1 = ln4, 2 = ln16
	double x1[] = { TMath::Log(1), TMath::Log(4), TMath::Log(16) }; 

	// y1 = p/Fe, y2 = He/Fe, y3 = O/Fe
/*	double y1[] = { 1.02913, 1.38459, 1.57016 };
	double y2[] = { 1.07604, 1.56753, 1.87619 };
	double y3[] = { 1.45146, 2.11678, 3.39675 };*/

	// y1 = p-Fe mixture, y2 = He-Fe mixture, y3 = O-Fe mixture
/*	double y1[] = { 1.02913, 1.07604, 1.45146 };
	double y2[] = { 1.38459, 1.56753, 2.11678 };
	double y3[] = { 1.57016, 1.87619, 3.39675 };*/

	// y1 = p/Fe, y2 = He/Fe, y3 = O/Fe
/*	double y1[] = { 45.82, 54.16, 57.51 };
	double y2[] = { 47.10, 57.46, 61.94 };
	double y3[] = { 55.43, 64.71, 73.60 };*/

	// y1 = p-Fe mixture, y2 = He-Fe mixture, y3 = O-Fe mixture
	double y1[] = { 45.82, 47.10, 55.43 };
	double y2[] = { 54.16, 57.46, 64.71 };
	double y3[] = { 57.51, 61.94, 73.60 };
	// y4 = p-O mixture, y5 = He-O mixture, y6 = p-Fe mixture
	double y4[] = { 36.99, 34.37, 32.16 };
	double y5[] = { 50.14, 48.02, 47.14};
	double y6[] = { 17.18, 14.38, 10.21 };

	nrp = 3;

	TGraph *gr1 = new TGraph(nrp, x1, y1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(2);
	gr1->SetName("p-Fe");
	gr1->GetXaxis()->SetRangeUser(-0.1,3.5);
//	gr1->GetYaxis()->SetRangeUser(-0.1,4.);
	gr1->GetYaxis()->SetRangeUser(0.,90.);
	gr1->GetXaxis()->SetLimits(-0.1,3.5);
//	gr1->GetYaxis()->SetLimits(-0.1,4.);
	gr1->GetYaxis()->SetLimits(0.,90.);
	gr1->GetXaxis()->SetTitle("lnA");
//	gr1->GetYaxis()->SetTitle("Fitting function parameter");
	gr1->GetYaxis()->SetTitle("Fitting function inclination (deg)");
	gr1->Draw("AP");

	TGraph *gr2 = new TGraph(nrp, x1, y2);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(3);
	gr2->SetName("He-Fe");
	gr2->Draw("SAME;P");

	TGraph *gr3 = new TGraph(nrp, x1, y3);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(4);
	gr3->SetName("O-Fe");
	gr3->Draw("SAME;P");

	TGraph *gr4 = new TGraph(nrp, x1, y4);
	gr4->SetMarkerStyle(21);
	gr4->SetMarkerColor(2);
	gr4->SetMarkerSize(0.9);
	gr4->SetName("p-O");
	gr4->Draw("SAME;P");

	TGraph *gr5 = new TGraph(nrp, x1, y5);
	gr5->SetMarkerStyle(21);
	gr5->SetMarkerColor(3);
	gr5->SetMarkerSize(0.9);
	gr5->SetName("He-O");
	gr5->Draw("SAME;P");

	TGraph *gr6 = new TGraph(nrp, x1, y6);
	gr6->SetMarkerStyle(22);
	gr6->SetMarkerColor(2);
	gr6->SetName("p-He");
	gr6->Draw("SAME;P");

	cout << endl << "Fitting test_fitting_finction_param.cpp -----------------------------" << endl;
	gr1->Fit("pol2");
	gr1->GetFunction("pol2")->SetLineWidth(1);
	gr1->GetFunction("pol2")->SetLineColor(1);
	gr2->Fit("pol2");
	gr2->GetFunction("pol2")->SetLineWidth(1);
	gr2->GetFunction("pol2")->SetLineColor(1);
	gr3->Fit("pol2");
	gr3->GetFunction("pol2")->SetLineWidth(1);
	gr3->GetFunction("pol2")->SetLineColor(1);
	gr4->Fit("pol2");
	gr4->GetFunction("pol2")->SetLineWidth(1);
	gr4->GetFunction("pol2")->SetLineColor(1);
	gr5->Fit("pol2");
	gr5->GetFunction("pol2")->SetLineWidth(1);
	gr5->GetFunction("pol2")->SetLineColor(1);
	gr6->Fit("pol2");
	gr6->GetFunction("pol2")->SetLineWidth(1);
	gr6->GetFunction("pol2")->SetLineColor(1);

	c1->SaveAs("fit_fitting_function_param.pdf");
}
