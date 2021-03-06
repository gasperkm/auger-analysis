{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	double x1[] = { 1.*TMath::Log(1) + 0.*TMath::Log(56), 1.*TMath::Log(4) + 0.*TMath::Log(56), 1.*TMath::Log(16) + 0.*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)}; 
	double x2[] = { 0.85*TMath::Log(1) + 0.15*TMath::Log(56), 0.85*TMath::Log(4) + 0.15*TMath::Log(56), 0.85*TMath::Log(16) + 0.15*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double x3[] = { 0.75*TMath::Log(1) + 0.25*TMath::Log(56), 0.75*TMath::Log(4) + 0.25*TMath::Log(56), 0.75*TMath::Log(16) + 0.25*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double x4[] = { 0.5*TMath::Log(1) + 0.5*TMath::Log(56), 0.5*TMath::Log(4) + 0.5*TMath::Log(56), 0.5*TMath::Log(16) + 0.5*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double x5[] = { 0.25*TMath::Log(1) + 0.75*TMath::Log(56), 0.25*TMath::Log(4) + 0.75*TMath::Log(56), 0.25*TMath::Log(16) + 0.75*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double x6[] = { 0.15*TMath::Log(1) + 0.85*TMath::Log(56), 0.15*TMath::Log(4) + 0.85*TMath::Log(56), 0.15*TMath::Log(16) + 0.85*TMath::Log(56), 0.*TMath::Log(1) + 1.*TMath::Log(56)};
	double xe[] = { 0., 0., 0., 0. };

	double y1[] =   { 1.08159, 1.03207, 0.58933, -0.00021 };
	double y1le[] = { 0.02015, 0.02302, 0.03416, 0.02557  };
	double y1he[] = { 0.01645, 0.02173, 0.03493, 0.02278  };
	double y2[] =   { 0.91210, 0.85012, 0.51078, -0.00021 };
	double y2le[] = { 0.01635, 0.01726, 0.03113, 0.02557  };
	double y2he[] = { 0.01431, 0.01752, 0.03018, 0.02278  };
	double y3[] =   { 0.80524, 0.74677, 0.45123, -0.00021 };
	double y3le[] = { 0.01507, 0.01661, 0.02837, 0.02557  };
	double y3he[] = { 0.01397, 0.01616, 0.02599, 0.02278  };
	double y4[] =   { 0.53687, 0.49362, 0.29176, -0.00021 };
	double y4le[] = { 0.01727, 0.01842, 0.02589, 0.02557  };
	double y4he[] = { 0.01385, 0.01640, 0.02481, 0.02278  };
	double y5[] =   { 0.27435, 0.26004, 0.15589, -0.00021 };
	double y5le[] = { 0.01687, 0.01744, 0.02240, 0.02557 };
	double y5he[] = { 0.01803, 0.01940, 0.02353, 0.02278 };
	double y6[] =   { 0.17189, 0.15778, 0.10094, -0.00021 };
	double y6le[] = { 0.01707, 0.01879, 0.02007, 0.02557 };
	double y6he[] = { 0.01817, 0.01953, 0.02152, 0.02278 };

	nrp = 4;

	TGraphAsymmErrors *gr1 = new TGraphAsymmErrors(nrp, x1, y1, xe, xe, y1le, y1he);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(1);
	gr1->SetName("pure_600");
	gr1->GetXaxis()->SetRangeUser(-0.1,4.5);
	gr1->GetYaxis()->SetRangeUser(-0.1,1.6);
	gr1->GetXaxis()->SetLimits(-0.1,4.5);
	gr1->GetYaxis()->SetLimits(-0.1,1.6);
	gr1->GetXaxis()->SetTitle("lnA");
	gr1->GetYaxis()->SetTitle("MVA fraction");
	gr1->Draw("AP");

	TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(nrp, x2, y2, xe, xe, y2le, y2he);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(2);
	gr2->SetName("85-15_Fe-bgd");
	gr2->Draw("SAME;P");

	TGraphAsymmErrors *gr3 = new TGraphAsymmErrors(nrp, x3, y3, xe, xe, y3le, y3he);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(3);
	gr3->SetName("75-25_Fe-bgd");
	gr3->Draw("SAME;P");

	TGraphAsymmErrors *gr4 = new TGraphAsymmErrors(nrp, x4, y4, xe, xe, y4le, y4he);
	gr4->SetMarkerStyle(20);
	gr4->SetMarkerColor(4);
	gr4->SetName("50-50_Fe-bgd");
	gr4->Draw("SAME;P");

	TGraphAsymmErrors *gr5 = new TGraphAsymmErrors(nrp, x5, y5, xe, xe, y5le, y5he);
	gr5->SetMarkerStyle(20);
	gr5->SetMarkerColor(6);
	gr5->SetName("25-75_Fe-bgd");
	gr5->Draw("SAME;P");

	TGraphAsymmErrors *gr6 = new TGraphAsymmErrors(nrp, x6, y6, xe, xe, y6le, y6he);
	gr6->SetMarkerStyle(20);
	gr6->SetMarkerColor(7);
	gr6->SetName("15-85_Fe-bgd");
	gr6->Draw("SAME;P");

	TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
	c2->SetGrid();

	// 11 = p-Fe mixture, 12 = He-Fe mixture, 13 = O-Fe mixture
	double x11[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(56), 0.85*TMath::Log(1)  + 0.15*TMath::Log(56), 0.75*TMath::Log(1)  + 0.25*TMath::Log(56), 0.5*TMath::Log(1)  + 0.5*TMath::Log(56), 0.25*TMath::Log(1)  + 0.75*TMath::Log(56), 0.15*TMath::Log(1)  + 0.85*TMath::Log(56), 0.*TMath::Log(1)  + 1.*TMath::Log(56) };
	double x12[] = { 1.*TMath::Log(4)  + 0.*TMath::Log(56), 0.85*TMath::Log(4)  + 0.15*TMath::Log(56), 0.75*TMath::Log(4)  + 0.25*TMath::Log(56), 0.5*TMath::Log(4)  + 0.5*TMath::Log(56), 0.25*TMath::Log(4)  + 0.75*TMath::Log(56), 0.15*TMath::Log(4)  + 0.85*TMath::Log(56), 0.*TMath::Log(4)  + 1.*TMath::Log(56) };
	double x13[] = { 1.*TMath::Log(16) + 0.*TMath::Log(56), 0.85*TMath::Log(16) + 0.15*TMath::Log(56), 0.75*TMath::Log(16) + 0.25*TMath::Log(56), 0.5*TMath::Log(16) + 0.5*TMath::Log(56), 0.25*TMath::Log(16) + 0.75*TMath::Log(56), 0.15*TMath::Log(16) + 0.85*TMath::Log(56), 0.*TMath::Log(16) + 1.*TMath::Log(56) };
	double xerr[] = { 0., 0., 0., 0. , 0., 0., 0. };

	double y11[] =   { 1.08159, 0.91210, 0.80524, 0.53687, 0.27435, 0.17189, -0.00021 };
	double y11le[] = { 0.02015, 0.01635, 0.01507, 0.01727, 0.01687, 0.01707, 0.02557  };
	double y11he[] = { 0.01645, 0.01431, 0.01397, 0.01385, 0.01803, 0.01817, 0.02278  };
	double y12[] =   { 1.03207, 0.85012, 0.74677, 0.49362, 0.26004, 0.15778, -0.00021 };
	double y12le[] = { 0.02302, 0.01726, 0.01661, 0.01842, 0.01744, 0.01879, 0.02557  };
	double y12he[] = { 0.02173, 0.01752, 0.01616, 0.01640, 0.01940, 0.01953, 0.02278  };
	double y13[] =   { 0.58933, 0.51078, 0.45123, 0.29176, 0.15589, 0.10094, -0.00021 };
	double y13le[] = { 0.03416, 0.03113, 0.02837, 0.02589, 0.02240, 0.02007, 0.02557  };
	double y13he[] = { 0.03493, 0.03018, 0.02599, 0.02481, 0.02353, 0.02152, 0.02278  };

	nrp = 7;

	TGraphAsymmErrors *gr11 = new TGraphAsymmErrors(nrp, x11, y11, xerr, xerr, y11le, y11he);
	gr11->SetMarkerStyle(20);
	gr11->SetMarkerColor(2);
	gr11->SetName("p-sig");
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
	gr12->SetName("He-sig");
	gr12->Draw("SAME;P");

	TGraphAsymmErrors *gr13 = new TGraphAsymmErrors(nrp, x13, y13, xerr, xerr, y13le, y13he);
	gr13->SetMarkerStyle(20);
	gr13->SetMarkerColor(4);
	gr13->SetName("O-sig");
	gr13->Draw("SAME;P");

	cout << endl << "Fitting test_He-Fe.cpp ----------------------------------------------" << endl;
	TF1 *fit1 = new TF1("fit1", "([0]-[1])*((x - TMath::Log(56))/(TMath::Log(1)  - TMath::Log(56))) + [1]", -0.1, 4.5);
	fit1->SetParameter(0, 0.1);
	fit1->SetParameter(1, 0.);
	gr11->Fit("fit1");
	gr11->GetFunction("fit1")->SetLineWidth(1);
	gr11->GetFunction("fit1")->SetLineColor(1);
	TF1 *fit1 = new TF1("fit1", "([0]-[1])*((x - TMath::Log(56))/(TMath::Log(4)  - TMath::Log(56))) + [1]", -0.1, 4.5);
	fit1->SetParameter(0, 0.1);
	fit1->SetParameter(1, 0.);
	gr12->Fit("fit1");
	gr12->GetFunction("fit1")->SetLineWidth(1);
	gr12->GetFunction("fit1")->SetLineColor(1);
	TF1 *fit1 = new TF1("fit1", "([0]-[1])*((x - TMath::Log(56))/(TMath::Log(16) - TMath::Log(56))) + [1]", -0.1, 4.5);
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

	c2->SaveAs("fit-lnA_He-Fe.pdf");
}
