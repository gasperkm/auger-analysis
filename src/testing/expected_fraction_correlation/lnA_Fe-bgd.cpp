{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp = 7;

	// Points are: 0 = 100% sig, 1 = 50% sig + 50% Fe, 2 = 75% sig + 25% Fe, 3 = 25% sig + 75% Fe
	// x1, x11, x21: sig = p
	// x2, x12, x22: sig = He
	// x3, x13, x23: sig = O

	// p/Fe treatment
	double x1[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(56), 0.85*TMath::Log(1)  + 0.15*TMath::Log(56), 0.75*TMath::Log(1)  + 0.25*TMath::Log(56), 0.5*TMath::Log(1)  + 0.5*TMath::Log(56), 0.25*TMath::Log(1)  + 0.75*TMath::Log(56), 0.15*TMath::Log(1)  + 0.85*TMath::Log(56), 0.*TMath::Log(1)  + 1.*TMath::Log(56) };
	double x2[] = { 1.*TMath::Log(4)  + 0.*TMath::Log(56), 0.85*TMath::Log(4)  + 0.15*TMath::Log(56), 0.75*TMath::Log(4)  + 0.25*TMath::Log(56), 0.5*TMath::Log(4)  + 0.5*TMath::Log(56), 0.25*TMath::Log(4)  + 0.75*TMath::Log(56), 0.15*TMath::Log(4)  + 0.85*TMath::Log(56), 0.*TMath::Log(4)  + 1.*TMath::Log(56) };
	double x3[] = { 1.*TMath::Log(16) + 0.*TMath::Log(56), 0.85*TMath::Log(16) + 0.15*TMath::Log(56), 0.75*TMath::Log(16) + 0.25*TMath::Log(56), 0.5*TMath::Log(16) + 0.5*TMath::Log(56), 0.25*TMath::Log(16) + 0.75*TMath::Log(56), 0.15*TMath::Log(16) + 0.85*TMath::Log(56), 0.*TMath::Log(16) + 1.*TMath::Log(56) };
	double xerr[] = { 0., 0., 0., 0. , 0., 0., 0. };

	double y1[] =   { 1.02785, 0.85003, 0.75360, 0.50140, 0.25431, 0.15725, -0.00639 };
	double y1le[] = { 0.01372, 0.01097, 0.01147, 0.00957, 0.01062, 0.01065, 0.01321  };
	double y1he[] = { 0.01267, 0.01331, 0.01102, 0.00892, 0.01169, 0.01251, 0.01307  };
	double y2[] =   { 0.91522, 0.77423, 0.67983, 0.45472, 0.23379, 0.13996, -0.00639 };
	double y2le[] = { 0.01642, 0.01373, 0.01299, 0.00957, 0.01100, 0.01098, 0.01321  };
	double y2he[] = { 0.01416, 0.01305, 0.01192, 0.01006, 0.01209, 0.01287, 0.01307  };
	double y3[] =   { 0.49253, 0.41458, 0.36582, 0.24379, 0.12598, 0.07859, -0.00639 };
	double y3le[] = { 0.02404, 0.02191, 0.02137, 0.01416, 0.01501, 0.01314, 0.01321  };
	double y3he[] = { 0.02343, 0.02283, 0.02001, 0.01806, 0.01500, 0.01420, 0.01307  };

	TGraphAsymmErrors *gr1 = new TGraphAsymmErrors(nrp, x1, y1, xerr, xerr, y1le, y1he);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(2);
	gr1->SetName("p-sig_p-Fe");
	gr1->GetXaxis()->SetRangeUser(-0.1,4.5);
	gr1->GetYaxis()->SetRangeUser(-0.1,1.6);
	gr1->GetXaxis()->SetLimits(-0.1,4.5);
	gr1->GetYaxis()->SetLimits(-0.1,1.6);
	gr1->GetXaxis()->SetTitle("lnA");
	gr1->GetYaxis()->SetTitle("MVA fraction");
	gr1->Draw("AP");

	TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(nrp, x2, y2, xerr, xerr, y2le, y2he);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(3);
	gr2->SetName("He-sig_p-Fe");
	gr2->Draw("SAME;P");

	TGraphAsymmErrors *gr3 = new TGraphAsymmErrors(nrp, x3, y3, xerr, xerr, y3le, y3he);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(4);
	gr3->SetName("O-sig_p-Fe");
	gr3->Draw("SAME;P");

	// He/Fe treatment
	double x11[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(56), 0.85*TMath::Log(1)  + 0.15*TMath::Log(56), 0.75*TMath::Log(1)  + 0.25*TMath::Log(56), 0.5*TMath::Log(1)  + 0.5*TMath::Log(56), 0.25*TMath::Log(1)  + 0.75*TMath::Log(56), 0.15*TMath::Log(1)  + 0.85*TMath::Log(56), 0.*TMath::Log(1)  + 1.*TMath::Log(56) };
	double x12[] = { 1.*TMath::Log(4)  + 0.*TMath::Log(56), 0.85*TMath::Log(4)  + 0.15*TMath::Log(56), 0.75*TMath::Log(4)  + 0.25*TMath::Log(56), 0.5*TMath::Log(4)  + 0.5*TMath::Log(56), 0.25*TMath::Log(4)  + 0.75*TMath::Log(56), 0.15*TMath::Log(4)  + 0.85*TMath::Log(56), 0.*TMath::Log(4)  + 1.*TMath::Log(56) };
	double x13[] = { 1.*TMath::Log(16) + 0.*TMath::Log(56), 0.85*TMath::Log(16) + 0.15*TMath::Log(56), 0.75*TMath::Log(16) + 0.25*TMath::Log(56), 0.5*TMath::Log(16) + 0.5*TMath::Log(56), 0.25*TMath::Log(16) + 0.75*TMath::Log(56), 0.15*TMath::Log(16) + 0.85*TMath::Log(56), 0.*TMath::Log(16) + 1.*TMath::Log(56) };

	double y11[] =   { 1.08159, 0.91210, 0.80524, 0.53687, 0.27435, 0.17189, -0.00021 };
	double y11le[] = { 0.02015, 0.01635, 0.01507, 0.01727, 0.01687, 0.01707, 0.02557  };
	double y11he[] = { 0.01645, 0.01431, 0.01397, 0.01385, 0.01803, 0.01817, 0.02278  };
	double y12[] =   { 1.03207, 0.85012, 0.74677, 0.49362, 0.26004, 0.15778, -0.00021 };
	double y12le[] = { 0.02302, 0.01726, 0.01661, 0.01842, 0.01744, 0.01879, 0.02557  };
	double y12he[] = { 0.02173, 0.01752, 0.01616, 0.01640, 0.01940, 0.01953, 0.02278  };
	double y13[] =   { 0.58933, 0.51078, 0.45123, 0.29176, 0.15589, 0.10094, -0.00021 };
	double y13le[] = { 0.03416, 0.03113, 0.02837, 0.02589, 0.02240, 0.02007, 0.02557  };
	double y13he[] = { 0.03493, 0.03018, 0.02599, 0.02481, 0.02353, 0.02152, 0.02278  };

	TGraphAsymmErrors *gr11 = new TGraphAsymmErrors(nrp, x11, y11, xerr, xerr, y11le, y11he);
	gr11->SetMarkerStyle(21);
	gr11->SetMarkerColor(2);
	gr11->SetMarkerSize(0.9);
	gr11->SetName("p-sig_He-Fe");
	gr11->Draw("SAME;P");

	TGraphAsymmErrors *gr12 = new TGraphAsymmErrors(nrp, x12, y12, xerr, xerr, y12le, y12he);
	gr12->SetMarkerStyle(21);
	gr12->SetMarkerColor(3);
	gr12->SetMarkerSize(0.9);
	gr12->SetName("He-sig_He-Fe");
	gr12->Draw("SAME;P");

	TGraphAsymmErrors *gr13 = new TGraphAsymmErrors(nrp, x13, y13, xerr, xerr, y13le, y13he);
	gr13->SetMarkerStyle(21);
	gr13->SetMarkerColor(4);
	gr13->SetMarkerSize(0.9);
	gr13->SetName("O-sig_He-Fe");
	gr13->Draw("SAME;P");

	// O/Fe treatment
	double x21[] = { 1.*TMath::Log(1)  + 0.*TMath::Log(56), 0.85*TMath::Log(1)  + 0.15*TMath::Log(56), 0.75*TMath::Log(1)  + 0.25*TMath::Log(56), 0.5*TMath::Log(1)  + 0.5*TMath::Log(56), 0.25*TMath::Log(1)  + 0.75*TMath::Log(56), 0.15*TMath::Log(1)  + 0.85*TMath::Log(56), 0.*TMath::Log(1)  + 1.*TMath::Log(56) };
	double x22[] = { 1.*TMath::Log(4)  + 0.*TMath::Log(56), 0.85*TMath::Log(4)  + 0.15*TMath::Log(56), 0.75*TMath::Log(4)  + 0.25*TMath::Log(56), 0.5*TMath::Log(4)  + 0.5*TMath::Log(56), 0.25*TMath::Log(4)  + 0.75*TMath::Log(56), 0.15*TMath::Log(4)  + 0.85*TMath::Log(56), 0.*TMath::Log(4)  + 1.*TMath::Log(56) };
	double x23[] = { 1.*TMath::Log(16) + 0.*TMath::Log(56), 0.85*TMath::Log(16) + 0.15*TMath::Log(56), 0.75*TMath::Log(16) + 0.25*TMath::Log(56), 0.5*TMath::Log(16) + 0.5*TMath::Log(56), 0.25*TMath::Log(16) + 0.75*TMath::Log(56), 0.15*TMath::Log(16) + 0.85*TMath::Log(56), 0.*TMath::Log(16) + 1.*TMath::Log(56) };

	double y21[] =   { 1.46308, 1.21762, 1.07554, 0.71115, 0.37499, 0.23156, -0.01271 };
	double y21le[] = { 0.13022, 0.08077, 0.07354, 0.06823, 0.06788, 0.07495, 0.11759  };
	double y21he[] = { 0.12181, 0.07537, 0.06796, 0.06735, 0.07327, 0.07912, 0.12145  };
	double y22[] =   { 1.39969, 1.16361, 1.02974, 0.67875, 0.36019, 0.22234, -0.01271 };
	double y22le[] = { 0.12279, 0.07694, 0.07031, 0.06990, 0.06745, 0.07576, 0.11759  };
	double y22he[] = { 0.11333, 0.07374, 0.06730, 0.06875, 0.07440, 0.08046, 0.12145  };
	double y23[] =   { 1.06299, 0.84556, 0.74695, 0.50417, 0.26296, 0.16603, -0.01271 };
	double y23le[] = { 0.11278, 0.07845, 0.07469, 0.08941, 0.07622, 0.08220, 0.11759  };
	double y23he[] = { 0.11121, 0.07774, 0.07160, 0.08353, 0.08563, 0.08573, 0.12145  };

	TGraphAsymmErrors *gr21 = new TGraphAsymmErrors(nrp, x21, y21, xerr, xerr, y21le, y21he);
	gr21->SetMarkerStyle(22);
	gr21->SetMarkerColor(2);
	gr21->SetName("p-sig_O-Fe");
	gr21->Draw("SAME;P");

	TGraphAsymmErrors *gr22 = new TGraphAsymmErrors(nrp, x22, y22, xerr, xerr, y22le, y22he);
	gr22->SetMarkerStyle(22);
	gr22->SetMarkerColor(3);
	gr22->SetName("He-sig_O-Fe");
	gr22->Draw("SAME;P");

	TGraphAsymmErrors *gr23 = new TGraphAsymmErrors(nrp, x23, y23, xerr, xerr, y23le, y23he);
	gr23->SetMarkerStyle(22);
	gr23->SetMarkerColor(4);
	gr23->SetName("O-sig_O-Fe");
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

/*	for(int i = 1; i <= 56; i++)
	{
           TLine *l1 = new TLine(TMath::Log(i), -0.1, TMath::Log(i), 1.6);
           l1->SetLineWidth(1);
           l1->SetLineColor(6);
           l1->SetLineStyle(9);
           l1->Draw("same");
	}*/

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

	c1->SaveAs("fit-lnA_Fe-bgd.pdf");
}
