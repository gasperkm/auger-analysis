// From results/analysis_20180522
//   Performing MVA analysis with pure mock sets, Fe as background
void mock_datasets_nonorm()
{
	TCanvas *c1 = new TCanvas("c1","Mock dataset results (no normalization)", 800, 600);
	c1->SetGrid();

	// Proton results (p/Fe, He/Fe, O/Fe)
	double yp[] =     {0.886257  , 0.911003  , 0.958799  };
	double yperrl[] = {0.00723943, 0.00644291, 0.00689994};
	double yperrh[] = {0.00729527, 0.00732277, 0.00465232};
	// Helium results (p/Fe, He/Fe, O/Fe)
	double yhe[] =     {0.807995  , 0.84435  , 0.921707  };
	double yheerrl[] = {0.0088853 , 0.0103276, 0.00927622};
	double yheerrh[] = {0.00852614, 0.0124228, 0.0092512 };
	// Oxygen results (p/Fe, He/Fe, O/Fe)
	double yo[] =     {0.479961 , 0.551149 , 0.723432 };
	double yoerrl[] = {0.0158467, 0.0196834, 0.0244409};
	double yoerrh[] = {0.015955 , 0.0192758, 0.0230105};
	// lnA values of p, He and O
	double x[] =    {TMath::Log(1), TMath::Log(4), TMath::Log(16)};
	double xerr[] = {0., 0., 0.};

	// All graphs
	TGraphAsymmErrors *grp = new TGraphAsymmErrors(3, x, yp, xerr, xerr, yperrl, yperrh);
	TGraphAsymmErrors *grhe = new TGraphAsymmErrors(3, x, yhe, xerr, xerr, yheerrl, yheerrh);
	TGraphAsymmErrors *gro = new TGraphAsymmErrors(3, x, yo, xerr, xerr, yoerrl, yoerrh);
	grp->SetMarkerStyle(20);
	grp->SetMarkerColor(1);
	grp->SetLineColor(1);
	grp->SetLineWidth(2);
	grhe->SetMarkerStyle(21);
	grhe->SetMarkerColor(2);
	grhe->SetLineColor(2);
	grhe->SetLineWidth(2);
	gro->SetMarkerStyle(22);
	gro->SetMarkerColor(4);
	gro->SetLineColor(4);
	gro->SetLineWidth(2);

	grp->GetXaxis()->SetRangeUser(-0.1, 3.1);
	grp->GetYaxis()->SetRangeUser(0.0,1.6);
	grp->Draw("AP");
	grhe->Draw("P;SAME");
	gro->Draw("P;SAME");

	TCanvas *c2 = new TCanvas("c2","Shifted mock dataset results (no normalization)", 800, 600);
	c2->SetGrid();


	// Shift all proton values to 1
	double shift;

	for(int i = 0; i < 3; i++)
	{
		shift = 1 - yp[i];
		cout << i << " Shift value: " << shift << endl;

		yp[i] += shift;
		yhe[i] += shift;
		yo[i] += shift;
	}

	// Shifted graphs
	TGraphAsymmErrors *grshp = new TGraphAsymmErrors(3, x, yp, xerr, xerr, yperrl, yperrh);
	TGraphAsymmErrors *grshhe = new TGraphAsymmErrors(3, x, yhe, xerr, xerr, yheerrl, yheerrh);
	TGraphAsymmErrors *grsho = new TGraphAsymmErrors(3, x, yo, xerr, xerr, yoerrl, yoerrh);
	grshp->SetMarkerStyle(20);
	grshp->SetMarkerColor(1);
	grshp->SetLineColor(1);
	grshp->SetLineWidth(2);
	grshhe->SetMarkerStyle(21);
	grshhe->SetMarkerColor(2);
	grshhe->SetLineColor(2);
	grshhe->SetLineWidth(2);
	grsho->SetMarkerStyle(22);
	grsho->SetMarkerColor(4);
	grsho->SetLineColor(4);
	grsho->SetLineWidth(2);

	grshp->GetXaxis()->SetRangeUser(-0.1, 3.1);
	grshp->GetYaxis()->SetRangeUser(0.0,1.1);
	grshp->Draw("AP");
	grshhe->Draw("P;SAME");
	grsho->Draw("P;SAME");
}
