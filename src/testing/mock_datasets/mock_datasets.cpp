// From results/analysis_20180522
//   Performing MVA analysis with pure mock sets, Fe as background
void mock_datasets()
{
	TCanvas *c1 = new TCanvas("c1","Mock dataset results", 800, 600);
	c1->SetGrid();

	// Proton results (p/Fe, He/Fe, O/Fe)
	double yp[] =     {0.996595 , 1.08183  , 1.48698  };
	double yperrl[] = {0.0130088, 0.0180084, 0.091135 };
	double yperrh[] = {0.0118478, 0.0190923, 0.0787212};
	// Helium results (p/Fe, He/Fe, O/Fe)
	double yhe[] =     {0.897526 , 0.989121 , 1.4074   };
	double yheerrl[] = {0.0140774, 0.0204516, 0.0874252};
	double yheerrh[] = {0.0128294, 0.0231252, 0.0770919};
	// Oxygen results (p/Fe, He/Fe, O/Fe)
	double yo[] =     {0.482553 , 0.581554 , 0.984212 };
	double yoerrl[] = {0.0214123, 0.0298651, 0.0763944};
	double yoerrh[] = {0.0212432, 0.0293038, 0.068967 };
	// lnA values of p, He and O
	double x[] =    {TMath::Log(1), TMath::Log(4), TMath::Log(16)};
	double xerr[] = {0., 0., 0.};

	// All graphs
	TGraphAsymmErrors *grp = new TGraphAsymmErrors(3, x, yp, xerr, xerr, yperrl, yperrh);
	grp->SetName("grp");
	TGraphAsymmErrors *grhe = new TGraphAsymmErrors(3, x, yhe, xerr, xerr, yheerrl, yheerrh);
	grhe->SetName("grhe");
	TGraphAsymmErrors *gro = new TGraphAsymmErrors(3, x, yo, xerr, xerr, yoerrl, yoerrh);
	gro->SetName("gro");

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

	TCanvas *c2 = new TCanvas("c2","Shifted mock dataset results", 800, 600);
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
	grshp->SetName("grshp");
	TGraphAsymmErrors *grshhe = new TGraphAsymmErrors(3, x, yhe, xerr, xerr, yheerrl, yheerrh);
	grshhe->SetName("grshhe");
	TGraphAsymmErrors *grsho = new TGraphAsymmErrors(3, x, yo, xerr, xerr, yoerrl, yoerrh);
	grsho->SetName("grsho");

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

	// Comparison of MVA fractions and expected fractions
	TCanvas *c3 = new TCanvas("c3","Comparison to expected fractions", 800, 600);
	c3->SetGrid();

	// Expected fractions
	double expfrac[] = {1-TMath::Log(1)/TMath::Log(56), 1-TMath::Log(4)/TMath::Log(56), 1-TMath::Log(16)/TMath::Log(56), 1-TMath::Log(56)/TMath::Log(56)};
	double expfracerr[] = {0., 0., 0., 0.};
	double pFe[4], HeFe[4], OFe[4];
	double pFeerrl[4], HeFeerrl[4], OFeerrl[4];
	double pFeerrh[4], HeFeerrh[4], OFeerrh[4];

	int cnt = 0;
	for(int i = 0; i < 4; i++)
	{
		if(i == 3)
		{
			pFe[i]      = 0.00228328;
			pFeerrl[i]  = 0.0135256 ;
			pFeerrh[i]  = 0.0129225 ;
			HeFe[i]     = 0.00181489;
			HeFeerrl[i] = 0.0250775 ;
			HeFeerrh[i] = 0.0244534 ;
			OFe[i]      = 0.0197649 ;
			OFeerrl[i]  = 0.0810437 ;
			OFeerrh[i]  = 0.0789351 ;
		}
		else if(i == 0)
		{
			pFe[i]      = yp[0];
			pFeerrl[i]  = yperrl[0];
			pFeerrh[i]  = yperrh[0];
			HeFe[i]     = yp[1];
			HeFeerrl[i] = yperrl[1];
			HeFeerrh[i] = yperrh[1];
			OFe[i]      = yp[2];
			OFeerrl[i]  = yperrl[2];
			OFeerrh[i]  = yperrh[2];
		}
		else if(i == 1)
		{
			pFe[i]      = yhe[0];
			pFeerrl[i]  = yheerrl[0];
			pFeerrh[i]  = yheerrh[0];
			HeFe[i]     = yhe[1];
			HeFeerrl[i] = yheerrl[1];
			HeFeerrh[i] = yheerrh[1];
			OFe[i]      = yhe[2];
			OFeerrl[i]  = yheerrl[2];
			OFeerrh[i]  = yheerrh[2];
		}
		else if(i == 2)
		{
			pFe[i]      = yo[0];
			pFeerrl[i]  = yoerrl[0];
			pFeerrh[i]  = yoerrh[0];
			HeFe[i]     = yo[1];
			HeFeerrl[i] = yoerrl[1];
			HeFeerrh[i] = yoerrh[1];
			OFe[i]      = yo[2];
			OFeerrl[i]  = yoerrl[2];
			OFeerrh[i]  = yoerrh[2];
		}
	}

	TGraphAsymmErrors *grcompp = new TGraphAsymmErrors(4, expfrac, pFe, expfracerr, expfracerr, pFeerrl, pFeerrh);
	grcompp->SetName("grcompp");
	TGraphAsymmErrors *grcomphe = new TGraphAsymmErrors(4, expfrac, HeFe, expfracerr, expfracerr, HeFeerrl, HeFeerrh);
	grcomphe->SetName("grcomphe");
	TGraphAsymmErrors *grcompo = new TGraphAsymmErrors(4, expfrac, OFe, expfracerr, expfracerr, OFeerrl, OFeerrh);
	grcompo->SetName("grcompo");

	grcompp->SetMarkerStyle(20);
	grcompp->SetMarkerColor(1);
	grcompp->SetLineColor(1);
	grcompp->SetLineWidth(2);
	grcomphe->SetMarkerStyle(21);
	grcomphe->SetMarkerColor(2);
	grcomphe->SetLineColor(2);
	grcomphe->SetLineWidth(2);
	grcompo->SetMarkerStyle(22);
	grcompo->SetMarkerColor(4);
	grcompo->SetLineColor(4);
	grcompo->SetLineWidth(2);

	grcompp->GetXaxis()->SetRangeUser(-0.1, 1.1);
	grcompp->GetYaxis()->SetRangeUser(-0.1, 1.1);
	grcompp->Draw("AP");
	grcomphe->Draw("P;SAME");
	grcompo->Draw("P;SAME");
}
