// Read published lnA results to add them to the plot (type: 0 = EPOS, 1 = QGSJETII, 2 = SIBYLL)
int ReadLnaResults(vector<float> *val, int type)
{
   ifstream infile;
   char ctemp[1024];
   float *ftemp;
   int nrp = 0;

   string stemp;
   if(type == 0)
      stemp = "lnA_moments_epos.txt";
   else if(type == 1)
      stemp = "lnA_moments_qgs.txt";
   else if(type == 2)
      stemp = "lnA_moments_sib.txt";
   else
      return -1;

   ftemp = new float[10];

   infile.open(stemp.c_str(), ifstream::in);

   if(infile.is_open())
   {
      infile.getline(ctemp, 1024, '\n');
      cout << string(ctemp) << endl;
      
      while(1)
      {
	 for(int i = 0; i < 10; i++)
            infile >> ftemp[i];

         infile.ignore(1, '\n');
         
         if(infile.eof())
         	break;

//	 if(ftemp[0] > 18.40)
//	 {
            val->push_back(ftemp[0]);
            val->push_back(ftemp[1]);
            val->push_back(ftemp[2]);
            cout << "Point " << nrp << " = " << ftemp[0] << ", " << ftemp[2] << " (" << ftemp[3] << ")" << endl;

	    nrp++;
//	 }
      }
   }

   infile.close();

   delete[] ftemp;

   return nrp;
}

void SetStyle()
{
   TStyle *basestyle = new TStyle("basestyle", "basestyle");
   // Set title and label font sizes (with precision 3, these are given in pixels)
   basestyle->SetTextFont(63);
   basestyle->SetTextSize(18);
   basestyle->SetLabelFont(63, "xyz");
   basestyle->SetLabelSize(24, "xyz");
   basestyle->SetTitleFont(63, "xyz");
   basestyle->SetTitleSize(24, "xyz");
   
   // Set option and statistics
   basestyle->SetOptStat(0);
   basestyle->SetPalette(1,0);
   basestyle->SetOptTitle(0);
   basestyle->SetStatFontSize(0.024);
   basestyle->SetStatBorderSize(1);
   basestyle->SetStatColor(kGray);
   basestyle->SetStatX(0.925);
   basestyle->SetStatY(0.925);
   basestyle->SetStatW(0.13);

   // Set canvas and pads
   basestyle->SetCanvasBorderMode(0);
   basestyle->SetFrameBorderMode(0);
   basestyle->SetCanvasColor(0);
   basestyle->SetPadTickX(1);
   basestyle->SetPadTickY(1);
   basestyle->SetCanvasDefX(100);
   basestyle->SetCanvasDefY(50);
   basestyle->SetCanvasDefW(900);
   basestyle->SetCanvasDefH(600);
   basestyle->SetPadBorderMode(0);
   basestyle->SetPadBottomMargin(0.1);
   basestyle->SetPadTopMargin(0.04);
   basestyle->SetPadLeftMargin(0.125);
   basestyle->SetPadRightMargin(0.04);
   basestyle->SetPadColor(0);
   basestyle->SetPadGridX(kTRUE);
   basestyle->SetPadGridY(kTRUE);

   // Label and title offsets
   basestyle->SetLabelOffset(0.015,"xyz");
   basestyle->SetTitleOffset(1.6, "x");
   basestyle->SetTitleOffset(1.9, "y");
   basestyle->SetEndErrorSize(8);

   gROOT->SetStyle("basestyle");
   gROOT->ForceStyle(1);
}

void SetAxisTitles(TGraph *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void SetAxisTitles(TGraphErrors *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void fraction_conv()
{
	int linestyle = 9;

	// Data X_max values from https://www.auger.unam.mx/AugerWiki/XmaxHeatIcrc2017/?action=AttachFile&do=get&target=Xmax_moments_icrc17_.txt
	vector<double> xbins;
	vector<double> dbins;
	vector<double> dbinsErr;
	double dtemp[10];
	char ctemp[1024];
	int itemp[3];

	int nbins = 0;

	ifstream ifile;
	ifile.open("Xmax_moments_icrc17_.txt", ifstream::in);
	if(ifile.is_open())
	{
		ifile.getline(ctemp, 1024);
		cout << string(ctemp) << endl;

		while(1)
		{
			for(int i = 0; i < 10; i++)
				ifile >> dtemp[i];

			ifile.ignore(1, '\n');

			if(ifile.eof())
				break;

			xbins.push_back(dtemp[0]);
			dbins.push_back(dtemp[2]);
			dbinsErr.push_back(dtemp[3]);
			cout << "Point " << nbins << " = " << dtemp[0] << ", " << dtemp[2] << " (" << dtemp[3] << ")" << endl;
			nbins++;
		}
	}
	else
	{
		cout << "File not found!" << endl;
		return;
	}

	ifile.close();

	// Simulations estimated from ICRC 2017 X_max plot)
	double xlimit[] = {xbins[0], xbins[xbins.size()-1]};
	double plimitEpos[] = {705., 841.};
	double felimitEpos[] = {600., 746.};
	double plimitQgs[] = {693., 821.};
	double felimitQgs[] = {590., 728.};

	// Calculating straight line: y = K*x + C
	double pKEpos = (plimitEpos[1] - plimitEpos[0])/(xlimit[1] - xlimit[0]);
	double feKEpos = (felimitEpos[1] - felimitEpos[0])/(xlimit[1] - xlimit[0]);
	double pKQgs = (plimitQgs[1] - plimitQgs[0])/(xlimit[1] - xlimit[0]);
	double feKQgs = (felimitQgs[1] - felimitQgs[0])/(xlimit[1] - xlimit[0]);

	double pCEpos = plimitEpos[0] - xlimit[0]*pKEpos;
	double feCEpos = felimitEpos[0] - xlimit[0]*feKEpos;
	double pCQgs = plimitQgs[0] - xlimit[0]*pKQgs;
	double feCQgs = felimitQgs[0] - xlimit[0]*feKQgs;

	cout << "pKEpos = " << pKEpos << ", pCEpos = " << pCEpos << endl;
	cout << "feKEpos = " << feKEpos << ", feCEpos = " << feCEpos << endl;
	cout << "pKQgs = " << pKQgs << ", pCQgs = " << pCQgs << endl;
	cout << "feKQgs = " << feKQgs << ", feCQgs = " << feCQgs << endl;

	vector<double> pbinsEpos;
	vector<double> febinsEpos;
	vector<double> pbinsQgs;
	vector<double> febinsQgs;

	SetStyle();

	TLegend *legend;
	int c_Legend = TColor::GetColor("#ffff66");
	int nrfigs;
	TCanvas *c1 = new TCanvas("c1","",1200,900);
	c1->SetGrid();
	TGraph *pgrEpos = new TGraph(nbins);
	TGraph *fegrEpos = new TGraph(nbins);
	TGraph *pgrQgs = new TGraph(nbins);
	TGraph *fegrQgs = new TGraph(nbins);
	TGraphErrors *dgr = new TGraphErrors(nbins);

	cout << "Xmax values:" << endl;
	for(int i = 0; i < nbins; i++)
	{
		pbinsEpos.push_back(pKEpos*xbins[i] + pCEpos);
		febinsEpos.push_back(feKEpos*xbins[i] + feCEpos);
		cout << "Epos values " << i << ": E = " << xbins[i] << ",\tp = " << pbinsEpos[i] << ",\tFe = " << febinsEpos[i] << ",\tdata = " << dbins[i] << " (" << dbinsErr[i] << ")" << endl;

		pbinsQgs.push_back(pKQgs*xbins[i] + pCQgs);
		febinsQgs.push_back(feKQgs*xbins[i] + feCQgs);
		cout << "Qgs values " << i << ": E = " << xbins[i] << ",\tp = " << pbinsQgs[i] << ",\tFe = " << febinsQgs[i] << ",\tdata = " << dbins[i] << " (" << dbinsErr[i] << ")" << endl;

		pgrEpos->SetPoint(i, xbins[i], pbinsEpos[i]);
		fegrEpos->SetPoint(i, xbins[i], febinsEpos[i]);

		pgrQgs->SetPoint(i, xbins[i], pbinsQgs[i]);
		fegrQgs->SetPoint(i, xbins[i], febinsQgs[i]);

		dgr->SetPoint(i, xbins[i], dbins[i]);
		dgr->SetPointError(i, 0, dbinsErr[i]);
	}

	SetAxisTitles(pgrEpos, "log(E/eV)", "<X_{max}> (g/cm^{2})");
//	pgrEpos->GetXaxis()->SetTitle("log(E/eV)");
//	pgrEpos->GetYaxis()->SetTitle("<X_{max}> (g/cm^{2})");
	pgrEpos->GetYaxis()->SetRangeUser(felimitQgs[0], plimitEpos[1]);
	pgrEpos->SetLineColor(2);
	pgrEpos->SetLineWidth(2);
	pgrEpos->Draw("AL");

	fegrEpos->SetLineColor(4);
	fegrEpos->SetLineWidth(2);
	fegrEpos->Draw("L;SAME");

	pgrQgs->SetLineColor(2);
	pgrQgs->SetLineWidth(2);
	pgrQgs->SetLineStyle(linestyle);
	pgrQgs->Draw("L;SAME");

	fegrQgs->SetLineColor(4);
	fegrQgs->SetLineWidth(2);
	fegrQgs->SetLineStyle(linestyle);
	fegrQgs->Draw("L;SAME");

	dgr->SetMarkerColor(1);
	dgr->SetMarkerStyle(22);
	dgr->SetLineColor(1);
	dgr->SetLineWidth(2);
	dgr->Draw("PL;SAME");

	nrfigs = 5;
	legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	legend->SetFillStyle(1001);
	legend->SetFillColor(c_Legend);
	legend->AddEntry(pgrEpos, "EPOS protons", "l");
	legend->AddEntry(fegrEpos, "EPOS iron", "l");
	legend->AddEntry(pgrQgs, "QGSJET-II protons", "l");
	legend->AddEntry(fegrQgs, "QGSJET-II iron", "l");
	legend->AddEntry(dgr, "Auger data", "pl");
	legend->SetBorderSize(1);
	legend->SetMargin(0.3);
	legend->Draw("same");

	c1->SaveAs("xmax_plot.pdf");

	pgrEpos->GetXaxis()->SetRangeUser(18.5,19.7);
	pgrEpos->GetYaxis()->SetRangeUser(650,plimitEpos[1]);
	c1->SaveAs("xmax_plot_zoom.pdf");

	TCanvas *c2 = new TCanvas("c2","",1200,900);
	c2->SetGrid();
	TGraphErrors *grfracEpos = new TGraphErrors(nbins);
	TGraphErrors *grfracQgs = new TGraphErrors(nbins);

	cout << "Fraction values:" << endl;
	for(int i = 0; i < nbins; i++)
	{
		dtemp[0] = (dbins[i] - febinsEpos[i])/(pbinsEpos[i] - febinsEpos[i]);
		dtemp[1] = ((dbins[i]+dbinsErr[i]) - febinsEpos[i])/(pbinsEpos[i] - febinsEpos[i]) - dtemp[0];
		cout << "Epos values: E = " << xbins[i] << ",\tfraction = " << dtemp[0] << " (" << dtemp[1] << ")" << endl;
		grfracEpos->SetPoint(i, xbins[i], dtemp[0]);
		grfracEpos->SetPointError(i, 0, dtemp[1]);

		dtemp[0] = (dbins[i] - febinsQgs[i])/(pbinsQgs[i] - febinsQgs[i]);
		dtemp[1] = ((dbins[i]+dbinsErr[i]) - febinsQgs[i])/(pbinsQgs[i] - febinsQgs[i]) - dtemp[0];
		cout << "Qgs values: E = " << xbins[i] << ",\tfraction = " << dtemp[0] << " (" << dtemp[1] << ")" << endl;

		grfracQgs->SetPoint(i, xbins[i], dtemp[0]);
		grfracQgs->SetPointError(i, 0, dtemp[1]);
	}

	SetAxisTitles(grfracEpos, "log(E/eV)", "Proton fraction");
//	grfracEpos->GetXaxis()->SetTitle("log(E/eV)");
//	grfracEpos->GetYaxis()->SetTitle("Proton fraction");
	grfracEpos->GetYaxis()->SetRangeUser(-0.1,1+0.25*1.1);
	grfracEpos->SetMarkerColor(1);
	grfracEpos->SetMarkerStyle(22);
	grfracEpos->SetLineColor(1);
	grfracEpos->SetLineWidth(2);
	grfracEpos->Draw("APL");

	grfracQgs->SetMarkerColor(1);
	grfracQgs->SetMarkerStyle(22);
	grfracQgs->SetLineColor(1);
	grfracQgs->SetLineWidth(2);
	grfracQgs->SetLineStyle(linestyle);
	grfracQgs->Draw("PL;SAME");

	nrfigs = 2;
	legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	legend->SetFillStyle(1001);
	legend->SetFillColor(c_Legend);
	legend->AddEntry(grfracEpos, "Auger data fraction (EPOS)", "pl");
	legend->AddEntry(grfracQgs, "Auger data fraction (QGSJET-II)", "pl");
	legend->SetBorderSize(1);
	legend->SetMargin(0.3);
	legend->Draw("same");

	c2->SaveAs("xmax_to_fraction_conversion.pdf");

	grfracEpos->GetXaxis()->SetRangeUser(18.5,19.7);
	c2->SaveAs("xmax_to_fraction_conversion_zoom.pdf");

	// Published lnA data values
	vector<float> returnVal;
	itemp[0] = ReadLnaResults(&returnVal, 0);  // EPOS
	cout << "Number of points = " << itemp[0] << endl;
	TGraphErrors *publnaEpos = new TGraphErrors(itemp[0]);

	for(int i = 0; i < itemp[0]; i++)
	{
		publnaEpos->SetPoint(i, returnVal[3*i], returnVal[3*i+1]);
		publnaEpos->SetPointError(i, 0, returnVal[3*i+2]);
	}

	returnVal.clear();
	itemp[0] = ReadLnaResults(&returnVal, 1);  // QGSJET
	cout << "Number of points = " << itemp[0] << endl;
	TGraphErrors *publnaQgs = new TGraphErrors(itemp[0]);

	for(int i = 0; i < itemp[0]; i++)
	{
		publnaQgs->SetPoint(i, returnVal[3*i], returnVal[3*i+1]);
		publnaQgs->SetPointError(i, 0, returnVal[3*i+2]);
	}

	// Plotting lnA stuff
	TCanvas *c3 = new TCanvas("c3","",1200,900);
	c3->SetGrid();
	TGraphErrors *grlnaEpos = new TGraphErrors(nbins);
	TGraphErrors *grlnaQgs = new TGraphErrors(nbins);

	cout << "lnA values:" << endl;
	for(int i = 0; i < nbins; i++)
	{
		dtemp[0] = (dbins[i] - febinsEpos[i])/(pbinsEpos[i] - febinsEpos[i]);
		dtemp[1] = ((dbins[i]+dbinsErr[i]) - febinsEpos[i])/(pbinsEpos[i] - febinsEpos[i]) - dtemp[0];

		dtemp[2] = (1.-dtemp[0])*TMath::Log(56);
		dtemp[3] = (1.-(dtemp[0]+dtemp[1]))*TMath::Log(56);
		dtemp[3] = TMath::Abs(dtemp[2] - dtemp[3]);

		cout << "Epos values: E = " << xbins[i] << ",\tlnA = " << dtemp[2] << " (" << dtemp[3] << ")" << endl;
		grlnaEpos->SetPoint(i, xbins[i], dtemp[2]);
		grlnaEpos->SetPointError(i, 0, dtemp[3]);

		dtemp[0] = (dbins[i] - febinsQgs[i])/(pbinsQgs[i] - febinsQgs[i]);
		dtemp[1] = ((dbins[i]+dbinsErr[i]) - febinsQgs[i])/(pbinsQgs[i] - febinsQgs[i]) - dtemp[0];

		dtemp[2] = (1.-dtemp[0])*TMath::Log(56);
		dtemp[3] = (1.-(dtemp[0]+dtemp[1]))*TMath::Log(56);
		dtemp[3] = TMath::Abs(dtemp[2] - dtemp[3]);

		cout << "Qgs values: E = " << xbins[i] << ",\tlnA = " << dtemp[2] << " (" << dtemp[3] << ")" << endl;
		grlnaQgs->SetPoint(i, xbins[i], dtemp[2]);
		grlnaQgs->SetPointError(i, 0, dtemp[3]);
	}

	SetAxisTitles(grlnaEpos, "log(E/eV)", "<lnA> of data events");
//	grlnaEpos->GetXaxis()->SetTitle("log(E/eV)");
//	grlnaEpos->GetYaxis()->SetTitle("<lnA> of data events");
	grlnaEpos->GetYaxis()->SetRangeUser(-0.7, 5.);
	grlnaEpos->SetMarkerColor(1);
	grlnaEpos->SetMarkerStyle(22);
	grlnaEpos->SetLineColor(1);
	grlnaEpos->SetLineWidth(2);
	grlnaEpos->Draw("APL");

	grlnaQgs->SetMarkerColor(1);
	grlnaQgs->SetMarkerStyle(22);
	grlnaQgs->SetLineColor(1);
	grlnaQgs->SetLineWidth(2);
	grlnaQgs->SetLineStyle(linestyle);
	grlnaQgs->Draw("PL;SAME");

	publnaEpos->SetMarkerColor(2);
	publnaEpos->SetMarkerStyle(20);
	publnaEpos->SetLineColor(2);
	publnaEpos->SetLineWidth(2);
	publnaEpos->Draw("PL;SAME");

	publnaQgs->SetMarkerColor(2);
	publnaQgs->SetMarkerStyle(20);
	publnaQgs->SetLineColor(2);
	publnaQgs->SetLineWidth(2);
	publnaQgs->SetLineStyle(linestyle);
	publnaQgs->Draw("PL;SAME");

	nrfigs = 4;
	legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	legend->SetFillStyle(1001);
	legend->SetFillColor(c_Legend);
	legend->AddEntry(publnaEpos, "Auger published data (EPOS)", "pl");
	legend->AddEntry(publnaQgs, "Auger published data (QGSJET-II)", "pl");
	legend->AddEntry(grlnaEpos, "Auger data (EPOS)", "pl");
	legend->AddEntry(grlnaQgs, "Auger data (QGSJET-II)", "pl");
	legend->SetBorderSize(1);
	legend->SetMargin(0.3);
	legend->Draw("same");

	TText *t = new TText();
	t->SetTextAlign(12);
	t->SetTextColor(28);
	t->SetTextSize(24);
	TLine *l = new TLine();
	l->SetLineWidth(2);
	l->SetLineStyle(7);
	l->SetLineColor(28);

	// Right side markings
	c3->Update();
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(16), "O");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");
	// Draw lines marking the particle types
	l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(16), gPad->GetUxmax(), TMath::Log(16));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

	c3->SaveAs("frac_to_lnA_conversion.pdf");

	grlnaEpos->GetXaxis()->SetRangeUser(18.5,19.7);
	grlnaEpos->Draw("APL");
	grlnaQgs->Draw("PL;SAME");
	publnaEpos->Draw("PL;SAME");
	publnaQgs->Draw("PL;SAME");
	legend->Draw("same");

	// Right side markings
	c3->Update();
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(16), "O");
	t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");
	// Draw lines marking the particle types
	l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(16), gPad->GetUxmax(), TMath::Log(16));
	l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

	c3->SaveAs("frac_to_lnA_conversion_zoom.pdf");
}
