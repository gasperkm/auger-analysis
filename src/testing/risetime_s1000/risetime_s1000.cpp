void SetColor(TTree *tr, int cur, int nrbins)
{
	int ci = 1738+cur;
	if(cur < nrbins-1.)
	{
		double colormix = (double)cur/(double)(nrbins-1.);
		TColor *color = new TColor(ci, colormix, 0, 1.-colormix, "", 1);
/*		tr->SetMarkerColor(ci);
		tr->SetMarkerSize(0.9);
		tr->SetMarkerStyle(20+cur);*/
		tr->SetLineColor(ci);
		tr->SetLineWidth(2);
	}
	else
	{
/*		tr->SetMarkerColor(3);
		tr->SetMarkerSize(0.9);
		tr->SetMarkerStyle(28);*/
		tr->SetLineColor(3);
		tr->SetLineWidth(2);
	}
}

void risetime_s1000()
{
	string inname = "temporary_event_file.root";
	string stemp[3];

	TFile *infile = TFile::Open(inname.c_str(), "READ");
	TList *tempkeyslist = (TList*)infile->GetListOfKeys();

	TTree *simptree = new TTree;
	TCanvas *c1 = new TCanvas("c1", "", 1200, 900);
	c1->SetLogy();

	for(int i = 0; i < infile->GetNkeys(); i++)
	{
		stemp[0] = string((tempkeyslist->At(i))->GetName());
		simptree = (TTree*)infile->Get(stemp[0].c_str());
		SetColor(simptree, i, infile->GetNkeys());

		if(i == 0)
			simptree->Draw("xmax");
		else
			simptree->Draw("xmax", "", "SAME");
	}

/*	stemp[1] = "output_xmax_" + stemp[0] + ".pdf";
	c1->SaveAs(stemp[1].c_str());*/
	c1->SaveAs("output_xmax.pdf");

	for(int i = 0; i < infile->GetNkeys(); i++)
	{
		stemp[0] = string((tempkeyslist->At(i))->GetName());
		simptree = (TTree*)infile->Get(stemp[0].c_str());
		SetColor(simptree, i, infile->GetNkeys());

		if(i == 0)
			simptree->Draw("shwsize");
		else
			simptree->Draw("shwsize", "", "SAME");
	}

/*	stemp[1] = "output_shwsize_" + stemp[0] + ".pdf";
	c1->SaveAs(stemp[1].c_str());*/
	c1->SaveAs("output_shwsize.pdf");

	for(int i = 0; i < infile->GetNkeys(); i++)
	{
		stemp[0] = string((tempkeyslist->At(i))->GetName());
		simptree = (TTree*)infile->Get(stemp[0].c_str());
		SetColor(simptree, i, infile->GetNkeys());

		if(i == 0)
			simptree->Draw("risetimerecalc");
		else
			simptree->Draw("risetimerecalc", "", "SAME");
	}

/*	stemp[1] = "output_risetimerecalc_" + stemp[0] + ".pdf";
	c1->SaveAs(stemp[1].c_str());*/
	c1->SaveAs("output_risetimerecalc.pdf");

	for(int i = 0; i < infile->GetNkeys(); i++)
	{
		stemp[0] = string((tempkeyslist->At(i))->GetName());
		simptree = (TTree*)infile->Get(stemp[0].c_str());
		SetColor(simptree, i, infile->GetNkeys());

		if(i == 0)
			simptree->Draw("shwsize*shwsize/risetimerecalc");
		else
			simptree->Draw("shwsize*shwsize/risetimerecalc", "", "SAME");
	}

/*	stemp[1] = "output_risetimerecalc_" + stemp[0] + ".pdf";
	c1->SaveAs(stemp[1].c_str());*/
	c1->SaveAs("output_shwsize-div-risetimerecalc.pdf");
}
