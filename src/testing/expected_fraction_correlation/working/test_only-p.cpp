{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	double x11[] = { 0.5, 0.75, 0.25 };
	double x12[] = { 0.65561, 0.82780, 0.48341 };
	double x13[] = { 0.82780, 0.91390, 0.74171 };

	double y11[] = { 0.50140, 0.75360, 0.25431 };
	double y12[] = { 0.74186, 0.87109, 0.61165 };
	double y13[] = { 0.95279, 0.98224, 0.92899 };

	nrp = 3;

	TGraph *gr11 = new TGraph(nrp, x11, y11);
	gr11->SetMarkerStyle(20);
	gr11->SetMarkerColor(1);
	gr11->SetName("p-Fe_p-Fe_data");
	gr11->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr11->GetYaxis()->SetRangeUser(-0.1,1.1);
	gr11->GetXaxis()->SetLimits(-0.1,1.1);
	gr11->GetYaxis()->SetLimits(-0.1,1.1);
	gr11->GetXaxis()->SetTitle("Expected proton fraction");
	gr11->GetYaxis()->SetTitle("Proton fraction");
	gr11->Draw("AP");

	TGraph *gr12 = new TGraph(nrp, x12, y12);
	gr12->SetMarkerStyle(20);
	gr12->SetMarkerColor(2);
	gr12->SetName("p-Fe_p-O_data");
	gr12->Draw("SAME;P");

	TGraph *gr13 = new TGraph(nrp, x13, y13);
	gr13->SetMarkerStyle(20);
	gr13->SetMarkerColor(4);
	gr13->SetName("p-Fe_p-He_data");
	gr13->Draw("SAME;P");

	cout << endl << "Fitting test_only-p.cpp ------------------------------------------" << endl;
	gr11->Fit("pol1");
	gr12->Fit("pol1");
	gr13->Fit("pol1");
}
