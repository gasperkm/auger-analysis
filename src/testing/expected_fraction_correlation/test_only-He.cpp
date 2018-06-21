{
	TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
	c1->SetGrid();
	int nrp;

	double x11[] = { 0.32780, 0.49171, 0.16390 };
	double x12[] = { 0.48341, 0.56951, 0.39732 };

	double y11[] = { 0.49362, 0.74677, 0.26004 };
	double y12[] = { 0.69518, 0.79732, 0.59112 };

	nrp = 3;

	TGraph *gr11 = new TGraph(nrp, x11, y11);
	gr11->SetMarkerStyle(20);
	gr11->SetMarkerColor(1);
	gr11->SetName("p-Fe_He-Fe_data");
	gr11->GetXaxis()->SetRangeUser(-0.1,1.1);
	gr11->GetYaxis()->SetRangeUser(-0.1,1.2);
	gr11->GetXaxis()->SetLimits(-0.1,1.1);
	gr11->GetYaxis()->SetLimits(-0.1,1.2);
	gr11->GetXaxis()->SetTitle("Expected proton fraction");
	gr11->GetYaxis()->SetTitle("Proton fraction");
	gr11->Draw("AP");

	TGraph *gr12 = new TGraph(nrp, x12, y12);
	gr12->SetMarkerStyle(20);
	gr12->SetMarkerColor(2);
	gr12->SetName("p-Fe_He-O_data");
	gr12->Draw("SAME;P");

	cout << endl << "Fitting test_only-He.cpp -----------------------------------------" << endl;
	gr11->Fit("pol1");
	gr12->Fit("pol1");
}
