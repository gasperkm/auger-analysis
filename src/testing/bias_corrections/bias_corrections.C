double FDXmaxCorr(double *x, double *par)
{
	double logE = TMath::Log10(x[0]);
	double a = 6.5/(TMath::Exp((logE - 18.23)/0.41) + 1.);
	double b = 0.93*(logE - 18.);
	return -a + 3.4 - b;
}

double OrigEnergy(double *x, double *par)
{
	return x[0];
}

double HECOEnergyCorr(double *x, double *par)
{
	double logE = TMath::Log10(x[0]);
	double a = TMath::Log(logE - 16.5);
	return x[0]*(1 + 0.04536243 - 0.08317193*a);
}

double HECOXmaxCorr(double *x, double *par)
{
	double logE;

	if(par[0] == 0)
		logE = TMath::Log10(x[0]);
	else
		logE = TMath::Log10(HECOEnergyCorr(x, par));

	double xmaxRecBias = -2.96976097 - 0.99199218*(logE - 17.80937342);
	double bLWBias = 6.5/(TMath::Exp((logE - 18.23)/0.41) + 1.);
	double xalignBias;
	if(logE < 17.55)
		xalignBias = 0.12 - 6.43*(logE - 17.55);
	else
		xalignBias = 0.12 - 0.27*(logE - 17.55);

	return xmaxRecBias + bLWBias + 0.5*xalignBias;
}

void bias_corrections()
{
	double inenergy = 1.e+19;

//	cout << "Enter the energy to calculate the energy and Xmax corrections: ";
//	cin >> inenergy;

	cout << "Calculating FD standard correction..." << endl;

	TCanvas *c1 = new TCanvas("c1", "FD standard Xmax correction", 800, 600);
	c1->SetGrid();
	c1->SetLogx();

//	TF1 *fdxmaxcor = new TF1("fdxmaxcor", "-6.5/(TMath::Exp((TMath::Log10(x) - 18.23)/0.41) + 1.) + 3.4 - 0.93*(TMath::Log10(x) - 18.)", 1.e+18, 1.e+20);
	TF1 *fdxmaxcor = new TF1("fdxmaxcor", FDXmaxCorr, TMath::Power(10,18.5), TMath::Power(10,20.), 0);

//	cout << "  Xmax correction is: " << FDXmaxCorr(inenergy) << endl;
	cout << "  Xmax correction is (E = 1.e+19): " << fdxmaxcor->Eval(inenergy) << endl;

	fdxmaxcor->Draw();

	cout << "Calculating HECO corrections..." << endl;

	TCanvas *c2 = new TCanvas("c2", "HECO energy correction", 800, 600);
	c2->SetGrid();
	c2->SetLogx();
	c2->SetLogy();

	TF1 *hecoenercor = new TF1("hecoenercor", HECOEnergyCorr, TMath::Power(10,18.5), TMath::Power(10,20.), 0);
	TF1 *origener = new TF1("origener", OrigEnergy, TMath::Power(10,18.5), TMath::Power(10,20.), 0);

	cout << "  Energy correction is (E = 1.e+19): " << hecoenercor->Eval(inenergy) << endl;

	hecoenercor->Draw();
	origener->SetLineColor(2);
	origener->Draw("SAME");

	TCanvas *c3 = new TCanvas("c3", "HECO Xmax correction", 800, 600);
	c3->SetGrid();
	c3->SetLogx();

	TF1 *hecoxmaxcor1 = new TF1("hecoxmaxcor1", HECOXmaxCorr, TMath::Power(10,18.5), TMath::Power(10,20.), 1);
	hecoxmaxcor1->SetParameter(0, 1);

	TF1 *hecoxmaxcor2 = new TF1("hecoxmaxcor2", HECOXmaxCorr, TMath::Power(10,18.5), TMath::Power(10,20.), 1);
	hecoxmaxcor2->SetParameter(0, 0);

	cout << "  Xmax correction is (E = 1.e+19): " << hecoxmaxcor1->Eval(inenergy) << ", No energy correction: " << hecoxmaxcor2->Eval(inenergy) << endl;

	hecoxmaxcor1->Draw();
	hecoxmaxcor2->SetLineColor(2);
	hecoxmaxcor2->Draw("SAME");
}
