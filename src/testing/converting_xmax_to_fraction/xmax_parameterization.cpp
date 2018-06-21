// Input parameters: 0 = E0, 1 = X0, 2 = D, 3 = eta, 4 = delta, 5 = A [energies are in units of log(E/eV)]
double ParamFunction(double *x, double *par)
{
	// Expanding Log10(E/E0) into Log10(E) - Log10(E0)
	double logTerm1 = x[0] - par[0];
	// Expanding Log10(E/(E0*A)) into Log10(E) - Log10(E0) - Log10(A)
	double logTerm2 = logTerm1 - TMath::Log10(par[5]);
	return par[1] + par[2]*logTerm2 + par[3]*TMath::Log(par[5]) + par[4]*TMath::Log(par[5])*logTerm1;
}

double fEParam(double *x, double *par)
{
	double logTerm = x[0] - par[0];
	return par[3] - par[2]/TMath::Log(10) + par[4]*logTerm;
}

void xmax_parametrization()
{
	double E0 = 19.;
	// Arranging of parameters: X0, D, eta, delta
	vector<double> sibParam, eposParam, qgsParam;
	sibParam.push_back(819.3);
	sibParam.push_back(57.4);
	sibParam.push_back(-0.56);
	sibParam.push_back(0.68);
	eposParam.push_back(806.0);
	eposParam.push_back(56.3);
	eposParam.push_back(0.35);
	eposParam.push_back(1.04);
	qgsParam.push_back(790.1);
	qgsParam.push_back(54.2);
	qgsParam.push_back(-0.42);
	qgsParam.push_back(0.69);

	double xlimits[] = {17.00, 20.00};
	double ylimits[] = {570., 870.};

	TCanvas *c1 = new TCanvas("c1", "Hadronic models", 820, 800);
	c1->SetGrid();

	// Proton Sibyll
	TF1 *protSib = new TF1("protSib", ParamFunction, xlimits[0], xlimits[1], 6);
	protSib->SetParameter(0, E0);
	for(int i = 0; i < sibParam.size(); i++)
		protSib->SetParameter(i+1, sibParam[i]);
	protSib->SetParameter(5, 1.);

	protSib->GetYaxis()->SetRangeUser(ylimits[0],ylimits[1]);
	protSib->SetLineColor(2);
	protSib->SetLineStyle(10);
	protSib->Draw();

	// Iron Sibyll
	TF1 *ironSib = new TF1("ironSib", ParamFunction, xlimits[0], xlimits[1], 6);
	ironSib->SetParameter(0, E0);
	for(int i = 0; i < sibParam.size(); i++)
		ironSib->SetParameter(i+1, sibParam[i]);
	ironSib->SetParameter(5, 56.);

	ironSib->SetLineColor(4);
	ironSib->SetLineStyle(10);
	ironSib->Draw("same");

	// Proton Epos
	TF1 *protEpos = new TF1("protEpos", ParamFunction, xlimits[0], xlimits[1], 6);
	protEpos->SetParameter(0, E0);
	for(int i = 0; i < eposParam.size(); i++)
		protEpos->SetParameter(i+1, eposParam[i]);
	protEpos->SetParameter(5, 1.);

	protEpos->SetLineColor(2);
	protEpos->SetLineStyle(1);
	protEpos->Draw("same");

	// Iron Epos
	TF1 *ironEpos = new TF1("ironEpos", ParamFunction, xlimits[0], xlimits[1], 6);
	ironEpos->SetParameter(0, E0);
	for(int i = 0; i < eposParam.size(); i++)
		ironEpos->SetParameter(i+1, eposParam[i]);
	ironEpos->SetParameter(5, 56.);

	ironEpos->SetLineColor(4);
	ironEpos->SetLineStyle(1);
	ironEpos->Draw("same");

	// Proton Qgsjet
	TF1 *protQgs = new TF1("protQgs", ParamFunction, xlimits[0], xlimits[1], 6);
	protQgs->SetParameter(0, E0);
	for(int i = 0; i < qgsParam.size(); i++)
		protQgs->SetParameter(i+1, qgsParam[i]);
	protQgs->SetParameter(5, 1.);

	protQgs->SetLineColor(2);
	protQgs->SetLineStyle(9);
	protQgs->Draw("same");

	// Iron Qgsjet
	TF1 *ironQgs = new TF1("ironQgs", ParamFunction, xlimits[0], xlimits[1], 6);
	ironQgs->SetParameter(0, E0);
	for(int i = 0; i < qgsParam.size(); i++)
		ironQgs->SetParameter(i+1, qgsParam[i]);
	ironQgs->SetParameter(5, 56.);

	ironQgs->SetLineColor(4);
	ironQgs->SetLineStyle(9);
	ironQgs->Draw("same");

	// Plotting the fE parameter used for lnA estimation
	TCanvas *c2 = new TCanvas("c2", "fE parameter", 820, 800);
	c2->SetGrid();

	TF1 *fESib = new TF1("fESib", fEParam, xlimits[0], xlimits[1], 6);
	fESib->SetParameter(0, E0);
	for(int i = 0; i < sibParam.size(); i++)
		fESib->SetParameter(i+1, sibParam[i]);

	fESib->GetYaxis()->SetRangeUser(-27., -23.);
	fESib->SetLineColor(1);
	fESib->SetLineStyle(10);
	fESib->Draw();

	TF1 *fEEpos = new TF1("fEEpos", fEParam, xlimits[0], xlimits[1], 6);
	fEEpos->SetParameter(0, E0);
	for(int i = 0; i < eposParam.size(); i++)
		fEEpos->SetParameter(i+1, eposParam[i]);

	fEEpos->SetLineColor(1);
	fEEpos->SetLineStyle(1);
	fEEpos->Draw("same");

	TF1 *fEQgs = new TF1("fEQgs", fEParam, xlimits[0], xlimits[1], 6);
	fEQgs->SetParameter(0, E0);
	for(int i = 0; i < qgsParam.size(); i++)
		fEQgs->SetParameter(i+1, qgsParam[i]);

	fEQgs->SetLineColor(1);
	fEQgs->SetLineStyle(9);
	fEQgs->Draw("same");
}
