{
	// Applying a composition of p, He, O and Fe
	int nrmix = 4;
	TMatrixD *Amat = new TMatrixD(nrmix,nrmix);
	TMatrixD *Bmat = new TMatrixD(nrmix,nrmix);

	TMatrixD *ymu = new TMatrixD(nrmix,nrmix);
	TMatrixD *yeps = new TMatrixD(nrmix,nrmix);
	// p/Fe treatment
	ymu(0,0) = 1.01183;
	ymu(0,1) = 0.91139;
	ymu(0,2) = 0.49110;
	ymu(0,3) = 0.;
	yeps(0,0) = -0.00167;
	yeps(0,1) = 0.00054;
	yeps(0,2) = -0.00056;
	yeps(0,3) = 0.;
	// He/Fe treatment
	ymu(1,0) = 1.07430;
	ymu(1,1) = 1.00622;
	ymu(1,2) = 0.59449;
	ymu(1,3) = 0.;
	yeps(1,0) = 0.00543;
	yeps(1,1) = 0.00132;
	yeps(1,2) = 0.00624;
	yeps(1,3) = 0.;
	// O/Fe treatment
	ymu(2,0) = 1.43579;
	ymu(2,1) = 1.37353;
	ymu(2,2) = 1.01299;
	ymu(2,3) = 0.;
	yeps(2,0) = 0.00847;
	yeps(2,1) = 0.00862;
	yeps(2,2) = 0.00307;
	yeps(2,3) = 0.;
	// Iron constraint (f_p + f_He + f_O + f_Fe = 1)
	ymu(3,0) = 1.;
	ymu(3,1) = 1.;
	ymu(3,2) = 1.;
	ymu(3,3) = 1.;
	yeps(3,0) = 0.;
	yeps(3,1) = 0.;
	yeps(3,2) = 0.;
	yeps(3,3) = 0.;

	for(int i = 0; i < nrmix; i++)
		for(int j = 0; j < nrmix; j++)
			Amat(i,j) = ymu(i,j) - yeps(i,j);
/*	// p/Fe treatment
	Amat(0,0) = 1.01183 - (-0.00167);
	Amat(0,1) = 0.91139 - 0.00054;
	Amat(0,2) = 0.49110 - (-0.00056);
	// He/Fe treatment
	Amat(1,0) = 1.07430 - 0.00543;
	Amat(1,1) = 1.00622 - 0.00132;
	Amat(1,2) = 0.59449 - 0.00624;
	// O/Fe treatment
	Amat(2,0) = 1.43579 - 0.00847;
	Amat(2,1) = 1.37353 - 0.00862;
	Amat(2,2) = 1.01299 - 0.00307;*/

	*Bmat = *Amat;

	cout << "Amat:" << endl;
	for(int i = 0; i < nrmix; i++)
	{
		cout << Amat(i,0) << " | " << Amat(i,1) << " | " << Amat(i,2) << " | " << Amat(i,3) << endl;
	}
	cout << endl;

	double det;
	TMatrixD *AmatInv;
	AmatInv = Bmat->Invert(&det);
	cout << "Determinant = " << det << endl << endl;

	cout << "Inverted Amat:" << endl;
	for(int i = 0; i < nrmix; i++)
	{
		cout << AmatInv(i,0) << " | " << AmatInv(i,1) << " | " << AmatInv(i,2) << " | " << AmatInv(i,3) << endl;
	}
	cout << endl;

	// Example with 50% p, 25% He and 25% Fe
	TMatrixD *Ymat = new TMatrixD(nrmix,1);
	// 50% p + 25% He + 25% Fe
	// p/Fe treatment
	Ymat(0,0) = 0.72908;
	// He/Fe treatment
	Ymat(1,0) = 0.79156;
	// O/Fe treatment
	Ymat(2,0) = 1.06776;
/*	// 50% p + 25% O + 25% Fe
	// p/Fe treatment
	Ymat(0,0) = 0.62489;
	// He/Fe treatment
	Ymat(1,0) = 0.69079;
	// O/Fe treatment
	Ymat(2,0) = 0.97350;*/
/*	// 50% He + 25% O + 25% Fe
	// p/Fe treatment
	Ymat(0,0) = 0.58639;
	// He/Fe treatment
	Ymat(1,0) = 0.65818;
	// O/Fe treatment
	Ymat(2,0) = 0.94258;*/
/*	// 75% p + 25% Fe
	// p/Fe treatment
	Ymat(0,0) = 0.75360;
	// He/Fe treatment
	Ymat(1,0) = 0.80524;
	// O/Fe treatment
	Ymat(2,0) = 1.07554;*/
/*	// Example data point (18.9-19.0, epos, hybrid, test_energyrange)
	// p/Fe treatment
	Ymat(0,0) = 0.55941;
	// He/Fe treatment
	Ymat(1,0) = 0.67013;
	// O/Fe treatment
	Ymat(2,0) = 1.01530;*/
	// Iron constraint (f_p + f_He + f_O + f_Fe = 1)
	Ymat(3,0) = 1.;

	cout << "Ymat (unmodified):" << endl;
	cout << Ymat(0,0) << " | " << Ymat(1,0) << " | " << Ymat(2,0) << " | " << Ymat(3,0) << endl << endl;

	for(int j = 0; j < nrmix; j++)
	{
		// p/Fe treatment
		Ymat(0,0) -= yeps(0,j);
		// He/Fe treatment
		Ymat(1,0) -= yeps(1,j);
		// O/Fe treatment
		Ymat(2,0) -= yeps(2,j);
	}

	cout << "Ymat:" << endl;
	cout << Ymat(0,0) << " | " << Ymat(1,0) << " | " << Ymat(2,0) << " | " << Ymat(3,0) << endl << endl;

	TMatrixD *Fmat = new TMatrixD(nrmix,1);
	Fmat->Mult(*AmatInv, *Ymat);

	cout << "Fmat:" << endl;
	cout << Fmat(0,0) << " | " << Fmat(1,0) << " | " << Fmat(2,0) << " | " << Fmat(3,0) << endl << endl;

	cout << "p/Fe treatment:" << endl;
	cout << "  y (proton) = " << Amat(0,0)*Fmat(0,0) + yeps(0,0) << endl;
	cout << "  y (helium) = " << Amat(0,1)*Fmat(1,0) + yeps(0,1) << endl;
	cout << "  y (oxygen) = " << Amat(0,2)*Fmat(2,0) + yeps(0,2) << endl;
	cout << "  y (sum)    = " << Amat(0,0)*Fmat(0,0) + yeps(0,0) + Amat(0,1)*Fmat(1,0) + yeps(0,1) + Amat(0,2)*Fmat(2,0) + yeps(0,2) << endl << endl;

	cout << "He/Fe treatment:" << endl;
	cout << "  y (proton) = " << Amat(1,0)*Fmat(0,0) + yeps(1,0) << endl;
	cout << "  y (helium) = " << Amat(1,1)*Fmat(1,0) + yeps(1,1) << endl;
	cout << "  y (oxygen) = " << Amat(1,2)*Fmat(2,0) + yeps(1,2) << endl;
	cout << "  y (sum)    = " << Amat(1,0)*Fmat(0,0) + yeps(1,0) + Amat(1,1)*Fmat(1,0) + yeps(1,1) + Amat(1,2)*Fmat(2,0) + yeps(1,2) << endl << endl;

	cout << "O/Fe treatment:" << endl;
	cout << "  y (proton) = " << Amat(2,0)*Fmat(0,0) + yeps(2,0) << endl;
	cout << "  y (helium) = " << Amat(2,1)*Fmat(1,0) + yeps(2,1) << endl;
	cout << "  y (oxygen) = " << Amat(2,2)*Fmat(2,0) + yeps(2,2) << endl;
	cout << "  y (sum)    = " << Amat(2,0)*Fmat(0,0) + yeps(2,0) + Amat(2,1)*Fmat(1,0) + yeps(2,1) + Amat(2,2)*Fmat(2,0) + yeps(2,2) << endl << endl;

	cout << "lnA = " << Fmat(0,0)*TMath::Log(1) + Fmat(1,0)*TMath::Log(4) + Fmat(2,0)*TMath::Log(16) + Fmat(3,0)*TMath::Log(56) << endl << endl;

	cout << "Composition:" << endl;
	cout << "  proton = " << Fmat(0,0)*100. << "%" << endl;
	cout << "  helium = " << Fmat(1,0)*100. << "%" << endl;
	cout << "  oxygen = " << Fmat(2,0)*100. << "%" << endl;
	cout << "  iron   = " << Fmat(3,0)*100. << "%" << endl;
}
