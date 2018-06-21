{
	// Applying a composition of p, He and Fe
	int nrmix = 2;
	TMatrixD *Amat = new TMatrixD(nrmix,nrmix);
	TMatrixD *Bmat = new TMatrixD(nrmix,nrmix);

	TMatrixD *ymu = new TMatrixD(nrmix,nrmix);
	TMatrixD *yeps = new TMatrixD(nrmix,nrmix);
	// p/Fe treatment
	ymu(0,0) = 1.01183;
	ymu(0,1) = 0.91139;
	yeps(0,0) = -0.00167;
	yeps(0,1) = 0.00054;
	// He/Fe treatment
	ymu(1,0) = 1.07430;
	ymu(1,1) = 1.00622;
	yeps(1,0) = 0.00543;
	yeps(1,1) = 0.00132;

	for(int i = 0; i < nrmix; i++)
		for(int j = 0; j < nrmix; j++)
			Amat(i,j) = ymu(i,j) - yeps(i,j);

	*Bmat = *Amat;

	cout << "Amat:" << endl;
	for(int i = 0; i < nrmix; i++)
	{
		cout << Amat(i,0) << " | " << Amat(i,1) << endl;
	}
	cout << endl;

	double det;
	TMatrixD *AmatInv;
	AmatInv = Bmat->Invert(&det);
	cout << "Determinant = " << det << endl << endl;

	cout << "Inverted Amat:" << endl;
	for(int i = 0; i < nrmix; i++)
	{
		cout << AmatInv(i,0) << " | " << AmatInv(i,1) << endl;
	}
	cout << endl;

	// Example with 50% p, 25% He and 25% Fe
	TMatrixD *Ymat = new TMatrixD(nrmix,1);
	// p/Fe treatment
	Ymat(0,0) = 0.72908;
	// He/Fe treatment
	Ymat(1,0) = 0.79156;
/*	// p/Fe treatment
	Ymat(0,0) = 0.71797;
	// He/Fe treatment
	Ymat(1,0) = 0.77597;*/
/*	// p/Fe treatment
	Ymat(0,0) = 0.74043;
	// He/Fe treatment
	Ymat(1,0) = 0.80613;*/

	cout << "Ymat (unmodified):" << endl;
	cout << Ymat(0,0) << " | " << Ymat(1,0) << endl << endl;

	for(int j = 0; j < nrmix; j++)
	{
		// p/Fe treatment
		Ymat(0,0) -= yeps(0,j);
		// He/Fe treatment
		Ymat(1,0) -= yeps(1,j);
	}

	cout << "Ymat:" << endl;
	cout << Ymat(0,0) << " | " << Ymat(1,0) << endl << endl;

	TMatrixD *Fmat = new TMatrixD(nrmix,1);
	Fmat->Mult(*AmatInv, *Ymat);

	cout << "Fmat:" << endl;
	cout << Fmat(0,0) << " | " << Fmat(1,0) << endl << endl;

	cout << "p/Fe treatment:" << endl;
	cout << "  y (proton) = " << Amat(0,0)*Fmat(0,0) + yeps(0,0) << endl;
	cout << "  y (helium) = " << Amat(0,1)*Fmat(1,0) + yeps(0,1) << endl;
	cout << "  y (sum)    = " << Amat(0,0)*Fmat(0,0) + yeps(0,0) + Amat(0,1)*Fmat(1,0) + yeps(0,1) << endl << endl;

	cout << "He/Fe treatment:" << endl;
	cout << "  y (proton) = " << Amat(1,0)*Fmat(0,0) + yeps(1,0) << endl;
	cout << "  y (helium) = " << Amat(1,1)*Fmat(1,0) + yeps(1,1) << endl;
	cout << "  y (sum)    = " << Amat(1,0)*Fmat(0,0) + yeps(1,0) + Amat(1,1)*Fmat(1,0) + yeps(1,1) << endl << endl;

	cout << "lnA = " << Fmat(0,0)*TMath::Log(1) + Fmat(1,0)*TMath::Log(4) + (1.-Fmat(0,0)-Fmat(1,0))*TMath::Log(56) << endl << endl;

	cout << "Composition:" << endl;
	cout << "  proton = " << Fmat(0,0)*100. << "%" << endl;
	cout << "  helium = " << Fmat(1,0)*100. << "%" << endl;
	cout << "  iron   = " << (1.-Fmat(0,0)-Fmat(1,0))*100. << "%" << endl;
}
