#include <vector>

void PrintMatrix(TMatrixD *mat, string matname)
{
	cout << "------" << endl;
	cout << matname << ":" << endl;
	for(int i = 0; i < mat->GetNrows(); i++)
	{
		for(int j = 0; j < mat->GetNcols(); j++)
		{
			cout << mat(i,j);

			if(j+1 != mat->GetNcols())
				cout << " | ";
			else
				cout << endl;
		}
	}
	cout << "------" << endl;
	cout << endl;
}

void lnA_inverting_matrix()
{
	// Select elements to include during the process (p, He, O, Fe)
	vector<int> includeElem;
        includeElem.push_back(1);
        includeElem.push_back(1);
        includeElem.push_back(1);
        includeElem.push_back(1);

	vector<string> includeName;
        includeName.push_back("proton");
        includeName.push_back("helium");
        includeName.push_back("oxygen");
        includeName.push_back("iron");

	vector<double> includeLna;
        includeLna.push_back(TMath::Log(1));
        includeLna.push_back(TMath::Log(4));
        includeLna.push_back(TMath::Log(16));
        includeLna.push_back(TMath::Log(56));

	int nrall = includeElem.size();
	int nrmix = 0;
	for(int i = 0; i < nrall; i++)
	{
		if(includeElem[i] == 1)
			nrmix++;
	}

	if(includeElem[3] == 0)
	{
		cout << "Error! Iron must always be selected for the analysis." << endl;
		return;
	}

	if(nrmix < 2)
	{
		cout << "Error! Analysis only supports a mixture of 3 elements or more." << endl;
		return;
	}

	TMatrixD *ymu = new TMatrixD(nrall,nrall);
	TMatrixD *yeps = new TMatrixD(nrall,nrall);
	TMatrixD *Amat = new TMatrixD(nrmix,nrmix);
	TMatrixD *Bmat = new TMatrixD(nrmix,nrmix);
	TMatrixD *AmatInv;
	TMatrixD *Ymat = new TMatrixD(nrmix,1);
	TMatrixD *Fmat = new TMatrixD(nrmix,1);

	cout << endl << "Explanation: Rows in matrices means different sig-bgd treatments. Columns in matrices means a different mock data mixture (p-Fe, He-Fe, O-Fe)." << endl << endl;

	// y_mu and y_epsilon values from fitting the p-Fe, He-Fe and O-Fe mixed mock dataset
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

	PrintMatrix(ymu, "y_mu");
	PrintMatrix(yeps, "y_epsilon");

	// Enter values from determining proton fraction with MVA analysis
	// 50% p + 25% He + 25% Fe
	// p/Fe treatment
	if(nrmix == 2)
	{
		Ymat(0,0) = 0.72908;
		Ymat(1,0) = 1.;
	}
	// He/Fe treatment
	if(nrmix == 3)
	{
		Ymat(0,0) = 0.72908;
		Ymat(1,0) = 0.79156;
		Ymat(2,0) = 1.;
	}
	// O/Fe treatment
	if(nrmix == 4)
	{
		Ymat(0,0) = 0.72908;
		Ymat(1,0) = 0.79156;
		Ymat(2,0) = 1.06776;
		Ymat(3,0) = 1.;
	}

	int jcnt;

	for(int i = 0; i < nrmix; i++) // go over rows (sig-bgd treatments) -> only take the first nrmix rows
	{
		jcnt = 0;
		for(int j = 0; j < nrall; j++)	// go over cols (mock data mixtures)
		{
			if( includeElem[j] == 1 )
				Amat(i,jcnt) = ymu(i,j) - yeps(i,j);

			if( includeElem[j] == 1 )
				jcnt++;
		}
	}

	*Bmat = *Amat;
	PrintMatrix(Amat, "Amat (y_mu - y_epsilon)");

	AmatInv = Bmat->Invert();
	PrintMatrix(AmatInv, "AmatInv");

	PrintMatrix(Ymat, "Ymat (MVA analysis results)");

	for(int i = 0; i < nrmix; i++)
	{
		jcnt = 0;
		for(int j = 0; j < nrall; j++)
		{
			Ymat(i,0) -= yeps(i,jcnt);

			if( includeElem[j] == 1 )
				jcnt++;
		}
	}

	PrintMatrix(Ymat, "Ymat");

	Fmat->Mult(*AmatInv, *Ymat);
	PrintMatrix(Fmat, "Fmat");

	double sum;

	for(int i = 0; i < nrmix; i++)
	{
		cout << includeName[i] << "/iron treatment:" << endl;

		sum = 0;
		jcnt = 0;
		for(int j = 0; j < nrall; j++)
		{
			if(includeElem[j] == 1)
			{
				cout << "  y (" << includeName[j] << ") = " << Amat(i,jcnt)*Fmat(jcnt,0) + yeps(i,j) << endl;
				sum += Amat(i,jcnt)*Fmat(jcnt,0) + yeps(i,j);

				jcnt++;
			}
		}
		cout << "  y (sum)    = " << sum << endl << endl;
	}

	double iron = 1.;

	cout << "lnA = ";
	sum = 0;
	jcnt = 0;
       	for(int i = 0; i < nrall; i++)
	{
		if(includeElem[i] == 1)
		{
			sum += Fmat(jcnt,0)*includeLna[i];
			iron -= Fmat(jcnt,0);

			jcnt++;
		}
	}
	
	cout << sum + iron*TMath::Log(56) << endl << endl;

	return;
	
/*	// Applying a composition of p, He, O and Fe
	int nrmix = 3;
	TMatrixD *Amat = new TMatrixD(nrmix,nrmix);
	TMatrixD *Bmat = new TMatrixD(nrmix,nrmix);

	TMatrixD *ymu = new TMatrixD(nrmix,nrmix);
	TMatrixD *yeps = new TMatrixD(nrmix,nrmix);
	// p/Fe treatment
	ymu(0,0) = 1.01183;
	ymu(0,1) = 0.91139;
	ymu(0,2) = 0.49110;
	yeps(0,0) = -0.00167;
	yeps(0,1) = 0.00054;
	yeps(0,2) = -0.00056;
	// He/Fe treatment
	ymu(1,0) = 1.07430;
	ymu(1,1) = 1.00622;
	ymu(1,2) = 0.59449;
	yeps(1,0) = 0.00543;
	yeps(1,1) = 0.00132;
	yeps(1,2) = 0.00624;
	// O/Fe treatment
	ymu(2,0) = 1.43579;
	ymu(2,1) = 1.37353;
	ymu(2,2) = 1.01299;
	yeps(2,0) = 0.00847;
	yeps(2,1) = 0.00862;
	yeps(2,2) = 0.00307;

	for(int i = 0; i < nrmix; i++)
		for(int j = 0; j < nrmix; j++)
			Amat(i,j) = ymu(i,j) - yeps(i,j);

	*Bmat = *Amat;

	cout << "Amat:" << endl;
	for(int i = 0; i < nrmix; i++)
	{
		cout << Amat(i,0) << " | " << Amat(i,1) << " | " << Amat(i,2) << endl;
	}
	cout << endl;

	double det;
	TMatrixD *AmatInv;
	AmatInv = Bmat->Invert(&det);
	cout << "Determinant = " << det << endl << endl;

	cout << "Inverted Amat:" << endl;
	for(int i = 0; i < nrmix; i++)
	{
		cout << AmatInv(i,0) << " | " << AmatInv(i,1) << " | " << AmatInv(i,2) << endl;
	}
	cout << endl;

	// Example with 50% p, 25% He and 25% Fe
	TMatrixD *Ymat = new TMatrixD(nrmix,1);
	// 50% p + 25% He + 25% Fe*/
/*	// p/Fe treatment
	Ymat(0,0) = 0.72908;
	// He/Fe treatment
	Ymat(1,0) = 0.79156;
	// O/Fe treatment
	Ymat(2,0) = 1.06776;*/
	// 50% p + 25% O + 25% Fe
/*	// p/Fe treatment
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
	// Example data point (18.9-19.0, epos, hybrid, test_energyrange)
/*	// p/Fe treatment
	Ymat(0,0) = 0.55941;
	// He/Fe treatment
	Ymat(1,0) = 0.67013;
	// O/Fe treatment
	Ymat(2,0) = 1.01530;

	cout << "Ymat (unmodified):" << endl;
	cout << Ymat(0,0) << " | " << Ymat(1,0) << " | " << Ymat(2,0) << endl << endl;

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
	cout << Ymat(0,0) << " | " << Ymat(1,0) << " | " << Ymat(2,0) << endl << endl;

	TMatrixD *Fmat = new TMatrixD(nrmix,1);
	Fmat->Mult(*AmatInv, *Ymat);

	cout << "Fmat:" << endl;
	cout << Fmat(0,0) << " | " << Fmat(1,0) << " | " << Fmat(2,0) << endl << endl;

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

	cout << "lnA = " << Fmat(0,0)*TMath::Log(1) + Fmat(1,0)*TMath::Log(4) + Fmat(2,0)*TMath::Log(16) + (1.-Fmat(0,0)-Fmat(1,0)-Fmat(2,0))*TMath::Log(56) << endl << endl;

	cout << "Composition:" << endl;
	cout << "  proton = " << Fmat(0,0)*100. << "%" << endl;
	cout << "  helium = " << Fmat(1,0)*100. << "%" << endl;
	cout << "  oxygen = " << Fmat(2,0)*100. << "%" << endl;
	cout << "  iron   = " << (1.-Fmat(0,0)-Fmat(1,0)-Fmat(2,0))*100. << "%" << endl;*/
}
