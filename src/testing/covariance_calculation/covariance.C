void covariance()
{
	double *event[4];
	double *eventerr[4];
	double *norm[4];
	double *normerr[4];
	double dtemp[3];

	for(int i = 0; i < 4; i++)
	{
		event[i] = new double[3];
		eventerr[i] = new double[3];
		norm[i] = new double[3];
		normerr[i] = new double[3];
	}

	double xmaxlimit[] = {676.088, 1096.93};
	double s1000limit[] = {73.3858, 472.542};
	double t12limit[] = {90.4633, 1189.9};

	event[0][0] = 806.927;
	event[1][0] = 779.257;
	event[2][0] = 760.044;
	event[3][0] = 761.580;
	eventerr[0][0] = 6.418;
	eventerr[1][0] = 18.1718;
	eventerr[2][0] = 25.9463;
	eventerr[3][0] = 9.54991;

	event[0][1] = 275.574;
	event[1][1] = 244.039;
	event[2][1] = 215.782;
	event[3][1] = 517.124;
	eventerr[0][1] = 9.18342;
	eventerr[1][1] = 5.62556;
	eventerr[2][1] = 6.47204;
	eventerr[3][1] = 14.8658;

	event[0][2] = 255.298;
	event[1][2] = 340.692;
	event[2][2] = 186.756;
	event[3][2] = 330.928;
	eventerr[0][2] = 10.0483;
	eventerr[1][2] = 83.4246;
	eventerr[2][2] = 7.75369;
	eventerr[3][2] = 4.86356;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			dtemp[0] = event[i][j] + eventerr[i][j];

			if(j == 0)
			{
				norm[i][j] = 2*(event[i][j] - xmaxlimit[0])/(xmaxlimit[1] - xmaxlimit[0]) - 1.;
				normerr[i][j] = 2*(dtemp[0] - xmaxlimit[0])/(xmaxlimit[1] - xmaxlimit[0]) - 1.;
				normerr[i][j] = normerr[i][j] - norm[i][j];
			}
			else if(j == 1)
			{
				norm[i][j] = 2*(event[i][j] - s1000limit[0])/(s1000limit[1] - s1000limit[0]) - 1.;
				normerr[i][j] = 2*(dtemp[0] - s1000limit[0])/(s1000limit[1] - s1000limit[0]) - 1.;
				normerr[i][j] = normerr[i][j] - norm[i][j];
			}
			else if(j == 2)
			{
				norm[i][j] = 2*(event[i][j] - t12limit[0])/(t12limit[1] - t12limit[0]) - 1.;
				normerr[i][j] = 2*(dtemp[0] - t12limit[0])/(t12limit[1] - t12limit[0]) - 1.;
				normerr[i][j] = normerr[i][j] - norm[i][j];
			}
		}

		cout << i << ": Norm    = " << norm[i][0] << "\t" << norm[i][1] << "\t" << norm[i][2] << endl;
		cout << i << ": NormErr = " << normerr[i][0] << "\t" << normerr[i][1] << "\t" << normerr[i][2] << endl;
	}

	double mean[3], meanerr[3];

	mean[0] = 0;
	mean[1] = 0;
	mean[2] = 0;
	meanerr[0] = 0;
	meanerr[1] = 0;
	meanerr[2] = 0;

	for(int i = 0; i < 4; i++)
	{
		mean[0] += norm[i][0];
		mean[1] += norm[i][1];
		mean[2] += norm[i][2];

		meanerr[0] += normerr[i][0];
		meanerr[1] += normerr[i][1];
		meanerr[2] += normerr[i][2];
	}

	mean[0] = mean[0]/4;
	mean[1] = mean[1]/4;
	mean[2] = mean[2]/4;

	meanerr[0] = meanerr[0]/4;
	meanerr[1] = meanerr[1]/4;
	meanerr[2] = meanerr[2]/4;

	cout << "Mean    = " << mean[0] << "\t" << mean[1] << "\t" << mean[2] << endl;
	cout << "MeanErr = " << meanerr[0] << "\t" << meanerr[1] << "\t" << meanerr[2] << endl;

	double *diff[4];
	double *differr[4];

	for(int i = 0; i < 4; i++)
	{
		diff[i] = new double[3];
		differr[i] = new double[3];
	}

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			diff[i][j] = norm[i][j] - mean[j];
			differr[i][j] = normerr[i][j] - meanerr[j];
		}

		cout << i << ": Diff    = " << diff[i][0] << "\t" << diff[i][1] << "\t" << diff[i][2] << endl;
		cout << i << ": DiffErr = " << differr[i][0] << "\t" << differr[i][1] << "\t" << differr[i][2] << endl;
	}

	TMatrixD *cov = new TMatrixD(3,3);
	TMatrixD *coverr = new TMatrixD(3,3);

	cout << "Matrix:" << endl;
	for(int k = 0; k < 3; k++)
	{
		for(int j = 0; j < 3; j++)
		{
			dtemp[0] = 0;
			for(int i = 0; i < 4; i++)
			{
				dtemp[1] = diff[i][j];
				dtemp[2] = diff[i][k];
				dtemp[0] += (dtemp[1]*dtemp[2]);
			}
			dtemp[0] = dtemp[0]/3.;
			(*cov)(j,k) = dtemp[0];
		}

		cout << (*cov)(0,k) << " | " << (*cov)(1,k) << " | " << (*cov)(2,k) << endl;
	}

	cout << "Matrix (err):" << endl;
	for(int k = 0; k < 3; k++)
	{
		for(int j = 0; j < 3; j++)
		{
			dtemp[0] = 0;
			for(int i = 0; i < 4; i++)
			{
				dtemp[1] = differr[i][j];
				dtemp[2] = differr[i][k];
				dtemp[0] += (dtemp[1]*dtemp[2]);
			}
			dtemp[0] = dtemp[0]/3.;
			(*coverr)(j,k) = dtemp[0];
		}

		cout << (*coverr)(0,k) << " | " << (*coverr)(1,k) << " | " << (*coverr)(2,k) << endl;
	}

	TMatrixDEigen *eigencov = new TMatrixDEigen((const TMatrixD)*cov);
	TMatrixD *dcov = new TMatrixD(3,3);
	(*dcov) = eigencov->GetEigenValues();

	cout << "Diagonal values: " << endl;
	dtemp[0] = 0;
	for(int j = 0; j < 3; j++)
	{
		dtemp[0] += TMath::Power((*dcov)(j,j), 2);
		cout << (*dcov)(j,j) << " ";
	}
	cout << endl;

	cout << "Final error = " << TMath::Sqrt(dtemp[0]) << endl;

	TMatrixDEigen *eigencoverr = new TMatrixDEigen((const TMatrixD)*coverr);
	TMatrixD *dcoverr = new TMatrixD(3,3);
	(*dcoverr) = eigencoverr->GetEigenValues();

	cout << "Diagonal values (err): " << endl;
	dtemp[0] = 0;
	for(int j = 0; j < 3; j++)
	{
		dtemp[0] += TMath::Power((*dcoverr)(j,j), 2);
		cout << (*dcoverr)(j,j) << " ";
	}
	cout << endl;

	cout << "Final error (err) = " << TMath::Sqrt(dtemp[0]) << endl;

	double sigma[3], sigmaerr[3];

	sigma[0] = 0;
	sigma[1] = 0;
	sigma[2] = 0;

	sigmaerr[0] = 0;
	sigmaerr[1] = 0;
	sigmaerr[2] = 0;

	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 4; i++)
		{
			dtemp[0] = norm[i][j] - mean[j];
			dtemp[0] = TMath::Power(dtemp[0], 2);
			sigma[j] += dtemp[0];

			dtemp[0] = normerr[i][j] - meanerr[j];
			dtemp[0] = TMath::Power(dtemp[0], 2);
			sigmaerr[j] += dtemp[0];
		}
	}

	sigma[0] = TMath::Sqrt(sigma[0]/3.);
	sigma[1] = TMath::Sqrt(sigma[1]/3.);
	sigma[2] = TMath::Sqrt(sigma[2]/3.);

	sigmaerr[0] = TMath::Sqrt(sigmaerr[0]/3.);
	sigmaerr[1] = TMath::Sqrt(sigmaerr[1]/3.);
	sigmaerr[2] = TMath::Sqrt(sigmaerr[2]/3.);

	cout << "Sigma    = " << sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << endl;
	cout << "SigmaErr = " << sigmaerr[0] << "\t" << sigmaerr[1] << "\t" << sigmaerr[2] << endl;
}
