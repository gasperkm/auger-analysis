void PrintMatrix(TMatrixD *mat, string matname)
{
   cout.precision(5);

   cout << "------" << endl;
   cout << matname << ":" << endl;
   for(int i = 0; i < mat->GetNrows(); i++)
   {
      for(int j = 0; j < mat->GetNcols(); j++)
      {
         cout << fixed << (*mat)(i,j);
         
         if(j+1 != mat->GetNcols())
            cout << " | ";
         else
            cout << endl;
      }
   }
   cout << "------" << endl;
   cout << endl;
}

void testing()
{
   TMatrixD *Amat = new TMatrixD(4, 4);
   TMatrixD *AmatInv = new TMatrixD(4, 4);
   (*Amat)(0,0) = 1.01350;
   (*Amat)(0,1) = 0.91085;
   (*Amat)(0,2) = 0.49166;
   (*Amat)(0,3) = 0.00000;
   (*Amat)(1,0) = 1.06887;
   (*Amat)(1,1) = 1.00490;
   (*Amat)(1,2) = 0.58825;
   (*Amat)(1,3) = 0.00000;
   (*Amat)(2,0) = 1.42732;
   (*Amat)(2,1) = 1.36491;
   (*Amat)(2,2) = 1.00992;
   (*Amat)(2,3) = 0.00000;
   (*Amat)(3,0) = 1.00000;
   (*Amat)(3,1) = 1.00000;
   (*Amat)(3,2) = 1.00000;
   (*Amat)(3,3) = 1.00000;

   PrintMatrix(Amat, "Amat");

   *AmatInv = *Amat;
   AmatInv->Invert();

   PrintMatrix(AmatInv, "AmatInv");

   TMatrixD *Ymat = new TMatrixD(4, 1);
   (*Ymat)(0,0) = 0.82498 - (-0.00056);
   (*Ymat)(1,0) = 0.89049 - 0.00433;
   (*Ymat)(2,0) = 1.19810 - 0.00672;
   (*Ymat)(3,0) = 1.00000;

   PrintMatrix(Ymat, "Ymat");

   TMatrixD *Fmat = new TMatrixD(4, 1);
   Fmat->Mult(*AmatInv, *Ymat);

   PrintMatrix(Fmat, "Fmat");
}
