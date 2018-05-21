{
   int side = -1;
   double elemVal;

   while(side < 2)
   {
      cout << "Select size of the matrix: ";
      cin >> side;
   }

   TMatrixD *m = new TMatrixD(side,side);

   cout << "Enter matrix elements (row,col):" << endl;
   for(int row = 0; row < side; row++)
   {
      for(int col = 0; col < side; col++)
      {
         cout << "(" << row << "," << col << ") = ";
	 cin >> elemVal;
	 m(row,col) = elemVal;
      }
   }
   cout << endl;
   for(int row = 0; row < side; row++)
   {
      for(int col = 0; col < side; col++)
         cout << m(row,col) << " ";
      cout << endl;
   }
   cout << "Determinant = " << m->Determinant() << endl << endl;

   TMatrixDEigen *eigenMat = new TMatrixDEigen((const TMatrixD)m);

   TMatrixD eigenVal = eigenMat->GetEigenValues();
   for(int row = 0; row < side; row++)
   {
      for(int col = 0; col < side; col++)
         cout << eigenVal(row,col) << " ";
      cout << endl;
   }
   cout << endl;

   TMatrixD &eigenVec = eigenMat->GetEigenVectors();
   for(int row = 0; row < side; row++)
   {
      for(int col = 0; col < side; col++)
         cout << eigenVec(row,col) << " ";
      cout << endl;
   }
   cout << endl;
}
