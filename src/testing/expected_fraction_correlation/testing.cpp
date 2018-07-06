/*#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>*/

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

   TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
   c1->SetGrid();

   double x1[] = { 18.55, 18.65, 18.75, 18.85, 18.95, 19.05, 19.15, 19.25, 19.35, 19.45, 19.60, 19.80 };
   double xerr[] = { 0., 0., 0., 0. , 0., 0., 0. };

   double y1[] =   { 1.571, 1.1748, 1.2328, 2.3929, 2.5323, 2.4814, 2.8472, 2.5339, 3.0567, 2.7454, 3.4476, 2.9892 };
   double y1le[] = { 0., 0.0775, 0., 0.0461, 0.0399, 0.055, 0.018, 0.0227, 0.0642, 0.0118, 0.1847, 0. };
   double y1he[] = { 0.2034, 0.0707, 0.2051, 0.0648, 0.0413, 0.0317, 0.0748, 0.0228, 0., 0.0611, 0., 0.3587 };

   TGraphAsymmErrors *gr1 = new TGraphAsymmErrors(12, x1, y1, xerr, xerr, y1le, y1he);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerColor(4);
   gr1->SetName("data");
   gr1->GetXaxis()->SetRangeUser(18.4,19.9);
   gr1->GetYaxis()->SetRangeUser(-0.1,4.5);
   gr1->GetXaxis()->SetLimits(18.4,19.9);
   gr1->GetYaxis()->SetLimits(-0.1,4.5);
   gr1->GetXaxis()->SetTitle("log(E/eV)");
   gr1->GetYaxis()->SetTitle("lnA");
   gr1->Draw("APL");

   c1->SaveAs("fit-lnA_data.pdf");

   TCanvas *c2 = new TCanvas("c2", "", 1200, 800);
   c2->SetGrid();

/*
1.5710	0.5834	0.4473	0
1.1748	0.0775	0.4631	1
1.2328	0.0623	0.9087	2
2.3929	0.5102	0.0648	3
2.5323	0.3666	0.0413	4
2.4814	0.4377	0.0317	5
2.8472	0.1796	0.0748	6
2.5339	0.0226	0.0228	7
3.0567	0.2033	0.0638	8
2.7454	0.0118	0.0611	9
3.4476	0.3425	0.0147	10
2.9892	0.0222	0.3858	11
*/

   double y2[] =   { 1.5710, 1.1748, 1.2328, 2.3929, 2.5323, 2.4814, 2.8472, 2.5339, 3.0567, 2.7454, 3.4476, 2.9892 };
   double y2le[] = { 0.5834, 0.0775, 0.0623, 0.5102, 0.3666, 0.4377, 0.1796, 0.0226, 0.2033, 0.0118, 0.3425, 0.0222 };
   double y2he[] = { 0.4473, 0.4631, 0.9087, 0.0648, 0.0413, 0.0317, 0.0748, 0.0228, 0.0638, 0.0611, 0.0147, 0.3858 };

   TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(12, x1, y2, xerr, xerr, y2le, y2he);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerColor(4);
   gr2->SetName("data-werr");
   gr2->GetXaxis()->SetRangeUser(18.4,19.9);
   gr2->GetYaxis()->SetRangeUser(-0.1,4.5);
   gr2->GetXaxis()->SetLimits(18.4,19.9);
   gr2->GetYaxis()->SetLimits(-0.1,4.5);
   gr2->GetXaxis()->SetTitle("log(E/eV)");
   gr2->GetYaxis()->SetTitle("lnA");
   gr2->Draw("APL");

   c2->SaveAs("fit-lnA_data-werr.pdf");

   TCanvas *c3 = new TCanvas("c3", "", 1200, 800);
   c3->SetGrid();

/*
1.4282	0.7753	0.3066
1.1132	0.3148	0.6389
0.8299	0.0558	0.8743
1.5154	0.8996	0.1655
1.7888	0.8544	0.0632
1.3756	0.6458	0.3221
0.7442	0.0329	0.9314
1.3548	0.4961	0.3453
1.0328	0.2711	0.6841
1.1250	0.3829	0.5854
1.0400	0.3400	0.6317
0.6953	0.0343	0.2800
*/

   double y3[] =   { 1.4282, 1.1132, 0.8299, 1.5154, 1.7888, 1.3756, 0.7442, 1.3548, 1.0328, 1.1250, 1.0400, 0.6953 };
   double y3le[] = { 0.7753, 0.3148, 0.0558, 0.8996, 0.8544, 0.6458, 0.0329, 0.4961, 0.2711, 0.3829, 0.3400, 0.0343 };
   double y3he[] = { 0.3066, 0.6389, 0.8743, 0.1655, 0.0632, 0.3221, 0.9314, 0.3453, 0.6841, 0.5854, 0.6317, 0.2800 };

   TGraphAsymmErrors *gr3 = new TGraphAsymmErrors(12, x1, y3, xerr, xerr, y3le, y3he);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerColor(4);
   gr3->SetName("mock_mixture");
   gr3->GetXaxis()->SetRangeUser(18.4,19.9);
   gr3->GetYaxis()->SetRangeUser(-0.1,4.5);
   gr3->GetXaxis()->SetLimits(18.4,19.9);
   gr3->GetYaxis()->SetLimits(-0.1,4.5);
   gr3->GetXaxis()->SetTitle("log(E/eV)");
   gr3->GetYaxis()->SetTitle("lnA");
   gr3->Draw("APL");

   c3->SaveAs("fit-lnA_50p-35He-15Fe.pdf");
}
