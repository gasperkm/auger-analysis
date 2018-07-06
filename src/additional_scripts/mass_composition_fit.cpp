#define _STANDALONE_ 1
#include "separate_functions.h"
#include "mass_composition_fit.h"
#include <iomanip>

using namespace std;

MassComposition::MassComposition(vector<string> *primVals, char *infile, int intype)
{
   nrall = 0;

   dtemp = new double[3];

   ptype = new PrimPart();

   includeElem = new vector<int>;
   includePart = new vector<int>;
   mvafraction = new vector<double>;
   fractionPoints = new vector<double>;

   midLna = new double;
   midFraction = new vector<double>;
   midComposition = new vector<double>;
/*   midFracNerr = new vector<double>;
   midFracPerr = new vector<double>;*/

   // Define the particles that will be used in the mix
   for(int i = 0; i < primVals->size(); i++)
      IncludePrimaryType(primVals->at(i));

   // Prepare y_mu and y_epsilon matrices with fit values from fraction/lnA lines
   PrepareFitValues();

   CalculateActivatedTypes();

   // Get the proton fraction values from the file
   ReadFromFile(infile, intype);
}

MassComposition::~MassComposition()
{
   delete ptype;
   delete includeElem;
   delete includePart;

   delete ymu;
   delete yeps;

   delete[] dtemp;

   delete midLna;
   delete midFraction;
   delete midComposition;
/*   delete midFracNerr;
   delete midFracPerr;*/

   delete mvafraction;
   delete fractionPoints;

   delete permuteMnp;
}

void MassComposition::PreparePermutations()
{
   // Make permutations of all mean, negative and positive values (for uncertainty calculation)
   permuteMnp = new vector<string>;
   string *stemp = new string;
   *stemp = "MNP";

   vector<int> *permIndex = new vector<int>(stemp->size());
//   cout << "nrtreat = " << nrtreat << ", nrall = " << nrall << ", nrmix = " << nrmix << endl;
   icnt = 0;
//   cout << "icnt = " << icnt << endl;
   jcnt = 0;
   startMnp = icnt;

   PermuteMNP(stemp, permIndex, icnt, jcnt);

//   cout << endl << "All permutations:" << endl;
   for(int i = 0; i < permuteMnp->size(); i++)
   {
      *stemp = permuteMnp->at(i);
      Reverse(stemp);
      permuteMnp->at(i) = *stemp;
//      cout << "i = " << i << ": " << permuteMnp->at(i) << endl;
   }
   cout << endl;

   nrMnp = permuteMnp->size();

   delete stemp;
   delete permIndex;
}

void MassComposition::PermuteMNP(string *st, vector<int> *permIndex, int depth, int count)
{
   string *stemp = new string;

   if(depth == st->size())
   {
      ++depth;
      *stemp = "";
      for(int i = startMnp; i < st->size(); ++i)
         *stemp += st->at(permIndex->at(i));
      permuteMnp->push_back(*stemp);

      return;
   }

   for(int i = 0; i < st->size(); ++i)
   {
      permIndex->at(depth) = i;
      PermuteMNP(st, permIndex, depth+1, count);
   }

   delete stemp;
}

int MassComposition::GetNrMnp()
{
   return nrMnp;
}

string MassComposition::GetMnpPermut(int mnp)
{
   return permuteMnp->at(mnp);
}

void MassComposition::ReadFromFile(char *infile, int intype)
{
   if(intype == 0)
   {
      ifstream *ifs = new ifstream;
      ifs->open(infile, ifstream::in);

      icnt = 0;
      if(ifs->is_open())
      {
         while(1)
         {
            *ifs >> dtemp[0] >> dtemp[1] >> dtemp[2];
            fractionPoints->push_back(dtemp[0]);
            fractionPoints->push_back(dtemp[1]);
            fractionPoints->push_back(dtemp[2]);

	    icnt++;

            ifs->ignore(1,' ');
            if(ifs->eof()) break;
         }
      }

      nrtreat = icnt;
      nrpoints = 1;

      ifs->close();
      delete ifs;
   }
   else if(intype == 1)
   {
      string *stemp = new string[3];
      vector<double> *normval = 0;
      TBranch *normbranch = 0;

      TFile *tfile = TFile::Open(infile, "READ");
      TTree *tree;
      for(int i = 0; i < tfile->GetNkeys(); i++)
      {
         stemp[0] = "fracTree_" + ptype->GetShortName(includePart->at(i)) + "-Fe_treatment";
         tfile->GetObject(stemp[0].c_str(), tree);
	 tree->SetBranchAddress("normFraction", &normval, &normbranch);

	 for(int j = 0; j < tree->GetEntries(); j++)
	 {
            tree->GetEntry(j);
	    normbranch->GetEntry(tree->LoadTree(j));

	    nrpoints = normval->size()/3.;

	    for(int k = 0; k < normval->size(); k++)
               fractionPoints->push_back(normval->at(k));
	 }
      }

      nrtreat = tfile->GetNkeys();

      tfile->Close();

      delete[] stemp;
   }
}

// Include a particle type into the mix (or all possible particle types)
void MassComposition::IncludePrimaryType(string type)
{
   includeElem->push_back(1);
   includePart->push_back(ptype->GetZ(&type));
   nrall++;
}

void MassComposition::ActivatePrimaryType(int type)
{
   if(type >= nrall)
      cout << "Error! Can't activate primary above included range (" << nrall-1 << ")" << endl;
   else
      includeElem->at(type) = 1;

   CalculateActivatedTypes();
}

void MassComposition::DeactivatePrimaryType(int type)
{
   if(nrmix < 3)
   {
      cout << "Error! At least 2 elements need to be activated." << endl;
      return;
   }

   if(type == nrall-1)
   {
      cout << "Error! Can't deactivate heaviest element in the mix, because it sets a constraint for sum of all individual element fractions." << endl;
      return;
   }

   if(type >= nrall)
   {
      cout << "Error! Can't deactivate primary above included range (" << nrall-1 << ")" << endl;
      return;
   }

   includeElem->at(type) = 0;
   CalculateActivatedTypes();
}

void MassComposition::CalculateActivatedTypes()
{
   nrmix = 0;
   for(int i = 0; i < nrall; i++)
   {
      if(includeElem->at(i) == 1)
         nrmix++;
   }

//   cout << "nrmix = " << nrmix << endl;
}

int MassComposition::CheckNrYvalues()
{
   if(mvafraction->size() < nrmix-1)
   {
      cout << "Not enough Y values in input file (" << mvafraction->size() << "/" << nrmix-1 << "). Removing heaviest element from analysis." << endl;
//      DeactivatePrimaryType(nrmix-1);
//      CalculateActivatedTypes();
      return 1;
   }
   else
      return -1;
}

void MassComposition::GetYfitValues(int yfits)
{
   char *ctemp = new char[1024];
   string *stemp = new string;
   vector<double> *tempvec = new vector<double>;
   tempvec->clear();

   ifstream *ifs = new ifstream;
   if(yfits == 0)
      *stemp = string(rootdir) + "/input/mass_composition_" + ptype->GetShortName(includePart->at(nrall-1)) + "-ymu.dat";
   else if(yfits == 1)
      *stemp = string(rootdir) + "/input/mass_composition_" + ptype->GetShortName(includePart->at(nrall-1)) + "-yeps.dat";

//   cout << "Yfit values file: " << *stemp << endl;
   ifs->open(stemp->c_str(), ifstream::in);

   if(ifs->is_open())
   {
      while(1)
      {
         if(ifs->peek() == '#')
            ifs->getline(ctemp, 1024, '\n');
	 else
	 {
            *ifs >> dtemp[0];
	    tempvec->push_back(dtemp[0]);

	    ifs->ignore(1,' ');
            if(ifs->eof()) break;
	 }
      }
   }
   else
      cout << "Error! Yfit values file " << *stemp << " does not exist." << endl;

   icnt = (int)TMath::Sqrt(tempvec->size());
   for(int i = 0; i < icnt; i++)
   {
      for(int j = 0; j < icnt; j++)
      {
         if(yfits == 0)
            (*ymu)(i,j) = tempvec->at(icnt*i + j);
	 else if(yfits == 1)
            (*yeps)(i,j) = tempvec->at(icnt*i + j);
      }
   }

   ifs->close();
   delete ifs;
   delete tempvec;
   delete[] ctemp;
   delete stemp;
}

void MassComposition::PrepareFitValues()
{
   // Explanation:
   //   Rows in matrices means different sig-bgd treatments (ex. p/Fe, He/Fe and O/Fe).
   //   Columns in matrices means a different mock data mixture (ex. p-Fe, He-Fe, O-Fe).
   //   y_mu is the proton fraction of the lightest element in the mock data mixture
   //   y_epsilon is the proton fraction of the heaviest element in the mock data mixture
   //   Values for iron are an additional constraint: sum(f_i) = 1
   ymu = new TMatrixD(nrall, nrall);
   yeps = new TMatrixD(nrall, nrall);

   GetYfitValues(0);
   GetYfitValues(1);

/*   // p/Fe treatment
   (*ymu)(0,0)  = 1.01183;
   (*ymu)(0,1)  = 0.91139;
   (*ymu)(0,2)  = 0.49110;
   (*ymu)(0,3)  = 0.;
   (*yeps)(0,0) = -0.00167;
   (*yeps)(0,1) = 0.00054;
   (*yeps)(0,2) = -0.00056;
   (*yeps)(0,3) = 0.;
   // He/Fe treatment
   (*ymu)(1,0)  = 1.07430;
   (*ymu)(1,1)  = 1.00622;
   (*ymu)(1,2)  = 0.59449;
   (*ymu)(1,3)  = 0.;
   (*yeps)(1,0) = 0.00543;
   (*yeps)(1,1) = 0.00132;
   (*yeps)(1,2) = 0.00624;
   (*yeps)(1,3) = 0.;
   // O/Fe treatment
   (*ymu)(2,0)  = 1.43579;
   (*ymu)(2,1)  = 1.37353;
   (*ymu)(2,2)  = 1.01299;
   (*ymu)(2,3)  = 0.;
   (*yeps)(2,0) = 0.00847;
   (*yeps)(2,1) = 0.00862;
   (*yeps)(2,2) = 0.00307;
   (*yeps)(2,3) = 0.;
   // Iron constraint (f_p + f_He + f_O + f_Fe = 1)
   (*ymu)(3,0)  = 1.;
   (*ymu)(3,1)  = 1.;
   (*ymu)(3,2)  = 1.;
   (*ymu)(3,3)  = 1.;
   (*yeps)(3,0) = 0.;
   (*yeps)(3,1) = 0.;
   (*yeps)(3,2) = 0.;
   (*yeps)(3,3) = 0.;*/

/*   PrintMatrix(ymu, "y_mu");
   PrintMatrix(yeps, "y_epsilon");*/
}

void MassComposition::PrintMixture()
{
   cout << "Elements in mixture: ";

   icnt = 0;
   for(int i = 0; i < nrall; i++)
   {
      if(includeElem->at(i) == 1)
      {
         if(i != nrall-1)
            cout << ptype->GetName(includePart->at(i)) << ", ";
	 else
            cout << ptype->GetName(includePart->at(i)) << endl;
         icnt++;
      }
   }

   cout << endl;
}

void MassComposition::PrepareAmatrix()
{
   // Explanation: A linear system of equations can be solved by using A*F = Y, where A is the matrix with (y_mu - y_epsilon) values, F is the matrix with individual element fractions and Y is the matrix of proton fraction values from MVA analysis reduced by y_epsilon values for all used mixtures

   Amat = new TMatrixD(nrmix, nrmix);
   ymuSub = new TMatrixD(nrmix, nrmix);
   yepsSub = new TMatrixD(nrmix, nrmix);

   icnt = 0;
   for(int i = 0; i < nrall; i++) // go over rows (sig-bgd treatments) -> only take the first nrmix-1 rows + final row (constraint)
   {
      jcnt = 0;
      for(int j = 0; j < nrall; j++)	// go over cols (mock data mixtures)
      {
         if( includeElem->at(j) == 1 )
	 {
            if( (i < nrmix-1) || (i == nrall-1) )
	    {
               (*Amat)(icnt,jcnt) = (*ymu)(i,j) - (*yeps)(i,j);
               (*ymuSub)(icnt,jcnt) = (*ymu)(i,j);
               (*yepsSub)(icnt,jcnt) = (*yeps)(i,j);
	    }
	 }
         
         if( includeElem->at(j) == 1 )
            jcnt++;
      }

      if(i < nrmix-1)
         icnt++;
   }

   PrintMatrix(ymuSub, "y_mu");
   PrintMatrix(yepsSub, "y_epsilon");
   PrintMatrix(Amat, "Amat");
}

void MassComposition::InvertAmatrix()
{
   // Explanation: To get individual element fractions, we need to invert the A matrix for calculation F = A^(-1)*Y
   
   AmatInv = new TMatrixD(nrmix, nrmix);

   *AmatInv = *Amat;
   AmatInv->Invert();

   PrintMatrix(AmatInv, "AmatInv");
}

void MassComposition::PrepareYmatrix()
{
   Ymat = new TMatrixD(nrmix, 1);
   Shiftmat = new TMatrixD(nrmix-1, 1);
   Deltamat = new TMatrixD(nrmix-1, nrmix-1);
   
   icnt = 0;
   for(int i = 0; i < nrmix-1; i++)
   {
      (*Ymat)(i, 0) = mvafraction->at(i);
      icnt++;
   }
   (*Ymat)(icnt, 0) = 1.;

   PrintMatrix(Ymat, "Ymat (MVA analysis results)");

   for(int i = 0; i < nrmix-1; i++)
   {
      dtemp[0] = 0;
      for(int j = 0; j < nrmix-1; j++)
         dtemp[0] += (*yepsSub)(i,j);

      (*Shiftmat)(i, 0) = dtemp[0]/(double)(nrmix-1);

      dtemp[1] = 0;
      for(int j = 0; j < nrmix-1; j++)
      {
         dtemp[0] = (*yepsSub)(i, j) - (*Shiftmat)(i, 0);
         (*Deltamat)(i, j) = dtemp[0];
	 dtemp[1] += dtemp[0];
      }

      dtemp[2] = dtemp[1] + (*Shiftmat)(i, 0);
      (*Ymat)(i, 0) -= dtemp[2];
   }

   PrintMatrix(Shiftmat, "Shiftmat");
   PrintMatrix(Deltamat, "Deltamat");
   PrintMatrix(Ymat, "Ymat (modified by y_epsilon)");
}

void MassComposition::CalculateFractions()
{
   // Explanation: Using F = A^(-1)*Y, calculate the individual element fractions
   Fmat = new TMatrixD(nrmix, 1);

   Fmat->Mult(*AmatInv, *Ymat);

   PrintMatrix(Fmat, "Fmat");
}

void MassComposition::SaveResults(int mnp)
{
   if((mnp >= 0) && (mnp < 27))
   {
      midFraction->clear();
      midComposition->clear();
   }
   else
      return;

   // Save fractions
   for(int i = 0; i < nrmix-1; i++)
   {
      cout << ptype->GetShortName(includePart->at(i)) << "/Fe treatment: " << endl;

      dtemp[0] = 0;

      jcnt = 0;
      for(int j = 0; j < nrall; j++)
      {
         if( (includeElem->at(j) == 1) && (j != nrall-1) )
         {
            dtemp[1] = (*Amat)(i,jcnt)*(*Fmat)(jcnt,0) + (*Deltamat)(i,jcnt);
            cout << "  y (" << ptype->GetName(includePart->at(j)) << ") = " << dtemp[1] << endl;
            dtemp[0] += dtemp[1];

            jcnt++;
         }
      }

      cout << "  y (shift)  = " << (*Shiftmat)(i, 0) << endl;
      dtemp[1] = dtemp[0] + (*Shiftmat)(i, 0);

      cout << "Saving " << permuteMnp->at(mnp) << " value (fraction)" << endl;
      midFraction->push_back(dtemp[1]);
/*      // Negative error
      else if(mnp == 1)
      {
         cout << "Saving negative error (fraction)" << endl;
         midFracNerr->push_back(dtemp[1]);
//         finalFracNerr->push_back(TMath::Abs(finalFraction->at(i) - dtemp[1]));
      }
      // Positive error
      else if(mnp == 2)
      {
         cout << "Saving positive error (fraction)" << endl;
         midFracPerr->push_back(dtemp[1]);
//         finalFracPerr->push_back(TMath::Abs(finalFraction->at(i) - dtemp[1]));
      }*/
   }

   // Save lnA
   cout << "lnA = ";
   dtemp[1] = 0;
   icnt = 0;
   for(int i = 0; i < nrall; i++)
   {
      if(includeElem->at(i) == 1)
      {
         dtemp[0] = (*Fmat)(icnt, 0)*TMath::Log(ptype->GetA(includePart->at(i)));
	 if(i != nrall-1)
	    cout << (*Fmat)(icnt, 0) << "*ln(" << ptype->GetA(includePart->at(i)) << ") + ";
	 else
	    cout << (*Fmat)(icnt, 0) << "*ln(" << ptype->GetA(includePart->at(i)) << ") = ";
	 dtemp[1] += dtemp[0];
         icnt++;
      }
   }

   cout << dtemp[1] << endl << endl;

   cout << "Saving " << permuteMnp->at(mnp) << " value (lnA)" << endl;
   *midLna = dtemp[1];
/*   // Negative error
   else if(mnp == 1)
   {
      cout << "Saving negative error (lnA)" << endl;
//      finalLna->push_back(TMath::Abs(finalLna->at(0) - dtemp[1]));
      *midLna = dtemp[1];
   }
   // Positive error
   else if(mnp == 2)
   {
      cout << "Saving positive error (lnA)" << endl;
//      finalLna->push_back(TMath::Abs(finalLna->at(0) - dtemp[1]));
      *midLna = dtemp[1];
   }*/

   // Print out individual element fractions (mass composition)
   cout << "Mass composition:" << endl;
   icnt = 0;
   cout.precision(2);
   for(int i = 0; i < nrall; i++)
   {
      if(includeElem->at(i) == 1)
      {
         cout << "  " << ptype->GetName(includePart->at(i)) << " =\t" << 100.*(*Fmat)(icnt, 0) << "%" <<  endl;
	 midComposition->push_back((*Fmat)(icnt, 0));
         icnt++;
      }
      else
      {
	 midComposition->push_back(-404);
      }
   }

   if(nrall == 3)
      midComposition->push_back(-404);
   if(nrall == 2)
      midComposition->push_back(-404);

   cout << "midComposition (" << midComposition->size() << "):" << endl;
   for(int i = 0; i < midComposition->size(); i++)
      cout << "  " << midComposition->at(i) << endl;
}

int MassComposition::CheckResults()
{
   icnt = 0;
   jcnt = 0;
   dtemp[1] = 1.;

   for(int i = 0; i < nrall; i++)
   {
      if(includeElem->at(i) == 1)
      {
	 dtemp[0] = (*Fmat)(icnt, 0);
         if(dtemp[0] < dtemp[1])
            dtemp[1] = dtemp[0];

         if(dtemp[0] <= 0.)
         {
	    // Trying to deactivate heaviest element
	    if(i == nrall-1)
            {
	       if(jcnt == 0)
	       {
                  cout << "Element fraction for " << ptype->GetName(includePart->at(i)) << " is negative. Removing it from next analysis." << endl;
                  DeactivatePrimaryType(i);
                  cout << endl << "Finished analysis (with errors)" << endl;
                  cout << "---------------------------------------------------" << endl << endl;
                  return 1;
	       }
	    }
	    // Deactivating lighter elements
	    else
	    {
               cout << "Element fraction for " << ptype->GetName(includePart->at(i)) << " is negative. Removing it from next analysis." << endl;
               DeactivatePrimaryType(i);
	       jcnt++;
	    }
         }

	 icnt++;
      }
   }

   if(dtemp[1] < 0)
   {
      cout << endl << "Rerunning analysis" << endl;
      cout << "---------------------------------------------------" << endl << endl;
      return -1;
   }
   else
   {
      cout << endl << "Finished analysis" << endl;
      cout << "---------------------------------------------------" << endl << endl;
      return 0;
   }
}

void MassComposition::DrawResults(int type, int mnp)
{
   // Draw results on a canvas for a p-Fe treatment
   RootStyle *mystyle = new RootStyle();
   mystyle->SetBaseStyle();
   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TF1 *func[3];
   string *stemp = new string[2];

   gStyle->SetEndErrorSize(5);

   // Plot the points of pure mixture (100%)
   TGraph *grpure;
   double *lnpure = new double[nrall];
   double *ypure = new double[nrall];

   lnpure[0] = TMath::Log(ptype->GetA(includePart->at(nrall-1)));
   ypure[0] = (*Shiftmat)(type, 0);

   // Draw mixture functions
   for(int j = 0; j < nrall-1; j++)
   {
      stemp[0] = ptype->GetShortName(includePart->at(j)) + "-" + ptype->GetShortName(includePart->at(nrall-1)) + "_func";
      stemp[1] = "((" + ToString((*ymu)(type, j), 5) + ") - (" + ToString((*yeps)(type, j), 5) + "))*((x - " + ToString(TMath::Log(ptype->GetA(includePart->at(nrall-1))), 5) + ")/(" + ToString(TMath::Log(ptype->GetA(includePart->at(j))), 5) + " - " + ToString(TMath::Log(ptype->GetA(includePart->at(nrall-1))), 5) + ")) + (" + ToString((*yeps)(type, j), 5) + ")";
//      cout << stemp[0] << endl;
//      cout << stemp[1] << endl;

      func[j] = new TF1(stemp[0].c_str(), stemp[1].c_str(), -0.1, 4.5);

      lnpure[j+1] = TMath::Log(ptype->GetA(includePart->at(j)));
      ypure[j+1] = func[j]->Eval(TMath::Log(ptype->GetA(includePart->at(j))));

      func[j]->SetLineWidth(1);
      func[j]->SetLineColor(1);
   }

   grpure = new TGraph(nrall, lnpure, ypure);
   mystyle->SetGraphColor(grpure, 1);
   mystyle->SetAxisTitles(grpure, "<lnA>", "MVA fraction");
   grpure->GetXaxis()->SetRangeUser(-0.1, 4.5);
   grpure->GetYaxis()->SetRangeUser(-0.1, 1.6);
   grpure->GetXaxis()->SetLimits(-0.1, 4.5);
   grpure->GetYaxis()->SetLimits(-0.1, 1.6);
   grpure->Draw("AP");

   for(int j = 0; j < nrall-1; j++)
      func[j]->Draw("SAME");

   // Draw lines where composition is supposed to end
//   if(mnp == 0)
      dtemp[0] = GetFinalFraction(type);
/*   else if(mnp == 1)
      dtemp[0] = GetFinalFracNerr(type);
   else if(mnp == 2)
      dtemp[0] = GetFinalFracPerr(type);*/

   TLine *l1 = new TLine(GetFinalLna(), -0.1, GetFinalLna(), 1.6);
   l1->SetLineWidth(1);
   l1->SetLineColor(6);
   l1->SetLineStyle(1);
   l1->Draw("SAME");
   TLine *l2 = new TLine(-0.1, dtemp[0], 4.5, dtemp[0]);
   l2->SetLineWidth(1);
   l2->SetLineColor(6);
   l2->SetLineStyle(1);
   l2->Draw("SAME");

   TLatex *lineText = new TLatex();
   lineText->SetTextAlign(13);
   lineText->SetTextColor(6);
   stemp[0] = ToString(GetFinalLna(), 4);
   lineText->DrawLatex(GetFinalLna()+0.02, 1.55, stemp[0].c_str());

   lineText->SetTextAlign(31);
   lineText->SetTextColor(6);
   stemp[0] = ToString(dtemp[0], 4);
   lineText->DrawLatex(4.4, dtemp[0]+0.01, stemp[0].c_str());

   lineText->SetTextAlign(11);
   lineText->SetTextColor(1);
   icnt = 0;
   for(int i = 0; i < nrall; i++)
   {
      if(includeElem->at(i) == 1)
      {
         stemp[0] = ptype->GetName(includePart->at(i)) + " = " + ToString(100.*(*Fmat)(icnt, 0), 2) + "%";
         lineText->DrawLatex(3.51, (1.45-icnt*0.05), stemp[0].c_str());
         icnt++;
      }
   }

   // Plot the mass composition of data
   TGraph *gr;
   double *lnval = new double[nrall];
   double *yval = new double[nrall];

   lnval[0] = TMath::Log(ptype->GetA(includePart->at(nrall-1)));
   yval[0] = (*Shiftmat)(type, 0);
//   cout << lnval[0] << "\t" << yval[0] << endl;

   dtemp[0] = lnval[0];
   dtemp[1] = yval[0];

   jcnt = 0;
   for(int j = 0; j < nrall; j++)
   {
      if( (includeElem->at(j) == 1) && (j != nrall-1) )
      {
         lnval[jcnt+1] = dtemp[0] - (*Fmat)(jcnt,0)*(TMath::Log(ptype->GetA(includePart->at(nrall-1))) - TMath::Log(ptype->GetA(includePart->at(j))));

	 stemp[0] = ToString(dtemp[0], 5) + " - " + ToString((*Fmat)(jcnt,0), 5) + "*(" + ToString(TMath::Log(ptype->GetA(includePart->at(nrall-1))), 5) + " - " + ToString(TMath::Log(ptype->GetA(includePart->at(j))), 5) + ")";
//	 cout << stemp[0] << endl;

	 dtemp[0] = lnval[jcnt+1];

	 yval[jcnt+1] = dtemp[1] + (*Amat)(type,jcnt)*(*Fmat)(jcnt,0) + (*Deltamat)(type,jcnt);
	 dtemp[1] = yval[jcnt+1];

//         cout << lnval[jcnt+1] << "\t" << yval[jcnt+1] << endl;

         jcnt++;
      }
   }

   gr = new TGraph(nrmix, lnval, yval);
   mystyle->SetGraphColor(gr, 0);
   gr->Draw("LP;SAME");
/*
   stemp[0] = "mkdir -p " + string(rootdir) + "/results/mass_composition_plots";
   system(stemp[0].c_str());

   stemp[0] = "rm -fr " + string(rootdir) + "/results/mass_composition_plots/*";
   system(stemp[0].c_str());
*/
   if(type == 0)
      stemp[0] = string(rootdir) + "/results/mass_composition_plots/lnA_fraction_sig-" + ptype->GetShortName(includePart->at(nrall-1)) + "_p-Fe-treatment";
   else if(type == 1)
      stemp[0] = string(rootdir) + "/results/mass_composition_plots/lnA_fraction_sig-" + ptype->GetShortName(includePart->at(nrall-1)) + "_He-Fe-treatment";
   else if(type == 2)
      stemp[0] = string(rootdir) + "/results/mass_composition_plots/lnA_fraction_sig-" + ptype->GetShortName(includePart->at(nrall-1)) + "_O-Fe-treatment";

/*   if(mnp == 0)
      stemp[0] += "_mean";
   else if(mnp == 1)
      stemp[0] += "_nerr";
   else if(mnp == 2)
      stemp[0] += "_perr";*/

   stemp[0] += "_" + permuteMnp->at(mnp);

   stemp[0] += ".pdf";
   cout << stemp[0] << endl;
   c1->SaveAs(stemp[0].c_str());

   for(int j = 0; j < nrall-1; j++)
      delete func[j];

   delete l1;
   delete l2;
   delete gr;
   delete lineText;
   delete mystyle;
   delete c1;
   delete[] stemp;

   delete[] lnval;
   delete[] yval;

   delete[] lnpure;
   delete[] ypure;
}

double MassComposition::GetFinalLna()
{
   return *midLna;
}

double MassComposition::GetFinalFraction(int type)
{
/*   cout << "nrfrac = " << midFraction->size() << endl;
   for(int i = 0; i < midFraction->size(); i++)
      cout << "i = " << i << ": " << midFraction->at(i) << " (" << midFraction->at(type) << ")" << endl;*/
   return midFraction->at(type);
}

double MassComposition::GetFinalComposition(int type)
{
   return midComposition->at(type);
}

/*double MassComposition::GetFinalFracNerr(int type)
{
   return midFracNerr->at(type);
}

double MassComposition::GetFinalFracPerr(int type)
{
   return midFracPerr->at(type);
}*/

int MassComposition::GetNrMix()
{
   return nrmix;
}

int MassComposition::GetNrTreat()
{
   return nrtreat;
}

int MassComposition::GetNrPoints()
{
   return nrpoints;
}

int MassComposition::SelectPoint(int point, int mnp)
{
   if( point >= nrpoints )
   {
      cout << "Error! Point (" << point << ") does not exist. Highest poins is " << nrpoints-1 <<  ". Rerun analysis and select the correct point." << endl;
      point = nrpoints-1;

      return -1;
   }

   int *itemp = new int[nrtreat];
   string *stemp = new string;

   *stemp = permuteMnp->at(mnp);

   cout.precision(5);
   cout << "Selecting point " << point << " and type " << mnp << " (" << *stemp << ")." << endl;

   mvafraction->clear();

   for(int i = 0; i < nrtreat; i++)
   {
      if( stemp->at(i) == 'M' )
      {
         itemp[i] = 0;
//	 cout << "i = " << i << ": mean" << endl;
      }
      else if( stemp->at(i) == 'N' )
      {
         itemp[i] = -1;
//	 cout << "i = " << i << ": neg" << endl;
      }
      else if( stemp->at(i) == 'P' )
      {
         itemp[i] = 1;
//	 cout << "i = " << i << ": pos" << endl;
      }
   }

   // Take values for treatments treatment
   for(int i = 0; i < nrtreat; i++)
   {
      // Mean value
      if(itemp[i] == 0)
         mvafraction->push_back(fractionPoints->at(3*(point+i*nrpoints)));
      // Negative error
      else if(itemp[i] == -1)
         mvafraction->push_back(fractionPoints->at(3*(point+i*nrpoints)) - fractionPoints->at(3*(point+i*nrpoints)+1));
      // Positive error
      else if(itemp[i] == 1)
         mvafraction->push_back(fractionPoints->at(3*(point+i*nrpoints)) + fractionPoints->at(3*(point+i*nrpoints)+2));

/*      // Mean value
      if(mnp == 0)
         mvafraction->push_back(fractionPoints->at(3*(point+i*nrpoints)));
      // Negative error
      else if(mnp == 1)
         mvafraction->push_back(fractionPoints->at(3*(point+i*nrpoints)) - fractionPoints->at(3*(point+i*nrpoints)+1));
      // Positive error
      else if(mnp == 2)
         mvafraction->push_back(fractionPoints->at(3*(point+i*nrpoints)) + fractionPoints->at(3*(point+i*nrpoints)+2));*/

      cout.precision(5);
      cout << "i = " << i << ": Using value = " << mvafraction->at(i) << " (" << fractionPoints->at(3*(point+i*nrpoints)) << ", " << fractionPoints->at(3*(point+i*nrpoints)+1) << ", " << fractionPoints->at(3*(point+i*nrpoints)+2) << ")" << endl;
   }

   delete[] itemp;
   delete stemp;

   return 0;
}

void MassComposition::PrintMatrix(TMatrixD *mat, string matname)
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

void MassComposition::Reverse(string *instring)
{
   string *stemp = new string;
   *stemp = "";

   for(int i = instring->size()-1; i >= 0; i--)
   {
      *stemp += instring->at(i);
   }

   *instring = *stemp;
   delete stemp;
}

void MassComposition::DeleteAnalysisMatrices()
{
   delete Amat;
   delete ymuSub;
   delete yepsSub;
   delete AmatInv;
   delete Ymat;
   delete Fmat;

   delete Shiftmat;
   delete Deltamat;
}

int main(int argc, char **argv)
{
   int invalid[27];
   for(int i = 0; i < 27; i++)
      invalid[i] = -1;
   int *itemp = new int[6];
   string *stemp = new string;
   double *dtemp = new double[3];

   int allnr = 4;

   vector<double> *finalLna = new vector<double>;
   vector<string> *finalPermut = new vector<string>;
   vector<double> *finalFraction = new vector<double>;
   vector<double> *finalFracNerr = new vector<double>;
   vector<double> *finalFracPerr = new vector<double>;
   vector<double> *finalComp = new vector<double>;

   int filetype;
   if(argc > 1)
   {
      *stemp = string(argv[1]);
      itemp[2] = CompareExtension(stemp, "root");
      if(itemp[2] == 1)
         filetype = 1;
      else
         filetype = 0;
   }

   cout << "Select the point to analyze: ";
   cin >> itemp[3];

   itemp[5] = 26;
   
   stemp[0] = "mkdir -p " + string(rootdir) + "/results/mass_composition_plots";
   system(stemp[0].c_str());

   stemp[0] = "rm -fr " + string(rootdir) + "/results/mass_composition_plots/*";
   system(stemp[0].c_str());

   // Go over all mean, negative and positive values
   for(int i = 0; i < 27; i++)
   {
      itemp[0] = 0;
      itemp[1] = 0;

      while(invalid[i] != 0)
      {
         vector<string> *types = new vector<string>;
         types->push_back("Proton");
         types->push_back("Helium");
         if(itemp[1] < 2)
            types->push_back("Oxygen");
         if( itemp[1] < 1 )
            types->push_back("Iron");

         MassComposition *massComp;

         if(argc > 1)
         {
            cout << "Using file " << argv[1] << " as input MVA fraction file." << endl;
            massComp = new MassComposition(types, argv[1], filetype);
	    massComp->PreparePermutations();
            if(massComp->SelectPoint(itemp[3], i) == -1)
	    {
               delete massComp;
               delete[] itemp;
               delete stemp;
               delete[] dtemp;
               delete types;

               delete finalLna;
               delete finalPermut;
               delete finalFraction;
               delete finalComp;
	       return 1;
	    }
            invalid[i] = massComp->CheckNrYvalues();
         }
         else
         {
            cout << "Error! No input file." << endl;
	    delete massComp;
            delete[] itemp;
            delete stemp;
            delete[] dtemp;
            delete types;

            delete finalLna;
            delete finalPermut;
            delete finalFraction;
            delete finalComp;
            return 1;
         }

         while(invalid[i] < 0)
         {
//            massComp->PreparePermutations();
            massComp->PrintMixture();
            massComp->PrepareAmatrix();
            massComp->InvertAmatrix();
            massComp->PrepareYmatrix();
            massComp->CalculateFractions();
	    massComp->SaveResults(i);

	    for(int j = 0; j < (massComp->GetNrMix()-1); j++)
               massComp->DrawResults(j, i);

            invalid[i] = massComp->CheckResults();

            massComp->DeleteAnalysisMatrices();
	 }

	 // Store results in a separate file
	 if(invalid[i] == 0)
	 {
            finalLna->push_back(massComp->GetFinalLna());
	    finalPermut->push_back(massComp->GetMnpPermut(i));

            for(int j = 0; j < massComp->GetNrTreat(); j++)
	    {
               if(j < massComp->GetNrMix()-1)
                  finalFraction->push_back(massComp->GetFinalFraction(j));
	       else
                  finalFraction->push_back(-1.);
/*	       if(i == 0)
                  finalFraction->push_back(massComp->GetFinalFraction(j));*/
/*	       else if(i == 1)
                  finalFracNerr->push_back(massComp->GetFinalFracNerr(j));
	       else if(i == 2)
                  finalFracPerr->push_back(massComp->GetFinalFracPerr(j));*/
	    }

	    for(int j = 0; j < allnr; j++)
               finalComp->push_back(massComp->GetFinalComposition(j));
	 }

         delete massComp;
         delete types;

         cout << "invalid[" << i << "] = " << invalid[i] << endl;
         if( invalid[i] == 1 )
            itemp[1]++;

         if( itemp[0] > 2 )
            break;

         itemp[0]++;

/*         cout << "Continue? ";
         cin >> itemp[4];*/
      }

/*         cout << "Continue? ";
         cin >> itemp[4];*/
   }

/*   cout.precision(4);
   cout << "Final lnA (" << finalLna->size() <<  "):" << endl;
   for(int i = 0; i < finalLna->size(); i++)
      cout << "i = " << i << ", permut = " << finalPermut->at(i) << ": " << finalLna->at(i) << endl;

   cout << "Final fraction (" << finalFraction->size() <<  "):" << endl;
   for(int i = 0; i < finalFraction->size(); i++)
      cout << "i = " << i << " (mean) = " << finalFraction->at(i) << endl;*/

   cout.precision(4);
   cout << endl << "Final results:" << endl;
   cout << "# nr\trun \tlnA \t\tp/Fe\tHe/Fe\tO/Fe\t\tp   \tHe  \tO   \tFe" << endl;
   for(int i = 0; i < finalLna->size(); i++)
   {
      cout << i << "\t" << finalPermut->at(i) << "\t" << finalLna->at(i);
      cout << "\t";

      for(int j = 0; j < allnr-1; j++)
         cout << "\t" << finalFraction->at((allnr-1)*i+j);

      cout << "\t";

      for(int j = 0; j < allnr; j++)
      {
	 if(finalComp->at(allnr*i+j) == -404)
            cout << "\t" << "----";
	 else
            cout << "\t" << finalComp->at(allnr*i+j);
      }

      cout << endl;
   }

   cout.precision(4);
   cout << endl << "Mean value and low/high errors on lnA:" << endl;

   for(int i = 0; i < finalLna->size(); i++)
   {
      if(finalPermut->at(i) == "MMM")
      {
         dtemp[0] = finalLna->at(i);
	 break;
      }
   }

   dtemp[1] = dtemp[0];
   dtemp[2] = dtemp[0];
   for(int i = 0; i < finalLna->size(); i++)
   {
      if(finalLna->at(i) < dtemp[1])
         dtemp[1] = finalLna->at(i);

      if(finalLna->at(i) > dtemp[2])
         dtemp[2] = finalLna->at(i);
   }

   cout << dtemp[0] << "\t" << dtemp[1] << "\t" << dtemp[2] << endl;
   cout << dtemp[0] << "\t" << TMath::Abs(dtemp[0] - dtemp[1]) << "\t" << TMath::Abs(dtemp[0] - dtemp[2]) << endl;

   delete[] itemp;
   delete stemp;
   delete[] dtemp;

   delete finalLna;
   delete finalPermut;
   delete finalFraction;
/*   delete finalFracNerr;
   delete finalFracPerr;*/

   return 0;
}
