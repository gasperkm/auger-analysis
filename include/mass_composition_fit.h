#ifndef _MASS_COMPOSITION_FIT_H_
#define _MASS_COMPOSITION_FIT_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"
#include "primary_type.h"

class MassComposition
{
private:
   int icnt, jcnt;

   // Vectors with all element information
   vector<int> *includeElem;
   vector<int> *includePart;
   vector<double> *mvafraction;

   vector<double> *fractionPoints;
   int nrpoints;
   int nrtreat;

   vector<string> *permuteMnp;
   int startMnp;
   int nrMnp;

   // Number of all possible elements and number of elements included in the mix
   int nrall;
   int nrmix;

   // Primary particle type information
   PrimPart *ptype;

   // Matrices for calculating the element fraction
   TMatrixD *ymu;
   TMatrixD *yeps;
   TMatrixD *Amat;
   TMatrixD *AmatInv;
   TMatrixD *Ymat;
   TMatrixD *Fmat;

   TMatrixD *ymuSub;
   TMatrixD *yepsSub;

   double *dtemp;

   TMatrixD *Shiftmat;
   TMatrixD *Deltamat;

   double *midLna;
   vector<double> *midFraction;
   vector<double> *midComposition;
/*   vector<double> *midFracNerr;
   vector<double> *midFracPerr;*/

public:
   MassComposition(vector<string> *primVals, char *infile, int intype);
   virtual ~MassComposition();

   void PreparePermutations();
   void PermuteMNP(string *st, vector<int> *permIndex, int depth, int count);
   int GetNrMnp();
   string GetMnpPermut(int mnp);
   void ReadFromFile(char *infile, int intype);
   void IncludePrimaryType(string type);
   void ActivatePrimaryType(int type);
   void DeactivatePrimaryType(int type);
   void CalculateActivatedTypes();
   int CheckNrYvalues();
   void GetYfitValues(int yfits);
   void PrepareFitValues();
   void PrintMixture();
   void PrepareAmatrix();
   void InvertAmatrix();
   void PrepareYmatrix();
   void CalculateFractions();
   void PrintResults();
   void SaveResults(int mnp);
   int CheckResults();
   void DrawResults(int type, int mnp);
   double GetFinalLna();
   double GetFinalFraction(int type);
   double GetFinalComposition(int type);
/*   double GetFinalFracNerr(int type);
   double GetFinalFracPerr(int type);*/
   int GetNrMix();
   int GetNrTreat();
   int GetNrPoints();
   int SelectPoint(int point, int mnp);
   void PrintMatrix(TMatrixD *mat, string matname);
   void Reverse(string *instring);
   void DeleteAnalysisMatrices();
};

#endif
