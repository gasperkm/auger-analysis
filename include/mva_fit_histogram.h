#ifndef _MVA_FIT_HISTOGRAM_H_
#define _MVA_FIT_HISTOGRAM_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"
#include "primary_type.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#define MAXINFILES 40

class MvaFitHist
{
private:

   int ret;

   string *stemp;
   int *itemp;
   double *dtemp;

   string *filename;
   int *nrelem;
   int *nrbins;
   int *nrkeys;
   double *norm;
   int *nentries;
   double *xlim;
   double *ylim;
   double *yrange, *yresidrange;
   bool *otherSettings;
   string *method;
   int *himodel;
   int *simProd;
   int *dataNum;
   string *dataName;
   double *mvalim;
   vector<int> *includePart;

   TH1F *result[MAXINFILES];
   TH1F *simhist[MAXINFILES];
   TH1F *sims[MAXINFILES];
   TH1F *data, *datacomp;
   TH1F *sim1, *sim2, *sim;
   TH1F *newdata, *newsig;

   vector<int> *trees;
   vector<string> *treeNames;

   TFile *ifile;

   RootStyle *mystyle;
   PrimPart *ptype;

   double *mvacut;
   int *nrsim;
   int *nrparam;
   bool *iszero;
   bool *hasmiddist;

   TObjArray *mc;
   TFractionFitter *fracFitter;
   TVirtualFitter *fracVirtFitter;

   // Composition and lnA values
   double *midLna;
   vector<double> *midComposition;
   vector<double> *midCompositionErr;
   double *chi2value;
   double *pvalue;
   int *ndfvalue;
   double *midEnergy;
   double *midStep;

   // Variables for chi2 test
   double *residVal;
   int *chiGood;

   // Plotting a number of files on a single plot
   double *xval;
   double *xerr;
   double *yval;
   double *yerrhi;
   double *yerrlo;
   int *compnrp;
   vector<double> *composition;

   TPad *pads[10];

public:
   MvaFitHist();
   virtual ~MvaFitHist();

   void SetPrimaries(vector<int> *primVals);
   void SetMethod(string *inMethod);
   void SetHImodel(int *inModel);
   void SetSimProduction(int *inProd);
   void SetData(int *inTree, string *inTreeName);
   void SetNrBins(int *inNrBins);
   void SetXaxisLimits(bool use, double *inLimits);
   void SetYaxisLimits(bool use, double *inLimits);
   void SetOtherSettings(bool *inConst);

   void ResetHistograms();

   double fcnConstrainedFunc(const double *par);
   double fcnFunc(const double *par);
   double fcn(const double *par);

   void SetInputFile(string *inFile);
   void PrepareHistograms(int run);
   int StartFitting();
   void PlotSumResiduals(double *sigFrac, double *sigFracErr, int status, double selStep);

   double GetChi2();
   double GetPvalue();
   int GetNdf();
   double GetStep();

   double GetFinalLna(int me);
   double GetFinalComposition(int type);
   double GetFinalCompositionErr(int type, int hilo);
   double GetFinalEnergy();

   void PrintoutFinalResults(int *nrp, vector<double> *lnas, vector<double> *chi2, vector<double> *pvalues, vector<double> *ndfs, vector<double> *comps, vector<double> *energies, vector<double> *steps);
   void PlotLnaComposition();

   int ReadLnaResultsFD(vector<float> *val);
   int ReadLnaResultsSD(vector<float> *val, int array);
   int ReadFracResultsFD(vector<float> *val, bool uselna);
};

#endif
