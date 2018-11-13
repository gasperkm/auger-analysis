#ifndef _MVA_FIT_HISTOGRAM_H_
#define _MVA_FIT_HISTOGRAM_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"
#include "primary_type.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class MvaFitHist
{
private:
   // Variables needed for the chi2 calculation
   int *nrsim;
   double *simnorm;
   TH1F *simhist[40];
   TH1F *sims[40];

   int nrsig, nrback, nrdata;
   double signorm, bgdnorm;
   TH1F *result[40];
   TH1F *simsig, *simback, *data, *datacomp;
   TH1F *sim1, *sim2, *sim;
   TH1F *newdata, *newsig, *newback;

   // Number of histogram bins, xrange limits and y ranges
   int nrbins, nrkeys;
   double xlim[2];
   double yrange[2], yresidrange[2];

   // Normalization factor and number of entries
   double *norm;
   int *nentries;

   double emptyhisterr;

   // Temporary variables for mid-calculations
   string *stemp;
   int *itemp;
   double *dtemp;
   double *dFitError;
   int *iFitError;

   // Plotting style information
   RootStyle *mystyle;

   // Variables for chi2 test
   double *residVal;
   int chiGood;

   // Variables for step and initial value during minimization
   double mvacut;
   double *minstep;
   double *minvar;
   int nrparam;
   int nrelem;
   int stepsel;

   // Procedure (0 = range of ratios, 1 = minimization)
   int fitproc;
   // Settings option (0 = keep same settings, 1 = change for each file)
   int setopt;
   // Step size of ratio change
   double ratioStep;
   // Selection type for trees
   int treeSelect;
   // Selected trees
   vector<int> trees;
   vector<string> treeNames;
   // Type of constraint (with all parameters or all-1)
   int constraint;
   int obsAnalysis;
   bool removeZero;

   // Composition and lnA values
   double *midLna;
   vector<double> *midComposition;
   vector<double> *midCompositionErr;
   double *chi2value;
   double *pvalue;
   int *ndfvalue;
   double *midEnergy;
   double *midStep;

   TObjArray *mc;
   TFractionFitter *fracFitter;
   TVirtualFitter *fracVirtFitter;

   PrimPart *ptype;

   int rangeTransform;
   bool *iszero;

   bool scaleHist;

//   int equateVal;
//   bool equateMcEvents;
   double fitrange[2];

   vector<int> *includePart;

   // Filename
   string filename;
   TFile *ifile;

public:
   MvaFitHist(vector<string> *primVals);
   virtual ~MvaFitHist();

   double fcnATan(const double *par);
   double fcnInvATan(const double *par);
   double fcnConstrainedFunc(const double *par);
   double fcnFunc(const double *par);
   double fcn(const double *par);

   void IncludePrimaryType(string type);
   void PrepareHistograms(int run, string *fname, int *proc, double *step);
   void PlotDistributions();
   int StartFitting();
   void PlotSumResiduals(double *sigFrac, double *sigFracErr, int status, double selStep);
   void ResetHistograms();

   int GetNrElem();
   double GetFinalLna(int me);
   double GetFinalComposition(int type);
   double GetFinalCompositionErr(int type, int hilo);
   double GetFinalEnergy();
   int GetNrBins();
   double GetChi2();
   double GetPvalue();
   int GetNdf();
   double GetStep();
   int GetAnalysisType();

   int GetElemType(int type);

   double PearsonChi2(TH1 *sim, TH1 *data, int chi2dist);
};

#endif
