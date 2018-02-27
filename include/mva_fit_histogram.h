#ifndef _MVA_FIT_HISTOGRAM_H_
#define _MVA_FIT_HISTOGRAM_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"

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
   TH1F *simsig, *simback, *data;
   TH1F *sim1, *sim2, *sim;

   // Number of histogram bins, xrange limits and y ranges
   int nrbins, nrkeys;
   double xlim[2];
   double yrange[2], yresidrange[2];

   // Normalization factor and number of entries
   double *norm;
   int *nentries;

   // Temporary variables for mid-calculations
   string *stemp;
   int *itemp;
   double *dtemp;

   // Plotting style information
   RootStyle *mystyle;

   // Variables for chi2 test
   double *residVal, chiVal, chiProb;
   int chiNdf, chiGood;

   // Variables for step and initial value during minimization
   double *minstep;
   double *minvar;
   int nrparam;

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

   // Filename
   string filename;
   TFile *f;
public:
   MvaFitHist();
   MvaFitHist(int bincount, double xlow, double xhigh);
   virtual ~MvaFitHist();

//   void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
   double fcn(const double *par);

   void PrepareHistograms(int run, string *fname, int *proc, double *step);
   void PlotDistributions();
   void StartFitting();
   void PlotSumResiduals(double *sigFrac, double *sigFracErr);
};

#endif
