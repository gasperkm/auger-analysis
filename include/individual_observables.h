#ifndef _INDIVIDUAL_OBSERVABLES_H_
#define _INDIVIDUAL_OBSERVABLES_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"

class IndividObs
{
private:
/*   // Variables needed for the chi2 calculation
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
   double appvar;
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
   TFile *f;*/
   int sigcount, backcount;
   int nrobs;
   string mvaprintout;

   vector<int> trees;
   vector<string> treeNames;
   double *mean;
   double *minrange;
   double *maxrange;

   double range[2];

   vector<double> obsvect[3];
   vector<double> obsvect_neg[3];
   vector<double> obsvect_pos[3];

   vector<string> observables;
   vector<float> obssel;
   vector<float> tempmin;
   vector<float> tempmax;
   vector<string> tempdesc;

   string filename;

   string *stemp;
   double *dtemp;
   int *itemp;
   float *ftemp;
   char *ctemp;

   double obscut[3];

   double negerror, poserror, meannegerror, meanposerror;

   float obsvars[3];

   TH1F *signalhist, *backgroundhist, *purS;

   int nbin_hist_high;

   int sign;

   int sigpurCutBest, backpurCutBest, sigbgdCutBest;
public:
   IndividObs(string infile, vector<int> inv);
   virtual ~IndividObs();

   int ReadObservables();

   void PrepareFolders();

   void AppendPrintout(string print);
   void GetObsVectors(int obs);
   void CalculateCut(int obs);
   void ApplyObsCut(int obs, int type);
   void PrintoutResults(int obs);

   int GetNrObservables();
};

#endif
