#ifndef _HISTOGRAMING_H_
#define _HISTOGRAMING_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"
#include "primary_type.h"
#include "observables.h"
#include "popups.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#define MAXINFILES 40

class ScatHist
{
private:

   int ret;
   int type; // type defines if this will be a histogram (0) or a scatter plot (1)

   string *stemp;
   int *itemp;
   double *dtemp;

   int *nrbins;
   int *nrobs;
   int *nrtrees;
   int *nrfiles;
   vector<string> *filename;
   string *method;
   vector<string> *obser;
   vector<int> *trees;
   vector<int> *treeColor;
   double *xlim;
   double *ylim;
   double *yrange;
   bool *otherSettings;

   TH1F *treeHist[3*MAXINFILES];
   TGraphAsymmErrors *treeAGraph[3*MAXINFILES];
   TLegend *legend;
   TLatex *uoflowtext;
   const int legendFill = 1001;
   const int c_MvaCut = TColor::GetColor("#ffff66");

   int *uflowcount, *oflowcount, *totcount[4]; // counting number of underflow and overflow events
   double *max;
   double *min;

   RootStyle *mystyle;

   TFile *ifile[3];

   TTree *getTree[3];

   float *obsvars;
   float *obsvars_neg;
   float *obsvars_pos;

   float *obsvars2;
   float *obsvars2_neg;
   float *obsvars2_pos;

   bool *isgood;

   double *mvalim;

   FSDialog *legendselectDialog;

/*   string *filename;
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
   vector<double> *composition;*/

public:
   ScatHist(int *inType);
   virtual ~ScatHist();

   void SetFilenames(vector<string> *inFiles);
   void SetMethod(string *inMethod);
   void SetObservables(vector<string> *inObs);
   void SetTrees(vector<int> *inTrees, vector<int> *inColor);
   void SetNrBins(int *inNrBins);
   void SetXaxisLimits(bool use, double *inLimits);
   void SetYaxisLimits(bool use, double *inLimits);
   void SetOtherSettings(bool *inConst);

   int StartPlotting(Observables *genObser);

   bool InvertAxis(string *obs1, string *obs2);

/*   void SetPrimaries(vector<int> *primVals);
   void SetMethod(string *inMethod);
   void SetHImodel(int *inModel);
   void SetSimProduction(int *inProd);
   void SetData(int *inTree);
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
   int ReadFracResultsFD(vector<float> *val, bool uselna);*/
};

#endif
