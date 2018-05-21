#ifndef _MVAEFFICIENCY_H_
#define _MVAEFFICIENCY_H_

#include "root_include.h"

class MvaEfficiency
{
private:
   TString methodName;
   TString methodTitle;
   TH1 *sig, *bgd;
   TH1 *origSigE, *origBgdE;
   TH1 *sigE, *bgdE;
   TH1 *purS, *sSig, *effpurS;
   TCanvas *effcanvas;
   TLatex *line1, *line2, *line3, *line4;
   double plotXmin, plotXmax;
   TGaxis *rightAxis;
   double maxSignificance, maxSignificanceErr;
   int maxLenTitle;
   int sigpurCutBest;
   int sigbgdCutBest;
   int maxbin;
   bool found[2];
   string plotlocation;

   TCanvas *c;
   TLegend *legend1, *legend2;
   TLine *effline;

   void SetResultHist();
   TString GetFormula();
   TString GetLatexFormula();
   void PrintResults();
public:
   MvaEfficiency(int sigN, int bgdN, string *plotdir);
   virtual ~MvaEfficiency();

   float fNSignal, fNBackground;
   TString fFormula;
   double optimalCut, sigbgdCut, sigpurCut;

   void RunMvaEfficiency(TString *inname);
   void SetPlotLimits();
   void ReadHistogram(TFile *file);
   void SetFormula(const TString& f);
   void UpdateSignificanceHist();
   void DrawHistogram();

   double GetHistValue(int type, int sigbgdpur);
};

#endif
