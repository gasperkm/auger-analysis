#ifndef _MVAEFFICIENCY_H_
#define _MVAEFFICIENCY_H_

#include "workstation.h"
#include "root_include.h"
#include <string>
/*#include "TH1.h"
#include "TROOT.h"
#include "TList.h"
#include "TIterator.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH2.h"
#include "TFormula.h"
#include "TFile.h"
#include "TApplication.h"
#include "TKey.h"
#include "TClass.h"
#include "TGaxis.h"
#include "TString.h"

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TBranchRef.h"
#include "TLeafI.h"
#include "TLeafF.h"
#include "TLeafC.h"

#include "TGWindow.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"*/

using namespace std;

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
   string plotdirectory;

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

   void WriteoutMethod();
};

#endif
