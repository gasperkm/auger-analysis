#ifndef _ROOT_STYLE_H_
#define _ROOT_STYLE_H_

#include "workstation.h"
#include "root_include.h"
#include "separate_functions.h"
#include <string>

using namespace std;

class RootStyle
{
private:
   int labelsize, labelfont, titlesize, titlefont, textsize, textfont;
   int c_SignalLine, c_SignalFill, c_BackgroundLine, c_BackgroundFill, c_DataLine, c_DataFill, c_DataNormLine, c_DataNormFill, c_ResidLine, c_ResidFill;
   int signalfill, backgroundfill, datafill, residfill;
   double legendBaseHeight;
   TStyle *basestyle;

   const double padTotDiffFactor = 0.003;
   const double padHeightDiffFactor = 0.004;
   const double padMarginDiffFactor = 0.014;
   const double xTitleFactor = 3.6;
   const double yTitleFactor = 2.25;
   const double xTitleSingleFactor = 2.4;
   const double yTitleSingleFactor = 2.85;
public:
   RootStyle();

   void SetBaseStyle();

   void SetAxisTitles(TF1 *plot, string xtitle, string ytitle);
   void SetAxisTitles(TH1 *plot, string xtitle, string ytitle);
   void SetAxisTitles(TH2 *plot, string xtitle, string ytitle);
   void SetAxisTitles(TGraph *plot, string xtitle, string ytitle);
   void SetAxisTitles(TGraph2D *plot, string xtitle, string ytitle, string ztitle);
   void SetAxisTitles(TGraphErrors *plot, string xtitle, string ytitle);
   void SetAxisTitles(TGraphAsymmErrors *plot, string xtitle, string ytitle);

   void SetHistColor(TH1 *plot, int sigbackdata);
   void SetFuncColor(TF1 *plot, int sigbackdata);
   void SetGraphColor(TGraph *plot, int sigbackdata);
   void SetGraphColor(TGraphErrors *plot, int sigbackdata);
   void SetGraphColor(TGraphAsymmErrors *plot, int sigbackdata);

   void SetHistColorSimple(TH1 *plot, int sigbackdata);

   double SetLegendHeight(int nrfigs);

   void SetColorScale(TH1 *plot, int cur, int nrscale);
   void SetColorScale(TGraph *plot, int cur, int nrscale);
   void SetColorScale(TGraphErrors *plot, int cur, int nrscale);
   void SetColorScale(TGraphAsymmErrors *plot, int cur, int nrscale);

   double GetPlotWidth(int size);
   double GetPlotHeight(int type, int size);

   void SetSinglePlot(int xsize, int ysize, TCanvas *inCanv);
   void SetMultiPlot(int xsize, int ysize, int nrpads, TCanvas *inCanv);
   void SetPaddedPlot(int nrpads, TCanvas *inCanv, TPad **inPads);

   double GetPaddedXoffset(int nrpads, TCanvas *inCanv);
   double GetPaddedYoffset(int nrpads, TCanvas *inCanv);

   double GetSingleXoffset(TCanvas *inCanv);
   double GetSingleYoffset(TCanvas *inCanv);
};

#endif
