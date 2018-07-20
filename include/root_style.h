#ifndef _ROOT_STYLE_H_
#define _ROOT_STYLE_H_

#include "workstation.h"
#include "root_include.h"
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
public:
   RootStyle();

   void SetBaseStyle();

   void SetAxisTitles(TF1 *plot, string xtitle, string ytitle);
   void SetAxisTitles(TH1 *plot, string xtitle, string ytitle);
   void SetAxisTitles(TH2 *plot, string xtitle, string ytitle);
   void SetAxisTitles(TGraph *plot, string xtitle, string ytitle);
   void SetAxisTitles(TGraphErrors *plot, string xtitle, string ytitle);
   void SetAxisTitles(TGraphAsymmErrors *plot, string xtitle, string ytitle);

   void SetHistColor(TH1 *plot, int sigbackdata);
   void SetGraphColor(TGraph *plot, int sigbackdata);
   void SetGraphColor(TGraphErrors *plot, int sigbackdata);
   void SetGraphColor(TGraphAsymmErrors *plot, int sigbackdata);

   void SetHistColorSimple(TH1 *plot, int sigbackdata);

   double SetLegendHeight(int nrfigs);

   void SetColorScale(TH1 *plot, int cur, int nrscale);
   void SetColorScale(TGraph *plot, int cur, int nrscale);
   void SetColorScale(TGraphErrors *plot, int cur, int nrscale);
   void SetColorScale(TGraphAsymmErrors *plot, int cur, int nrscale);
};

#endif
