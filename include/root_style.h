#ifndef _ROOT_STYLE_H_
#define _ROOT_STYLE_H_

#include "workstation.h"
#include "root_include.h"

class RootStyle
{
private:
   int labelsize, labelfont, titlesize, titlefont, textsize, textfont;
   int c_SignalLine, c_SignalFill, c_BackgroundLine, c_BackgroundFill, c_DataLine, c_DataFill, c_DataNormLine, c_DataNormFill, c_ResidLine, c_ResidFill;
   int signalfill, backgroundfill, datafill, residfill;
   TStyle *basestyle;
public:
   RootStyle();

   void SetBaseStyle();

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
};

#endif
