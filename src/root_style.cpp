#include "root_style.h"
#include <iostream>
#include <string>

using namespace std;

RootStyle::RootStyle()
{
   labelsize = 24;
   labelfont = 63;
   titlesize = 24;
   titlefont = 63;
   textsize = 18;
   textfont = 63;

   c_SignalLine = TColor::GetColor("#0000ee");
   c_SignalFill = TColor::GetColor("#7d99d1");
   c_BackgroundLine = TColor::GetColor("#ff0000");
   c_BackgroundFill = TColor::GetColor("#ff0000");
   c_DataLine = TColor::GetColor("#000000");
   c_DataFill = TColor::GetColor("#808080");
   c_DataNormLine = TColor::GetColor("#00bb00");
   c_DataNormFill = TColor::GetColor("#00bb00");
   c_ResidLine = TColor::GetColor("#000000");
   c_ResidFill = TColor::GetColor("#808080");

   signalfill = 1001;
   backgroundfill = 3554;
   datafill = 3002;
   residfill = 3001;

   legendBaseHeight = 0.033;

   basestyle = new TStyle("basestyle", "basestyle");
}

// Setting the base style settings
void RootStyle::SetBaseStyle()
{
   // Set title and label font sizes (with precision 3, these are given in pixels)
   basestyle->SetTextFont(textfont);
   basestyle->SetTextSize(textsize);
   basestyle->SetLabelFont(labelfont, "xyz");
   basestyle->SetLabelSize(labelsize, "xyz");
   basestyle->SetTitleFont(titlefont, "xyz");
   basestyle->SetTitleSize(titlesize, "xyz");
   
   // Set option and statistics
   basestyle->SetOptStat(0);
   basestyle->SetPalette(1,0);
   basestyle->SetOptTitle(0);
   basestyle->SetStatFontSize(0.024);
   basestyle->SetStatBorderSize(1);
   basestyle->SetStatColor(kGray);
   basestyle->SetStatX(0.925);
   basestyle->SetStatY(0.925);
   basestyle->SetStatW(0.13);

   // Set canvas and pads
   basestyle->SetCanvasBorderMode(0);
   basestyle->SetFrameBorderMode(0);
   basestyle->SetCanvasColor(0);
   basestyle->SetPadTickX(1);
   basestyle->SetPadTickY(1);
   basestyle->SetCanvasDefX(100);
   basestyle->SetCanvasDefY(50);
   basestyle->SetCanvasDefW(900);
   basestyle->SetCanvasDefH(600);
   basestyle->SetPadBorderMode(0);
   basestyle->SetPadBottomMargin(0.1);
   basestyle->SetPadTopMargin(0.04);
   basestyle->SetPadLeftMargin(0.125);
   basestyle->SetPadRightMargin(0.04);
   basestyle->SetPadColor(0);
   basestyle->SetPadGridX(kTRUE);
   basestyle->SetPadGridY(kTRUE);

   // Label and title offsets
   basestyle->SetLabelOffset(0.015,"xyz");
   basestyle->SetTitleOffset(1.6, "x");
   basestyle->SetTitleOffset(1.9, "y");
   basestyle->SetEndErrorSize(10);

   gROOT->SetStyle("basestyle");
   gROOT->ForceStyle(1);
}

// Setting the axis titles
void RootStyle::SetAxisTitles(TF1 *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

// Setting the axis titles
void RootStyle::SetAxisTitles(TH1 *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void RootStyle::SetAxisTitles(TH2 *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void RootStyle::SetAxisTitles(TGraph *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void RootStyle::SetAxisTitles(TGraphErrors *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void RootStyle::SetAxisTitles(TGraphAsymmErrors *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

// Setting colors for histograms
void RootStyle::SetHistColor(TH1 *plot, int sigbackdata)
{
   if(sigbackdata == 0)
   {
      plot->SetLineColor(c_BackgroundLine);
      plot->SetLineWidth(2);
      plot->SetFillColor(c_BackgroundFill);
      plot->SetFillStyle(backgroundfill);
   }
   else if(sigbackdata == 1)
   {
      plot->SetLineColor(c_SignalLine);
      plot->SetLineWidth(2);
      plot->SetFillColorAlpha(c_SignalFill, 0.75);
      plot->SetFillStyle(signalfill);
   }
   else if(sigbackdata == 2)
   {
      plot->SetLineColor(c_DataLine);
      plot->SetLineWidth(2);
      plot->SetFillColor(c_DataFill);
      plot->SetFillStyle(datafill);
   }
   else if(sigbackdata == 3)
   {
      plot->SetLineColor(c_ResidLine);
      plot->SetLineWidth(2);
      plot->SetFillColor(c_ResidFill);
      plot->SetFillStyle(residfill);
   }
}

void RootStyle::SetHistColorSimple(TH1 *plot, int sigbackdata)
{
   if(sigbackdata == 0)
   {
      plot->SetLineColor(c_BackgroundLine);
      plot->SetLineWidth(2);
   }
   else if(sigbackdata == 1)
   {
      plot->SetLineColor(c_SignalLine);
      plot->SetLineWidth(2);
   }
   else if(sigbackdata == 2)
   {
      plot->SetLineColor(c_DataLine);
      plot->SetLineWidth(2);
   }
   else if(sigbackdata == 3)
   {
      plot->SetLineColor(c_ResidLine);
      plot->SetLineWidth(2);
   }
}

// Setting colors for graphs
void RootStyle::SetGraphColor(TGraph *plot, int sigbackdata)
{
   if(sigbackdata == 0)
   {
      plot->SetLineColorAlpha(c_BackgroundLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_BackgroundLine, 0.75);
      plot->SetMarkerStyle(21);
      plot->SetMarkerSize(1.4);
   }
   else if(sigbackdata == 1)
   {
      plot->SetLineColorAlpha(c_SignalLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_SignalLine, 0.75);
      plot->SetMarkerStyle(20);
      plot->SetMarkerSize(1.4);
   }
   else if(sigbackdata == 2)
   {
      plot->SetLineColorAlpha(c_DataLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_DataLine, 0.75);
      plot->SetMarkerStyle(22);
      plot->SetMarkerSize(1.6);
   }
   else if(sigbackdata == 3)
   {
      plot->SetLineColorAlpha(c_DataNormLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_DataNormLine, 0.75);
      plot->SetMarkerStyle(23);
      plot->SetMarkerSize(1.6);
   }
}

void RootStyle::SetGraphColor(TGraphErrors *plot, int sigbackdata)
{
   if(sigbackdata == 0)
   {
      plot->SetLineColorAlpha(c_BackgroundLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_BackgroundLine, 0.75);
      plot->SetMarkerStyle(21);
      plot->SetMarkerSize(1.4);
   }
   else if(sigbackdata == 1)
   {
      plot->SetLineColorAlpha(c_SignalLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_SignalLine, 0.75);
      plot->SetMarkerStyle(20);
      plot->SetMarkerSize(1.4);
   }
   else if(sigbackdata == 2)
   {
      plot->SetLineColorAlpha(c_DataLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_DataLine, 0.75);
      plot->SetMarkerStyle(22);
      plot->SetMarkerSize(1.6);
   }
   else if(sigbackdata == 3)
   {
      plot->SetLineColorAlpha(c_DataNormLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_DataNormLine, 0.75);
      plot->SetMarkerStyle(23);
      plot->SetMarkerSize(1.6);
   }
}

void RootStyle::SetGraphColor(TGraphAsymmErrors *plot, int sigbackdata)
{
   if(sigbackdata == 0)
   {
      plot->SetLineColorAlpha(c_BackgroundLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_BackgroundLine, 0.75);
      plot->SetMarkerStyle(21);
      plot->SetMarkerSize(1.4);
   }
   else if(sigbackdata == 1)
   {
      plot->SetLineColorAlpha(c_SignalLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_SignalLine, 0.75);
      plot->SetMarkerStyle(20);
      plot->SetMarkerSize(1.4);
   }
   else if(sigbackdata == 2)
   {
      plot->SetLineColorAlpha(c_DataLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_DataLine, 0.75);
      plot->SetMarkerStyle(22);
      plot->SetMarkerSize(1.6);
   }
   else if(sigbackdata == 3)
   {
      plot->SetLineColorAlpha(c_DataNormLine, 0.75);
      plot->SetLineWidth(2);
      plot->SetMarkerColorAlpha(c_DataNormLine, 0.75);
      plot->SetMarkerStyle(23);
      plot->SetMarkerSize(1.6);
   }
}

double RootStyle::SetLegendHeight(int nrfigs)
{
   return nrfigs*legendBaseHeight;
}
