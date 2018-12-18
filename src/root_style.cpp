#include "root_style.h"
#include <iostream>

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
   basestyle->SetTitleOffset(1.9, "z");
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

void RootStyle::SetAxisTitles(TGraph2D *plot, string xtitle, string ytitle, string ztitle)
{
   plot->SetTitle("");
   plot->GetHistogram()->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetHistogram()->GetXaxis()->CenterTitle();
   plot->GetHistogram()->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetHistogram()->GetYaxis()->CenterTitle();
   plot->GetHistogram()->GetZaxis()->SetTitle(ztitle.c_str());
   plot->GetHistogram()->GetZaxis()->CenterTitle();
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

// Setting colors for functions
void RootStyle::SetFuncColor(TF1 *plot, int sigbackdata)
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

void RootStyle::SetColorScale(TH1 *plot, int cur, int nrscale)
{
   int *ci = new int;
   *ci = 1738+cur;
   double *colormix = new double;
   *colormix = (double)cur/(double)(nrscale-1.);
   TColor *color = new TColor(*ci, *colormix, 0, 1.-(*colormix), "", 1);

   plot->SetMarkerColor(*ci);
   plot->SetLineColor(*ci);

   delete ci;
   delete colormix;
/*   gr->SetMarkerSize(0.9);
   gr->SetMarkerStyle(20+cur);
   gr->SetLineWidth(2);*/
}

void RootStyle::SetColorScale(TGraph *plot, int cur, int nrscale)
{
   int *ci = new int;
   *ci = 1738+cur;
   double *colormix = new double;
   *colormix = (double)cur/(double)(nrscale-1.);
   TColor *color = new TColor(*ci, *colormix, 0, 1.-(*colormix), "", 1);

   plot->SetMarkerColor(*ci);
   plot->SetLineColor(*ci);

   delete ci;
   delete colormix;
}

void RootStyle::SetColorScale(TGraphErrors *plot, int cur, int nrscale)
{
   int *ci = new int;
   *ci = 1738+cur;
   double *colormix = new double;
   *colormix = (double)cur/(double)(nrscale-1.);
   TColor *color = new TColor(*ci, *colormix, 0, 1.-(*colormix), "", 1);

   plot->SetMarkerColor(*ci);
   plot->SetLineColor(*ci);

   delete ci;
   delete colormix;
}

void RootStyle::SetColorScale(TGraphAsymmErrors *plot, int cur, int nrscale)
{
   int *ci = new int;
   *ci = 1738+cur;
   double *colormix = new double;
   *colormix = (double)cur/(double)(nrscale-1.);
   TColor *color = new TColor(*ci, *colormix, 0, 1.-(*colormix), "", 1);

   plot->SetMarkerColor(*ci);
   plot->SetLineColor(*ci);

   delete ci;
   delete colormix;
}

double RootStyle::GetPlotWidth(int size)
{
   // Narrow
   if(size == -1)
      return 1000;
   // Medium/Default
   else if(size == 0)
      return 1200;
   // High
   else if(size == 1)
      return 1400;
}
double RootStyle::GetPlotHeight(int type, int size)
{
   // Single plots
   if(type == 0)
   {
      // Very narrow
      if(size == -2)
         return 400;
      // Narrow
      else if(size == -1)
         return 600;
      // Medium/Default
      else if(size == 0)
         return 800;
      // High
      else if(size == 1)
         return 1000;
      // Very high
      else if(size == 2)
         return 1000;
   }
   // Multipad plots
   else if(type == 1)
   {
      // Very narrow
      if(size == -2)
         return 200;
      // Narrow
      else if(size == -1)
         return 300;
      // Medium/Default
      else if(size == 0)
         return 400;
      // High
      else if(size == 1)
         return 500;
      // Very high
      else if(size == 2)
         return 600;
   }
}

void RootStyle::SetSinglePlot(int xsize, int ysize, TCanvas *inCanv)
{
   inCanv->SetCanvasSize(GetPlotWidth(xsize), GetPlotHeight(0, ysize));
   inCanv->Modified();
   inCanv->Update();
}

void RootStyle::SetMultiPlot(int xsize, int ysize, int nrpads, TCanvas *inCanv)
{
   inCanv->SetCanvasSize(GetPlotWidth(xsize), nrpads*GetPlotHeight(1, ysize));
   inCanv->Modified();
   inCanv->Update();
}

void RootStyle::SetPaddedPlot(int nrpads, TCanvas *inCanv, TPad **inPads)
{
   double *dtemp = new double[4];
   string *stemp = new string;

   inCanv->cd();

   dtemp[1] = 1.;
   for(int i = 0; i < nrpads; i++)
   {
      dtemp[2] = 1./nrpads - padTotDiffFactor;
//      cout << "1/nrpads - padTotDiffFactor = " << dtemp[2] << endl;
      if(i != nrpads-1)
         dtemp[3] = dtemp[2] - padHeightDiffFactor;
      else
         dtemp[3] = dtemp[2] + (nrpads-1)*padHeightDiffFactor;

//      cout << "Each pad height = " << dtemp[3] << endl;

      dtemp[0] = dtemp[1] - dtemp[3];

//      cout << "Pad limits = " << dtemp[0] << ", " << dtemp[1] << endl;

      *stemp = "pad" + ToString(i+1);
      inPads[i] = new TPad(stemp->c_str(), "", 0.005, dtemp[0], 0.995, dtemp[1]);

      if(i == nrpads-1)
         inPads[i]->SetBottomMargin(0.1 + padMarginDiffFactor*nrpads);

//      cout << " IN: Pad " << i << " bottom margin = " << inPads[i]->GetBottomMargin() << endl;

      dtemp[1] -= dtemp[3];

      inPads[i]->Draw();
   }

   inCanv->Modified();
   inCanv->Update();

   delete stemp;
   delete[] dtemp;
}

double RootStyle::GetPaddedXoffset(int nrpads, TCanvas *inCanv)
{
//   return xTitleFactor*(inCanv->GetWindowHeight())*nrpads/(inCanv->GetWindowWidth());
   return xTitleFactor*(inCanv->GetWindowHeight())/(inCanv->GetWindowWidth());
}

double RootStyle::GetPaddedYoffset(int nrpads, TCanvas *inCanv)
{
//   return yTitleFactor*(inCanv->GetWindowHeight())*nrpads/(inCanv->GetWindowWidth());
   return yTitleFactor*(inCanv->GetWindowHeight())/(inCanv->GetWindowWidth());
}

double RootStyle::GetSingleXoffset(TCanvas *inCanv)
{
   return xTitleSingleFactor*(inCanv->GetWh())/(inCanv->GetWw());
}

double RootStyle::GetSingleYoffset(TCanvas *inCanv)
{
   return yTitleSingleFactor*(inCanv->GetWh())/(inCanv->GetWw());
}
