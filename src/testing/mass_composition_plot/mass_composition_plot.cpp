using namespace std;

void SetColor(TGraphAsymmErrors *gr, int cur, int nrbins)
{
   int ci = 1738+cur;
   double colormix = (double)cur/(double)(nrbins-1.);
   TColor *color = new TColor(ci, colormix, 0, 1.-colormix, "", 1);

   gr->SetMarkerColor(ci);
   gr->SetMarkerSize(0.9);
   gr->SetMarkerStyle(20+cur);
   gr->SetLineColor(ci);
   gr->SetLineWidth(2);
}

void SetAxisTitles(TGraphAsymmErrors *plot, string xtitle, string ytitle)
{
   plot->SetTitle("");
   plot->GetXaxis()->SetTitle(xtitle.c_str());
   plot->GetXaxis()->CenterTitle();
   plot->GetYaxis()->SetTitle(ytitle.c_str());
   plot->GetYaxis()->CenterTitle();
}

void mass_composition_plot(char *filename)
{
   int labelsize = 22;
   int labelfont = 63;
   int titlesize = 22;
   int titlefont = 63;
   int textsize = 18;
   int textfont = 63;

   int c_SignalLine = TColor::GetColor("#0000ee");
   int c_SignalFill = TColor::GetColor("#7d99d1");
   int c_BackgroundLine = TColor::GetColor("#ff0000");
   int c_BackgroundFill = TColor::GetColor("#ff0000");
   int c_DataLine = TColor::GetColor("#000000");
   int c_DataFill = TColor::GetColor("#808080");
   int c_ResidLine = TColor::GetColor("#000000");
   int c_ResidFill = TColor::GetColor("#808080");

   int signalfill = 1001;
   int backgroundfill = 3554;
   int datafill = 3002;
   int residfill = 3001;

   TStyle *basestyle = new TStyle("basestyle", "basestyle");
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

   gROOT->SetStyle("basestyle");
   gROOT->ForceStyle(1);

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   gStyle->SetEndErrorSize(6);
   TGraphAsymmErrors *grComp[4];// = new TGraphAsymmErrors(nrbins, xbin[0], ybinSig[0], xbin[1], xbin[2], ybinSig[1], ybinSig[2]);

   int nrcomp = -1;
   int nrpoints = -1;
   float *xbin;
   float *xbinError;
   float *ybinComp[4];
   float *ybinCompError[4];

   ifstream ifile;
   ifile.open(filename);

   int count = 0;

   ifile >> nrpoints;

   xbin = new float[nrpoints];
   xbinError = new float[nrpoints];
   for(int i = 0; i < 4; i++)
   {
      ybinComp[i] = new float[nrpoints];
      ybinCompError[i] = new float[nrpoints];
   }

   while(!ifile.eof())
   {
      ifile >> nrcomp;
//      cout << nrcomp << endl;
      ifile >> xbin[count] >> xbinError[count];
      xbinError[count] = 0.00001;//xbinError[count]/1.05;
//      cout << xbin[count] << "\t" << xbinError[count] << endl;

      for(int i = 0; i < nrcomp; i++)
      {
         ifile >> ybinComp[i][count] >> ybinCompError[i][count];
	 ybinComp[i][count] = ybinComp[i][count]*100.;
	 ybinCompError[i][count] = ybinCompError[i][count]*100.;
//         cout << ybinComp[i][count] << "\t" << ybinCompError[i][count] << endl;
      }

      count++;
   }

   float limit[2] = {-1, -1};
   char comma;

   cout << "Set the y-axis range (comma separated, use -1 for either number to keep original range): ";
   cin >> limit[0] >> comma >> limit[1];
   if(limit[0] == -1)
      limit[0] = -10.;
   if(limit[1] == -1)
      limit[1] = 110.;

   cout << "Y-axis range is: " << limit[0] << ", " << limit[1] << endl;

   for(int i = 0; i < nrcomp; i++)
   {
      grComp[i] = new TGraphAsymmErrors(nrpoints, xbin, ybinComp[i], xbinError, xbinError, ybinCompError[i], ybinCompError[i]);
      SetColor(grComp[i], i, nrcomp);
      grComp[i]->GetYaxis()->SetRangeUser(limit[0],limit[1]);

      if(i == 0)
      {
	 SetAxisTitles(grComp[i], "FD energy [log(E/eV)]", "Signal fraction of data events");
         grComp[i]->Draw("ALP");
      }
      else
         grComp[i]->Draw("LP;SAME");
   }

   c1->SaveAs("output.pdf");

   ifile.close();
}
