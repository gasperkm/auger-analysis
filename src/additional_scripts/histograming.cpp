#define _STANDALONE_ 1
#include "workstation.h"
#include "separate_functions.h"

using namespace std;

int main(int argc, char **argv)
{
   string filename;
   if(argc > 1)
      filename = string(argv[1]);
   else
   {
      cerr << "Error! No input file supplied. Rerun program and add input file as argument." << endl;
      return 1;
   }
   
//   gSystem->Load("libTree.so");

   //-----------------------------------------
   ifstream ifs;
   string *stemp = new string[5];
   int *itemp = new int[2];
   float *ftemp = new float[2];
   char *ctemp = new char[1024];

   // Use observables text file to determine which observables are saved in output file and what their ranges are
   stemp[0] = string(rootdir) + "/input/observables.txt";
   ifs.open(stemp[0].c_str(), ifstream::in);
   cout << "Observables file: " << stemp[0] << endl;

   vector<string> observables;
   vector<float> obssel;
   vector<float> tempmin;
   vector<float> tempmax;
   vector<string> tempdesc;

   if(ifs.is_open())
   {
      // Empty vectors that will hold:
      // - observable names (observables)
      // - observable initial selection (obssel) - not really needed here
      // - minimum range value (tempmin)
      // - maximum range value (tempmax)
      // - observable description to be used on plot axes (tempdesc)
      if(!observables.empty())
         observables.erase(observables.begin(), observables.end());
      if(!obssel.empty())
         obssel.erase(obssel.begin(), obssel.end());
      if(!tempmin.empty())
         tempmin.erase(tempmin.begin(), tempmin.end());
      if(!tempmax.empty())
         tempmax.erase(tempmax.begin(), tempmax.end());
      if(!tempdesc.empty())
         tempdesc.erase(tempdesc.begin(), tempdesc.end());

      // Read from the observables.txt file
      while(1)
      {
	 if(ifs.peek() == '#')
	 {
            ifs.getline(ctemp, 1024);
            cout << "Line: " << ctemp << endl;
	 }
	 else
	 {
            ifs >> stemp[1] >> itemp[0] >> ftemp[0] >> ftemp[1];
            if(ifs.eof()) break;
            observables.push_back(stemp[1]);
	    obssel.push_back((bool)itemp[0]);
	    tempmin.push_back(ftemp[0]);
	    tempmax.push_back(ftemp[1]);

            ifs.getline(ctemp, 1024);
	    stemp[2] = string(ctemp);
            tempdesc.push_back(stemp[2]);

            cout << "Values: " << stemp[1] << "\t" << itemp[0] << "\t" << ftemp[0] << "\t" << ftemp[1] << "\t" << stemp[2] << endl;
	 }
      }
   }
   else
   {
      cerr << "Observables list error: Can't open the observables list file " << stemp[0] << ". Please make sure it exists. If not, create it. Rerun program." << endl;
      return 1;
   }

   // Save the number of observables into itemp[1]
   itemp[1] = observables.size();
   cout << "Number of observables = " << itemp[1] << endl;

   ifs.close();

   //-----------------------------------------
   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/histograms";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/histograms/hist_*";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/histograms/risetime-comparison_*";
   system(stemp[0].c_str());

   // Prepare colors for signal, background and MVA cut line
   int c_SignalLine     = TColor::GetColor("#0000ee");
   int c_SignalFill     = TColor::GetColor("#7d99d1");
   int c_AllLine        = TColor::GetColor("#ff0000");
   int c_AllFill        = TColor::GetColor("#ff0000");
   int c_MvaCut         = TColor::GetColor("#ffff66");

   // Prepare plotting objects (legend, line, histograms, graphs,...) - for a maximum of 50 observables
   TLegend *legend;
   TLine *line[3];
   TH1F *treeHist[50];
   TH1F *riseHist;
   TH1F *riserecalcHist;
   TLatex uoflowtext;

   // Prepare other observables
   int legendFill = 1001;		// filling style for legend
   float *xhistlimit;			// array for x axis limits
   xhistlimit = new float[3];
   double *dtemp;			// temporary deouble number
   dtemp = new double[2];

   float *max;				// maximum value in the histogram (in order to scale y axis so that histogram is below legend)
   max = new float[50];
   float *maxobs;
   maxobs = new float[4];
   int *riseobs;
   riseobs = new int[2];
   int *uflowcount, *oflowcount, *totcount; // counting number of underflow and overflow events
   uflowcount = new int[50];
   oflowcount = new int[50];
   totcount = new int[50];

   bool isangle[2] = {false, false};	// check if observable is an angle (needs to be converted to degrees from radians) 
   bool isgood[2] = {true, true};	// check if observable has an invalid value -1

   float *obsvars;			// array for variable values that we read from output file
   obsvars = new float[150];

   float *minimumval, *maximumval;	// holders for minimum and maximum values of observables
   minimumval = new float[50];
   maximumval = new float[50];
   for(int i = 0; i < itemp[1]; i++)
   {
      minimumval[i] = 1.e+30;
      maximumval[i] = -1.e+30;
   }

   float completemin = 1.e+30, completemax = -1.e+30;	// minimum and maximum values of observables from all trees

   // Prepare canvases and clear any ROOT default statistics
   gStyle->SetOptStat(0);

   gStyle->SetCanvasBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetPalette(1,0);
   gStyle->SetOptTitle(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetStatFontSize(0.024);
   gStyle->SetStatBorderSize(1);
   gStyle->SetStatColor(kGray);
   gStyle->SetStatX(0.925);
   gStyle->SetStatY(0.925);
   gStyle->SetStatW(0.13);
   gStyle->SetTextFont(132);
   gStyle->SetTextSize(0.08);
   gStyle->SetLabelSize(0.03,"xyz");
   gStyle->SetLabelOffset(0.01,"xyz");
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetCanvasDefX(100);
   gStyle->SetCanvasDefY(50);
   gStyle->SetCanvasDefW(900);
   gStyle->SetCanvasDefH(600);
   gStyle->SetPadBottomMargin(0.1);
   gStyle->SetPadTopMargin(0.04);
   gStyle->SetPadLeftMargin(0.125);
   gStyle->SetPadRightMargin(0.04);

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);

   TCanvas *c2 = new TCanvas("c2","",1200,900);
   c2->SetGrid();
   c2->SetRightMargin(0.05);
   c2->SetTopMargin(0.05);

   // Open output file to read saved values
   TFile *ifile = new TFile(filename.c_str(), "READ");
   TTree *getTree;
   TList *tempkeyslist = (TList*)ifile->GetListOfKeys();

   // Loop over all trees
   for(int k = 0; k < ifile->GetNkeys(); k++)
   {
      // Zero minimum and maximum holders for observables
      for(int i = 0; i < itemp[1]; i++)
      {
         minimumval[i] = 1.e+30;
         maximumval[i] = -1.e+30;
	 uflowcount[i] = 0;
	 oflowcount[i] = 0;
	 totcount[i] = 0;

	 if(observables[i] == "risetime")
	    riseobs[0] = i;

	 if(observables[i] == "risetimerecalc")
	    riseobs[1] = i;
      }

      // Select the trees
      stemp[1] = string((tempkeyslist->At(k))->GetName());
      cout << endl << stemp[1] << " has been selected for evaluation." << endl;
      getTree = (TTree*)ifile->Get(stemp[1].c_str());
      cout << "Number of events in the tree = " << getTree->GetEntries() << endl;

      // Set addresses for reading observables
      for(int i = 0; i < itemp[1]; i++)
      {
         cout << "Setting observable " << observables[i] << endl;
         getTree->SetBranchAddress((observables[i]).c_str(), &obsvars[3*i]); // mean
         cout << "Mean variable: obsvars[3*" << i << "] = obsvars[" << 3*i << "]" << endl;
	 stemp[2] = observables[i] + "_neg";
         getTree->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+1]); // neg error
         cout << "Negerror variable: obsvars[3*" << i << "+1] = obsvars[" << 3*i+1 << "]" << endl;
	 stemp[2] = observables[i] + "_pos";
         getTree->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+2]); // pos error
         cout << "Poserror variable: obsvars[3*" << i << "+2] = obsvars[" << 3*i+2 << "]" << endl;
      }

      // Prepare plots to save values into (histograms need to able to rebin)
      if(tempmin[riseobs[0]] < tempmin[riseobs[1]])
         maxobs[0] = tempmin[riseobs[0]];
      else
         maxobs[0] = tempmin[riseobs[1]];
      
      if(tempmax[riseobs[0]] > tempmax[riseobs[1]])
         maxobs[1] = tempmax[riseobs[0]];
      else
         maxobs[1] = tempmax[riseobs[1]];
      
      maxobs[2] = (maxobs[1]-maxobs[0])*0.05;
      maxobs[0] = maxobs[0] - maxobs[2];
      maxobs[1] = maxobs[1] + maxobs[2];

      for(int i = 0; i < itemp[1]; i++)
      {
	 xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
         xhistlimit[0] = tempmin[i] - xhistlimit[2];
	 xhistlimit[1] = tempmax[i] + xhistlimit[2];

	 // Prepare histograms for all observables
	 stemp[2] = "treeHist" + ToString(i);
         treeHist[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
         treeHist[i]->SetBit(TH1::kCanRebin);

	 // Special histograms for risetime and recalculated risetime
	 if(observables[i] == "risetime")
	 {
            riseHist = new TH1F("riseHist", observables[i].c_str(), 100, maxobs[0], maxobs[1]);
            riseHist->SetBit(TH1::kCanRebin);
	 }
	 if(observables[i] == "risetimerecalc")
	 {
            riserecalcHist = new TH1F("riserecalcHist", observables[i].c_str(), 100, maxobs[0], maxobs[1]);
            riserecalcHist->SetBit(TH1::kCanRebin);
	 }

	 // Set maximum value to zero, and underflow and overflow counters to zero
         max[i] = 0.0;
         uflowcount[i] = 0;
         oflowcount[i] = 0;
      }

      // Loop over all entries
      for(int ievt = 0; ievt < getTree->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTree->GetEntry(ievt);

	 // Loop over all observables
         for(int i = 0; i < itemp[1]; i++)
         {
            // Check if observable is an angle in radians
            if( (observables[i] == "zenithSD") || (observables[i] == "azimuthSD") || (observables[i] == "zenithFD") || (observables[i] == "azimuthFD") || (observables[i] == "latitudeSD") || (observables[i] == "longitudeSD") || (observables[i] == "latitudeFD") || (observables[i] == "longitudeFD") )
               isangle[0] = true;
	    else
               isangle[0] = false;

            // Check if observable is invalid and has value -1
	    if(obsvars[3*i] == -1)
               isgood[0] = false;
	    else
               isgood[0] = true;

	    // Check the mean observable value in order to get minimum and maximum values of all trees
	    if(obsvars[3*i] < minimumval[i])
	    {
	       if(isgood[0])
	       {
                  if(isangle[0])
                     minimumval[i] = RadToDeg(obsvars[3*i]);
	          else
                     minimumval[i] = obsvars[3*i];
	       }
	    }
	    if(obsvars[3*i] > maximumval[i])
	    {
	       if(isgood[0])
	       {
                  if(isangle[0])
                     maximumval[i] = RadToDeg(obsvars[3*i]);
	          else
                     maximumval[i] = obsvars[3*i];
	       }
	    }

	    // Save observable values into the histograms
	    if(isgood[0])
	    {
	       xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
	       xhistlimit[1] = tempmax[i] + xhistlimit[2];

	       if(obsvars[3*i] > xhistlimit[1])
                  oflowcount[i]++;
	       else if(obsvars[3*i] < xhistlimit[0])
                  uflowcount[i]++;
               else
	       {
	          if(i == riseobs[0])
                     riseHist->Fill(obsvars[3*i]);

	          if(i == riseobs[1])
                     riserecalcHist->Fill(obsvars[3*i]);

                  if(isangle[0])
                     treeHist[i]->Fill(RadToDeg(obsvars[3*i]));
                  else
                     treeHist[i]->Fill(obsvars[3*i]);
	       }
                  
               totcount[i]++;
	    }
         } // Loop over all observables

//         cout << endl;
      } // Loop over all entries

      // Write the valid event information
      cout << "Valid events (" << getTree->GetTitle() << "):" << endl;
      for(int i = 0; i < itemp[1]; i++)
      {
         cout << " - " << observables[i] << " = " << totcount[i] << " (UFLOW=" << uflowcount[i] << ", OFLOW=" << oflowcount[i] << ")" << endl;
      }
      cout << endl;

      // Plot all observables
      for(int i = 0; i < itemp[1]; i++)
      {
         c1->cd();

         // Check maximum value in a histogram
         if((treeHist[i]->GetMaximum()) > max[i])
            max[i] = treeHist[i]->GetMaximum();

	 // Set line and fill attributes for histograms
         treeHist[i]->SetLineColor(c_SignalLine);
         treeHist[i]->SetLineWidth(2);
         treeHist[i]->SetFillColor(c_SignalFill);
         treeHist[i]->SetFillStyle(1001);

	 // Draw signal and background histograms
         treeHist[i]->Draw();

	 // Set axis ranges
	 xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
         xhistlimit[0] = tempmin[i] - xhistlimit[2];
	 xhistlimit[1] = tempmax[i] + xhistlimit[2];
         treeHist[i]->GetYaxis()->SetRangeUser(0.,max[i]*1.1);
         treeHist[i]->SetMaximum(max[i]*1.1);
         treeHist[i]->GetXaxis()->SetRange(xhistlimit[0], xhistlimit[1]);
         treeHist[i]->GetXaxis()->SetRangeUser(xhistlimit[0], xhistlimit[1]);

	 // Setup axis (title, labels)
         treeHist[i]->GetYaxis()->SetTitle("Number of events");
         treeHist[i]->GetXaxis()->SetTitle(tempdesc[i].c_str());
         treeHist[i]->GetXaxis()->SetTitleOffset(1.2);
         treeHist[i]->GetXaxis()->CenterTitle(kTRUE);
         treeHist[i]->GetXaxis()->SetLabelSize(0.028);
         treeHist[i]->GetXaxis()->SetLabelOffset(0.015);
         treeHist[i]->GetYaxis()->SetTitleOffset(1.3);
         treeHist[i]->GetYaxis()->CenterTitle(kTRUE);
         treeHist[i]->GetYaxis()->SetLabelSize(0.028);
         treeHist[i]->GetYaxis()->SetLabelOffset(0.015);
         treeHist[i]->SetTitle("");

         // Write down the underflow and overflow information
         uoflowtext.SetTextSize(0.017);
         uoflowtext.SetTextAlign(11);
         stemp[0] = "UFLOW = " + ToString(uflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[0], max[i]*1.105, stemp[0].c_str());
         uoflowtext.SetTextAlign(31);
         stemp[0] = "OFLOW = " + ToString(oflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[1], max[i]*1.105, stemp[0].c_str());

	 // Prepare save name for plots
         stemp[1] = string((tempkeyslist->At(k))->GetName());
         stemp[3] = RemoveFilename(&filename) + "/histograms/hist_" + stemp[1] + "_" + observables[i] + ".pdf";

	 // Save plots as PDF
         c1->SaveAs(stemp[3].c_str());
      } // Plot all observables

      // Plot risetime versus risetimerecalc
      c2->cd();
      maxobs[3] = 0.0;

      // Check maximum value in a histogram
      if((riseHist->GetMaximum()) > maxobs[3])
         maxobs[3] = riseHist->GetMaximum();
      if((riserecalcHist->GetMaximum()) > maxobs[3])
         maxobs[3] = riserecalcHist->GetMaximum();

      // Set line and fill attributes for histograms
      riserecalcHist->SetLineColor(c_SignalLine);
      riserecalcHist->SetLineWidth(2);
      riserecalcHist->SetFillColor(c_SignalFill);
      riserecalcHist->SetFillStyle(1001);
      riseHist->SetLineColor(c_AllLine);
      riseHist->SetLineWidth(2);
      riseHist->SetFillColor(c_AllFill);
      riseHist->SetFillStyle(3554);

      // Draw signal and background histograms
      riserecalcHist->Draw();
      riseHist->Draw("same");

      // Set axis ranges
      riserecalcHist->GetYaxis()->SetRangeUser(0.,maxobs[3]*1.2);
      riseHist->GetYaxis()->SetRangeUser(0.,maxobs[3]*1.2);
      riserecalcHist->SetMaximum(maxobs[3]*1.2);
      riseHist->SetMaximum(maxobs[3]*1.2);
      riserecalcHist->GetXaxis()->SetRange(maxobs[0], maxobs[1]);
      riserecalcHist->GetXaxis()->SetRangeUser(maxobs[0], maxobs[1]);

      // Setup axis (title, labels)
      riserecalcHist->GetYaxis()->SetTitle("Number of events");
      riserecalcHist->GetXaxis()->SetTitle(tempdesc[riseobs[1]].c_str());
      riserecalcHist->GetXaxis()->SetTitleOffset(1.2);
      riserecalcHist->GetXaxis()->CenterTitle(kTRUE);
      riserecalcHist->GetXaxis()->SetLabelSize(0.028);
      riserecalcHist->GetXaxis()->SetLabelOffset(0.015);
      riserecalcHist->GetYaxis()->SetTitleOffset(1.3);
      riserecalcHist->GetYaxis()->CenterTitle(kTRUE);
      riserecalcHist->GetYaxis()->SetLabelSize(0.028);
      riserecalcHist->GetYaxis()->SetLabelOffset(0.015);
      riserecalcHist->SetTitle("");

      // Draw a legend with info on signal and background events
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.12, gPad->GetLeftMargin()+.45, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);
      stemp[0] = "Offline risetime (" + ToString(totcount[riseobs[0]]) + ")";
      legend->AddEntry(riseHist, stemp[0].c_str(),"f");
      stemp[0] = "Recalculated risetime (" + ToString(totcount[riseobs[1]]) + ")";
      legend->AddEntry(riserecalcHist, stemp[0].c_str(),"f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      // Write down the underflow and overflow information
      uoflowtext.SetTextSize(0.017);
      uoflowtext.SetTextAlign(11);
      stemp[0] = "UFLOW = " + ToString(uflowcount[riseobs[0]]+uflowcount[riseobs[1]]);
      uoflowtext.DrawLatex(maxobs[0], maxobs[3]*1.205, stemp[0].c_str());
      uoflowtext.SetTextAlign(31);
      stemp[0] = "OFLOW = " + ToString(oflowcount[riseobs[0]]+oflowcount[riseobs[1]]);
      uoflowtext.DrawLatex(maxobs[1], maxobs[3]*1.205, stemp[0].c_str());

      // Prepare save name for plot
      stemp[1] = string((tempkeyslist->At(k))->GetName());
      stemp[3] = RemoveFilename(&filename) + "/histograms/risetime-comparison_" + stemp[1] + ".pdf";

      // Save plot as PDF
      c2->SaveAs(stemp[3].c_str());

      delete legend;

      // Delete all plot objects
      for(int i = 0; i < itemp[1]; i++)
      {
         delete treeHist[i];
      }
      delete riseHist;
      delete riserecalcHist;
   } // Loop over all trees

   ifile->Close();

   delete ifile;

   //-----------------------------------------

   delete[] stemp;
   delete[] itemp;
   delete[] ftemp;
   delete[] ctemp;
   delete c1;
   delete c2;
   delete[] dtemp;
   delete[] xhistlimit;
   delete[] max;
   delete[] obsvars;
   delete[] uflowcount;
   delete[] oflowcount;
   delete[] totcount;
   delete[] minimumval;
   delete[] maximumval;

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
