/*#define _STANDALONE_ 1*/
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "histograming.h"

using namespace std;

ScatHist::ScatHist(int *inType)
{
   stemp = new string[5];
   itemp = new int[4];
   dtemp = new double[4];

   xlim = new double[2];
   ylim = new double[2];

   uflowcount = new int[3*MAXINFILES];
   oflowcount = new int[3*MAXINFILES];
   totcount[0] = new int[3*MAXINFILES];
   totcount[1] = new int[3*MAXINFILES];
   totcount[2] = new int[3*MAXINFILES];
   totcount[3] = new int[3*MAXINFILES];

   type = *inType;

   mystyle = new RootStyle();

   isgood = new bool[3];

   uoflowtext = new TLatex();
}

ScatHist::~ScatHist()
{
   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   delete[] xlim;
   delete[] ylim;

   delete[] uflowcount;
   delete[] oflowcount;
   delete[] totcount[0];
   delete[] totcount[1];
   delete[] totcount[2];
   delete[] totcount[3];

   delete filename;
   delete obser;
   delete nrobs;
   delete trees;
   delete treeColor;
   delete nrtrees;
   delete nrfiles;

   delete nrbins;
   delete[] otherSettings;

   delete mystyle;

   delete[] isgood;

   delete uoflowtext;

   delete[] mvalim;
}

void ScatHist::SetFilenames(vector<string> *inFiles)
{
   filename = new vector<string>;
   for(int i = 0; i < inFiles->size(); i++)
      filename->push_back(inFiles->at(i));

   nrfiles = new int;
   *nrfiles = inFiles->size();
}

void ScatHist::SetMethod(string *inMethod)
{
   method = new string;
   mvalim = new double[2];
   mvalim[0] = GetMethodMin(*inMethod);
   mvalim[1] = GetMethodMax(*inMethod);
}

void ScatHist::SetObservables(vector<string> *inObs)
{
   obser = new vector<string>;
   for(int i = 0; i < inObs->size(); i++)
      obser->push_back(inObs->at(i));

   nrobs = new int;
   *nrobs = inObs->size();
}

void ScatHist::SetTrees(vector<int> *inTrees, vector<int> *inColor)
{
   trees = new vector<int>;
   treeColor = new vector<int>;

   nrtrees = new int;
   *nrtrees = 0;

   for(int i = 0; i < inTrees->size(); i++)
   {
      if(inTrees->at(i) != 0)
      {
         trees->push_back(inTrees->at(i));
	 (*nrtrees)++;

	 // Setting the signal color (1)
	 if(inColor->at(i) == 0)
	    treeColor->push_back(1);
	 // Setting the background color (0)
	 else if(inColor->at(i) == 1)
	    treeColor->push_back(0);
	 // Setting any other color
	 else
	    treeColor->push_back(inColor->at(i));
      }
   }
}

void ScatHist::SetNrBins(int *inNrBins)
{
   nrbins = new int;
   *nrbins = *inNrBins; 
}

void ScatHist::SetXaxisLimits(bool use, double *inLimits)
{
   if(use)
   {
      xlim[0] = inLimits[0]; 
      xlim[1] = inLimits[1]; 

      cout << "Custom X-axis limits: " << xlim[0] << ", " << xlim[1] << endl;
   }
}

void ScatHist::SetYaxisLimits(bool use, double *inLimits)
{
   if(use)
   {
      ylim[0] = inLimits[0]; 
      ylim[1] = inLimits[1]; 

      cout << "Custom Y-axis limits: " << ylim[0] << ", " << ylim[1] << endl;
   }
}

void ScatHist::SetOtherSettings(bool *inConst)
{
   otherSettings = new bool[4];
   otherSettings[0] = inConst[0];

   // Changing mean value for uncertainties (only for histograms)
   if(type == 0)
   {
      otherSettings[1] = inConst[1];
      otherSettings[2] = inConst[2];
   }
}

// Plot histograms or scatter plots
int ScatHist::StartPlotting(Observables *genObser)
{
   mystyle->SetBaseStyle();
   TCanvas *c1;
   gStyle->SetEndErrorSize(0);
//   TLatex chiText;

   // Comparison of different observables on the same plot
   if(otherSettings[0])
   {
      cout << "nrfiles = " << *nrfiles << endl;
      cout << "nrobs = " << *nrobs << endl;
      cout << "nrtrees = " << *nrtrees << endl;

      // Enough trees are needed to set colors
      if(*nrtrees < *nrfiles)
         return -3;
      if((*nrfiles == 1) && (*nrobs != *nrtrees))
         return -3;

      // Use two observables from a single file or one observable from two files to compare (always only one tree)
      if(*nrfiles > 3)
         return -2;
      else if((*nrobs == 2) && (*nrfiles == 1))
         cout << "Plotting two different observables from a single file." << endl;
      else if((*nrobs == 1) && (*nrfiles >= 2))
         cout << "Plotting the same observable from two or three files." << endl;
      else
      {
//         if((*nrobs > 1) && (*nrfiles >= 2))
            return -1;
/*	 else if((*nrfiles > 1) && (*nrobs == 2))
            return -1;*/
      }

      // Only works for histograms (at the moment)
      if(type == 1)
         return 0;

      for(int i = 0; i < *nrfiles; i++)
      {
         stemp[0] = "mkdir -p " + RemoveFilename(&(filename->at(i))) + "/histograms";
         system(stemp[0].c_str());
      }
     
      c1 = new TCanvas("c1","",1200,900);
      obsvars = new float[3*MAXINFILES];
      obsvars_neg = new float[3*MAXINFILES];
      obsvars_pos = new float[3*MAXINFILES];
     
      max = new double[*nrobs];
      min = new double[*nrobs];
     
      for(int i = 0; i < *nrfiles; i++)
      {
         cout << "Running through a new input file (" << i << "): " << filename->at(i) << endl;
         ifile[i] = new TFile((filename->at(i)).c_str(), "READ");
         itemp[0] = GetRootKeys(ifile[i], "TreeS");
      }
      
      // This variable is now purely for running through number of files or number of observables
      if(*nrfiles >= 2)
         *nrtrees = *nrfiles;
      else
         *nrtrees = *nrobs;
      cout << "nrtrees = " << *nrtrees << endl;

      cout << "Selected observables for plotting histograms are (" << *nrobs << "):" << endl;
      for(int j = 0; j < (*nrobs); j++)
         cout << "- " << obser->at(j) << endl;
      
      // Zero holders for underflow, overflow and total events (from each observable/file + all together)
      for(int j = 0; j < (*nrtrees); j++)
         totcount[0][j] = 0;
      totcount[0][*nrtrees] = 0;
     
      // Only the first tree is taken
      cout << "Selected trees for plotting are (" << *nrfiles << "):" << endl;
      stemp[0] = "TreeS" + ToString(trees->at(0));
      if(*nrfiles >= 2)
      {
         for(int i = 0; i < *nrfiles; i++)
	 {
            getTree[i] = (TTree*)ifile[i]->Get(stemp[0].c_str());
            cout << "- name = " << getTree[i]->GetName() << ", title = " << getTree[i]->GetTitle() << ", nr. events = " << getTree[i]->GetEntries() << endl;
	 }
      }
      else
      {
         getTree[0] = (TTree*)ifile[0]->Get(stemp[0].c_str());
         getTree[1] = (TTree*)ifile[0]->Get(stemp[0].c_str());
         cout << "- name = " << getTree[0]->GetName() << ", title = " << getTree[0]->GetTitle() << ", nr. events = " << getTree[0]->GetEntries() << endl;
      }

      // obsvars, treeHist and treeAGraph are setup, so that the iterator is = j
     
      // Set addresses for reading observables
      for(int j = 0; j < (*nrtrees); j++)
      {
         if(*nrfiles >= 2)
            stemp[0] = obser->at(0);
	 else
            stemp[0] = obser->at(j);

         cout << "Setting branch addresses for: " << stemp[0];
         getTree[j]->SetBranchAddress(stemp[0].c_str(), &obsvars[j]);
         if(stemp[0] != "MVA")
         {
            stemp[1] = stemp[0] + "_neg";
            cout << ", " << stemp[1];
            getTree[j]->SetBranchAddress(stemp[1].c_str(), &obsvars_neg[j]);
            stemp[1] = stemp[0] + "_pos";
            cout << ", " << stemp[1];
            getTree[j]->SetBranchAddress(stemp[1].c_str(), &obsvars_pos[j]);
         }
         cout << endl;
      }

      // Do the same for all runs
      for(int j = 0; j < (*nrtrees); j++)
      {
         stemp[0] = "treeHist" + ToString(j);
         if(*nrfiles >= 2)
            stemp[1] = obser->at(0);
	 else
            stemp[1] = obser->at(j);

         if(stemp[1] == "MVA")
         {
            dtemp[0] = mvalim[0];
            dtemp[1] = mvalim[1];
         }
         else
         {
            SetRangeSpecial(dtemp, genObser, getTree[j], &stemp[1]);
/*            dtemp[0] = genObser->GetMin(stemp[1]);
            dtemp[1] = genObser->GetMax(stemp[1]);*/

	    if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
	    {
               dtemp[0] = dtemp[0]/1.e+18;
               dtemp[1] = dtemp[1]/1.e+18;
	    }
         }

         treeHist[j] = new TH1F(stemp[0].c_str(), stemp[1].c_str(), *nrbins, dtemp[0], dtemp[1]);
#if ROOTVER == 5
         treeHist[j]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
         treeHist[j]->SetCanExtend(TH1::kXaxis);
#endif
         
         // Set underflow and overflow counters to 0
         max[0] = 0;
         uflowcount[0] = 0;
         oflowcount[0] = 0;
      }

      // Loop over all runs
      for(int j = 0; j < (*nrtrees); j++)
      {
         // Loop over all entries in a tree
         for(int ievt = 0; ievt < getTree[j]->GetEntries(); ievt++)
         {
            getTree[j]->GetEntry(ievt);

            if(*nrfiles >= 2)
               stemp[1] = obser->at(0);
	    else
               stemp[1] = obser->at(j);

            // Check if observable is invalid and has value -1
            if(stemp[1] == "MVA")
               isgood[0] = true;
            else
            {
               if(obsvars[j] == -1)
                  isgood[0] = false;
               else
                  isgood[0] = true;
            }
   
            // Save observable values into the histograms
            if(isgood[0])
            {
               // Take mean
               dtemp[2] = obsvars[j];

               if(stemp[1] == "MVA")
               {
                  dtemp[0] = mvalim[0];
                  dtemp[1] = mvalim[1];
               }
               else
               {
                  SetRangeSpecial(dtemp, genObser, getTree[j], &stemp[1]);
/*                  dtemp[0] = genObser->GetMin(stemp[1]);
                  dtemp[1] = genObser->GetMax(stemp[1]);*/

                  // Take mean + negative uncertainties
                  if(otherSettings[1])
                     dtemp[2] -= obsvars_neg[j];
                  // Take mean + positive uncertainties
                  if(otherSettings[2])
                     dtemp[2] += obsvars_pos[j];
               }

               if(dtemp[2] > dtemp[1])
                  oflowcount[0]++;
               else if(dtemp[2] < dtemp[0])
                  uflowcount[0]++;
               else
	       {
	          if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
                     dtemp[2] = dtemp[2]/1.e+18;

                  treeHist[j]->Fill(dtemp[2]);
//                  cout << "Filling value = " << dtemp[2] << endl;
	       }
   
               totcount[0][j]++;
               totcount[0][*nrtrees]++;
            }
         } // Loop over all entries
      } // Loop over all runs
   
      // Write the valid event information
      cout << "Valid events:" << endl;
      for(int j = 0; j < (*nrtrees); j++)
      {
         if(*nrfiles >= 2)
            stemp[1] = obser->at(0);
	 else
            stemp[1] = obser->at(j);

         cout << " - " << stemp[1] << " = " << totcount[0][j] << " (UFLOW=" << uflowcount[0] << ", OFLOW=" << oflowcount[0] << ")" << endl;
      }
      cout << endl;

      // Plot histograms
      c1->cd();
      mystyle->SetSinglePlot(0, -1, c1);

      // Prepare logarithmic Y-axis for energy
      if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
         c1->SetLogy(true);
      else
         c1->SetLogy(false);

      min[0] = 1;

      for(int j = 0; j < (*nrtrees); j++)
      {
         if(*nrfiles >= 2)
            stemp[1] = obser->at(0);
	 else
            stemp[1] = obser->at(j);

         if(stemp[1] == "MVA")
         {
            dtemp[0] = mvalim[0];
            dtemp[1] = mvalim[1];
            stemp[2] = "MVA variable";
         }
         else
         {
            SetRangeSpecial(dtemp, genObser, getTree[0], &stemp[1]);
/*            dtemp[0] = genObser->GetMin(stemp[1]);
            dtemp[1] = genObser->GetMax(stemp[1]);*/
            stemp[2] = genObser->GetLabel(stemp[1]);

            if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
	    {
               dtemp[0] = dtemp[0]/1.e+18;
               dtemp[1] = dtemp[1]/1.e+18;
	       ReplacePart(&stemp[2], "(eV)", "(EeV)");
	    }
	    cout << "Label = " << stemp[2] << endl;
         }

         // Scale all to get a distribution
         treeHist[j]->Scale(1./totcount[0][j]);
   
         // Check maximum values in a histogram
         if(treeHist[j]->GetMaximum() > max[0])
            max[0] = treeHist[j]->GetMaximum();

	 // Set minimum, if we have a logarithmic axis
         if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
	 {
            if(0.5/totcount[0][j] < min[0])
               min[0] = 0.5/totcount[0][j];
	 }
	 else
            min[0] = 0;

	 cout << "Minimum, " << j << " = " << min[0] << endl;
	 cout << "Maximum, " << j << " = " << max[0] << endl;
   
         // Plot the histogram
         if(j == 0)
            treeHist[j]->Draw();
         else
            treeHist[j]->Draw("same");
      }
   
      for(int j = 0; j < (*nrtrees); j++)
      {
         // Setup other plotting options
         mystyle->SetHistColor((TH1*)treeHist[j], treeColor->at(j));
   
         treeHist[j]->GetXaxis()->SetRange(dtemp[0], dtemp[1]);
         treeHist[j]->GetXaxis()->SetRangeUser(dtemp[0], dtemp[1]);
         if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
	 {
            treeHist[j]->GetYaxis()->SetRangeUser(min[0], TMath::Power(10, 1.2*max[0]));
            treeHist[j]->SetMaximum(TMath::Power(10, 1.2*max[0]));
	 }
	 else
	 {
            treeHist[j]->GetYaxis()->SetRangeUser(min[0], 1.2*max[0]);
            treeHist[j]->SetMaximum(1.2*max[0]);
	 }
         mystyle->SetAxisTitles((TH1*)treeHist[j], stemp[2], "Number of events (normalized)");
         treeHist[j]->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
         treeHist[j]->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
      }
   
      cout << "Preparing a legend" << endl;
      // Draw a legend with into on the number of event of all trees
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(*nrtrees)), gPad->GetLeftMargin()+.50, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);

      string *slabels = new string[*nrtrees];
      for(int j = 0; j < (*nrtrees); j++)
      {
         if(*nrfiles >= 2)
	 {
	    slabels[j] = "Set legend text for file " + ToString(j+1) + ": ";
	    itemp[0] = 0;
	 }
	 else
	 {
	    slabels[j] = "Set legend text for observable " + ToString(j+1) + ": ";
	    itemp[0] = j;
	 }

	 if(treeHist[j]->GetMean() >= 10000.)
            stemp[j] = string(getTree[j]->GetTitle()) + ", " + obser->at(itemp[0]) + " [N = " + ToString(totcount[0][j]) + ", mean = " + ToSciString(treeHist[j]->GetMean(), 2) + "]";
	 else
            stemp[j] = string(getTree[j]->GetTitle()) + ", " + obser->at(itemp[0]) + " [N = " + ToString(totcount[0][j]) + ", mean = " + ToString(treeHist[j]->GetMean(), 2) + "]";

/*         if(*nrfiles >= 2)
	 {
	    slabels[j] = "Set legend text for file " + ToString(j+1) + ": ";
	    if(treeHist[j]->GetMean() >= 10000.)
               stemp[j] = string(getTree[j]->GetTitle()) + ", " + obser->at(0) + " [N = " + ToString(totcount[0][j]) + ", mean = " + ToSciString(treeHist[j]->GetMean(), 2) + "]";
	    else
               stemp[j] = string(getTree[j]->GetTitle()) + ", " + obser->at(0) + " [N = " + ToString(totcount[0][j]) + ", mean = " + ToString(treeHist[j]->GetMean(), 2) + "]";
	 }
	 else
	 {
	    slabels[j] = "Set legend text for observable " + ToString(j+1) + ": ";
            stemp[j] = string(getTree[j]->GetTitle()) + ", " + obser->at(j) + " [N = " + ToString(totcount[0][j]) + "]";
	 }*/
      }
      itemp[0] = 1110;

      legendselectDialog = new FSDialog(wxT("Legend title selection"), wxSize(600,150+30*(*nrtrees)), "Select titles for both observables/files that will appear in the legend of the observable/file comparison plot.\n", slabels, *nrtrees, stemp, &itemp[0]);
      if(legendselectDialog->ShowModal() == wxID_OK)
      {
         for(int j = 0; j < (*nrtrees); j++)
	 {
/*            stemp[j+1] = legendselectDialog->GetFSValue(j);

	    if(treeHist[j]->GetMean() >= 10000.)
               stemp[0] = stemp[j+1] + " [N = " + ToString(totcount[0][j]) + ", mean = " + ToSciString(treeHist[j]->GetMean(), 2) + "]";
            else
               stemp[0] = stemp[j+1] + " [N = " + ToString(totcount[0][j]) + ", mean = " + ToString(treeHist[j]->GetMean(), 2) + "]";
            legend->AddEntry(treeHist[j], stemp[0].c_str(), "f");*/

            stemp[j+1] = legendselectDialog->GetFSValue(j);
            legend->AddEntry(treeHist[j], stemp[j+1].c_str(), "f");
	 }
      }
      legendselectDialog->Destroy();
      delete legendselectDialog;
      delete[] slabels;

      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      stemp[1] = obser->at(0);

      // Write down the underflow and overflow information
      uoflowtext->SetTextSize(0.017);
      uoflowtext->SetTextAlign(11);
      stemp[0] = "UFLOW = " + ToString(uflowcount[0]);
      if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
         uoflowtext->DrawLatex(dtemp[0], TMath::Power(10, max[0]*1.205), stemp[0].c_str());
      else
         uoflowtext->DrawLatex(dtemp[0], max[0]*1.205, stemp[0].c_str());
      uoflowtext->SetTextAlign(31);
      stemp[0] = "OFLOW = " + ToString(oflowcount[0]);
      if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
         uoflowtext->DrawLatex(dtemp[1], TMath::Power(10, max[0]*1.205), stemp[0].c_str());
      else
         uoflowtext->DrawLatex(dtemp[1], max[0]*1.205, stemp[0].c_str());

      stemp[2] = "comparehist_TreeS" + ToString(trees->at(0));
      if(*nrfiles >= 2)
      {
	 for(int i = 0; i < *nrfiles; i++)
	 {
            stemp[0] = RemoveFilename(&(filename->at(i))) + "/histograms/" + stemp[2] + "_" + obser->at(0) + ".pdf";
            c1->SaveAs(stemp[0].c_str());
            cout << "Saving histogram to file (" << i << "): " << stemp[0] << endl;
	 }
      }
      else
      {
         stemp[0] = RemoveFilename(&(filename->at(0))) + "/histograms/" + stemp[2] + "_" + obser->at(0) + "-" + obser->at(1) + ".pdf";
         c1->SaveAs(stemp[0].c_str());
         cout << "Saving histogram to file: " << stemp[0] << endl;
      }

      for(int j = 0; j < (*nrtrees); j++)
         delete treeHist[j];
   
      delete[] max;
      delete[] obsvars;
      delete[] obsvars_neg;
      delete[] obsvars_pos;

      if(*nrfiles >= 2)
      {
	 for(int i = 0; i < *nrfiles; i++)
            delete ifile[i];
      }
      else
         delete ifile[0];
      delete c1;
   }
   // Don't mix observables on the same plot
   else
   {
      // Run over all selected files
      for(int i = 0; i < filename->size(); i++)
      {
         // Create histogram plot directory, if it doesn't exist already
         stemp[0] = "mkdir -p " + RemoveFilename(&(filename->at(i))) + "/histograms";
         system(stemp[0].c_str());
     
         c1 = new TCanvas("c1","",1200,900);
         obsvars = new float[3*MAXINFILES];
         obsvars_neg = new float[3*MAXINFILES];
         obsvars_pos = new float[3*MAXINFILES];

	 if(type == 1)
	    yrange = new double[2*(*nrobs)];
     
         max = new double[*nrobs];
         min = new double[*nrobs];
     
         cout << "Running through a new input file: " << filename->at(i) << endl;
         ifile[0] = new TFile((filename->at(i)).c_str(), "READ");
         itemp[0] = GetRootKeys(ifile[0], "TreeS");
     
	 if(type == 0)
	 {
            cout << "Selected observables for plotting histograms are (" << *nrobs << "):" << endl;
            for(int j = 0; j < (*nrobs); j++)
               cout << "- " << obser->at(j) << endl;
	 }
	 else if(type == 1)
	 {
            cout << "Selected observables for plotting scatter plots are (" << *nrobs << "):" << endl;
            itemp[1] = 0;
            for(int j = 0; j < (*nrobs); j++)
	    {
               for(int j2 = 0; j2 < (*nrobs); j2++)
	       {
		  if(j2 > j)
		  {
                     cout << "- " << obser->at(j) << " (" << j << ") vs. " << obser->at(j2) << " (" << j2 << "), plot = " << itemp[1] << endl;
		     itemp[1]++;
		  }
	       }
	    }
	 }
     
         // Zero holders for underflow, overflow and total events (from each tree + all together)
         for(int j = 0; j < (*nrobs); j++)
         {
            for(int k = 0; k < (*nrtrees); k++)
               totcount[k][j] = 0;
            totcount[*nrtrees][j] = 0;
         }
     
         // Choose the trees that were selected
         cout << "Selected trees for plotting are (" << *nrtrees << "):" << endl;
         for(int k = 0; k < (*nrtrees); k++)
         {
            stemp[0] = "TreeS" + ToString(trees->at(k));
            getTree[k] = (TTree*)ifile[0]->Get(stemp[0].c_str());
            cout << "- name = " << getTree[k]->GetName() << ", title = " << getTree[k]->GetTitle() << ", nr. events = " << getTree[k]->GetEntries() << endl;
         }
     
         // obsvars, treeHist and treeAGraph are setup, so that the iterator is = j*(*nrtrees)+k
     
         // Set addresses for reading observables
         for(int j = 0; j < (*nrobs); j++)
         {
     	    stemp[0] = obser->at(j);

	    cout << "Setting branch addresses for: " << stemp[0];
            // Do the same for all trees
            for(int k = 0; k < (*nrtrees); k++)
            {
               getTree[k]->SetBranchAddress(stemp[0].c_str(), &obsvars[j*(*nrtrees)+k]);
               if(stemp[0] != "MVA")
	       {
                  stemp[1] = stemp[0] + "_neg";
		  if(k == 0)
	             cout << ", " << stemp[1];
                  getTree[k]->SetBranchAddress(stemp[1].c_str(), &obsvars_neg[j*(*nrtrees)+k]);
                  stemp[1] = stemp[0] + "_pos";
		  if(k == 0)
	             cout << ", " << stemp[1];
                  getTree[k]->SetBranchAddress(stemp[1].c_str(), &obsvars_pos[j*(*nrtrees)+k]);
	       }
            }
	    cout << endl;
         }
     
         itemp[1] = 0;
         // Prepare histograms for all observables
         for(int j = 0; j < (*nrobs); j++)
         {
            // Do the same for all trees
            for(int k = 0; k < (*nrtrees); k++)
            {
               // Setting up the tree histograms
               if(type == 0)
	       {
                  stemp[0] = "treeHist" + ToString(j*(*nrtrees)+k);
     	          stemp[1] = obser->at(j);

		  if(stemp[1] == "MVA")
                  {
                     dtemp[0] = mvalim[0];
		     dtemp[1] = mvalim[1];
		  }
		  else
		  {
                     SetRangeSpecial(dtemp, genObser, getTree[k], &stemp[1]);
/*                     dtemp[0] = genObser->GetMin(stemp[1]);
		     dtemp[1] = genObser->GetMax(stemp[1]);*/

//		     cout << "Set range = " << dtemp[0] << "\t" << dtemp[1] << endl;

	             if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
                     {
                        dtemp[0] = dtemp[0]/1.e+18;
                        dtemp[1] = dtemp[1]/1.e+18;
		     }
		  }

                  treeHist[j*(*nrtrees)+k] = new TH1F(stemp[0].c_str(), stemp[1].c_str(), *nrbins, dtemp[0], dtemp[1]);
#if ROOTVER == 5
                  treeHist[j*(*nrtrees)+k]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
                  treeHist[j*(*nrtrees)+k]->SetCanExtend(TH1::kXaxis);
#endif
	       }
               // Setting up the tree scatter plots
	       else if(type == 1)
	       {
                  // Prepare histograms for all observables
                  for(int j2 = 0; j2 < (*nrobs); j2++)
                  {
                     if(j2 > j)
		     {
                        stemp[0] = "treeAGraph" + ToString(itemp[1]);
			cout << "i = " << i << ", j = " << j << ", k = " << k << ", j2 = " << j2 << ", Name = " << stemp[0] << endl;
                        treeAGraph[itemp[1]] = new TGraphAsymmErrors();
		        treeAGraph[itemp[1]]->SetName(stemp[0].c_str());
			itemp[1]++;
		     }
		  }
	       }
            }
   
            // Set underflow and overflow counters to 0
            max[j] = 0;
            uflowcount[j] = 0;
            oflowcount[j] = 0;

	    if(type == 1)
	    {
	       yrange[2*j] = 1.e+18;
	       yrange[2*j+1] = -1.e+18;
	    }
         }

	 itemp[1] = -1;
   
	 // Populate histograms
         if(type == 0)
         {
            // Loop over all trees
            for(int k = 0; k < (*nrtrees); k++)
            {
               // Loop over all entries in a tree
               for(int ievt = 0; ievt < getTree[k]->GetEntries(); ievt++)
               {
                  getTree[k]->GetEntry(ievt);
   
                  // Loop over all observables
                  for(int j = 0; j < (*nrobs); j++)
                  {
                     stemp[1] = obser->at(j);

                     // Check if observable is invalid and has value -1
		     if(stemp[1] == "MVA")
		     {
                        isgood[0] = true;
		     }
		     else
		     {
                        if(obsvars[j*(*nrtrees)+k] == -1)
                           isgood[0] = false;
                        else
                           isgood[0] = true;
		     }
   
                     // Save observable values into the histograms
                     if(isgood[0])
                     {
		        // Take mean
                        dtemp[2] = obsvars[j*(*nrtrees)+k];

		        if(stemp[1] == "MVA")
		        {
                           dtemp[0] = mvalim[0];
		           dtemp[1] = mvalim[1];
		        }
		        else
		        {
                           SetRangeSpecial(dtemp, genObser, getTree[k], &stemp[1]);
/*                           dtemp[0] = genObser->GetMin(stemp[1]);
		           dtemp[1] = genObser->GetMax(stemp[1]);*/

		           // Take mean + negative uncertainties
                           if(otherSettings[1])
                              dtemp[2] -= obsvars_neg[j*(*nrtrees)+k];
		           // Take mean + positive uncertainties
                           if(otherSettings[2])
                              dtemp[2] += obsvars_pos[j*(*nrtrees)+k];
		        }

                        if(dtemp[2] > dtemp[1])
                           oflowcount[j]++;
                        else if(dtemp[2] < dtemp[0])
                           uflowcount[j]++;
                        else
                        {
	                   if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
                              dtemp[2] = dtemp[2]/1.e+18;

                           treeHist[j*(*nrtrees)+k]->Fill(dtemp[2]);
//                           cout << "Filling value = " << dtemp[2] << endl;
			}
   
                        totcount[k][j]++;
                        totcount[*nrtrees][j]++;
                     }
		  } // Loop over all observables
               } // Loop over all entries
            } // Loop over all trees
         }
	 // Populate scatter plots
         else if(type == 1)
         {
            for(int j = 0; j < (*nrobs); j++)
	    {
               stemp[1] = obser->at(j);

               for(int j2 = 0; j2 < (*nrobs); j2++)
               {
                  stemp[2] = obser->at(j2);

                  if(j2 > j)
                  {
                     for(int k = 0; k < (*nrtrees); k++)
	             {
                        itemp[2] = 0;
                        for(int ievt = 0; ievt < getTree[k]->GetEntries(); ievt++)
	                {
                           getTree[k]->GetEntry(ievt);
         
                           // Check if observable is invalid and has value -1
                           if(stemp[1] == "MVA")
                              isgood[0] = true;
                           else
                           {
                              if(obsvars[j*(*nrtrees)+k] == -1)
                                 isgood[0] = false;
                              else
                                 isgood[0] = true;
                           }
         
                           // Check if observable 2 is invalid and has value -1
                           if(stemp[2] == "MVA")
                              isgood[1] = true;
                           else
                           {
                              if(obsvars[j2*(*nrtrees)+k] == -1)
                                 isgood[1] = false;
                              else
                                 isgood[1] = true;
                           }
         
                           // Save observable values into the scatter plots
                           if(isgood[0] && isgood[1])
                           {
                              if(itemp[2] == 0)
                              {
                                 itemp[1]++;
			         cout << "Using graph = " << itemp[1] << " (i = " << i << ", j = " << j << ", k = " << k << ", j2 = " << j2 << ")" << endl;
                              }

                              // Take means
                              dtemp[2] = obsvars[j*(*nrtrees)+k];
                              dtemp[3] = obsvars[j2*(*nrtrees)+k];
         
                              // Get minimum and maximum values for observables and set plotting points
                              if(stemp[1] == "MVA")
                              {
                                 if(dtemp[2] < yrange[2*j])
                                    yrange[2*j] = dtemp[2];
                                 if(dtemp[2] > yrange[2*j+1])
                                    yrange[2*j+1] = dtemp[2];

                                 if(dtemp[3] < yrange[2*j2])
                                    yrange[2*j2] = dtemp[3] - obsvars_neg[j2*(*nrtrees)+k];
                                 if(dtemp[3] > yrange[2*j2+1])
                                    yrange[2*j2+1] = dtemp[3] + obsvars_pos[j2*(*nrtrees)+k];
         
                                 treeAGraph[itemp[1]]->SetPoint(itemp[2], dtemp[2], dtemp[3]);
                                 treeAGraph[itemp[1]]->SetPointError(itemp[2], 0., 0., obsvars_neg[j2*(*nrtrees)+k], obsvars_pos[j2*(*nrtrees)+k]);
                                 cout << "Setting point " << itemp[2] << " (" << dtemp[2] << ", " << dtemp[3] << ") to graph " << itemp[1] << endl;
                              }
                              else if(stemp[2] == "MVA")
                              {
                                 if(dtemp[2] < yrange[2*j])
                                    yrange[2*j] = dtemp[2] - obsvars_neg[j*(*nrtrees)+k];
                                 if(dtemp[2] > yrange[2*j+1])
                                    yrange[2*j+1] = dtemp[2] + obsvars_pos[j*(*nrtrees)+k];

                                 if(dtemp[3] < yrange[2*j2])
                                    yrange[2*j2] = dtemp[3];
                                 if(dtemp[3] > yrange[2*j2+1])
                                    yrange[2*j2+1] = dtemp[3];
         
                                 treeAGraph[itemp[1]]->SetPoint(itemp[2], dtemp[3], dtemp[2]);
                                 treeAGraph[itemp[1]]->SetPointError(itemp[2], 0., 0., obsvars_neg[j*(*nrtrees)+k], obsvars_pos[j*(*nrtrees)+k]);
                                 cout << "Setting point " << itemp[2] << " (" << dtemp[3] << ", " << dtemp[2] << ") to graph " << itemp[1] << endl;
                              }
                              else
                              {
                                 if(dtemp[2] < yrange[2*j])
                                    yrange[2*j] = dtemp[2] - obsvars_neg[j*(*nrtrees)+k];
                                 if(dtemp[2] > yrange[2*j+1])
                                    yrange[2*j+1] = dtemp[2] + obsvars_pos[j*(*nrtrees)+k];
         
                                 if(dtemp[3] < yrange[2*j2])
                                    yrange[2*j2] = dtemp[3] - obsvars_neg[j2*(*nrtrees)+k];
                                 if(dtemp[3] > yrange[2*j2+1])
                                    yrange[2*j2+1] = dtemp[3] + obsvars_pos[j2*(*nrtrees)+k];
			      
				 if(InvertAxis(&stemp[1], &stemp[2]))
				 {
                                    treeAGraph[itemp[1]]->SetPoint(itemp[2], dtemp[3], dtemp[2]);
                                    treeAGraph[itemp[1]]->SetPointError(itemp[2], obsvars_neg[j2*(*nrtrees)+k], obsvars_pos[j2*(*nrtrees)+k], obsvars_neg[j*(*nrtrees)+k], obsvars_pos[j*(*nrtrees)+k]);
                                    cout << "Setting point " << itemp[2] << " (" << dtemp[3] << ", " << dtemp[2] << ") to graph " << itemp[1] << endl;
				 }
				 else
				 {
                                    treeAGraph[itemp[1]]->SetPoint(itemp[2], dtemp[2], dtemp[3]);
                                    treeAGraph[itemp[1]]->SetPointError(itemp[2], obsvars_neg[j*(*nrtrees)+k], obsvars_pos[j*(*nrtrees)+k], obsvars_neg[j2*(*nrtrees)+k], obsvars_pos[j2*(*nrtrees)+k]);
                                    cout << "Setting point " << itemp[2] << " (" << dtemp[2] << ", " << dtemp[3] << ") to graph " << itemp[1] << endl;
				 }
                              }
         
			      if((j+1 == j2) && (j2+1 == (*nrobs)))
			      {
                                 totcount[k][j]++;
                                 totcount[*nrtrees][j]++;
                                 totcount[k][j2]++;
                                 totcount[*nrtrees][j2]++;
			      }
			      else if(j+1 == j2)
			      {
                                 totcount[k][j]++;
                                 totcount[*nrtrees][j]++;
			      }
         
                              itemp[2]++;
                           }
	                } // Loop over all entries
                     } // Loop over all trees
                  } // Check j2 > j
               } // Loop over all observables (2)
            } // Loop over all observables
         }
   
         // Write the valid event information
         for(int k = 0; k < (*nrtrees); k++)
         {
            cout << "Valid events (" << getTree[k]->GetTitle() << "):" << endl;
            for(int j = 0; j < (*nrobs); j++)
	    {
               if(type == 0)
                  cout << " - " << obser->at(j) << " = " << totcount[k][j] << " (UFLOW=" << uflowcount[j] << ", OFLOW=" << oflowcount[j] << ")" << endl;
               else if(type == 1)
                  cout << " - " << obser->at(j) << " = " << totcount[k][j] << " (yrange min = " << yrange[2*j] << ", yrange max = " << yrange[2*j+1] << ")" << endl;
	    }
            cout << endl;
         }
   
	 itemp[1] = 0;
         // Plot all observables (with all selected trees in the same plot)
         // Plot histograms
         if(type == 0)
	 {
            for(int j = 0; j < (*nrobs); j++)
            {
               c1->cd();
               mystyle->SetSinglePlot(0, -1, c1);
               stemp[1] = obser->at(j);

	       // Prepare logarithmic Y-axis for energy
	       if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
                  c1->SetLogy(true);
	       else
                  c1->SetLogy(false);

               min[j] = 1;

               if(stemp[1] == "MVA")
               {
                  dtemp[0] = mvalim[0];
                  dtemp[1] = mvalim[1];
	          stemp[2] = "MVA variable";
               }
               else
               {
                  SetRangeSpecial(dtemp, genObser, getTree[0], &stemp[1]);
/*                  dtemp[0] = genObser->GetMin(stemp[1]);
                  dtemp[1] = genObser->GetMax(stemp[1]);*/
                  stemp[2] = genObser->GetLabel(stemp[1]);

	          if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
		  {
                     dtemp[0] = dtemp[0]/1.e+18;
                     dtemp[1] = dtemp[1]/1.e+18;
		     ReplacePart(&stemp[2], "(eV)", "(EeV)");
		  }
               }

               // Go over all trees
               for(int k = 0; k < (*nrtrees); k++)
               {
                  // Scale all to get a distribution
                  treeHist[j*(*nrtrees)+k]->Scale(1./totcount[k][j]);
   
                  // Check maximum values in a histogram
                  if(treeHist[j*(*nrtrees)+k]->GetMaximum() > max[j])
                     max[j] = treeHist[j*(*nrtrees)+k]->GetMaximum();

		  // Set minimum, if we have logarithmic y-axis
	          if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
		  {
                     if(0.5/totcount[k][j] < min[j])
                        min[j] = 0.5/totcount[k][j];
		  }
		  else
                     min[j] = 0;

		  cout << "Minimum " << j << ", " << k << " = " << min[j] << endl;
		  cout << "Maximum " << j << ", " << k << " = " << max[j] << endl;
   
                  // Plot the histogram
                  if(k == 0)
                     treeHist[j*(*nrtrees)+k]->Draw();
                  else
                     treeHist[j*(*nrtrees)+k]->Draw("same");
               }
   
               // Setup other plotting options
               for(int k = 0; k < (*nrtrees); k++)
               {
                  mystyle->SetHistColor((TH1*)treeHist[j*(*nrtrees)+k], treeColor->at(k));
   
                  treeHist[j*(*nrtrees)+k]->GetXaxis()->SetRange(dtemp[0], dtemp[1]);
                  treeHist[j*(*nrtrees)+k]->GetXaxis()->SetRangeUser(dtemp[0], dtemp[1]);
	          if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
		  {
                     treeHist[j*(*nrtrees)+k]->GetYaxis()->SetRangeUser(min[j], TMath::Power(10, 1.2*max[j]));
                     treeHist[j*(*nrtrees)+k]->SetMaximum(TMath::Power(10, 1.2*max[j]));
		  }
                  else
		  {
                     treeHist[j*(*nrtrees)+k]->GetYaxis()->SetRangeUser(min[j], 1.2*max[j]);
                     treeHist[j*(*nrtrees)+k]->SetMaximum(1.2*max[j]);
		  }
                  mystyle->SetAxisTitles((TH1*)treeHist[j*(*nrtrees)+k], stemp[2], "Number of events (normalized)");
                  treeHist[j*(*nrtrees)+k]->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
                  treeHist[j*(*nrtrees)+k]->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
               }
   
               // Draw a legend with into on the number of event of all trees
               legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(*nrtrees)), gPad->GetLeftMargin()+.50, 1-gPad->GetTopMargin());
               legend->SetFillStyle(legendFill);
               legend->SetFillColor(c_MvaCut);

/*	       // Ask the user, if he wishes to change the legend titles
               string *slabels = new string[*nrtrees];
               for(int k = 0; k < (*nrtrees); k++)
               {
                  slabels[k] = "Set legend text for tree " + ToString(k+1) + ": ";
                  if(treeHist[j*(*nrtrees)+k]->GetMean() >= 10000.)
                     stemp[k] = string(getTree[k]->GetTitle()) + ", " + obser->at(j) + " [N = " + ToString(totcount[k][j]) + ", mean = " + ToSciString(treeHist[j*(*nrtrees)+k]->GetMean(), 2) + "]";
                  else
                     stemp[k] = string(getTree[k]->GetTitle()) + ", " + obser->at(j) + " [N = " + ToString(totcount[k][j]) + ", mean = " + ToString(treeHist[j*(*nrtrees)+k]->GetMean(), 2) + "]";
               }
               itemp[0] = 1110;

               legendselectDialog = new FSDialog(wxT("Legend title selection"), wxSize(600,150+30*(*nrtrees)), "Select titles for both observables/files that will appear in the legend of the observable/file comparison plot.\n", slabels, *nrtrees, stemp, &itemp[0]);
               if(legendselectDialog->ShowModal() == wxID_OK)
               {
                  for(int k = 0; k < (*nrtrees); k++)
		  {
                     stemp[k+1] = legendselectDialog->GetFSValue(k);
                     legend->AddEntry(treeHist[j*(*nrtrees)+k], stemp[k+1].c_str(), "f");
		  }
               }
               legendselectDialog->Destroy();
               delete legendselectDialog;
               delete[] slabels;*/

               for(int k = 0; k < (*nrtrees); k++)
               {
	          if(treeHist[j*(*nrtrees)+k]->GetMean() >= 10000.)
                     stemp[0] = string(getTree[k]->GetTitle()) + " [N = " + ToString(totcount[k][j]) + ", mean = " + ToSciString(treeHist[j*(*nrtrees)+k]->GetMean(), 2) + "]";
                  else
                     stemp[0] = string(getTree[k]->GetTitle()) + " [N = " + ToString(totcount[k][j]) + ", mean = " + ToString(treeHist[j*(*nrtrees)+k]->GetMean(), 2) + "]";
                  legend->AddEntry(treeHist[j*(*nrtrees)+k], stemp[0].c_str(), "f");
               }
               legend->SetBorderSize(1);
               legend->SetMargin(0.3);
               legend->Draw("same");

               // Write down the underflow and overflow information
               uoflowtext->SetTextSize(0.017);
               uoflowtext->SetTextAlign(11);
               stemp[0] = "UFLOW = " + ToString(uflowcount[j]);
	       if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
                  uoflowtext->DrawLatex(dtemp[0], TMath::Power(10, 1.205*max[j]), stemp[0].c_str());
	       else
                  uoflowtext->DrawLatex(dtemp[0], 1.205*max[j], stemp[0].c_str());
               uoflowtext->SetTextAlign(31);
               stemp[0] = "OFLOW = " + ToString(oflowcount[j]);
	       if((stemp[1] == "energySD") || (stemp[1] == "energyFD"))
                  uoflowtext->DrawLatex(dtemp[1], TMath::Power(10, 1.205*max[j]), stemp[0].c_str());
	       else
                  uoflowtext->DrawLatex(dtemp[1], 1.205*max[j], stemp[0].c_str());

               stemp[2] = "mergedhist";
               for(int k = 0; k < (*nrtrees); k++)
                  stemp[2] += "_TreeS" + ToString(trees->at(k));
               stemp[0] = RemoveFilename(&(filename->at(i))) + "/histograms/" + stemp[2] + "_" + stemp[1] + ".pdf";
   
               c1->SaveAs(stemp[0].c_str());
               cout << "Saving histogram to file: " << stemp[0] << endl;
	    }
	 }
	 // Plot scatter plots
	 else if(type == 1)
	 {
            // Go over all observables
            for(int j = 0; j < (*nrobs); j++)
            {
               stemp[1] = obser->at(j);

               // Go over all observables (2)
               for(int j2 = 0; j2 < (*nrobs); j2++)
               {
                  stemp[2] = obser->at(j2);
                  if(j2 > j)
                  {
                     c1->cd();
                     mystyle->SetSinglePlot(0, 0, c1);
                     // Go over all trees
                     for(int k = 0; k < (*nrtrees); k++)
                     {
                        // Plot the scatter plots
                        if(k == 0)
                           treeAGraph[itemp[1]]->Draw("ap");
                        else
                           treeAGraph[itemp[1]+k]->Draw("p;same");
		     }
   
                     // Setup other plotting options
                     for(int k = 0; k < (*nrtrees); k++)
                     {
                        mystyle->SetGraphColor(treeAGraph[itemp[1]+k], treeColor->at(k));
         
			if(InvertAxis(&stemp[1], &stemp[2]))
			{
                           dtemp[0] = yrange[2*j];
                           dtemp[1] = yrange[2*j+1] - yrange[2*j];
			}
			else
			{
                           dtemp[0] = yrange[2*j2];
                           dtemp[1] = yrange[2*j2+1] - yrange[2*j2];
			}

                        treeAGraph[itemp[1]]->GetYaxis()->SetRange(-0.05*dtemp[1]+dtemp[0], 1.20*dtemp[1]+dtemp[0]);
                        treeAGraph[itemp[1]+k]->GetYaxis()->SetRangeUser(-0.05*dtemp[1]+dtemp[0], 1.20*dtemp[1]+dtemp[0]);
                        treeAGraph[itemp[1]+k]->SetMarkerStyle(1);
                        treeAGraph[itemp[1]+k]->SetLineWidth(1);
                        stemp[3] = genObser->GetLabel(stemp[1]);
                        stemp[4] = genObser->GetLabel(stemp[2]);

			if(InvertAxis(&stemp[1], &stemp[2]))
                           mystyle->SetAxisTitles(treeAGraph[itemp[1]+k], stemp[4], stemp[3]);
			else
                           mystyle->SetAxisTitles(treeAGraph[itemp[1]+k], stemp[3], stemp[4]);
                        treeAGraph[itemp[1]+k]->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c1));
                        treeAGraph[itemp[1]+k]->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c1));
                     }
   
                     // Draw a legend with into on the number of event of all trees
                     legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(*nrtrees)), gPad->GetLeftMargin()+.50, 1-gPad->GetTopMargin());
                     legend->SetFillStyle(legendFill);
                     legend->SetFillColor(c_MvaCut);
                     for(int k = 0; k < (*nrtrees); k++)
                     {
                        stemp[0] = string(getTree[k]->GetTitle()) + " [N = " + ToString(totcount[k][j]) + /*", mean = " + ToString(treeAGraph[itemp[1]]->GetMean(), 2) +*/ "]";
                        legend->AddEntry(treeAGraph[itemp[1]+k], stemp[0].c_str(), "l");
                     }
                     legend->SetBorderSize(1);
                     legend->SetMargin(0.3);
                     legend->Draw("same");

                     stemp[3] = "mergedscat";
                     for(int k = 0; k < (*nrtrees); k++)
                        stemp[3] += "_TreeS" + ToString(trees->at(k));

		     if(InvertAxis(&stemp[1], &stemp[2]))
                        stemp[0] = RemoveFilename(&(filename->at(i))) + "/histograms/" + stemp[3] + "_" + stemp[2] + "-" + stemp[1] + ".pdf";
		     else
                        stemp[0] = RemoveFilename(&(filename->at(i))) + "/histograms/" + stemp[3] + "_" + stemp[1] + "-" + stemp[2] + ".pdf";
   
                     c1->SaveAs(stemp[0].c_str());
                     cout << "Saving scatter plot to file: " << stemp[0] << endl;

		     itemp[1] += (*nrtrees);
	          }
	       }
	    }
         }
   
         itemp[1] = 0;
         for(int j = 0; j < (*nrobs); j++)
	 {
            for(int k = 0; k < (*nrtrees); k++)
	    {
               if(type == 0)
                  delete treeHist[j*(*nrtrees)+k];
               else if(type == 1)
	       {
                  for(int j2 = 0; j2 < (*nrobs); j2++)
	          {
                     if(j2 > j)
		     {
                        delete treeAGraph[itemp[1]];
			itemp[1]++;
		     }
		  }
	       }
	    }
	 }
   
         delete[] max;
         delete[] obsvars;
         delete[] obsvars_neg;
         delete[] obsvars_pos;

	 if(type == 1)
	 {
	    delete[] yrange;
	 }

         delete ifile[0];
         delete c1;
      }
   }

   return 0;
}

bool ScatHist::InvertAxis(string *obs1, string *obs2)
{
   if(*obs2 == "MVA")
      return true;

   bool *btemp = new bool[2];

   if(*obs1 == "energySD")
      btemp[0] = true;
   else if(*obs1 == "energyFD")
      btemp[0] = true;
   else if(*obs1 == "zenithSD")
      btemp[0] = true;
   else if(*obs1 == "zenithFD")
      btemp[0] = true;
   else
      btemp[0] = false;

   if(*obs2 == "energySD")
      btemp[1] = true;
   else if(*obs2 == "energyFD")
      btemp[1] = true;
   else if(*obs2 == "zenithSD")
      btemp[1] = true;
   else if(*obs2 == "zenithFD")
      btemp[1] = true;
   else
      btemp[1] = false;

   if(!btemp[0] && btemp[1])
   {
      delete[] btemp;
      return true;
   }
   else
   {
      delete[] btemp;
      return false;
   }
}

void ScatHist::SetRangeSpecial(double *range, Observables *genObser, TTree *tempTree, string *obsname)
{
   double *tempen = new double;

   *tempen = tempTree->GetMaximum("energyFD");
//   cout << "Energy maximum = " << *tempen << ", " << TMath::Log10(*tempen) << endl;
   *tempen = TMath::Log10(*tempen);

   if( (*obsname == "shwsize") || (*obsname == "deltas38") )
   {
      if(*tempen < 19.05)
      {
         if(*obsname == "shwsize")
         {
            range[0] = genObser->GetMin(*obsname);
            range[1] = genObser->GetMax(*obsname)/6.;
         }
      
         if(*obsname == "deltas38")
         {
            range[0] = genObser->GetMin(*obsname)/4.;
            range[1] = genObser->GetMax(*obsname)/3.;
         }

//         cout << "Set range (<19.0) = " << range[0] << "\t" << range[1] << endl;
      }
      else if(*tempen < 19.55)
      {
         if(*obsname == "shwsize")
         {
            range[0] = genObser->GetMin(*obsname);
            range[1] = genObser->GetMax(*obsname)/2.3;
         }
      
         if(*obsname == "deltas38")
         {
            range[0] = genObser->GetMin(*obsname)/2.;
            range[1] = genObser->GetMax(*obsname)/2.;
         }

//         cout << "Set range (<19.5) = " << range[0] << "\t" << range[1] << endl;
      }
      else
      {
         range[0] = genObser->GetMin(*obsname);
         range[1] = genObser->GetMax(*obsname);

//         cout << "Set range (>19.5) = " << range[0] << "\t" << range[1] << endl;
      }
   }
   else
   {
      range[0] = genObser->GetMin(*obsname);
      range[1] = genObser->GetMax(*obsname);

//      cout << "Set range (all) = " << range[0] << "\t" << range[1] << endl;
   }

   delete tempen;
}
