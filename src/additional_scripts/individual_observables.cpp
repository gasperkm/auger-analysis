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

   gSystem->Load("libTree.so");

   //-----------------------------------------
   ifstream ifs;
   string *stemp = new string[5];
   double *dtemp = new double[2];
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

   // Save the number of observables (+MVA) into itemp[1]
   itemp[1] = observables.size();
   cout << "Number of observables = " << itemp[1] << endl;

   ifs.close();

   //-----------------------------------------
   // Select observable to check
   int obscheck = -1;
   cerr << "Possible choices of observables:" << endl;
   for(int i = 0; i < observables.size(); i++)
      cerr << " - " << observables[i] << endl;
   cerr << endl << "Select the observable that you wish to check: ";
   cin >> stemp[3];

   for(int i = 0; i < observables.size(); i++)
   {
      if(observables[i] == stemp[3])
         obscheck = i;
   }

   if(obscheck == -1)
   {
      cerr << "Error! Invalid observable selected, rerun program." << endl;
      return 1;
   }

   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/individual_plots";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/individual_plots/purS_" + observables[obscheck] + "_*";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/individual_plots/apply_" + observables[obscheck] + "_*";
   system(stemp[0].c_str());

   TLegend *legend;
   int legendFill = 1001;

   TLatex uoflowtext;

   // Array for variable values of the observable we selected
   float *obsvars;
   obsvars = new float[3];

   // Vectors for observables from all trees (at most 10)
   vector<double> obsvect[10];
   vector<double> obsvect_neg[10];
   vector<double> obsvect_pos[10];
   // Average errors for observables from the Data tree
   double negerror = 0;
   double poserror = 0;
   // Information for all trees on title, mean value and observable range
   vector<string> treetitle;
   double mean[10];
   double minrange[10];
   double maxrange[10];
   double completerange[3] = {1.e+30, -1.e+30, -1.};

   // Erase all previous vectors values and mean value
   for(int k = 0; k < 10; k++)
   {
      if(!obsvect[k].empty())
         obsvect[k].erase(obsvect[k].begin(), obsvect[k].end());
      if(!obsvect_neg[k].empty())
         obsvect_neg[k].erase(obsvect_neg[k].begin(), obsvect_neg[k].end());
      if(!obsvect_pos[k].empty())
         obsvect_pos[k].erase(obsvect_pos[k].begin(), obsvect_pos[k].end());
      mean[k] = 0;
      minrange[k] = 1.e+30;
      maxrange[k] = -1.e+30;
   }
   if(!treetitle.empty())
      treetitle.erase(treetitle.begin(), treetitle.end());

   // Loop over all trees
   const int nbin_hist_high = 10000;
   int sign = 0;
   double range[2];
   bool found[3] = {false, false, false};
   int sigpurCutBest;
   int backpurCutBest;
   int sigbgdCutBest;
   int comcut = -1;
   double selectcut = -1;

   // Open output file to read saved values
   TFile *ifile = new TFile(filename.c_str(), "READ");
   TTree *getTree;
   TList *tempkeyslist = (TList*)ifile->GetListOfKeys();

   for(int k = 0; k < ifile->GetNkeys()/2; k++)
   {
      negerror = 0;
      poserror = 0;

      // Select the trees (taking only the newest save with MVA variable values)
      stemp[1] = string((tempkeyslist->At(2*k))->GetName()) + ";2";
      cout << endl << stemp[1] << " has been selected for evaluation." << endl;
      getTree = (TTree*)ifile->Get(stemp[1].c_str());
      cout << "Number of events in the tree = " << getTree->GetEntries() << endl;

      // Save the tree title (perticle type)
      treetitle.push_back(string(getTree->GetTitle()));

      getTree->SetBranchAddress((observables[obscheck]).c_str(), &obsvars[0]); // mean
      stemp[2] = observables[obscheck] + "_neg";
      getTree->SetBranchAddress(stemp[2].c_str(), &obsvars[1]); // neg error
      stemp[2] = observables[obscheck] + "_pos";
      getTree->SetBranchAddress(stemp[2].c_str(), &obsvars[2]); // pos error

      // Loop over all entries
      for(int ievt = 0; ievt < getTree->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTree->GetEntry(ievt);
	 obsvect[k].push_back(obsvars[0]);
	 obsvect_neg[k].push_back(obsvars[1]);
	 obsvect_pos[k].push_back(obsvars[2]);

	 if(treetitle[k] == "Data")
	 {
	    negerror += obsvars[1];
	    poserror += obsvars[2];
	 }

	 // If not a data tree, determine the overall range of all trees
/*	 if(treetitle[k] != "Data")
	 {*/
            if(minrange[k] > obsvect[k][ievt])
               minrange[k] = obsvect[k][ievt];
   
            if(maxrange[k] < obsvect[k][ievt])
               maxrange[k] = obsvect[k][ievt];
//	 }

	 // Calculate the mean value of a tree
         mean[k] += obsvect[k][ievt];

//	 cout << "  " << ievt << ":\t" << obsvect[k][ievt] << "\t(" << obsvect_neg[k][ievt] << ",\t" << obsvect_pos[k][ievt] << ")" << endl;
      }

      if(treetitle[k] != "Data")
      {
         cout << "Nr. points: " << nbin_hist_high << endl;
         cout << "Minimum: " << minrange[k] << endl;
         cout << "Maximum: " << maxrange[k] << endl;
      }
      else
      {
         negerror = negerror/obsvect_neg[k].size();
         poserror = poserror/obsvect_pos[k].size();
         cout << "Negative error: " << negerror << endl;
         cout << "Positive error: " << poserror << endl;
      }

      mean[k] = mean[k]/obsvect[k].size();
      cout << "Mean: " << mean[k] << endl;

      if(completerange[0] > minrange[k])
         completerange[0] = minrange[k];

      if(completerange[1] < maxrange[k])
         completerange[1] = maxrange[k];

      cout << "Complete minimum: " << completerange[0] << endl;
      cout << "Complete maximum: " << completerange[1] << endl;
   }

   // Select a custom range for the histograms
   int customrange = -1;
   cerr << endl << "Select custom x-axis ranges for histograms (use -1 for automatic range = " << completerange[0] << "," << completerange[1] << ")? ";
   cin >> customrange;

   if(customrange != -1)
   {
      cerr << "  - Currently selected minimum range for observable " << observables[obscheck] << " (" << completerange[0] << "). Select custom minimum range: ";
      cin >> completerange[0];
      cerr << "  - Currently selected maximum range for observable " << observables[obscheck] << " (" << completerange[1] << "). Select custom maximum range: ";
      cin >> completerange[1];
   }

   completerange[2] = (completerange[1]-completerange[0])*0.05;
   completerange[0] = completerange[0] - completerange[2];
   completerange[1] = completerange[1] + completerange[2];

   cout << endl;

   // Prepare colors for signal, background and observable cut line
   int c_SignalLine     = TColor::GetColor("#0000ee");
   int c_SignalFill     = TColor::GetColor("#7d99d1");
   int c_AllLine        = TColor::GetColor("#ff0000");
   int c_AllFill        = TColor::GetColor("#ff0000");
   int c_MvaCut         = TColor::GetColor("#ffff66");

   // Prepare canvas and histograms for efficiency/purity
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

   TH1F *signalhist;
   TH1F *backgroundhist;
   TH1F *purS;

   // Prepare canvas and histograms for applying the cut
   TCanvas *c2 = new TCanvas("c2","",1200,900);
   c2->SetGrid();
   c2->SetRightMargin(0.05);
   c2->SetTopMargin(0.05);

   int applycount = 0;
   TH1F *basesig[50];
   TH1F *baseback[50];
   TLine *line[3];
   int *sigcount, *backcount;		// counting for signal and background events
   double *max;
   sigcount = new int[50];
   backcount = new int[50];
   max = new double[50];

   // Go through pairs of simulation trees and determine the best cut
   for(int k = 0; k < ifile->GetNkeys()/2; k++) // SIGNAL
   {
      for(int j = 0; j < ifile->GetNkeys()/2; j++) // BACKGROUND
      {
         // Only run, if not a data tree, not equal trees and only once per pair
         if( (treetitle[k] != "Data") && (treetitle[j] != "Data") && (treetitle[j] != treetitle[k]) && (treetitle[j] < treetitle[k]) )
	 {
            cout << "Performing combination: " << treetitle[k] << " and " << treetitle[j] << endl;
            cerr << "Performing combination: " << treetitle[k] << " and " << treetitle[j] << endl;

	    // Zero variables
            found[0] = false;
            found[1] = false;
            found[2] = false;
	    
	    // Determine the range of both trees
	    range[0] = (minrange[k] < minrange[j]) ? minrange[k] : minrange[j];
	    range[1] = (maxrange[k] > maxrange[j]) ? maxrange[k] : maxrange[j];

	    // Create the histograms that will hold efficiencies
	    stemp[0] = treetitle[k] + ", " + observables[obscheck];
            signalhist = new TH1F("signalhist", stemp[0].c_str(), nbin_hist_high, range[0], range[1]);
	    stemp[0] = treetitle[j] + ", " + observables[obscheck];
	    backgroundhist = new TH1F("backgroundhist", stemp[0].c_str(), nbin_hist_high, range[0], range[1]);

            cout << "Signal mean value = " << mean[k] << endl;
            cout << "Background mean value = " << mean[j] << endl;

	    // Determine the sign (depends on the mean values)
	    sign = (mean[0] > mean[1]) ? 1 : -1;

	    // Go through all signal events and create the signal efficiency histogram
            for(int ievt = 0; ievt < obsvect[k].size(); ievt++)
            {
               itemp[2] = ((obsvect[k][ievt] - range[0])/(range[1] - range[0]))*nbin_hist_high + 1;
               if(sign > 0 && itemp[2] > nbin_hist_high) continue;
               if(sign < 0 && itemp[2] < 1             ) continue;
               if(sign > 0 && itemp[2] < 1             ) itemp[2] = 1;
               if(sign < 0 && itemp[2] > nbin_hist_high) itemp[2] = nbin_hist_high;

	       if(sign > 0)
	          for(int ibin = 1; ibin <= itemp[2]; ibin++)
	             signalhist->AddBinContent(ibin, 1);
	       else if(sign < 0)
	          for(int ibin = itemp[2]+1; ibin <= nbin_hist_high; ibin++)
	             signalhist->AddBinContent(ibin, 1);
               else
	       {
		  cout << "Error! Can't determine the sign (signal)." << endl;
                  return 1;
	       }
            }

	    // Go through all background events and create the background efficiency histogram
            for(int ievt = 0; ievt < obsvect[j].size(); ievt++)
            {
               itemp[2] = ((obsvect[j][ievt] - range[0])/(range[1] - range[0]))*nbin_hist_high + 1;
               if(sign > 0 && itemp[2] > nbin_hist_high) continue;
               if(sign < 0 && itemp[2] < 1             ) continue;
               if(sign > 0 && itemp[2] < 1             ) itemp[2] = 1;
               if(sign < 0 && itemp[2] > nbin_hist_high) itemp[2] = nbin_hist_high;

	       if(sign > 0)
	          for(int ibin = 1; ibin <= itemp[2]; ibin++)
	             backgroundhist->AddBinContent(ibin, 1);
	       else if(sign < 0)
	          for(int ibin = itemp[2]+1; ibin <= nbin_hist_high; ibin++)
	             backgroundhist->AddBinContent(ibin, 1);
               else
	       {
		  cout << "Error! Can't determine the sign (background)." << endl;
                  return 1;
	       }
            }

	    // Scale both histograms to maximum value of 1
            signalhist->Scale(1.0/(signalhist->GetMaximum()));
            backgroundhist->Scale(1.0/(backgroundhist->GetMaximum()));

	    // Printout signal histogram values
	    cout << "Signal histogram:" << endl;
	    for(int i = 1; i <= nbin_hist_high; i++)
	       cout << "  " << i << "\t" << signalhist->GetBinContent(i) << endl;
	    cout << endl;

	    // Printout background histogram values
	    cout << "Background histogram:" << endl;
	    for(int i = 1; i <= nbin_hist_high; i++)
	       cout << "  " << i << "\t" << backgroundhist->GetBinContent(i) << endl;
	    cout << endl;

	    // Create the signal purity histogram
	    purS = new TH1F("purS", "Purity of signal", nbin_hist_high, range[0], range[1]);

	    cout << "Purity histogram:" << endl;
	    for(int i = 1; i <= nbin_hist_high; i++)
	    {
               ftemp[0] = (signalhist->GetBinContent(i))*obsvect[k].size();
               ftemp[1] = (backgroundhist->GetBinContent(i))*obsvect[j].size();
	       ftemp[2] = (ftemp[0]+ftemp[1] == 0) ? 0 : ftemp[0]/(ftemp[0]+ftemp[1]);
	       cout << "  " << i << "\t" << ftemp[0] << "\t" << ftemp[1] << "\t" << ftemp[2] << endl;
	       purS->SetBinContent(i, ftemp[2]);
	    }

	    c1->cd();
	    signalhist->SetLineWidth(2);
	    signalhist->SetLineColor(4);
	    signalhist->SetLineStyle(1);
	    signalhist->GetXaxis()->SetTitle(tempdesc[obscheck].c_str());
	    signalhist->GetYaxis()->SetTitle("Efficiency/Purity");
            signalhist->GetXaxis()->SetTitleOffset(1.2);
            signalhist->GetXaxis()->CenterTitle(kTRUE);
            signalhist->GetXaxis()->SetLabelSize(0.028);
            signalhist->GetXaxis()->SetLabelOffset(0.015);
            signalhist->GetYaxis()->SetTitleOffset(1.3);
            signalhist->GetYaxis()->CenterTitle(kTRUE);
            signalhist->GetYaxis()->SetLabelSize(0.028);
            signalhist->GetYaxis()->SetLabelOffset(0.015);
            signalhist->GetXaxis()->SetRange(range[0], range[1]);
            signalhist->GetXaxis()->SetRangeUser(range[0], range[1]);
            signalhist->SetTitle("");
	    signalhist->Draw();

	    backgroundhist->SetLineWidth(2);
	    backgroundhist->SetLineColor(2);
	    backgroundhist->SetLineStyle(1);
	    backgroundhist->Draw("same");

	    purS->SetLineWidth(2);
	    purS->SetLineColor(4);
	    purS->SetLineStyle(7);
	    purS->Draw("same");

            // Draw a legend with info on signal and background efficiency
	    if(sign > 0)
               legend = new TLegend(gPad->GetLeftMargin(), gPad->GetBottomMargin(), gPad->GetLeftMargin()+.2, gPad->GetBottomMargin()+0.12);
	    else if(sign < 0)
               legend = new TLegend(1.-gPad->GetRightMargin()-.2, gPad->GetBottomMargin(), 1.-gPad->GetRightMargin(), gPad->GetBottomMargin()+0.12);
	    else
            {
               cout << "Error! Can't determine the sign (drawing a legend)." << endl;
               return 1;
            }
            legend->SetFillStyle(legendFill);
            legend->AddEntry(signalhist, "Signal efficiency", "l");
            legend->AddEntry(backgroundhist, "Background efficiency", "l");
            legend->AddEntry(purS, "Signal purity", "l");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");

            // Determine the best signal to background cut
	    if(sign > 0)
	    {
	       for(int i = 1; i <= nbin_hist_high; i++)
	       {
                  if( (purS->GetBinContent(i) >= signalhist->GetBinContent(i)) && (!found[0]) && (purS->GetBinContent(i) > 0) )
                  {
                     sigpurCutBest = i;
                     found[0] = true;
                  }

                  if( (purS->GetBinContent(i) >= backgroundhist->GetBinContent(i)) && (!found[1]) && (purS->GetBinContent(i) > 0) )
                  {
                     backpurCutBest = i;
                     found[1] = true;
                  }

                  if( (1.-backgroundhist->GetBinContent(i) >= signalhist->GetBinContent(i)) && (!found[2]) && (purS->GetBinContent(i) > 0) )
                  {
                     sigbgdCutBest = i;
                     found[2] = true;
                  }

                  if(found[0] && found[1] && found[2]) break;
	       }

               if(purS->GetBinContent(sigpurCutBest) > signalhist->GetBinContent(sigpurCutBest))
                  sigpurCutBest = sigpurCutBest - 1;

               if(purS->GetBinContent(backpurCutBest) > backgroundhist->GetBinContent(backpurCutBest))
                  backpurCutBest = backpurCutBest - 1;
          
               if(1.-backgroundhist->GetBinContent(sigbgdCutBest) > signalhist->GetBinContent(sigbgdCutBest))
                  sigbgdCutBest = sigbgdCutBest - 1;
	    }
	    else if(sign < 0)
	    {
	       for(int i = nbin_hist_high; i >= 1; i--)
	       {
                  if( (purS->GetBinContent(i) >= signalhist->GetBinContent(i)) && (!found[0]) )
                  {
                     sigpurCutBest = i;
                     found[0] = true;
                  }

                  if( (purS->GetBinContent(i) >= backgroundhist->GetBinContent(i)) && (!found[1]) )
                  {
                     backpurCutBest = i;
                     found[1] = true;
                  }

                  if( (1.-backgroundhist->GetBinContent(i) >= signalhist->GetBinContent(i)) && (!found[2]) )
                  {
                     sigbgdCutBest = i;
                     found[2] = true;
                  }

                  if(found[0] && found[1] && found[2]) break;
	       }

               if(purS->GetBinContent(sigpurCutBest) > signalhist->GetBinContent(sigpurCutBest))
                  sigpurCutBest = sigpurCutBest - 1;

               if(purS->GetBinContent(backpurCutBest) > backgroundhist->GetBinContent(backpurCutBest))
                  backpurCutBest = backpurCutBest - 1;
          
               if(1.-backgroundhist->GetBinContent(sigbgdCutBest) > signalhist->GetBinContent(sigbgdCutBest))
                  sigbgdCutBest = sigbgdCutBest - 1;
	    }
	    else
            {
               cout << "Error! Can't determine the sign (determining best cut)." << endl;
               return 1;
            }

	    // Inform about the best cut
            cout << "Best sig cut (" << sigpurCutBest << ", " << signalhist->GetXaxis()->GetBinCenter(sigpurCutBest) << "), sig = " << (signalhist->GetBinContent(sigpurCutBest))*100. << ", back = " << (1.-backgroundhist->GetBinContent(sigpurCutBest))*100. << ", purity = " << (purS->GetBinContent(sigpurCutBest))*100. << endl;
            cout << "Best back cut (" << backpurCutBest << ", " << backgroundhist->GetXaxis()->GetBinCenter(backpurCutBest) << "), sig = " << (signalhist->GetBinContent(backpurCutBest))*100. << ", back = " << (1.-backgroundhist->GetBinContent(backpurCutBest))*100. << ", purity = " << (purS->GetBinContent(backpurCutBest))*100. << endl;
            cout << "Equal cut (" << sigbgdCutBest << ", " << signalhist->GetXaxis()->GetBinCenter(sigbgdCutBest) << "), sig = " << (signalhist->GetBinContent(sigbgdCutBest))*100. << ", back = " << (1.-backgroundhist->GetBinContent(sigbgdCutBest))*100. << ", purity = " << (purS->GetBinContent(sigbgdCutBest))*100. << endl;
            cerr << "Best sig cut (" << sigpurCutBest << ", " << signalhist->GetXaxis()->GetBinCenter(sigpurCutBest) << "), sig = " << (signalhist->GetBinContent(sigpurCutBest))*100. << ", back = " << (1.-backgroundhist->GetBinContent(sigpurCutBest))*100. << ", purity = " << (purS->GetBinContent(sigpurCutBest))*100. << endl;
            cerr << "Best back cut (" << backpurCutBest << ", " << backgroundhist->GetXaxis()->GetBinCenter(backpurCutBest) << "), sig = " << (signalhist->GetBinContent(backpurCutBest))*100. << ", back = " << (1.-backgroundhist->GetBinContent(backpurCutBest))*100. << ", purity = " << (purS->GetBinContent(backpurCutBest))*100. << endl;
            cerr << "Equal cut (" << sigbgdCutBest << ", " << signalhist->GetXaxis()->GetBinCenter(sigbgdCutBest) << "), sig = " << (signalhist->GetBinContent(sigbgdCutBest))*100. << ", back = " << (1.-backgroundhist->GetBinContent(sigbgdCutBest))*100. << ", purity = " << (purS->GetBinContent(sigbgdCutBest))*100. << endl;
	    cerr << endl;

            // Select cut to apply to data
	    if(comcut == -1)
            {
               cerr << "Select Best sig cut (0), Best back cut (1), Equal cut (2) or a custom cut (3)? ";
               cin >> comcut;
	    }

	    if(comcut == 0)
	    {
               selectcut = signalhist->GetXaxis()->GetBinCenter(sigpurCutBest);
	       cerr << "Cut " << selectcut << " selected (combination " << treetitle[k] << " vs. " << treetitle[j] << ")" << endl;
	       cerr << "Cut errors are " << selectcut-negerror << " and " << selectcut+poserror << endl;
	    }
            else if(comcut == 1)
	    {
               selectcut = backgroundhist->GetXaxis()->GetBinCenter(backpurCutBest);
	       cerr << "Cut " << selectcut << " selected (combination " << treetitle[k] << " vs. " << treetitle[j] << ")" << endl;
	       cerr << "Cut errors are " << selectcut-negerror << " and " << selectcut+poserror << endl;
	    }
            else if(comcut == 2)
	    {
               selectcut = signalhist->GetXaxis()->GetBinCenter(sigbgdCutBest);
	       cerr << "Cut " << selectcut << " selected (combination " << treetitle[k] << " vs. " << treetitle[j] << ")" << endl;
	       cerr << "Cut errors are " << selectcut-negerror << " and " << selectcut+poserror << endl;
	    }
	    else
	    {
               cerr << "Select the cut to apply to simulations and data (combination " << treetitle[k] << " vs. " << treetitle[j] << "): ";
               cin >> selectcut;
	       cerr << "Cut " << selectcut << " selected (combination " << treetitle[k] << " vs. " << treetitle[j] << ")" << endl;
	       cerr << "Cut errors are " << selectcut-negerror << " and " << selectcut+poserror << endl;
	    }

            signalhist->GetXaxis()->SetRange(completerange[0], completerange[1]);
            signalhist->GetXaxis()->SetRangeUser(completerange[0], completerange[1]);
	    c1->Update();

    	    // Add observable cut line histograms
            line[0] = new TLine(selectcut, 0., selectcut, signalhist->GetMaximum());
            line[0]->SetLineWidth(2);
            line[0]->SetLineStyle(1);
            line[0]->SetLineColor(kOrange+2);
            line[0]->Draw("same");

            line[1] = new TLine(selectcut-negerror, 0., selectcut-negerror, signalhist->GetMaximum());
            line[1]->SetLineWidth(2);
            line[1]->SetLineStyle(7);
            line[1]->SetLineColor(kOrange+2);
            line[1]->Draw("same");

            line[2] = new TLine(selectcut+poserror, 0., selectcut+poserror, signalhist->GetMaximum());
            line[2]->SetLineWidth(2);
            line[2]->SetLineStyle(7);
            line[2]->SetLineColor(kOrange+2);
            line[2]->Draw("same");

	    stemp[0] = RemoveFilename(&filename) + "/individual_plots/purS_" + observables[obscheck] + "_" + treetitle[k] + "-" + treetitle[j] + ".pdf";
	    c1->SaveAs(stemp[0].c_str());

	    delete signalhist;
	    delete backgroundhist;
	    delete purS;

            for(int i = 0; i < ifile->GetNkeys()/2; i++) // ALL TREES
	    {
	       // Prepare signal and background histograms
	       stemp[0] = "basesig" + ToString(applycount);
	       stemp[2] = treetitle[i] + ", " + observables[obscheck];
               basesig[applycount] = new TH1F(stemp[0].c_str(), stemp[2].c_str(), 100, completerange[0], completerange[1]);
               basesig[applycount]->SetBit(TH1::kCanRebin);
               stemp[0] = "baseback" + ToString(applycount);
	       stemp[2] = treetitle[i] + ", " + observables[obscheck];
               baseback[applycount] = new TH1F(stemp[0].c_str(), stemp[2].c_str(), 100, completerange[0], completerange[1]);
               baseback[applycount]->SetBit(TH1::kCanRebin);

               sigcount[applycount] = 0;
               backcount[applycount] = 0;
	       max[applycount] = 0;

	       // Apply cut to current tree
	       for(int ievt = 0; ievt < obsvect[i].size(); ievt++)
	       {
		  // Value is above the histogram range
		  if(obsvect[i][ievt] > completerange[1])
		  {
//	             cerr << "Overflow bin (" << obsvect[i][ievt] << ")!" << endl;
                     baseback[applycount]->SetBinContent(0, baseback[applycount]->GetBinContent(0)+1);

		     if(sign > 0)
		        sigcount[applycount]++;
		     else if(sign < 0)
		        backcount[applycount]++;
		  }
		  // Value is below the histogram range
		  else if(obsvect[i][ievt] < completerange[0])
		  {
//	             cerr << "Underflow bin (" << obsvect[i][ievt] << ")!" << endl;
                     basesig[applycount]->SetBinContent(0, basesig[applycount]->GetBinContent(0)+1);

		     if(sign > 0)
		        backcount[applycount]++;
		     else if(sign < 0)
		        sigcount[applycount]++;
		  }
		  // Value goes into the histogram
		  else
		  {
		     if(sign > 0)
		     {
		        // Event above MVA cut (signal)
                        if(obsvect[i][ievt] >= selectcut)
		        {
                           basesig[applycount]->Fill(obsvect[i][ievt]);
		           sigcount[applycount]++;
		        }
		        // Event above MVA cut (background)
		        else
		        {
                           baseback[applycount]->Fill(obsvect[i][ievt]);
		           backcount[applycount]++;
		        }
		     }
		     else if(sign < 0)
		     {
		        // Event above MVA cut (background)
                        if(obsvect[i][ievt] >= selectcut)
		        {
                           baseback[applycount]->Fill(obsvect[i][ievt]);
		           backcount[applycount]++;
		        }
		        // Event above MVA cut (signal)
		        else
		        {
                           basesig[applycount]->Fill(obsvect[i][ievt]);
		           sigcount[applycount]++;
		        }
		     }
	             else
                     {
                        cout << "Error! Can't determine the sign (applying to simulations and data)." << endl;
                        return 1;
                     }
		  }
	       }

               // Write the signal vs. background count information
	       dtemp[0] = 100.*(double)backcount[applycount]/((double)sigcount[applycount] + (double)backcount[applycount]);
	       dtemp[1] = 100.*(double)sigcount[applycount]/((double)sigcount[applycount] + (double)backcount[applycount]);
               cout << "Signal vs. background (" << treetitle[i] << "):" << endl;
               cout << " - " << observables[obscheck] << " = " << sigcount[applycount] << " vs. " << backcount[applycount] << " (" << ToString(dtemp[1],3) << "% vs. " << ToString(dtemp[0],3) << "%)" << endl;
	       cout << " - Number of underflow events = " << basesig[applycount]->GetBinContent(0) << endl;
	       cout << " - Number of overflow events = " << baseback[applycount]->GetBinContent(0) << endl;
               cout << endl;
               cerr << "Signal vs. background (" << treetitle[i] << "):" << endl;
               cerr << " - " << observables[obscheck] << " = " << sigcount[applycount] << " vs. " << backcount[applycount] << " (" << ToString(dtemp[1],3) << "% vs. " << ToString(dtemp[0],3) << "%)" << endl;
	       cerr << " - Number of underflow events = " << basesig[applycount]->GetBinContent(0) << endl;
	       cerr << " - Number of overflow events = " << baseback[applycount]->GetBinContent(0) << endl;
               cerr << endl;

               // Check maximum value in a histogram
               if((basesig[applycount]->GetMaximum()) > max[applycount])
                  max[applycount] = basesig[applycount]->GetMaximum();
               if((baseback[applycount]->GetMaximum()) > max[applycount])
                  max[applycount] = baseback[applycount]->GetMaximum();

	       // Set line and fill attributes for histograms
               basesig[applycount]->SetLineColor(c_SignalLine);
               basesig[applycount]->SetLineWidth(2);
               basesig[applycount]->SetFillColor(c_SignalFill);
               basesig[applycount]->SetFillStyle(1001);
               baseback[applycount]->SetLineColor(c_AllLine);
               baseback[applycount]->SetLineWidth(2);
               baseback[applycount]->SetFillColor(c_AllFill);
               baseback[applycount]->SetFillStyle(3554);

	       // Draw signal and background histograms
	       c2->cd();
               basesig[applycount]->Draw();
               baseback[applycount]->Draw("same");

               // Set axis ranges
               basesig[applycount]->GetYaxis()->SetRangeUser(0.,max[applycount]*1.2);
               basesig[applycount]->SetMaximum(max[applycount]*1.2);
               basesig[applycount]->GetXaxis()->SetRange(completerange[0], completerange[1]);
               basesig[applycount]->GetXaxis()->SetRangeUser(completerange[0], completerange[1]);

	       // Setup axis (title, labels)
               basesig[applycount]->GetYaxis()->SetTitle("Number of events");
               basesig[applycount]->GetXaxis()->SetTitle(tempdesc[obscheck].c_str());
               basesig[applycount]->GetXaxis()->SetTitleOffset(1.2);
               basesig[applycount]->GetXaxis()->CenterTitle(kTRUE);
               basesig[applycount]->GetXaxis()->SetLabelSize(0.028);
               basesig[applycount]->GetXaxis()->SetLabelOffset(0.015);
               basesig[applycount]->GetYaxis()->SetTitleOffset(1.3);
               basesig[applycount]->GetYaxis()->CenterTitle(kTRUE);
               basesig[applycount]->GetYaxis()->SetLabelSize(0.028);
               basesig[applycount]->GetYaxis()->SetLabelOffset(0.015);
               basesig[applycount]->SetTitle("");

    	       // Add observable cut line histograms
               line[0] = new TLine(selectcut, 0., selectcut, basesig[applycount]->GetMaximum());
               line[0]->SetLineWidth(2);
               line[0]->SetLineStyle(1);
               line[0]->SetLineColor(kOrange+2);
               line[0]->Draw("same");

               line[1] = new TLine(selectcut-negerror, 0., selectcut-negerror, basesig[applycount]->GetMaximum());
               line[1]->SetLineWidth(2);
               line[1]->SetLineStyle(7);
               line[1]->SetLineColor(kOrange+2);
               line[1]->Draw("same");

               line[2] = new TLine(selectcut+poserror, 0., selectcut+poserror, basesig[applycount]->GetMaximum());
               line[2]->SetLineWidth(2);
               line[2]->SetLineStyle(7);
               line[2]->SetLineColor(kOrange+2);
               line[2]->Draw("same");

               // Draw a legend with info on signal and background events
               legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.12, gPad->GetLeftMargin()+.45, 1-gPad->GetTopMargin());
               legend->SetFillStyle(legendFill);
               legend->SetFillColor(c_MvaCut);
	       dtemp[0] = 100.*(double)backcount[applycount]/((double)sigcount[applycount] + (double)backcount[applycount]);
               stemp[0] = "Cut background events (" + ToString(backcount[applycount]) + " = " + ToString(dtemp[0], 2) + "%)";
               legend->AddEntry(baseback[applycount],stemp[0].c_str(),"f");
	       dtemp[0] = 100.*(double)sigcount[applycount]/((double)sigcount[applycount] + (double)backcount[applycount]);
               stemp[0] = "Cut signal events (" + ToString(sigcount[applycount]) + " = " + ToString(dtemp[0], 2) + "%)";
               legend->AddEntry(basesig[applycount],stemp[0].c_str(),"f");
               legend->SetBorderSize(1);
               legend->SetMargin(0.3);
               legend->Draw("same");

	       // Write down the underflow and overflow information
	       uoflowtext.SetTextSize(0.017);
	       uoflowtext.SetTextAlign(11);
               stemp[0] = "UFLOW = " + ToString(basesig[applycount]->GetBinContent(0));
	       uoflowtext.DrawLatex(completerange[0], max[applycount]*1.205, stemp[0].c_str());
	       uoflowtext.SetTextAlign(31);
               stemp[0] = "OFLOW = " + ToString(baseback[applycount]->GetBinContent(0));
	       uoflowtext.DrawLatex(completerange[1], max[applycount]*1.205, stemp[0].c_str());

	       // Save plot as PDF
               stemp[0] = RemoveFilename(&filename) + "/individual_plots/apply_" + observables[obscheck] + "_cut-" + treetitle[k] + "-" + treetitle[j] + "_" + treetitle[i] + ".pdf";
               c2->SaveAs(stemp[0].c_str());

               applycount++;
	    }
	 }
      }
   }
   //-----------------------------------------
   
   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;
   delete[] ftemp;
   delete[] ctemp;
   delete[] obsvars;

   cerr << "Analysis program finished correctly." << endl;
   
   return 0;
}
