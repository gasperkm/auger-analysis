#define _STANDALONE_ 1
#include "workstation.h"
#include "separate_functions.h"
#include "root_style.h"

using namespace std;

int main(int argc, char **argv)
{
   string filename;
   if(argc > 1)
      filename = string(argv[1]);
   else
   {
      cerr << "Error! No input file supplied. Rerun program and add input file as argument (mvatree_file.root)." << endl;
      return 1;
   }
   
//   gSystem->Load("libTree.so");

   //-----------------------------------------
   ifstream ifs;
   string *stemp = new string[5];
   int *itemp = new int[3];
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
   // Give option to plot all separately or plot two trees together TODO
   int doublePlot = 0;
   cout << endl << "Would you like to plot trees separately (0) or two trees together (1) or three trees together (2)? ";
   cin >> itemp[0];
   doublePlot = itemp[0];

   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/histograms";
   system(stemp[0].c_str());
   if(doublePlot == 2)
   {
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/histograms/triphist_*";
      system(stemp[0].c_str());
   }
   else if(doublePlot == 1)
   {
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/histograms/dualhist_*";
      system(stemp[0].c_str());
   }
   else if(doublePlot == 0)
   {
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/histograms/hist_*";
      system(stemp[0].c_str());
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/histograms/risetime-comparison_*";
      system(stemp[0].c_str());
   }
   else
   {
      cout << "Error with plotting selection!" << endl;
      return 0;
   }

   // Prepare colors for MVA cut line
   int c_MvaCut         = TColor::GetColor("#ffff66");

   // Prepare plotting objects (legend, line, histograms, graphs,...) - for a maximum of 50 observables
   TLegend *legend;
   int nrfigs;  // number of plots on same graph at the same time
   TLine *line[3];
   TH1F *treeHist[50], *treeHistSig[50], *treeHistBack[50], *treeHistData[50];
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
   int *uflowcount, *oflowcount, *totcount, *totcountSig, *totcountBack, *totcountData; // counting number of underflow and overflow events
   uflowcount = new int[50];
   oflowcount = new int[50];
   totcount = new int[50];
   totcountSig = new int[50];
   totcountBack = new int[50];
   totcountData = new int[50];

   bool isangle[2] = {false, false};	// check if observable is an angle (needs to be converted to degrees from radians) 
   bool isgood[2] = {true, true};	// check if observable has an invalid value -1
   bool iszenith[2] = {false, false};	// check if observable is a zenith angle

   float *obsvars;			// array for variable values that we read from output file
   obsvars = new float[300];
   int backShift = 100;	// background shift for saving values
   int dataShift = 200;	// data shift for saving values

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
   RootStyle *mystyle = new RootStyle();
   mystyle->SetBaseStyle();

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TCanvas *c2 = new TCanvas("c2","",1200,900);

   // Open output file to read saved values
   TFile *ifile = new TFile(filename.c_str(), "READ");
   TTree *getTree, *getTreeDual[2], *getTreeTrip[3];
   TList *tempkeyslist = (TList*)ifile->GetListOfKeys();

   // Plot three histograms together
   if(doublePlot == 2)
   {
      string sigTree = "empty";
      string backTree = "empty";
      string dataTree = "empty";
      cout << endl << "Select signal, background and data trees to plot on the same graph. Available trees:" << endl;
      for(int k = 0; k < ifile->GetNkeys(); k++)
      {
         stemp[1] = string((tempkeyslist->At(k))->GetName());
         cout << "  " << k+1 << " = " << stemp[1] << endl;
      }
      cout << "Choose option for signal tree: ";
      cin >> itemp[0];
      if(itemp[0] != -1)
         sigTree = string((tempkeyslist->At(itemp[0]-1))->GetName());

      cout << "Choose option for background tree: ";
      cin >> itemp[0];
      if(itemp[0] != -1)
         backTree = string((tempkeyslist->At(itemp[0]-1))->GetName());

      cout << "Choose option for data tree: ";
      cin >> itemp[0];
      if(itemp[0] != -1)
         dataTree = string((tempkeyslist->At(itemp[0]-1))->GetName());

      cout << "Signal, background and data trees are: " << sigTree << ", " << backTree << ", " << dataTree << endl << endl;

      // Zero minimum and maximum holders for observables
      for(int i = 0; i < itemp[1]; i++)
      {
         minimumval[i] = 1.e+30;
         maximumval[i] = -1.e+30;
         uflowcount[i] = 0;
         oflowcount[i] = 0;
         totcountSig[i] = 0;
         totcountBack[i] = 0;
         totcountData[i] = 0;

/*         if(observables[i] == "risetime")
            riseobs[0] = i;

         if(observables[i] == "risetimerecalc")
            riseobs[1] = i;*/
      }

      // Run over signal tree observables
      cout << endl << sigTree << " has been selected for evaluation." << endl;
      getTreeTrip[0] = (TTree*)ifile->Get(sigTree.c_str());
      cout << "Number of events in the signal tree = " << getTreeTrip[0]->GetEntries() << endl;
      getTreeTrip[1] = (TTree*)ifile->Get(backTree.c_str());
      cout << "Number of events in the background tree = " << getTreeTrip[1]->GetEntries() << endl;
      getTreeTrip[2] = (TTree*)ifile->Get(dataTree.c_str());
      cout << "Number of events in the background tree = " << getTreeTrip[2]->GetEntries() << endl;

      // Set addresses for reading observables
      for(int i = 0; i < itemp[1]; i++)
      {
         cout << "Signal: Setting observable " << observables[i] << endl;
         getTreeTrip[0]->SetBranchAddress((observables[i]).c_str(), &obsvars[3*i]); // mean
         cout << "Mean variable: obsvars[3*" << i << "] = obsvars[" << 3*i << "]" << endl;
         stemp[2] = observables[i] + "_neg";
         getTreeTrip[0]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+1]); // neg error
         cout << "Negerror variable: obsvars[3*" << i << "+1] = obsvars[" << 3*i+1 << "]" << endl;
         stemp[2] = observables[i] + "_pos";
         getTreeTrip[0]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+2]); // pos error
         cout << "Poserror variable: obsvars[3*" << i << "+2] = obsvars[" << 3*i+2 << "]" << endl;

         cout << "Background: Setting observable " << observables[i] << endl;
         getTreeTrip[1]->SetBranchAddress((observables[i]).c_str(), &obsvars[3*i+backShift]); // mean
         cout << "Mean variable: obsvars[3*" << i << "+backShift] = obsvars[" << 3*i+backShift << "]" << endl;
         stemp[2] = observables[i] + "_neg";
         getTreeTrip[1]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+1+backShift]); // neg error
         cout << "Negerror variable: obsvars[3*" << i << "+1+backShift] = obsvars[" << 3*i+1+backShift << "]" << endl;
         stemp[2] = observables[i] + "_pos";
         getTreeTrip[1]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+2+backShift]); // pos error
         cout << "Poserror variable: obsvars[3*" << i << "+2+backShift] = obsvars[" << 3*i+2+backShift << "]" << endl;

         cout << "Data: Setting observable " << observables[i] << endl;
         getTreeTrip[2]->SetBranchAddress((observables[i]).c_str(), &obsvars[3*i+dataShift]); // mean
         cout << "Mean variable: obsvars[3*" << i << "+dataShift] = obsvars[" << 3*i+dataShift << "]" << endl;
         stemp[2] = observables[i] + "_neg";
         getTreeTrip[2]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+1+dataShift]); // neg error
         cout << "Negerror variable: obsvars[3*" << i << "+1+dataShift] = obsvars[" << 3*i+1+dataShift << "]" << endl;
         stemp[2] = observables[i] + "_pos";
         getTreeTrip[2]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+2+dataShift]); // pos error
         cout << "Poserror variable: obsvars[3*" << i << "+2+dataShift] = obsvars[" << 3*i+2+dataShift << "]" << endl;
      }

      // Prepare plots to save values into (histograms need to able to rebin)
      maxobs[2] = (maxobs[1]-maxobs[0])*0.05;
      maxobs[0] = maxobs[0] - maxobs[2];
      maxobs[1] = maxobs[1] + maxobs[2];

      for(int i = 0; i < itemp[1]; i++)
      {
         xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
         xhistlimit[0] = tempmin[i] - xhistlimit[2];
         xhistlimit[1] = tempmax[i] + xhistlimit[2];

         // Prepare histograms for all observables - Signal
         stemp[2] = "treeHistSig" + ToString(i);
         treeHistSig[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
#if ROOTVER == 5
         treeHistSig[i]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
         treeHistSig[i]->SetCanExtend(TH1::kXaxis);
#endif

         // Prepare histograms for all observables - Background
         stemp[2] = "treeHistBack" + ToString(i);
         treeHistBack[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
#if ROOTVER == 5
         treeHistBack[i]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
         treeHistBack[i]->SetCanExtend(TH1::kXaxis);
#endif

         // Prepare histograms for all observables - Data
         stemp[2] = "treeHistData" + ToString(i);
         treeHistData[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
#if ROOTVER == 5
         treeHistData[i]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
         treeHistData[i]->SetCanExtend(TH1::kXaxis);
#endif

         // Set maximum value to zero, and underflow and overflow counters to zero
         max[i] = 0.0;
         uflowcount[i] = 0;
         oflowcount[i] = 0;
      }

      // Loop over all entries (signal)
      for(int ievt = 0; ievt < getTreeTrip[0]->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTreeTrip[0]->GetEntry(ievt);

         // Loop over all observables
         for(int i = 0; i < itemp[1]; i++)
         {
            // Check if observable is an angle in radians
            if( (observables[i] == "zenithSD") || (observables[i] == "azimuthSD") || (observables[i] == "zenithFD") || (observables[i] == "azimuthFD") || (observables[i] == "latitudeSD") || (observables[i] == "longitudeSD") || (observables[i] == "latitudeFD") || (observables[i] == "longitudeFD") )
	    {
               if( (observables[i] == "zenithSD") || (observables[i] == "zenithFD") )
                  iszenith[0] = true;
	       else
                  iszenith[0] = false;
               isangle[0] = true;
	    }
            else
               isangle[0] = false;

            // Check if observable is invalid and has value -1
            if(obsvars[3*i] == -1)
               isgood[0] = false;
            else
               isgood[0] = true;

            // Check the mean observable value in order to get minimum and maximum values of all trees
	    if(isangle[0])
	    {
               if(iszenith[0])
	          ftemp[0] = obsvars[3*i];
//	          ftemp[0] = SecTheta(obsvars[3*i]);
//	          ftemp[0] = RadToDeg(obsvars[3*i]);
	       else
	          ftemp[0] = RadToDeg(obsvars[3*i]);
	    }
	    else
	       ftemp[0] = obsvars[3*i];

            if(isgood[0])
            {
               if(ftemp[0] < minimumval[i])
                  minimumval[i] = ftemp[0];
               if(ftemp[0] > maximumval[i])
                  maximumval[i] = ftemp[0];
	    }

            // Save observable values into the histograms
            if(isgood[0])
            {
               xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] = tempmax[i] + xhistlimit[2];

               if(ftemp[0] > xhistlimit[1])
                  oflowcount[i]++;
               else if(ftemp[0] < xhistlimit[0])
                  uflowcount[i]++;
               else
                  treeHistSig[i]->Fill(ftemp[0]);
                  
               totcountSig[i]++;
            }
         } // Loop over all observables
      } // Loop over all entries (signal)

      // Loop over all entries (background)
      for(int ievt = 0; ievt < getTreeTrip[1]->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTreeTrip[1]->GetEntry(ievt);

         // Loop over all observables
         for(int i = 0; i < itemp[1]; i++)
         {
            // Check if observable is an angle in radians
            if( (observables[i] == "zenithSD") || (observables[i] == "azimuthSD") || (observables[i] == "zenithFD") || (observables[i] == "azimuthFD") || (observables[i] == "latitudeSD") || (observables[i] == "longitudeSD") || (observables[i] == "latitudeFD") || (observables[i] == "longitudeFD") )
	    {
               if( (observables[i] == "zenithSD") || (observables[i] == "zenithFD") )
                  iszenith[0] = true;
	       else
                  iszenith[0] = false;
               isangle[0] = true;
	    }
            else
               isangle[0] = false;

            // Check if observable is invalid and has value -1
            if(obsvars[3*i+backShift] == -1)
               isgood[0] = false;
            else
               isgood[0] = true;

            // Check the mean observable value in order to get minimum and maximum values of all trees
	    if(isangle[0])
	    {
               if(iszenith[0])
	          ftemp[0] = obsvars[3*i+backShift];
//	          ftemp[0] = SecTheta(obsvars[3*i+backShift]);
//	          ftemp[0] = RadToDeg(obsvars[3*i+backShift]);
	       else
	          ftemp[0] = RadToDeg(obsvars[3*i+backShift]);
	    }
	    else
	       ftemp[0] = obsvars[3*i+backShift];

            if(isgood[0])
            {
               if(ftemp[0] < minimumval[i])
                  minimumval[i] = ftemp[0];
               if(ftemp[0] > maximumval[i])
                  maximumval[i] = ftemp[0];
            }

            // Save observable values into the histograms
            if(isgood[0])
            {
               xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] = tempmax[i] + xhistlimit[2];

               if(ftemp[0] > xhistlimit[1])
                  oflowcount[i]++;
               else if(ftemp[0] < xhistlimit[0])
                  uflowcount[i]++;
               else
               {
                  if(i == riseobs[0])
                     riseHist->Fill(ftemp[0]);

                  if(i == riseobs[1])
                     riserecalcHist->Fill(ftemp[0]);

                  treeHistBack[i]->Fill(ftemp[0]);
               }
                  
               totcountBack[i]++;
            }
         } // Loop over all observables
      } // Loop over all entries (background)

      // Loop over all entries (data)
      for(int ievt = 0; ievt < getTreeTrip[2]->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTreeTrip[2]->GetEntry(ievt);

         // Loop over all observables
         for(int i = 0; i < itemp[1]; i++)
         {
            // Check if observable is an angle in radians
            if( (observables[i] == "zenithSD") || (observables[i] == "azimuthSD") || (observables[i] == "zenithFD") || (observables[i] == "azimuthFD") || (observables[i] == "latitudeSD") || (observables[i] == "longitudeSD") || (observables[i] == "latitudeFD") || (observables[i] == "longitudeFD") )
	    {
               if( (observables[i] == "zenithSD") || (observables[i] == "zenithFD") )
                  iszenith[0] = true;
	       else
                  iszenith[0] = false;
               isangle[0] = true;
	    }
            else
               isangle[0] = false;

            // Check if observable is invalid and has value -1
            if(obsvars[3*i+dataShift] == -1)
               isgood[0] = false;
            else
               isgood[0] = true;

            // Check the mean observable value in order to get minimum and maximum values of all trees
	    if(isangle[0])
	    {
               if(iszenith[0])
	          ftemp[0] = obsvars[3*i+dataShift];
//	          ftemp[0] = SecTheta(obsvars[3*i+dataShift]);
//	          ftemp[0] = RadToDeg(obsvars[3*i+dataShift]);
	       else
	          ftemp[0] = RadToDeg(obsvars[3*i+dataShift]);
	    }
	    else
	       ftemp[0] = obsvars[3*i+dataShift];

            if(isgood[0])
            {
               if(ftemp[0] < minimumval[i])
                  minimumval[i] = ftemp[0];
               if(ftemp[0] > maximumval[i])
                  maximumval[i] = ftemp[0];
            }

            // Save observable values into the histograms
            if(isgood[0])
            {
               xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] = tempmax[i] + xhistlimit[2];

               if(ftemp[0] > xhistlimit[1])
                  oflowcount[i]++;
               else if(ftemp[0] < xhistlimit[0])
                  uflowcount[i]++;
               else
               {
                  if(i == riseobs[0])
                     riseHist->Fill(ftemp[0]);

                  if(i == riseobs[1])
                     riserecalcHist->Fill(ftemp[0]);

                  treeHistData[i]->Fill(ftemp[0]);
               }
                  
               totcountData[i]++;
            }
         } // Loop over all observables
      } // Loop over all entries (data)

      // Write the valid event information
      cout << "Valid events (" << getTreeTrip[0]->GetTitle() << ", " << getTreeTrip[1]->GetTitle() << " and " << getTreeTrip[2]->GetTitle() << "):" << endl;
      for(int i = 0; i < itemp[1]; i++)
      {
         cout << " - " << observables[i] << " = " << totcountSig[i]+totcountBack[i]+totcountData[i] << " (UFLOW=" << uflowcount[i] << ", OFLOW=" << oflowcount[i] << ")" << endl;
      }
      cout << endl;

      // Find energy observable (FD)
      for(int i = 0; i < itemp[1]; i++)
      {
         if(observables[i] == "energyFD")
            itemp[2] = i;
      }

      // Plot all observables
      for(int i = 0; i < itemp[1]; i++)
      {
         c1->cd();

	 treeHistSig[i]->Scale(1./totcountSig[i]);
	 treeHistBack[i]->Scale(1./totcountBack[i]);
	 treeHistData[i]->Scale(1./totcountData[i]);

         // Check maximum value in a histogram
         if((treeHistSig[i]->GetMaximum()) > max[i])
            max[i] = treeHistSig[i]->GetMaximum();
         if((treeHistBack[i]->GetMaximum()) > max[i])
            max[i] = treeHistBack[i]->GetMaximum();
         if((treeHistData[i]->GetMaximum()) > max[i])
            max[i] = treeHistData[i]->GetMaximum();

         // Set line and fill attributes for histograms
         mystyle->SetHistColor((TH1*)treeHistBack[i], 0);
         mystyle->SetHistColor((TH1*)treeHistSig[i], 1);
         mystyle->SetHistColor((TH1*)treeHistData[i], 2);

         // Draw signal and background histograms
         treeHistSig[i]->Draw();
         treeHistBack[i]->Draw("SAME");
         treeHistData[i]->Draw("SAME");

/*	 cout << "Observable name = " << observables[i] << ", observable min = " << minimumval[i] << ", observable max = " << maximumval[i] << endl;

         // Set axis ranges
	 if( (observables[i] == "energyFD") )
	 {
            xhistlimit[2] = (maximumval[i]-minimumval[i])*0.05;
            xhistlimit[0] = minimumval[i] - xhistlimit[2];
            xhistlimit[1] = maximumval[i] + xhistlimit[2];
            cout << "  limits = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	 }
	 else if(observables[i] == "shwsize")
	 {
            xhistlimit[1] = TMath::Exp(-38.+2.238*TMath::Log10((minimumval[itemp[2]]+maximumval[itemp[2]])/2.));
            xhistlimit[2] = (xhistlimit[1]-tempmin[i])*0.05;
            xhistlimit[0] = tempmin[i] - xhistlimit[2];
            xhistlimit[1] += xhistlimit[2];
            cout << "  limits = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	 }
	 else
	 {*/
            xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
            xhistlimit[0] = tempmin[i] - xhistlimit[2];
            xhistlimit[1] = tempmax[i] + xhistlimit[2];
/*            cout << "  overwritten = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	 }*/
         treeHistSig[i]->GetYaxis()->SetRangeUser(0.,max[i]+0.25*max[i]);
         treeHistSig[i]->SetMaximum(max[i]*1.2);
         treeHistSig[i]->GetXaxis()->SetRange(xhistlimit[0], xhistlimit[1]);
         treeHistSig[i]->GetXaxis()->SetRangeUser(xhistlimit[0], xhistlimit[1]);

         // Setup axis (title, labels)
         mystyle->SetAxisTitles((TH1*)treeHistSig[i], tempdesc[i], "Number of events (normalized)");

         // Draw a legend with info on signal and background events
	 nrfigs = 3;
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
         legend->SetFillStyle(legendFill);
         legend->SetFillColor(c_MvaCut);
         stemp[0] = "Pure signal (" + ToString(totcountSig[i]) + ")";
         legend->AddEntry(treeHistSig[i], stemp[0].c_str(),"f");
         stemp[0] = "Pure background (" + ToString(totcountBack[i]) + ")";
         legend->AddEntry(treeHistBack[i], stemp[0].c_str(),"f");
         stemp[0] = "Data (" + ToString(totcountData[i]) + ")";
         legend->AddEntry(treeHistData[i], stemp[0].c_str(),"f");
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");

         // Write down the underflow and overflow information
         uoflowtext.SetTextSize(0.017);
         uoflowtext.SetTextAlign(11);
         stemp[0] = "UFLOW = " + ToString(uflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[0], max[i]*1.205, stemp[0].c_str());
         uoflowtext.SetTextAlign(31);
         stemp[0] = "OFLOW = " + ToString(oflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[1], max[i]*1.205, stemp[0].c_str());

         // Prepare save name for plots
         stemp[1] = sigTree + "_" + backTree + "_" + dataTree;
         stemp[3] = RemoveFilename(&filename) + "/histograms/triphist_" + stemp[1] + "_" + observables[i] + ".pdf";

         // Save plots as PDF
         c1->SaveAs(stemp[3].c_str());
      } // Plot all observables
   }
   // Plot two histograms together
   else if(doublePlot == 1)
   {
      // Decide which tree will be "signal" and which "background"
      string sigTree = "empty";
      string backTree = "empty";
      cout << endl << "Select signal and background trees to plot on the same graph. Available trees:" << endl;
      for(int k = 0; k < ifile->GetNkeys(); k++)
      {
         stemp[1] = string((tempkeyslist->At(k))->GetName());
         cout << "  " << k+1 << " = " << stemp[1] << endl;
      }
      cout << "Choose option for signal tree: ";
      cin >> itemp[0];
      if(itemp[0] != -1)
         sigTree = string((tempkeyslist->At(itemp[0]-1))->GetName());

      cout << "Choose option for background tree: ";
      cin >> itemp[0];
      if(itemp[0] != -1)
         backTree = string((tempkeyslist->At(itemp[0]-1))->GetName());

      cout << "Signal and background trees are: " << sigTree << ", " << backTree << endl << endl;

      // Zero minimum and maximum holders for observables
      for(int i = 0; i < itemp[1]; i++)
      {
         minimumval[i] = 1.e+30;
         maximumval[i] = -1.e+30;
         uflowcount[i] = 0;
         oflowcount[i] = 0;
         totcountSig[i] = 0;
         totcountBack[i] = 0;

/*         if(observables[i] == "risetime")
            riseobs[0] = i;

         if(observables[i] == "risetimerecalc")
            riseobs[1] = i;*/
      }

      // Run over signal tree observables
      cout << endl << sigTree << " has been selected for evaluation." << endl;
      getTreeDual[0] = (TTree*)ifile->Get(sigTree.c_str());
      cout << "Number of events in the signal tree = " << getTreeDual[0]->GetEntries() << endl;
      getTreeDual[1] = (TTree*)ifile->Get(backTree.c_str());
      cout << "Number of events in the background tree = " << getTreeDual[1]->GetEntries() << endl;

      // Set addresses for reading observables
      for(int i = 0; i < itemp[1]; i++)
      {
         cout << "Signal: Setting observable " << observables[i] << endl;
         getTreeDual[0]->SetBranchAddress((observables[i]).c_str(), &obsvars[3*i]); // mean
         cout << "Mean variable: obsvars[3*" << i << "] = obsvars[" << 3*i << "]" << endl;
         stemp[2] = observables[i] + "_neg";
         getTreeDual[0]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+1]); // neg error
         cout << "Negerror variable: obsvars[3*" << i << "+1] = obsvars[" << 3*i+1 << "]" << endl;
         stemp[2] = observables[i] + "_pos";
         getTreeDual[0]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+2]); // pos error
         cout << "Poserror variable: obsvars[3*" << i << "+2] = obsvars[" << 3*i+2 << "]" << endl;

         cout << "Background: Setting observable " << observables[i] << endl;
         getTreeDual[1]->SetBranchAddress((observables[i]).c_str(), &obsvars[3*i+backShift]); // mean
         cout << "Mean variable: obsvars[3*" << i << "+backShift] = obsvars[" << 3*i+backShift << "]" << endl;
         stemp[2] = observables[i] + "_neg";
         getTreeDual[1]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+1+backShift]); // neg error
         cout << "Negerror variable: obsvars[3*" << i << "+1+backShift] = obsvars[" << 3*i+1+backShift << "]" << endl;
         stemp[2] = observables[i] + "_pos";
         getTreeDual[1]->SetBranchAddress(stemp[2].c_str(), &obsvars[3*i+2+backShift]); // pos error
         cout << "Poserror variable: obsvars[3*" << i << "+2+backShift] = obsvars[" << 3*i+2+backShift << "]" << endl;
      }

      // Prepare plots to save values into (histograms need to able to rebin)
      maxobs[2] = (maxobs[1]-maxobs[0])*0.05;
      maxobs[0] = maxobs[0] - maxobs[2];
      maxobs[1] = maxobs[1] + maxobs[2];

      for(int i = 0; i < itemp[1]; i++)
      {
         xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
         xhistlimit[0] = tempmin[i] - xhistlimit[2];
         xhistlimit[1] = tempmax[i] + xhistlimit[2];

         // Prepare histograms for all observables - Signal
         stemp[2] = "treeHistSig" + ToString(i);
         treeHistSig[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
#if ROOTVER == 5
         treeHistSig[i]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
         treeHistSig[i]->SetCanExtend(TH1::kXaxis);
#endif

         // Prepare histograms for all observables - Background
         stemp[2] = "treeHistBack" + ToString(i);
         treeHistBack[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
#if ROOTVER == 5
         treeHistBack[i]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
         treeHistBack[i]->SetCanExtend(TH1::kXaxis);
#endif

         // Set maximum value to zero, and underflow and overflow counters to zero
         max[i] = 0.0;
         uflowcount[i] = 0;
         oflowcount[i] = 0;
      }

      // Loop over all entries (signal)
      for(int ievt = 0; ievt < getTreeDual[0]->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTreeDual[0]->GetEntry(ievt);

         // Loop over all observables
         for(int i = 0; i < itemp[1]; i++)
         {
            // Check if observable is an angle in radians
            if( (observables[i] == "zenithSD") || (observables[i] == "azimuthSD") || (observables[i] == "zenithFD") || (observables[i] == "azimuthFD") || (observables[i] == "latitudeSD") || (observables[i] == "longitudeSD") || (observables[i] == "latitudeFD") || (observables[i] == "longitudeFD") )
	    {
               if( (observables[i] == "zenithSD") || (observables[i] == "zenithFD") )
                  iszenith[0] = true;
	       else
                  iszenith[0] = false;
               isangle[0] = true;
	    }
            else
               isangle[0] = false;

            // Check if observable is invalid and has value -1
            if(obsvars[3*i] == -1)
               isgood[0] = false;
            else
               isgood[0] = true;

            // Check the mean observable value in order to get minimum and maximum values of all trees
	    if(isangle[0])
	    {
               if(iszenith[0])
	          ftemp[0] = obsvars[3*i];
//	          ftemp[0] = SecTheta(obsvars[3*i]);
//	          ftemp[0] = RadToDeg(obsvars[3*i]);
	       else
	          ftemp[0] = RadToDeg(obsvars[3*i]);
	    }
	    else
	       ftemp[0] = obsvars[3*i];

            if(isgood[0])
            {
               if(ftemp[0] < minimumval[i])
                  minimumval[i] = ftemp[0];
               if(ftemp[0] > maximumval[i])
                  maximumval[i] = ftemp[0];
            }

            // Save observable values into the histograms
            if(isgood[0])
            {
               xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] = tempmax[i] + xhistlimit[2];

               if(ftemp[0] > xhistlimit[1])
                  oflowcount[i]++;
               else if(ftemp[0] < xhistlimit[0])
                  uflowcount[i]++;
               else
                  treeHistSig[i]->Fill(ftemp[0]);
                  
               totcountSig[i]++;
            }
         } // Loop over all observables
      } // Loop over all entries (signal)

      // Loop over all entries (background)
      for(int ievt = 0; ievt < getTreeDual[1]->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTreeDual[1]->GetEntry(ievt);

         // Loop over all observables
         for(int i = 0; i < itemp[1]; i++)
         {
            // Check if observable is an angle in radians
            if( (observables[i] == "zenithSD") || (observables[i] == "azimuthSD") || (observables[i] == "zenithFD") || (observables[i] == "azimuthFD") || (observables[i] == "latitudeSD") || (observables[i] == "longitudeSD") || (observables[i] == "latitudeFD") || (observables[i] == "longitudeFD") )
	    {
               if( (observables[i] == "zenithSD") || (observables[i] == "zenithFD") )
                  iszenith[0] = true;
	       else
                  iszenith[0] = false;
               isangle[0] = true;
	    }
            else
               isangle[0] = false;

            // Check if observable is invalid and has value -1
            if(obsvars[3*i+backShift] == -1)
               isgood[0] = false;
            else
               isgood[0] = true;

            // Check the mean observable value in order to get minimum and maximum values of all trees
	    if(isangle[0])
	    {
               if(iszenith[0])
	          ftemp[0] = obsvars[3*i+backShift];
//	          ftemp[0] = SecTheta(obsvars[3*i+backShift]);
//	          ftemp[0] = RadToDeg(obsvars[3*i+backShift]);
	       else
	          ftemp[0] = RadToDeg(obsvars[3*i+backShift]);
	    }
	    else
	       ftemp[0] = obsvars[3*i+backShift];

            if(isgood[0])
            {
               if(ftemp[0] < minimumval[i])
                  minimumval[i] = ftemp[0];
               if(ftemp[0] > maximumval[i])
                  maximumval[i] = ftemp[0];
            }

            // Save observable values into the histograms
            if(isgood[0])
            {
               xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] = tempmax[i] + xhistlimit[2];

               if(ftemp[0] > xhistlimit[1])
                  oflowcount[i]++;
               else if(ftemp[0] < xhistlimit[0])
                  uflowcount[i]++;
               else
               {
                  if(i == riseobs[0])
                     riseHist->Fill(ftemp[0]);

                  if(i == riseobs[1])
                     riserecalcHist->Fill(ftemp[0]);

                  treeHistBack[i]->Fill(ftemp[0]);
               }
                  
               totcountBack[i]++;
            }
         } // Loop over all observables
      } // Loop over all entries (background)

      // Write the valid event information
      cout << "Valid events (" << getTreeDual[0]->GetTitle() << " and " << getTreeDual[1]->GetTitle() << "):" << endl;
      for(int i = 0; i < itemp[1]; i++)
      {
         cout << " - " << observables[i] << " = " << totcountSig[i]+totcountBack[i] << " (UFLOW=" << uflowcount[i] << ", OFLOW=" << oflowcount[i] << ")" << endl;
      }
      cout << endl;

      // Plot all observables
      for(int i = 0; i < itemp[1]; i++)
      {
         c1->cd();

         // Check maximum value in a histogram
         if((treeHistSig[i]->GetMaximum()) > max[i])
            max[i] = treeHistSig[i]->GetMaximum();
         if((treeHistBack[i]->GetMaximum()) > max[i])
            max[i] = treeHistBack[i]->GetMaximum();

         // Set line and fill attributes for histograms
         mystyle->SetHistColor((TH1*)treeHistBack[i], 0);
         mystyle->SetHistColor((TH1*)treeHistSig[i], 1);

         // Draw signal and background histograms
         treeHistSig[i]->Draw();
         treeHistBack[i]->Draw("SAME");

/*	 cout << "Observable name = " << observables[i] << ", observable min = " << minimumval[i] << ", observable max = " << maximumval[i] << endl;

         // Set axis ranges
	 if( (observables[i] == "energyFD") )
	 {
            xhistlimit[2] = (maximumval[i]-minimumval[i])*0.05;
            xhistlimit[0] = minimumval[i] - xhistlimit[2];
            xhistlimit[1] = maximumval[i] + xhistlimit[2];
            cout << "  limits = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	 }
	 else if(observables[i] == "shwsize")
	 {
            xhistlimit[1] = TMath::Exp(-38.+2.238*TMath::Log10((minimumval[itemp[2]]+maximumval[itemp[2]])/2.));
            xhistlimit[2] = (xhistlimit[1]-tempmin[i])*0.05;
            xhistlimit[0] = tempmin[i] - xhistlimit[2];
            xhistlimit[1] += xhistlimit[2];
            cout << "  limits = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	 }
	 else
	 {*/
            xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
            xhistlimit[0] = tempmin[i] - xhistlimit[2];
            xhistlimit[1] = tempmax[i] + xhistlimit[2];
/*            cout << "  overwritten = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	 }*/
         treeHistSig[i]->GetYaxis()->SetRangeUser(0.,max[i]+0.25*max[i]);
         treeHistSig[i]->SetMaximum(max[i]*1.1);
         treeHistSig[i]->GetXaxis()->SetRange(xhistlimit[0], xhistlimit[1]);
         treeHistSig[i]->GetXaxis()->SetRangeUser(xhistlimit[0], xhistlimit[1]);

         // Setup axis (title, labels)
         mystyle->SetAxisTitles((TH1*)treeHistSig[i], tempdesc[i], "Number of events");

         // Draw a legend with info on signal and background events
	 nrfigs = 2;
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
         legend->SetFillStyle(legendFill);
         legend->SetFillColor(c_MvaCut);
         stemp[0] = "Pure signal (" + ToString(totcountSig[i]) + ")";
         legend->AddEntry(treeHistSig[i], stemp[0].c_str(),"f");
         stemp[0] = "Pure background (" + ToString(totcountBack[i]) + ")";
         legend->AddEntry(treeHistBack[i], stemp[0].c_str(),"f");
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");

         // Write down the underflow and overflow information
         uoflowtext.SetTextSize(0.017);
         uoflowtext.SetTextAlign(11);
         stemp[0] = "UFLOW = " + ToString(uflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[0], max[i]*1.105, stemp[0].c_str());
         uoflowtext.SetTextAlign(31);
         stemp[0] = "OFLOW = " + ToString(oflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[1], max[i]*1.105, stemp[0].c_str());

         // Prepare save name for plots
         stemp[1] = sigTree + "_" + backTree;
         stemp[3] = RemoveFilename(&filename) + "/histograms/dualhist_" + stemp[1] + "_" + observables[i] + ".pdf";

         // Save plots as PDF
         c1->SaveAs(stemp[3].c_str());
      } // Plot all observables

      // Run over signal tree observables

/*      // Zero minimum and maximum holders for observables
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
      
      // Run over background tree observables
      cout << endl << backTree << " has been selected for evaluation." << endl;
      getTree = (TTree*)ifile->Get(backTree.c_str());
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
      // Run over background tree observables*/
   }
   // Plot all trees separately
   else
   {
      // Give option to plot data with a different color
      string dataTree = "empty";
      cout << endl << "Select data tree to plot with different color (use -1, if not needed). Available trees:" << endl;
      for(int k = 0; k < ifile->GetNkeys(); k++)
      {
         stemp[1] = string((tempkeyslist->At(k))->GetName());
         cout << "  " << k+1 << " = " << stemp[1] << endl;
      }
      cout << "Choose option: ";
      cin >> itemp[0];
      cout << endl;

      if(itemp[0] != -1)
         dataTree = string((tempkeyslist->At(itemp[0]-1))->GetName());

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
/*	    cout << "Observable name = " << observables[i] << ", observable min = " << minimumval[i] << ", observable max = " << maximumval[i] << endl;

            // Set axis ranges
	    if( (observables[i] == "energyFD") )
	    {
               xhistlimit[2] = (maximumval[i]-minimumval[i])*0.05;
               xhistlimit[0] = minimumval[i] - xhistlimit[2];
               xhistlimit[1] = maximumval[i] + xhistlimit[2];
               cout << "  limits = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	    }
	    else if(observables[i] == "shwsize")
	    {
               xhistlimit[1] = TMath::Exp(-38.+2.238*TMath::Log10((minimumval[itemp[2]]+maximumval[itemp[2]])/2.));
               xhistlimit[2] = (xhistlimit[1]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] += xhistlimit[2];
               cout << "  limits = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	    }
	    else
	    {*/
               xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
               xhistlimit[0] = tempmin[i] - xhistlimit[2];
               xhistlimit[1] = tempmax[i] + xhistlimit[2];
/*               cout << "  overwritten = " << xhistlimit[0] << ", " << xhistlimit[1] << endl;
	    }*/

            // Prepare histograms for all observables
            stemp[2] = "treeHist" + ToString(i);
            treeHist[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
#if ROOTVER == 5
            treeHist[i]->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
            treeHist[i]->SetCanExtend(TH1::kXaxis);
#endif

            // Special histograms for risetime and recalculated risetime
            if(observables[i] == "risetime")
            {
               riseHist = new TH1F("riseHist", observables[i].c_str(), 100, maxobs[0], maxobs[1]);
#if ROOTVER == 5
               riseHist->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
               riseHist->SetCanExtend(TH1::kXaxis);
#endif
            }
            if(observables[i] == "risetimerecalc")
            {
               riserecalcHist = new TH1F("riserecalcHist", observables[i].c_str(), 100, maxobs[0], maxobs[1]);
#if ROOTVER == 5
               riserecalcHist->SetBit(TH1::kCanRebin);
#elif ROOTVER == 6
               riserecalcHist->SetCanExtend(TH1::kXaxis);
#endif
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
	       {
                  if( (observables[i] == "zenithSD") || (observables[i] == "zenithFD") )
                     iszenith[0] = true;
	          else
                     iszenith[0] = false;
                  isangle[0] = true;
	       }
               else
                  isangle[0] = false;

               // Check if observable is invalid and has value -1
               if(obsvars[3*i] == -1)
                  isgood[0] = false;
               else
                  isgood[0] = true;

               // Check the mean observable value in order to get minimum and maximum values of all trees
	       if(isangle[0])
	       {
                  if(iszenith[0])
	             ftemp[0] = obsvars[3*i];
//	             ftemp[0] = SecTheta(obsvars[3*i]);
//	             ftemp[0] = RadToDeg(obsvars[3*i]);
	          else
	             ftemp[0] = RadToDeg(obsvars[3*i]);
	       }
	       else
	          ftemp[0] = obsvars[3*i];

               if(isgood[0])
               {
                  if(ftemp[0] < minimumval[i])
                     minimumval[i] = ftemp[0];
                  if(ftemp[0] > maximumval[i])
                     maximumval[i] = ftemp[0];
               }

               // Save observable values into the histograms
               if(isgood[0])
               {
                  xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
                  xhistlimit[0] = tempmin[i] - xhistlimit[2];
                  xhistlimit[1] = tempmax[i] + xhistlimit[2];

                  if(ftemp[0] > xhistlimit[1])
                     oflowcount[i]++;
                  else if(ftemp[0] < xhistlimit[0])
                     uflowcount[i]++;
                  else
                  {
                     if(i == riseobs[0])
                        riseHist->Fill(ftemp[0]);

                     if(i == riseobs[1])
                        riserecalcHist->Fill(ftemp[0]);

                     treeHist[i]->Fill(ftemp[0]);
                  }
                     
                  totcount[i]++;
               }
            } // Loop over all observables

//            cout << endl;
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
            if(stemp[1] == dataTree)
               mystyle->SetHistColor((TH1*)treeHist[i], 2);
            else
               mystyle->SetHistColor((TH1*)treeHist[i], 1);

            // Draw signal and background histograms
            treeHist[i]->Draw();

            // Set axis ranges
            xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
            xhistlimit[0] = tempmin[i] - xhistlimit[2];
            xhistlimit[1] = tempmax[i] + xhistlimit[2];
            treeHist[i]->GetYaxis()->SetRangeUser(0.,max[i]+0.25*max[i]);
            treeHist[i]->SetMaximum(max[i]*1.1);
            treeHist[i]->GetXaxis()->SetRange(xhistlimit[0], xhistlimit[1]);
            treeHist[i]->GetXaxis()->SetRangeUser(xhistlimit[0], xhistlimit[1]);

            // Setup axis (title, labels)
            mystyle->SetAxisTitles((TH1*)treeHist[i], tempdesc[i], "Number of events");

            // Draw a legend with info on total number of events
	    nrfigs = 1;
            legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
            legend->SetFillStyle(legendFill);
            legend->SetFillColor(c_MvaCut);
            stemp[0] = "Total (" + ToString(totcount[i]) + ")";
            legend->AddEntry(treeHist[i], stemp[0].c_str(),"f");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");

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
         mystyle->SetHistColor((TH1*)riserecalcHist, 1);
         mystyle->SetHistColor((TH1*)riseHist, 0);

         // Draw signal and background histograms
         riserecalcHist->Draw();
         riseHist->Draw("same");

         // Set axis ranges
         riserecalcHist->GetYaxis()->SetRangeUser(0.,maxobs[3]+0.25*maxobs[3]);
         riseHist->GetYaxis()->SetRangeUser(0.,maxobs[3]+0.25*maxobs[3]);
         riserecalcHist->SetMaximum(maxobs[3]*1.2);
         riseHist->SetMaximum(maxobs[3]*1.2);
         riserecalcHist->GetXaxis()->SetRange(maxobs[0], maxobs[1]);
         riserecalcHist->GetXaxis()->SetRangeUser(maxobs[0], maxobs[1]);

         // Setup axis (title, labels)
         mystyle->SetAxisTitles((TH1*)riserecalcHist, tempdesc[riseobs[1]], "Number of events");

         // Draw a legend with info on signal and background events
	 nrfigs = 2;
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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
   }

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
