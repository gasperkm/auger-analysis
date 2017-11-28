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

   // Add information for the MVA variable
   observables.push_back("MVA");
   obssel.push_back(1);
   tempmin.push_back(-0.5);
   tempmax.push_back(1.5);
   tempdesc.push_back("MVA variable");

   // Save the number of observables (+MVA) into itemp[1]
   itemp[1] = observables.size();
   cout << "Number of observables (+MVA) = " << itemp[1] << endl;

   ifs.close();

   //-----------------------------------------
   // Prepare MVA cut values (any can be used)
   float mvacut[4];
   int mvacutapply;
   mvacut[1] = -1;
   mvacut[2] = -1;
   cerr << "Set the MVA cut value: ";
   cin >> mvacut[0];
   cerr << "Set negative MVA cut value (set to -1 if not needed): ";
   cin >> mvacut[1];
   cerr << "Set positive MVA cut value (set to -1 if not needed): ";
   cin >> mvacut[2];

   // Plots use the mean, negative or positive MVA value as splitting element
   cerr << "Perform cut on mean value (0), negative value (-1) or positive value (1): ";
   cin >> mvacutapply;
   if( (mvacutapply != 0) && (mvacutapply != -1) && (mvacutapply != 1) )
   {
      cerr << "Error! Wrong cut selected, rerun program." << endl;
      return 1;
   }

   if(mvacutapply == 0)
      mvacut[3] = mvacut[0];
   else if((mvacutapply == -1) && (mvacut[1] != -1))
      mvacut[3] = mvacut[1];
   else if((mvacutapply == 1) && (mvacut[2] != -1))
      mvacut[3] = mvacut[2];
   else
   {
      cerr << "Error! Wrong cut selected, rerun program." << endl;
      return 1;
   }
   
   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/created_plots";
   system(stemp[0].c_str());
   if(mvacutapply == 0)
   {
      stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/created_plots/mean";
      system(stemp[0].c_str());
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/created_plots/mean/mva_analysis* " + RemoveFilename(&filename) + "/created_plots/mean/scatter* " + RemoveFilename(&filename) + "/created_plots/mean/skymap_*";
      system(stemp[0].c_str());
   }
   else if(mvacutapply == -1)
   {
      stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/created_plots/negerror";
      system(stemp[0].c_str());
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/created_plots/negerror/mva_analysis* " + RemoveFilename(&filename) + "/created_plots/negerror/scatter* " + RemoveFilename(&filename) + "/created_plots/negerror/skymap_*";
      system(stemp[0].c_str());
   }
   else if(mvacutapply == 1)
   {
      stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/created_plots/poserror";
      system(stemp[0].c_str());
      stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/created_plots/poserror/mva_analysis* " + RemoveFilename(&filename) + "/created_plots/poserror/scatter* " + RemoveFilename(&filename) + "/created_plots/poserror/skymap_*";
      system(stemp[0].c_str());
   }

   // Prepare colors for signal, background and MVA cut line
   int c_SignalLine     = TColor::GetColor("#0000ee");
   int c_SignalFill     = TColor::GetColor("#7d99d1");
   int c_AllLine        = TColor::GetColor("#ff0000");
   int c_AllFill        = TColor::GetColor("#ff0000");
   int c_MvaCut         = TColor::GetColor("#ffff66");

   // Prepare plotting objects (legend, line, histograms, graphs,...) - for a maximum of 50 observables
   TLegend *legend;
   TLatex uoflowtext;
   TLine *line[3];
   TH1F *basesig[50];
   TH1F *baseback[50];
   TH2F *galacticaitoffsig;
   TH2F *galacticaitoffback;
   TH2F *horizontalaitoffsig;
   TH2F *horizontalaitoffback;
   TGraph *scatsig[50];
   TGraph *scatback[50];
//   TGraph *aitoffsig;

   // Prepare other observables
   int legendFill = 1001;		// filling style for legend
   float *xhistlimit;			// array for x axis limits
   xhistlimit = new float[3];
   int *galacticNr;			// vector positions for galactic longitude and latitude
   galacticNr = new int[2];
/*   float *galacticVal;			// vector for values of galactic longitude and latitude
   galacticVal = new float[2];
   float *projectVal;			// vector for values of galactic longitude and latitude (projected in Aitoff)
   projectVal = new float[2];*/
   int *horizontalNr;			// vector positions for horizontal azimuth and zenith
   horizontalNr = new int[2];
/*   float *horizontalVal;		// vector for values of horizontal azimuth and zenith
   horizontalVal = new float[2];*/
//   int cnt;				// counter
   int scatcnt;				// counter for scatter plots
   int scatcntsig[50];
   int scatcntback[50];
   double *dtemp;			// temporary deouble number
   dtemp = new double[2];

   float *max;				// maximum value in the histogram (in order to scale y axis so that histogram is below legend)
   max = new float[50];
   int *sigcount, *backcount;		// counting for signal and background events
   sigcount = new int[50];
   backcount = new int[50];

   float *obsvars;			// array for variable values that we read from output file
   obsvars = new float[150];

   int mvanumber = 3*(itemp[1]-1);	// position in array, where MVA variable is located

   int *uflowcount, *oflowcount; // counting number of underflow and overflow events
   uflowcount = new int[50];
   oflowcount = new int[50];

   bool isangle[2] = {false, false};	// check if observable is an angle (needs to be converted to degrees from radians) 
   bool isgood[2] = {true, true};	// check if observable has an invalid value -1

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
   for(int k = 0; k < ifile->GetNkeys()/2; k++)
   {
      // Zero minimum and maximum holders for observables
      for(int i = 0; i < itemp[1]; i++)
      {
         minimumval[i] = 1.e+30;
         maximumval[i] = -1.e+30;
	 uflowcount[i] = 0;
	 oflowcount[i] = 0;
      }

      // Select the trees (taking only the newest save with MVA variable values)
      stemp[1] = string((tempkeyslist->At(2*k))->GetName()) + ";2";
      cout << endl << stemp[1] << " has been selected for evaluation." << endl;
      getTree = (TTree*)ifile->Get(stemp[1].c_str());
      cout << "Number of events in the tree = " << getTree->GetEntries() << endl;

      // Set addresses for reading observables
      for(int i = 0; i < itemp[1]; i++)
      {
         // MVA variable
	 if(i == itemp[1]-1)
	 {
            cout << "Setting MVA observable " << observables[i] << endl;
	    stemp[2] = observables[i];
            getTree->SetBranchAddress(stemp[2].c_str(), &obsvars[mvanumber]); // MVA
            cout << "Mean variable: obsvars[3*" << i << "] = obsvars[" << mvanumber << "]" << endl;
	 }
	 // All other observables (mean, neg and pos values)
	 else
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
      }

      scatcnt = 0;

      // Prepare plots to save values into (histograms need to able to rebin)
      for(int i = 0; i < itemp[1]; i++)
      {
	 xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
         xhistlimit[0] = tempmin[i] - xhistlimit[2];
	 xhistlimit[1] = tempmax[i] + xhistlimit[2];

	 // Determine which observables are galactic latitude and longitude
	 if(observables[i] == "longitudeFD")
            galacticNr[0] = i;
	 if(observables[i] == "latitudeFD")
            galacticNr[1] = i;
	 // Determine which observables are horizontal azimuth and zenith
	 if(observables[i] == "azimuthFD")
            horizontalNr[0] = i;
	 if(observables[i] == "zenithFD")
            horizontalNr[1] = i;

	 // Prepare histograms for all observables
	 stemp[2] = "basesig" + ToString(i);
         basesig[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
         basesig[i]->SetBit(TH1::kCanRebin);
         stemp[2] = "baseback" + ToString(i);
         baseback[i] = new TH1F(stemp[2].c_str(), observables[i].c_str(), 100, xhistlimit[0], xhistlimit[1]);
         baseback[i]->SetBit(TH1::kCanRebin);

	 // Set maximum value to zero, and signal and background counters to zero
         max[i] = 0.0;
         sigcount[i] = 0;
         backcount[i] = 0;
         uflowcount[i] = 0;
         oflowcount[i] = 0;

	 // Prepare scatter plots for all combinations of observables
	 for(int j = 0; j < itemp[1]; j++)
	 {
            if( (i < j) && (obssel[i] == 1) && (obssel[j] == 1) )
	    {
	       cout << "Preparing scatter plot " << scatcnt << " (" << i << ", " << j << ")" << endl;
               scatsig[scatcnt] = new TGraph();
               scatback[scatcnt] = new TGraph();
	       scatcnt++;
	    }
	 }
      }

      // Prepare 2D histograms for aitoff projection of galactic longitude and latitude (sky maps)
      galacticaitoffsig = new TH2F("galacticaitoffsig", "Sky map (signal)", 360, (tempmin[galacticNr[0]]-180.), (tempmax[galacticNr[0]]-180.), 180, tempmin[galacticNr[1]], tempmax[galacticNr[1]]);
      galacticaitoffback = new TH2F("galacticaitoffback", "Sky map (background)", 360, (tempmin[galacticNr[0]]-180.), (tempmax[galacticNr[0]]-180.), 180, tempmin[galacticNr[1]], tempmax[galacticNr[1]]);
      // Prepare 2D histograms for aitoff projection of horizontal azimuth and zenith (sky maps)
      horizontalaitoffsig = new TH2F("horizontalaitoffsig", "Sky map (signal)", 360, (tempmin[horizontalNr[0]]-180.), (tempmax[horizontalNr[0]]-180.), 180, tempmin[horizontalNr[1]], tempmax[horizontalNr[1]]);
      horizontalaitoffback = new TH2F("horizontalaitoffback", "Sky map (background)", 360, (tempmin[horizontalNr[0]]-180.), (tempmax[horizontalNr[0]]-180.), 180, tempmin[horizontalNr[1]], tempmax[horizontalNr[1]]);

//      aitoffsig = new TGraph();

      for(int i = 0; i < 50; i++)
      {
         scatcntsig[i] = 1;
         scatcntback[i] = 1;
      }

      // Loop over all entries
      for(int ievt = 0; ievt < getTree->GetEntries(); ievt++)
      {
         // Get an entry from the tree
         getTree->GetEntry(ievt);
//         cnt = 0;
         scatcnt = 0;

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

	    // Write out information for all observables, cut and event number
/*            if(i == itemp[1]-1)
	    {
               cout << "  MVA Event = " << ievt << ": value = " << obsvars[mvanumber] << " (cut = " << mvacut[3] << ") , observable = " << observables[i] << endl;
	    }
	    else
	    {
               if(isangle[0])
                  cout << "  Event = " << ievt << ": values = " << ToDeg(obsvars[3*i]) << ", " << ToDeg(obsvars[3*i+1]) << ", " << ToDeg(obsvars[3*i+2]) << ", observable = " << observables[i] << endl;
               else
                  cout << "  Event = " << ievt << ": values = " << obsvars[3*i] << ", " << obsvars[3*i+1] << ", " << obsvars[3*i+2] << ", observable = " << observables[i] << endl;
	    }*/

	    xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
            xhistlimit[0] = tempmin[i] - xhistlimit[2];
	    xhistlimit[1] = tempmax[i] + xhistlimit[2];

	    // Count observables, if above MVA cut (signal) and save them to signal plots
            if(obsvars[mvanumber] >= mvacut[3])
            {
	       if(isgood[0])
	       {
	          if(obsvars[3*i] > xhistlimit[1])
                     oflowcount[i]++;
	          else if(obsvars[3*i] < xhistlimit[0])
                     uflowcount[i]++;
                  else
	          {
                     if(isangle[0])
                        basesig[i]->Fill(RadToDeg(obsvars[3*i]));
	             else
                        basesig[i]->Fill(obsvars[3*i]);
		  }
                  sigcount[i]++;
	       }
	    }
	    // Count observables, if below MVA cut (background) and save them to background plots
	    else
	    {
	       if(isgood[0])
	       {
	          if(obsvars[3*i] > xhistlimit[1])
                     oflowcount[i]++;
	          else if(obsvars[3*i] < xhistlimit[0])
                     uflowcount[i]++;
                  else
	          {
                     if(isangle[0])
                        baseback[i]->Fill(RadToDeg(obsvars[3*i]));
	             else
                        baseback[i]->Fill(obsvars[3*i]);
		  }
	          backcount[i]++;
	       }
	    }

	    // Loop over all observables (scatter plots)
            for(int j = 0; j < itemp[1]; j++)
	    {
               if( (i < j) && (obssel[i] == 1) && (obssel[j] == 1) )
	       {
                  // Check if observable is an angle in radians
                  if( (observables[j] == "zenithSD") || (observables[j] == "azimuthSD") || (observables[j] == "zenithFD") || (observables[j] == "azimuthFD") || (observables[j] == "latitudeSD") || (observables[j] == "longitudeSD") || (observables[j] == "latitudeFD") || (observables[j] == "longitudeFD") )
                     isangle[1] = true;
	          else
                     isangle[1] = false;

                  // Check if observable is invalid and has value -1
	          if(obsvars[3*j] == -1)
                     isgood[1] = false;
	          else
                     isgood[1] = true;

	          // Count observables, if above MVA cut (signal) and save them to signal scatter plot
                  if(obsvars[mvanumber] >= mvacut[3])
                  {
	             if(isgood[0] && isgood[1])
	             {
			if(isangle[0] && isangle[1])
                           scatsig[scatcnt]->SetPoint(scatcntsig[scatcnt], RadToDeg(obsvars[3*i]), RadToDeg(obsvars[3*j]));
                        if(isangle[0])
                           scatsig[scatcnt]->SetPoint(scatcntsig[scatcnt], RadToDeg(obsvars[3*i]), obsvars[3*j]);
			else if(isangle[1])
                           scatsig[scatcnt]->SetPoint(scatcntsig[scatcnt], obsvars[3*i], RadToDeg(obsvars[3*j]));
	                else
                           scatsig[scatcnt]->SetPoint(scatcntsig[scatcnt], obsvars[3*i], obsvars[3*j]);

			scatcntsig[scatcnt]++;
	             }
	          }
	          // Count observables, if below MVA cut (background) and save them to background scatter plot
	          else
	          {
	             if(isgood[0] && isgood[1])
	             {
			if(isangle[0] && isangle[1])
                           scatback[scatcnt]->SetPoint(scatcntback[scatcnt], RadToDeg(obsvars[3*i]), RadToDeg(obsvars[3*j]));
                        if(isangle[0])
                           scatback[scatcnt]->SetPoint(scatcntback[scatcnt], RadToDeg(obsvars[3*i]), obsvars[3*j]);
			else if(isangle[1])
                           scatback[scatcnt]->SetPoint(scatcntback[scatcnt], obsvars[3*i], RadToDeg(obsvars[3*j]));
	                else
                           scatback[scatcnt]->SetPoint(scatcntback[scatcnt], obsvars[3*i], obsvars[3*j]);

			scatcntback[scatcnt]++;
	             }
	          }

	          if(isgood[0] && isgood[1])
		     scatcnt++;
	       }
	    } // Loop over all observables (scatter plots)
         } // Loop over all observables

	 // Fill 2D histogram sky maps
         if(obsvars[mvanumber] >= mvacut[3])
	 {
	    galacticaitoffsig->Fill( (RadToDeg(obsvars[3*galacticNr[0]])-180.), RadToDeg(obsvars[3*galacticNr[1]]));
	    horizontalaitoffsig->Fill( (RadToDeg(obsvars[3*horizontalNr[0]])-180.), RadToDeg(obsvars[3*horizontalNr[1]]));

//            GalactProject(obsvars[3*galacticNr[0]], obsvars[3*galacticNr[1]], galacticVal);
//	    aitoffsig->SetPoint(ievt, galacticVal[0], galacticVal[1]);
	 }
	 else
	 {
	    galacticaitoffback->Fill( (RadToDeg(obsvars[3*galacticNr[0]])-180.), RadToDeg(obsvars[3*galacticNr[1]]));
	    horizontalaitoffback->Fill( (RadToDeg(obsvars[3*horizontalNr[0]])-180.), RadToDeg(obsvars[3*horizontalNr[1]]));
	 }

//         cout << endl;
      } // Loop over all entries

      // Write the signal vs. background count information
      cout << "Signal vs. background (" << getTree->GetTitle() << "):" << endl;
      for(int i = 0; i < itemp[1]; i++)
      {
         if(i == itemp[1]-1)
            cout << " - MVA = " << sigcount[i] << " vs. " << backcount[i] << endl;
	 else
            cout << " - " << observables[i] << " = " << sigcount[i] << " vs. " << backcount[i] << endl;
      }
      cout << endl;

      // Plot all observables
      for(int i = 0; i < itemp[1]; i++)
      {
         c1->cd();

         // Check maximum value in a histogram
         if((basesig[i]->GetMaximum()) > max[i])
            max[i] = basesig[i]->GetMaximum();
         if((baseback[i]->GetMaximum()) > max[i])
            max[i] = baseback[i]->GetMaximum();

	 // Set line and fill attributes for histograms
         basesig[i]->SetLineColor(c_SignalLine);
         basesig[i]->SetLineWidth(2);
         basesig[i]->SetFillColor(c_SignalFill);
         basesig[i]->SetFillStyle(1001);
         baseback[i]->SetLineColor(c_AllLine);
         baseback[i]->SetLineWidth(2);
         baseback[i]->SetFillColor(c_AllFill);
         baseback[i]->SetFillStyle(3554);

	 // Draw signal and background histograms
         basesig[i]->Draw();
         baseback[i]->Draw("same");

	 // Set axis ranges
	 xhistlimit[2] = (tempmax[i]-tempmin[i])*0.05;
         xhistlimit[0] = tempmin[i] - xhistlimit[2];
	 xhistlimit[1] = tempmax[i] + xhistlimit[2];
/*	 if(i == itemp[1]-1)
            basesig[i]->GetXaxis()->SetRangeUser(-0.5, 1.5);*/
         basesig[i]->GetYaxis()->SetRangeUser(0.,max[i]*1.2);
         basesig[i]->SetMaximum(max[i]*1.2);
         basesig[i]->GetXaxis()->SetRange(xhistlimit[0], xhistlimit[1]);
         basesig[i]->GetXaxis()->SetRangeUser(xhistlimit[0], xhistlimit[1]);

	 // Setup axis (title, labels)
         basesig[i]->GetYaxis()->SetTitle("Number of events");
         basesig[i]->GetXaxis()->SetTitle(tempdesc[i].c_str());
         basesig[i]->GetXaxis()->SetTitleOffset(1.2);
         basesig[i]->GetXaxis()->CenterTitle(kTRUE);
         basesig[i]->GetXaxis()->SetLabelSize(0.028);
         basesig[i]->GetXaxis()->SetLabelOffset(0.015);
         basesig[i]->GetYaxis()->SetTitleOffset(1.3);
         basesig[i]->GetYaxis()->CenterTitle(kTRUE);
         basesig[i]->GetYaxis()->SetLabelSize(0.028);
         basesig[i]->GetYaxis()->SetLabelOffset(0.015);
         basesig[i]->SetTitle("");

         // Write down the underflow and overflow information
         uoflowtext.SetTextSize(0.017);
         uoflowtext.SetTextAlign(11);
         stemp[0] = "UFLOW = " + ToString(uflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[0], max[i]*1.205, stemp[0].c_str());
         uoflowtext.SetTextAlign(31);
         stemp[0] = "OFLOW = " + ToString(oflowcount[i]);
         uoflowtext.DrawLatex(xhistlimit[1], max[i]*1.205, stemp[0].c_str());

	 // Add MVA cut lines to MVA variable histograms (mean, neg and pos)
	 if(i == itemp[1]-1)
	 {
            line[0] = new TLine(mvacut[0], 0., mvacut[0], basesig[i]->GetMaximum());
            line[0]->SetLineWidth(2);
            line[0]->SetLineStyle(1);
            line[0]->SetLineColor(kOrange+2);
            line[0]->Draw("same");

	    if(mvacut[1] != -1)
	    {
               line[1] = new TLine(mvacut[1], 0., mvacut[1], basesig[i]->GetMaximum());
               line[1]->SetLineWidth(2);
               line[1]->SetLineStyle(7);
               line[1]->SetLineColor(kOrange+2);
               line[1]->Draw("same");
	    }

	    if(mvacut[2] != -1)
	    {
               line[2] = new TLine(mvacut[2], 0., mvacut[2], basesig[i]->GetMaximum());
               line[2]->SetLineWidth(2);
               line[2]->SetLineStyle(7);
               line[2]->SetLineColor(kOrange+2);
               line[2]->Draw("same");
	    }
	 }

         // Draw a legend with info on signal and background events
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.12, gPad->GetLeftMargin()+.45, 1-gPad->GetTopMargin());
         legend->SetFillStyle(legendFill);
         legend->SetFillColor(c_MvaCut);
	 dtemp[0] = 100.*(double)backcount[i]/((double)sigcount[i] + (double)backcount[i]);
         stemp[0] = "MVA cut background events (" + ToString(backcount[i]) + " = " + ToString(dtemp[0], 2) + "%)";
         legend->AddEntry(baseback[i],stemp[0].c_str(),"f");
	 dtemp[0] = 100.*(double)sigcount[i]/((double)sigcount[i] + (double)backcount[i]);
         stemp[0] = "MVA cut signal events (" + ToString(sigcount[i]) + " = " + ToString(dtemp[0], 2) + "%)";
         legend->AddEntry(basesig[i],stemp[0].c_str(),"f");
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");

	 // Prepare save name for plots
         stemp[1] = string((tempkeyslist->At(2*k))->GetName());
         if(mvacutapply == 0)
            stemp[3] = RemoveFilename(&filename) + "/created_plots/mean/mva_analysis_" + stemp[1] + "_" + observables[i] + ".pdf";
	 else if(mvacutapply == -1)
            stemp[3] = RemoveFilename(&filename) + "/created_plots/negerror/mva_analysis_" + stemp[1] + "_" + observables[i] + ".pdf";
	 else if(mvacutapply == 1)
            stemp[3] = RemoveFilename(&filename) + "/created_plots/poserror/mva_analysis_" + stemp[1] + "_" + observables[i] + ".pdf";
	 else
	 {
            cerr << "Error! Wrong cut selected, rerun program." << endl;
            return 1;
	 }

	 // Save plots as PDF
         c1->SaveAs(stemp[3].c_str());

         delete legend;
      } // Plot all observables

      scatcnt = 0;
      // Plot all scatter plots
      for(int i = 0; i < itemp[1]; i++)
      {
         for(int j = 0; j < itemp[1]; j++)
	 {
            if( (i < j) && (obssel[i] == 1) && (obssel[j] == 1) )
	    {
               // Prepare scatter plots style
	       scatsig[scatcnt]->SetMarkerStyle(20);
	       scatsig[scatcnt]->SetMarkerColor(c_SignalLine);
	       scatsig[scatcnt]->SetMarkerSize(0.7);
	       scatback[scatcnt]->SetMarkerStyle(21);
	       scatback[scatcnt]->SetMarkerColor(c_AllLine);
	       scatback[scatcnt]->SetMarkerSize(0.7);

	       // Set axis ranges
/*	       if(i == itemp[1]-1)
                  scatsig[scatcnt]->GetXaxis()->SetRangeUser(-0.5, 1.5);
	       else*/
		  dtemp[0] = (maximumval[i] - minimumval[i])*0.05;
                  scatsig[scatcnt]->GetXaxis()->SetRangeUser(minimumval[i] - dtemp[0], maximumval[i] + dtemp[0]);
/*	       if(j == itemp[1]-1)
                  scatsig[scatcnt]->GetYaxis()->SetRangeUser(-0.5, 1.5);
	       else*/
		  dtemp[1] = (maximumval[j] - minimumval[j])*0.05;
                  scatsig[scatcnt]->GetYaxis()->SetRangeUser(minimumval[j] - dtemp[1], maximumval[j] + dtemp[1]);

               // Setup axis (title, labels)
               scatsig[scatcnt]->GetYaxis()->SetTitle(tempdesc[j].c_str());
               scatsig[scatcnt]->GetXaxis()->SetTitle(tempdesc[i].c_str());
               scatsig[scatcnt]->GetXaxis()->SetTitleOffset(1.2);
               scatsig[scatcnt]->GetXaxis()->CenterTitle(kTRUE);
               scatsig[scatcnt]->GetXaxis()->SetLabelSize(0.028);
               scatsig[scatcnt]->GetXaxis()->SetLabelOffset(0.015);
               scatsig[scatcnt]->GetYaxis()->SetTitleOffset(1.3);
               scatsig[scatcnt]->GetYaxis()->CenterTitle(kTRUE);
               scatsig[scatcnt]->GetYaxis()->SetLabelSize(0.028);
               scatsig[scatcnt]->GetYaxis()->SetLabelOffset(0.015);
               scatsig[scatcnt]->SetTitle("");

               // Draw signal and background histograms
               scatsig[scatcnt]->Draw("AP");
               scatback[scatcnt]->Draw("P;SAME");

/*               // Draw a legend with info on signal and background events
               legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.12, gPad->GetLeftMargin()+.28, 1-gPad->GetTopMargin());
               legend->SetFillStyle(legendFill);
               legend->SetFillColor(c_MvaCut);
               stemp[0] = "MVA cut background events (" + ToString(backcount[i]) + ")";
               legend->AddEntry(scatback[scatcnt],stemp[0].c_str(),"p");
               stemp[0] = "MVA cut signal events (" + ToString(sigcount[i]) + ")";
               legend->AddEntry(scatsig[scatcnt],stemp[0].c_str(),"p");
               legend->SetBorderSize(1);
               legend->SetMargin(0.3);
               legend->Draw("SAME");*/

               // Prepare save name for scatter plots
               stemp[1] = string((tempkeyslist->At(2*k))->GetName());
               if(mvacutapply == 0)
                  stemp[3] = RemoveFilename(&filename) + "/created_plots/mean/scatter_" + stemp[1] + "_" + observables[i] + "-" + observables[j] + ".pdf";
               else if(mvacutapply == -1)
                  stemp[3] = RemoveFilename(&filename) + "/created_plots/negerror/scatter_" + stemp[1] + "_" + observables[i] + "-" + observables[j] + ".pdf";
               else if(mvacutapply == 1)
                  stemp[3] = RemoveFilename(&filename) + "/created_plots/poserror/scatter_" + stemp[1] + "_" + observables[i] + "-" + observables[j] + ".pdf";
               else
               {
                  cerr << "Error! Wrong cut selected, rerun program." << endl;
                  return 1;
               }
              
               // Save plots as PDF
               c1->SaveAs(stemp[3].c_str());
              
//               delete legend;

	       scatcnt++;
	    }
	 }
      } // Plot all scatter plots

      // Prepare save name for sky maps (galactic)
      stemp[1] = string((tempkeyslist->At(2*k))->GetName());
      if(mvacutapply == 0)
      {
         stemp[3] = RemoveFilename(&filename) + "/created_plots/mean/skymap_galactic_" + stemp[1] + "_signal.png";
         stemp[4] = RemoveFilename(&filename) + "/created_plots/mean/skymap_galactic_" + stemp[1] + "_background.png";
      }
      else if(mvacutapply == -1)
      {
         stemp[3] = RemoveFilename(&filename) + "/created_plots/negerror/skymap_galactic_" + stemp[1] + "_signal.png";
         stemp[4] = RemoveFilename(&filename) + "/created_plots/negerror/skymap_galactic_" + stemp[1] + "_background.png";
      }
      else if(mvacutapply == 1)
      {
         stemp[3] = RemoveFilename(&filename) + "/created_plots/poserror/skymap_galactic_" + stemp[1] + "_signal.png";
         stemp[4] = RemoveFilename(&filename) + "/created_plots/poserror/skymap_galactic_" + stemp[1] + "_background.png";
      }
      else
      {
         cerr << "Error! Wrong cut selected, rerun program." << endl;
         return 1;
      }

      // Plot 2D histogram sky maps (galactic)
      c2->cd();
      galacticaitoffsig->GetYaxis()->SetTitle("Galactic lattitude (deg)");
      galacticaitoffsig->GetXaxis()->SetTitle("Galactic longitude (deg)");
      galacticaitoffsig->GetXaxis()->SetTitleOffset(1.2);
      galacticaitoffsig->GetXaxis()->CenterTitle(kTRUE);
      galacticaitoffsig->GetXaxis()->SetLabelSize(0.028);
      galacticaitoffsig->GetXaxis()->SetLabelOffset(0.015);
      galacticaitoffsig->GetYaxis()->SetTitleOffset(1.3);
      galacticaitoffsig->GetYaxis()->CenterTitle(kTRUE);
      galacticaitoffsig->GetYaxis()->SetLabelSize(0.028);
      galacticaitoffsig->GetYaxis()->SetLabelOffset(0.015);
      galacticaitoffsig->SetTitle("Sky map (signal)");
      galacticaitoffsig->Draw("AITOFF");
      c2->SaveAs(stemp[3].c_str());
      galacticaitoffback->GetYaxis()->SetTitle("Galactic lattitude (deg)");
      galacticaitoffback->GetXaxis()->SetTitle("Galactic longitude (deg)");
      galacticaitoffback->GetXaxis()->SetTitleOffset(1.2);
      galacticaitoffback->GetXaxis()->CenterTitle(kTRUE);
      galacticaitoffback->GetXaxis()->SetLabelSize(0.028);
      galacticaitoffback->GetXaxis()->SetLabelOffset(0.015);
      galacticaitoffback->GetYaxis()->SetTitleOffset(1.3);
      galacticaitoffback->GetYaxis()->CenterTitle(kTRUE);
      galacticaitoffback->GetYaxis()->SetLabelSize(0.028);
      galacticaitoffback->GetYaxis()->SetLabelOffset(0.015);
      galacticaitoffback->SetTitle("Sky map (background)");
      galacticaitoffback->Draw("AITOFF");
      c2->SaveAs(stemp[4].c_str());

      // Prepare save name for sky maps (horizontal)
      stemp[1] = string((tempkeyslist->At(2*k))->GetName());
      if(mvacutapply == 0)
      {
         stemp[3] = RemoveFilename(&filename) + "/created_plots/mean/skymap_horizontal_" + stemp[1] + "_signal.png";
         stemp[4] = RemoveFilename(&filename) + "/created_plots/mean/skymap_horizontal_" + stemp[1] + "_background.png";
      }
      else if(mvacutapply == -1)
      {
         stemp[3] = RemoveFilename(&filename) + "/created_plots/negerror/skymap_horizontal_" + stemp[1] + "_signal.png";
         stemp[4] = RemoveFilename(&filename) + "/created_plots/negerror/skymap_horizontal_" + stemp[1] + "_background.png";
      }
      else if(mvacutapply == 1)
      {
         stemp[3] = RemoveFilename(&filename) + "/created_plots/poserror/skymap_horizontal_" + stemp[1] + "_signal.png";
         stemp[4] = RemoveFilename(&filename) + "/created_plots/poserror/skymap_horizontal_" + stemp[1] + "_background.png";
      }
      else
      {
         cerr << "Error! Wrong cut selected, rerun program." << endl;
         return 1;
      }

      // Plot 2D histogram sky maps (horizontal)
      horizontalaitoffsig->GetYaxis()->SetTitle("FD Azimuth (deg)");
      horizontalaitoffsig->GetXaxis()->SetTitle("FD Zenith (deg)");
      horizontalaitoffsig->GetXaxis()->SetTitleOffset(1.2);
      horizontalaitoffsig->GetXaxis()->CenterTitle(kTRUE);
      horizontalaitoffsig->GetXaxis()->SetLabelSize(0.028);
      horizontalaitoffsig->GetXaxis()->SetLabelOffset(0.015);
      horizontalaitoffsig->GetYaxis()->SetTitleOffset(1.3);
      horizontalaitoffsig->GetYaxis()->CenterTitle(kTRUE);
      horizontalaitoffsig->GetYaxis()->SetLabelSize(0.028);
      horizontalaitoffsig->GetYaxis()->SetLabelOffset(0.015);
      horizontalaitoffsig->SetTitle("Sky map (signal)");
      horizontalaitoffsig->Draw("AITOFF");
      c2->SaveAs(stemp[3].c_str());
      horizontalaitoffback->GetYaxis()->SetTitle("FD Azimuth (deg)");
      horizontalaitoffback->GetXaxis()->SetTitle("FD Zenith (deg)");
      horizontalaitoffback->GetXaxis()->SetTitleOffset(1.2);
      horizontalaitoffback->GetXaxis()->CenterTitle(kTRUE);
      horizontalaitoffback->GetXaxis()->SetLabelSize(0.028);
      horizontalaitoffback->GetXaxis()->SetLabelOffset(0.015);
      horizontalaitoffback->GetYaxis()->SetTitleOffset(1.3);
      horizontalaitoffback->GetYaxis()->CenterTitle(kTRUE);
      horizontalaitoffback->GetYaxis()->SetLabelSize(0.028);
      horizontalaitoffback->GetYaxis()->SetLabelOffset(0.015);
      horizontalaitoffback->SetTitle("Sky map (background)");
      horizontalaitoffback->Draw("AITOFF");
      c2->SaveAs(stemp[4].c_str());

/*      stemp[3] = "created_plots/skymap_aitoff_" + stemp[1] + "_signal.pdf";
      aitoffsig->SetMarkerStyle(20);
      aitoffsig->SetMarkerSize(0.7);
      aitoffsig->SetMarkerColor(2);
      aitoffsig->Draw("AP");
      c2->SaveAs(stemp[3].c_str());*/

      // Delete all plot objects
      for(int i = 0; i < itemp[1]; i++)
      {
         delete basesig[i];
         delete baseback[i];
      }
      delete galacticaitoffsig;
      delete galacticaitoffback;
      delete horizontalaitoffsig;
      delete horizontalaitoffback;
//      delete aitoffsig;
   } // Loop over all trees

   ifile->Close();

/*   // Minimum and maximum values over all observables and trees
   cout << endl << "Observables minimum and maximum values:" << endl;
   for(int i = 0; i < itemp[1]; i++)
   {
      cout << observables[i] << ":\t" << minimumval[i] << "\t" << maximumval[i] << endl;
   }*/

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
   delete[] galacticNr;
//   delete[] galacticVal;
//   delete[] projectVal;
   delete[] horizontalNr;
//   delete[] horizontalVal;
   delete[] max;
   delete[] obsvars;
   delete[] sigcount;
   delete[] backcount;
   delete[] minimumval;
   delete[] maximumval;

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
