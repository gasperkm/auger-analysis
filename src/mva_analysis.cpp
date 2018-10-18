#include "frame.h"
#include "separate_functions.h"
#include "adst_mva.h"
#include "popups.h"
#include "combine.h"
#include "mvaefficiency.h"
#include <algorithm>

// Add observables to the MVA analysis
int MyFrame::MvaNoteObservables(int count)
{
   string *stemp;

   // Write out the number of different selected observables
   stemp = new string;
   if(DBGSIG > 0)
      cout << "# MvaNoteObservables    #: " << "Number of different selected observables: " << count << endl;
   *stemp = string(rootdir) + "/results/observables_nr.dat";
   ofstream *fobser = new ofstream;
   fobser->open(stemp->c_str(), ofstream::out | ofstream::trunc);
   if(fobser->is_open())
      *fobser << count << endl;
   fobser->close();

   delete fobser;
   delete stemp;
   return count;
}

// Write out the MVA trees to contain only one of the eyes (best eye inside selection cuts or average of all eyes)
int MyFrame::MvaTreeFile(string *infilename, string *outfilename, int *nrEvents, int curenbin)
{
   string *stemp = new string[2];

   // Open the input file
   TFile *ifile = TFile::Open(infilename->c_str(), "READ");
   TTree *signalTempTree[nrkeys];

   // Cancel MVA analysis in case signal and background trees are the same
   if( (signalSelect->widgetCB)->GetSelection() == (backgroundSelect->widgetCB)->GetSelection() )
   {
      delete[] stemp;
      ifile->Close();
      return -1;
   }

   // Determine the DeltaS38 and Risetime Delta observables and implement their fits depending on the selected data tree
   if(multipleEnergyBins)
      stemp[0] = "mkdir -p " + string(*currentAnalysisDir) + "/../delta_conversion";
   else
      stemp[0] = "mkdir -p " + string(*currentAnalysisDir) + "/delta_conversion";

   if( (multipleEnergyBins && (curenbin == 0)) || (!multipleEnergyBins) )
   {
      ret = system(stemp[0].c_str());
      ret = SetDeltas(0, (dataSelect->widgetCB)->GetSelection()+1, ifile, true);
      ret = SetDeltas(1, (dataSelect->widgetCB)->GetSelection()+1, ifile, true);
      for(int j = 1; j <= nrkeys; j++)
      {
         stemp[0] = "TreeS" + ToString(j);
         stemp[1] = string(ifile->GetKey(stemp[0].c_str())->GetTitle());
         stemp[1] = RemovePath(&stemp[1]);

         if( string((dataSelect->widgetCB)->GetStringSelection()) == stemp[1] )
            j++;
         else
         {
            ret = SetDeltas(0, j, ifile, false);
            ret = SetDeltas(1, j, ifile, false);
         }
      }
   }
   
   // Open the output file and run the MvaSetTrees
   stemp[0] = "rm -fr " + *outfilename;
   system(stemp[0].c_str());
   TFile *ofile = TFile::Open(outfilename->c_str(), "RECREATE");
   for(int j = 1; j <= nrkeys; j++)
   {
      stemp[0] = "TreeS" + ToString(j);
      stemp[1] = string(ifile->GetKey(stemp[0].c_str())->GetTitle());
      stemp[1] = RemovePath(&stemp[1]);

      if(DBGSIG > 0)
      {
         if( string((signalSelect->widgetCB)->GetStringSelection()) == stemp[1] )
            cout << "# MvaTreeFile           #: " << "Using signal tree: " << stemp[1] << endl;
         else if( string((backgroundSelect->widgetCB)->GetStringSelection()) == stemp[1] )
            cout << "# MvaTreeFile           #: " << "Using background tree: " << stemp[1] << endl;
         else if( string((dataSelect->widgetCB)->GetStringSelection()) == stemp[1] )
            cout << "# MvaTreeFile           #: " << "Using data tree: " << stemp[1] << endl;
         else
            cout << "# MvaTreeFile           #: " << "Using other tree: " << stemp[1] << endl;
      }

      signalTempTree[j-1] = new TTree(stemp[0].c_str(), stemp[1].c_str());
      nrEvents[j-1] = MvaSetTrees(j, ifile, signalTempTree[j-1]);
      if(nrEvents[j-1] == -1)
      {
	 for(int i = 0; i < j-1; i++)
            delete signalTempTree[i];
         delete[] stemp;
         ofile->Close();
         ifile->Close();
         return -1;
      }
      signalTempTree[j-1]->Write();

      delete signalTempTree[j-1];
   }

   ofile->Close();
   ifile->Close();
   
   delete[] stemp;
   return 0;
}

// Set MVA trees and add them to the Factory
int MyFrame::MvaSetTrees(int type, TFile *ifile, TTree *outtree)
{
   string *stemp;
   int *itemp;
   float *ftemp;

   stemp = new string[3];
   itemp = new int[4];
   ftemp = new float[11];

   // SD observables for cut
   if(selcuttype == 0)
   {
      ret = Find(observables, "energySD");
      if(ret == -1)
      {
         AlertPopup("No SD energy observable found", "No SD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energySD.");
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      ret = Find(observables, "zenithSD");
      if(ret == -1)
      {
         AlertPopup("No SD zenith angle observable found", "No SD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithSD.");
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }
   }
   // FD observables for cut
   else if(selcuttype == 1)
   {
      ret = Find(observables, "energyFD");
      if(ret == -1)
      {
         AlertPopup("No FD energy observable found", "No FD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energyFD.");
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      ret = Find(observables, "zenithFD");
      if(ret == -1)
      {
         AlertPopup("No FD zenith angle observable found", "No FD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithFD.");
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }
   }

   ret = Find(observables, "risetimerecalc");
   if(ret == -1)
   {
      AlertPopup("No risetime observable found", "No risetime observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable risetimerecalc.");
      delete[] stemp;
      delete[] itemp;
      delete[] ftemp;
      return -1;
   }

   // Prepare the output tree and observables for them
   float *outobs, *outobs_neg, *outobs_pos;
   outobs = new float[nrobs];
   outobs_neg = new float[nrobs];
   outobs_pos = new float[nrobs];

   stemp[0] = "TreeS" + ToString(type);
   if(DBGSIG > 0)
      cout << "# MvaSetTrees           #: " << stemp[0] << endl;
   for(int i = 0; i < nrobs; i++)
   {
      stemp[1] = observables[i];
      stemp[2] = observables[i] + "/F";
      outtree->Branch(stemp[1].c_str(), &outobs[i], stemp[2].c_str());
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Mean out tree: " << stemp[1] << ", " << stemp[2] << endl;

      stemp[1] = observables[i] + "_neg";
      stemp[2] = observables[i] + "_neg/F";
      outtree->Branch(stemp[1].c_str(), &outobs_neg[i], stemp[2].c_str());
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Neg out tree: " << stemp[1] << ", " << stemp[2] << endl;

      stemp[1] = observables[i] + "_pos";
      stemp[2] = observables[i] + "_pos/F";
      outtree->Branch(stemp[1].c_str(), &outobs_pos[i], stemp[2].c_str());
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Pos out tree: " << stemp[1] << ", " << stemp[2] << endl;
   }

   // Prepare the input tree (temporary tree) and observables to read them
   if(DBGSIG > 1)
      cout << "# MvaSetTrees           #: " << "Tree name = " << stemp[0] << endl;
   TTree *tempTree = new TTree;
   tempTree = (TTree*)ifile->Get(stemp[0].c_str());

   Observables *invalues = new Observables(observables);
   Observables *invalues_neg = new Observables(observables);
   Observables *invalues_pos = new Observables(observables);

   // Read S38 and DeltaS38 and Risetime Delta fitting results
   vector<float> *S38fitresults;
   float *deltaS38fitresults;
   ifstream *fitFile;

   S38fitresults = new vector<float>;
   fitFile = new ifstream;
   if(multipleEnergyBins)
      stemp[0] = string(*currentAnalysisDir) + "/../delta_conversion/sectheta_s1000_fits_" + ToString(type) + ".txt";
   else
      stemp[0] = string(*currentAnalysisDir) + "/delta_conversion/sectheta_s1000_fits_" + ToString(type) + ".txt";
   fitFile->open(stemp[0].c_str(), ifstream::in );
   if(fitFile->is_open())
   {
      while(1)
      {
         ftemp[10] = 0.;
         for(int i = 0; i < 10; i++)
            *fitFile >> ftemp[i];
         
         if(fitFile->eof())
            break;

         if(ftemp[10] == ftemp[0])
            break;
         else
         {
            for(int i = 0; i < 10; i++)
               S38fitresults->push_back(ftemp[i]);
            ftemp[10] = ftemp[0];
         }
      }
   }
   fitFile->close();
   delete fitFile;

   deltaS38fitresults = new float[4];
   fitFile = new ifstream;
   if(multipleEnergyBins)
      stemp[0] = string(*currentAnalysisDir) + "/../delta_conversion/s38_data_fit.txt";
   else
      stemp[0] = string(*currentAnalysisDir) + "/delta_conversion/s38_data_fit.txt";
   fitFile->open(stemp[0].c_str(), ifstream::in );
   if(fitFile->is_open())
   {
      *fitFile >> ftemp[0] >> ftemp[1];
      for(int i = 0; i < 4; i++)
         *fitFile >> deltaS38fitresults[i];
   }
   fitFile->close();
   delete fitFile;

   cout << "Fit results (conversion S1000 to S38):" << endl;
   for(int i = 0; i < S38fitresults->size()/10; i++)
   {
      cout << i << ":"; 
      for(int j = 0; j < 10; j++)
         cout << "\t" << S38fitresults->at(10*i+j);
      cout << endl;
   }
   cout << endl;

   cout << "Fit results (conversion S38 to DeltaS38):" << endl;
   for(int i = 0; i < 2; i++)
   {
      cout << i/2 << ":\t" << deltaS38fitresults[2*i] << "\t" << deltaS38fitresults[2*i+1] << endl;
   }
   cout << endl;

   for(int i = 0; i < 3; i++)
   {
      stationDistance[i].clear();
      stationRisetime[i].clear();
   }
   stationHSat.clear();

   bdist = 0;
   bdistneg = 0;
   bdistpos = 0;
   brise = 0;
   briseneg = 0;
   brisepos = 0;
   bsat = 0;

   // Set input values from root files for all observables
   for(int i = 0; i < nrobs; i++)
   {
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Setting mean variable: " << invalues->GetName(i) << endl;
      // For systematics estimation, apply uncertainties to mean value
      tempTree->SetBranchAddress((invalues->GetName(i)).c_str(), &(invalues->obsstruct[i].value));
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Setting neg variable: " << invalues_neg->GetName(i) << endl;
      tempTree->SetBranchAddress((invalues_neg->GetName(i) + "_neg").c_str(), &(invalues_neg->obsstruct[i].value));
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Setting pos variable: " << invalues_pos->GetName(i) << endl;
      tempTree->SetBranchAddress((invalues_pos->GetName(i) + "_pos").c_str(), &(invalues_pos->obsstruct[i].value));
   }

   tempTree->SetBranchAddress("stationdistance", &stationDistance[0], &bdist);
   tempTree->SetBranchAddress("stationdistance_neg", &stationDistance[1], &bdistneg);
   tempTree->SetBranchAddress("stationdistance_pos", &stationDistance[2], &bdistpos);
   tempTree->SetBranchAddress("stationrisetime", &stationRisetime[0], &brise);
   tempTree->SetBranchAddress("stationrisetime_neg", &stationRisetime[1], &briseneg);
   tempTree->SetBranchAddress("stationrisetime_pos", &stationRisetime[2], &brisepos);
   tempTree->SetBranchAddress("stationhighsat", &stationHSat, &bsat);

   itemp[0] = 0;
   itemp[1] = 0;
   itemp[2] = 0;
   itemp[3] = 0;

   for(int j = 0; j < tempTree->GetEntries(); j++)
   {
      // Read all observables
      tempTree->GetEntry(j);

      // Read all vectors for SD stations
      tentry = tempTree->LoadTree(j);
      bdist->GetEntry(tentry);
      bdistneg->GetEntry(tentry);
      bdistpos->GetEntry(tentry);
      brise->GetEntry(tentry);
      briseneg->GetEntry(tentry);
      brisepos->GetEntry(tentry);
      bsat->GetEntry(tentry);

      // Apply any Energy and Xmax corrections to the data file (only on data)
      ret = Find(observables, "xmax");
      if(ret == -1)
      {
         AlertPopup("No xmax observable found", "No xmax observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable xmax.");
         delete tempTree;
         delete[] deltaS38fitresults;
         delete S38fitresults;
         delete[] outobs;
         delete[] outobs_neg;
         delete[] outobs_pos;
         delete invalues;
         delete invalues_neg;
         delete invalues_pos;
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }
      else
      {
         // Apply Xmax corrections to Auger FD standard data
         if((specialMva->widgetChBox[3])->IsChecked())
            invalues->ApplyCorrectionFD();

         // Apply Xmax and energy corrections to Auger HECO data
         if((specialMva->widgetChBox[4])->IsChecked())
         {
            invalues->ApplyCorrectionHECO();
            invalues_neg->ApplyCorrectionHECOErrors(invalues, -1);
            invalues_pos->ApplyCorrectionHECOErrors(invalues, 1);
         }

	 // Apply systematic estimation by adding statistical uncertainties to mean value
         stemp[0] = "TreeS" + ToString(type);
         stemp[1] = string(ifile->GetKey(stemp[0].c_str())->GetTitle());
         stemp[1] = RemovePath(&stemp[1]);
         if( string((dataSelect->widgetCB)->GetStringSelection()) == stemp[1] )	// set to == if changing data and to != if changing simulations
	 {
	    if((specialMva->widgetChBox[5])->IsChecked())
               invalues->ApplyUncertainty(invalues_neg, invalues_pos, 0);

	    if((specialMva->widgetChBox[6])->IsChecked())
               invalues->ApplyUncertainty(invalues_neg, invalues_pos, 1);
	 }
      }

      // Convert S1000 to S38 and then to DeltaS38
      ret = Find(observables, "deltas38");
      if(ret != -1)
      {
         invalues->ConvertToS38(ret, invalues_neg, invalues_pos, S38fitresults);
         invalues->ConvertToDeltaS38(ret, invalues_neg, invalues_pos, deltaS38fitresults);
      }

      // Convert risetimes to Risetime Deltas
      ret = Find(observables, "deltarisetime");
      if(ret != -1)
         invalues->ConvertToDelta(ret, invalues_neg, invalues_pos, S38fitresults);

      // Change zenith angle (theta) to sec(theta) value for MVA analysis
      ret = Find(observables, "zenithSD");
      if(ret != -1)
         invalues->SetupZenith(ret, invalues_neg, invalues_pos);

      ret = Find(observables, "zenithFD");
      if(ret != -1)
         invalues->SetupZenith(ret, invalues_neg, invalues_pos);

      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << "Event = " << j << endl;

      ret = IsInsideCuts(invalues, invalues_neg, invalues_pos, false, 0);

      bool *sdobservable = new bool;

      // Only write out events that are inside the selected cuts
      if(ret == 0)
      {
         itemp[0]++;

         for(int i = 0; i < nrobs; i++)
	 {
            outobs[i] = invalues->GetValue(i);
            outobs_neg[i] = invalues_neg->GetValue(i);
            outobs_pos[i] = invalues_pos->GetValue(i);
	 }
         outtree->Fill();
      }

      delete sdobservable;
   }

   cout << "# MvaSetTrees           #: " << "Number of events inside the cuts = " << itemp[0] << endl;
   cout << "# MvaSetTrees           #: " << "Number of events with multiple eyes inside the cuts = " << itemp[1] << endl;

   ret = itemp[0];

   delete[] deltaS38fitresults;
   delete S38fitresults;

   delete tempTree;
   delete[] outobs;
   delete[] outobs_neg;
   delete[] outobs_pos;
   delete invalues;
   delete invalues_neg;
   delete invalues_pos;
   delete[] stemp;
   delete[] itemp;
   delete[] ftemp;

   return ret;
}

// Set DeltaS38 and Risetime Delta observables from inputs
int MyFrame::SetDeltas(int s38rise, int type, TFile *ifile, bool isdata)
{
   cout << "# SetDeltas             #: " << "Running DeltaS38 and Delta implementation." << endl;

   string *stemp;
   int *itemp;
   float *ftemp;

   stemp = new string[3];
   itemp = new int[4];
   ftemp = new float[9];

   TF1 *fitfuncMid; 
   TGraphAsymmErrors *fitgraphMid;
   TF1 *fitfunc; 
   TGraphAsymmErrors *fitgraph;
   TCanvas *canvtemp;
   TF1 *tempfunc;
   RootStyle *mystyle = new RootStyle();
   mystyle->SetBaseStyle();

   stemp[0] = "TreeS" + ToString(type);
   TTree *tempTree = new TTree;
   tempTree = (TTree*)ifile->Get(stemp[0].c_str());

   if(s38rise == 0)
   {
      cout << "# SetDeltas             #: " << "Setting DeltaS38 (from tree " << stemp[0] << ")" << endl;

      if(multipleEnergyBins)
         stemp[1] = "rm -fr " + string(*currentAnalysisDir) + "/../delta_conversion/sectheta_s1000_fits_" + ToString(type) + ".txt";
      else
         stemp[1] = "rm -fr " + string(*currentAnalysisDir) + "/delta_conversion/sectheta_s1000_fits_" + ToString(type) + ".txt";
      ret = system(stemp[1].c_str());

      if(isdata)
      {
         if(multipleEnergyBins)
            stemp[1] = "rm -fr " + string(*currentAnalysisDir) + "/../delta_conversion/s38_data_fit.txt";
	 else
            stemp[1] = "rm -fr " + string(*currentAnalysisDir) + "/delta_conversion/s38_data_fit.txt";
         ret = system(stemp[1].c_str());
      }

      // Check if all needed observables are present
      ret = Find(observables, "energyFD");
      if(ret == -1)
      {
         AlertPopup("No FD energy observable found", "No FD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energyFD.");
	 delete tempTree;
	 delete mystyle;
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      ret = Find(observables, "zenithFD");
      if(ret == -1)
      {
         AlertPopup("No FD zenith angle observable found", "No FD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithFD.");
	 delete tempTree;
	 delete mystyle;
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      ret = Find(observables, "shwsize");
      if(ret == -1)
      {
         AlertPopup("No S1000 observable found", "No S1000 observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable shwsize.");
	 delete tempTree;
	 delete mystyle;
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      vector<double> *tempvector = new vector<double>;
      vector<double> *energy = new vector<double>;
      vector<double> *zenith = new vector<double>;
      vector<double> *shwsize = new vector<double>;

      Observables *invalues = new Observables(observables);
      Observables *invalues_neg = new Observables(observables);
      Observables *invalues_pos = new Observables(observables);

      int *selectedBin = new int;
      *selectedBin = (cutZenithBins->widgetCB)->GetSelection();

      float *fitparam = new float[4];
      float *fitparamErr = new float[4];

      for(int i = 0; i < nrobs; i++)
      {
         if(DBGSIG > 1)
            cout << "# SetDeltas             #: " << i << ": Setting mean variable: " << invalues->GetName(i) << endl;
         tempTree->SetBranchAddress((invalues->GetName(i)).c_str(), &(invalues->obsstruct[i].value));
         if(DBGSIG > 1)
            cout << "# SetDeltas             #: " << i << ": Setting neg variable: " << invalues_neg->GetName(i) << endl;
         tempTree->SetBranchAddress((invalues_neg->GetName(i) + "_neg").c_str(), &(invalues_neg->obsstruct[i].value));
         if(DBGSIG > 1)
            cout << "# SetDeltas             #: " << i << ": Setting pos variable: " << invalues_pos->GetName(i) << endl;
         tempTree->SetBranchAddress((invalues_pos->GetName(i) + "_pos").c_str(), &(invalues_pos->obsstruct[i].value));
      }

      itemp[0] = 0;
      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

	 if( (TMath::Log10(invalues->GetValue("energyFD")) > 18.5) && (TMath::Log10(invalues->GetValue("energyFD")) <= 20.0) )
	 {
            if( (invalues->GetValue("zenithFD") > InvSecTheta(zcutBins[*selectedBin],false)) && (invalues->GetValue("zenithFD") <= InvSecTheta(zcutBins[*selectedBin+1],false)) )
            {
	       if( (invalues->GetValue("energyFD") != -1) && (invalues->GetValue("zenithFD") != -1) && (invalues->GetValue("shwsize") != -1) )
	       {
	          tempvector->push_back(invalues->GetValue("energyFD"));

	          energy->push_back(invalues->GetValue("energyFD"));
	          energy->push_back(invalues_neg->GetValue("energyFD"));
	          energy->push_back(invalues_pos->GetValue("energyFD"));

	          zenith->push_back(invalues->GetValue("zenithFD"));
	          zenith->push_back(invalues_neg->GetValue("zenithFD"));
	          zenith->push_back(invalues_pos->GetValue("zenithFD"));

	          shwsize->push_back(invalues->GetValue("shwsize"));
	          shwsize->push_back(invalues_neg->GetValue("shwsize"));
	          shwsize->push_back(invalues_pos->GetValue("shwsize"));

                  cout << itemp[0] << " (" << j << "): energyFD = " << invalues->GetValue("energyFD") << ", " << invalues_neg->GetValue("energyFD") << ", " << invalues_pos->GetValue("energyFD") << ", zenithFD = " << invalues->GetValue("zenithFD") << ", " << invalues_neg->GetValue("zenithFD") << ", " << invalues_pos->GetValue("zenithFD") << ", shwsize = " << invalues->GetValue("shwsize") << ", " << invalues_neg->GetValue("shwsize") << ", " << invalues_pos->GetValue("shwsize") << endl;
	          itemp[0]++;
	       }
	    }
	 }
      }

      // Determine energy binning for calculating the S38 fit
      vector<double> *energybins = new vector<double>;
      itemp[0] = 0;
      if(isdata)
      {
         for(int i = 0; i < (19.5-18.5)/0.1; i++)
         {
            energybins->push_back(18.5 + 0.1*i);
            energybins->push_back(18.5 + 0.1*(i+1));
            itemp[0]++;
         }
         energybins->push_back(19.5);
         energybins->push_back(20.0);
         itemp[0]++;
      }
      else
      {
         for(int i = 0; i < (20.0-18.5)/0.1; i++)
         {
            energybins->push_back(18.5 + 0.1*i);
            energybins->push_back(18.5 + 0.1*(i+1));
            itemp[0]++;
         }
      }

      bool *goodrec = new bool;
      vector<double> *energyVect = new vector<double>;
      vector<double> *zenithVect = new vector<double>;
      vector<double> *seczenithVect = new vector<double>;
      vector<double> *shwsizeVect = new vector<double>;
      vector<double> *s38Vect = new vector<double>;

      canvtemp = new TCanvas("canvtemp","",1200,800);

      if(isdata)
      {
         fitgraph = new TGraphAsymmErrors();
         itemp[1] = 0;
      }

      for(int iEn = 0; iEn < itemp[0]; iEn++)
      {
	 // Select the energy bin to use
	 ftemp[0] = energybins->at(2*iEn);
	 ftemp[1] = energybins->at(2*iEn+1);
         cout << "Chosen energy bin = " << ftemp[0] << ", " << ftemp[1] << " (nr events = " << SelectionPass(tempvector, TMath::Power(10,ftemp[0]), TMath::Power(10,ftemp[1])) << ")" << endl;

	 // Clear vectors that will hold values for the selected bin
	 energyVect->clear();
	 zenithVect->clear();
	 seczenithVect->clear();
	 shwsizeVect->clear();
	 s38Vect->clear();

	 // Select only the values for the selected bin
	 for(int i = 0; i < shwsize->size()/3; i++)
	 {
            *goodrec = true;

	    if(TMath::Log10(energy->at(3*i)) <= ftemp[0])
               *goodrec = false;
	    if(TMath::Log10(energy->at(3*i)) > ftemp[1])
               *goodrec = false;

	    if(*goodrec)
	    {
	       energyVect->push_back(energy->at(3*i));
	       energyVect->push_back(energy->at(3*i+1));
	       energyVect->push_back(energy->at(3*i+2));
	       zenithVect->push_back(zenith->at(3*i));
	       zenithVect->push_back(zenith->at(3*i+1));
	       zenithVect->push_back(zenith->at(3*i+2));
	       shwsizeVect->push_back(shwsize->at(3*i));
	       shwsizeVect->push_back(shwsize->at(3*i+1));
	       shwsizeVect->push_back(shwsize->at(3*i+2));
	       // Convert theta to sec(theta):
	       //   sec(theta) = 1/cos(theta)
	       //   dsec(theta) = sin(theta)*dtheta/(cos(theta))^2
	       seczenithVect->push_back(SecTheta(zenith->at(3*i),false));
	       seczenithVect->push_back((TMath::Sin(zenith->at(3*i))*zenith->at(3*i+1))/TMath::Power(TMath::Cos(zenith->at(3*i)),2));
	       seczenithVect->push_back((TMath::Sin(zenith->at(3*i))*zenith->at(3*i+2))/TMath::Power(TMath::Cos(zenith->at(3*i)),2));
	    }
	 }

	 // Fill up a graph with sec(theta) vs S1000 values in order to fit with the fCIC function
	 fitgraphMid = new TGraphAsymmErrors();
         ftemp[2] = 0.;
         for(int i = 0; i < zenithVect->size()/3; i++)
	 {
            fitgraphMid->SetPoint(i, seczenithVect->at(3*i), shwsizeVect->at(3*i));
            fitgraphMid->SetPointError(i, seczenithVect->at(3*i+1), seczenithVect->at(3*i+2), shwsizeVect->at(3*i+1), shwsizeVect->at(3*i+2));

	    ftemp[2] += shwsizeVect->at(3*i);
	 }

	 ftemp[2] = ftemp[2]/(zenithVect->size()/3.);
	 cout << "Mean value of S1000 = " << ftemp[2] << endl;
	 fitgraphMid->Draw("AP");

         fitfuncMid = new TF1("fitfuncMid", "[0]+[1]*(TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2))+[2]*TMath::Power((TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2)),2)+[3]*TMath::Power((TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2)),3)", 1., 2.);
	 fitfuncMid->SetParameters(ftemp[2],ftemp[2]*1.,-ftemp[2]*1.5,-ftemp[2]*1.3);
	 fitgraphMid->Fit("fitfuncMid");

	 for(int i = 0; i < 4; i++)
	 {
            if(i == 0)
	    {
               fitparam[i] = fitfuncMid->GetParameter(i);
               fitparamErr[i] = fitfuncMid->GetParError(i);
	    }
	    else
	    {
               fitparam[i] = (fitfuncMid->GetParameter(i))/fitparam[0];
               fitparamErr[i] = fitfuncMid->GetParError(i);
	       fitparamErr[i] = TMath::Sqrt(TMath::Power((fitfuncMid->GetParError(i))/fitparam[0],2) + TMath::Power(-((fitfuncMid->GetParameter(i))*fitparamErr[0])/TMath::Power(fitparam[0],2),2));
	    }
	 }

	 cout << endl << "Fitting parameters of fCIC are (chi2/ndf = " << fitfuncMid->GetChisquare() << "/" << fitfuncMid->GetNDF() << " = " << (fitfuncMid->GetChisquare())/(fitfuncMid->GetNDF()) << "):" << endl;
	 cout << "- S1000 at 38 deg = " << fitfuncMid->Eval(SecTheta(38.,true)) << endl;
	 cout << "- S = " << fitparam[0] << " ± " << fitparamErr[0] << endl;
	 cout << "- a = " << fitparam[1] << " ± " << fitparamErr[1] << endl;
	 cout << "- b = " << fitparam[2] << " ± " << fitparamErr[2] << endl;
	 cout << "- c = " << fitparam[3] << " ± " << fitparamErr[3] << endl;
	 cout << endl;

	 WriteoutS38Fits(type, 0, ftemp[0], ftemp[1], 4, fitparam, fitparamErr);

	 if(isdata)
	 {
	    cout << "Calculating S38:" << endl;
	    CalculateS38(shwsizeVect, zenithVect, fitparam, fitparamErr, s38Vect);

	    for(int i = 0; i < zenithVect->size()/3; i++)
	    {
               fitgraph->SetPoint(itemp[1], energyVect->at(3*i)/1.e+18, s38Vect->at(3*i));
               fitgraph->SetPointError(itemp[1], energyVect->at(3*i+1)/1.e+18, energyVect->at(3*i+2)/1.e+18, s38Vect->at(3*i+1), s38Vect->at(3*i+2));
               itemp[1]++;
	    }
	 }

	 delete fitfuncMid;
	 delete fitgraphMid;
      }

      if(isdata)
      {
         fitgraph->Draw("AP");
         fitfunc = new TF1("fitfunc", "TMath::Power(x/[0],1./[1])", 2., 100.);
         fitfunc->SetParameters(0.1,1.);
         fitgraph->Fit("fitfunc", "0");
	 tempfunc = (TF1*)fitgraph->GetFunction("fitfunc");

         for(int i = 0; i < 2; i++)
         {
            fitparam[i] = fitfunc->GetParameter(i);
            fitparamErr[i] = fitfunc->GetParError(i);
         }
         
         cout << endl << "Fitting parameters for linear dependence between FD energy and S38:" << endl;
         cout << "- A = " << fitparam[0] << " ± " << fitparamErr[0] << endl;
         cout << "- B = " << fitparam[1] << " ± " << fitparamErr[1] << endl;
         cout << endl;

	 PrintS38Fit(fitgraph, tempfunc, fitparam, fitparamErr, mystyle);
         WriteoutS38Fits(type, 1, 18.5, 20.0, 2, fitparam, fitparamErr);

         delete fitfunc;
         delete fitgraph;
      }

      delete canvtemp;

      delete[] fitparam;
      delete[] fitparamErr;

      delete energyVect;
      delete zenithVect;
      delete seczenithVect;
      delete shwsizeVect;
      delete s38Vect;

      delete goodrec;
      delete selectedBin;
      delete energybins;

      delete invalues;
      delete invalues_neg;
      delete invalues_pos;

      delete tempvector;
      delete energy;
      delete zenith;
      delete shwsize;
   }
   else if(s38rise == 1)
   {
      cout << "# SetDeltas             #: " << "Setting Risetime Delta (from tree " << stemp[0] << ")" << endl;

      // Check if all needed observables are present
      ret = Find(observables, "energySD");
      if(ret == -1)
      {
         AlertPopup("No SD energy observable found", "No SD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energySD.");
	 delete tempTree;
	 delete mystyle;
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      ret = Find(observables, "zenithSD");
      if(ret == -1)
      {
         AlertPopup("No SD zenith angle observable found", "No SD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithSD.");
	 delete tempTree;
	 delete mystyle;
         delete[] stemp;
         delete[] itemp;
         delete[] ftemp;
         return -1;
      }

      Observables *invalues = new Observables(observables);
      Observables *invalues_neg = new Observables(observables);
      Observables *invalues_pos = new Observables(observables);

      bdist = 0;
      bdistneg = 0;
      bdistpos = 0;
      brise = 0;
      briseneg = 0;
      brisepos = 0;
      bsat = 0;

      for(int i = 0; i < 3; i++)
      {
         stationDistance[i].clear();
         stationRisetime[i].clear();
      }
      stationHSat.clear();

      int *selectedBin = new int;
      *selectedBin = (cutZenithBins->widgetCB)->GetSelection();

      float *fitparam = new float[3];
      float *fitparamErr = new float[3];

      for(int i = 0; i < nrobs; i++)
      {
         if(DBGSIG > 1)
            cout << "# SetDeltas             #: " << i << ": Setting mean variable: " << invalues->GetName(i) << endl;
         tempTree->SetBranchAddress((invalues->GetName(i)).c_str(), &(invalues->obsstruct[i].value));
         if(DBGSIG > 1)
            cout << "# SetDeltas             #: " << i << ": Setting neg variable: " << invalues_neg->GetName(i) << endl;
         tempTree->SetBranchAddress((invalues_neg->GetName(i) + "_neg").c_str(), &(invalues_neg->obsstruct[i].value));
         if(DBGSIG > 1)
            cout << "# SetDeltas             #: " << i << ": Setting pos variable: " << invalues_pos->GetName(i) << endl;
         tempTree->SetBranchAddress((invalues_pos->GetName(i) + "_pos").c_str(), &(invalues_pos->obsstruct[i].value));
      }

      tempTree->SetBranchAddress("stationdistance", &stationDistance[0], &bdist);
      tempTree->SetBranchAddress("stationdistance_neg", &stationDistance[1], &bdistneg);
      tempTree->SetBranchAddress("stationdistance_pos", &stationDistance[2], &bdistpos);
      tempTree->SetBranchAddress("stationrisetime", &stationRisetime[0], &brise);
      tempTree->SetBranchAddress("stationrisetime_neg", &stationRisetime[1], &briseneg);
      tempTree->SetBranchAddress("stationrisetime_pos", &stationRisetime[2], &brisepos);
      tempTree->SetBranchAddress("stationhighsat", &stationHSat, &bsat);

      itemp[0] = 0;
      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

         // Read all vectors for SD stations
         tentry = tempTree->LoadTree(j);
         bdist->GetEntry(tentry);
         bdistneg->GetEntry(tentry);
         bdistpos->GetEntry(tentry);
         brise->GetEntry(tentry);
         briseneg->GetEntry(tentry);
         brisepos->GetEntry(tentry);
         bsat->GetEntry(tentry);
      }

      delete[] fitparam;
      delete[] fitparamErr;

      delete invalues;
      delete invalues_neg;
      delete invalues_pos;
   }
   else
   {
      delete tempTree;
      delete mystyle;
      delete[] stemp;
      delete[] itemp;
      delete[] ftemp;
      return -1;
   }

   delete tempTree;
   delete mystyle;
   delete[] stemp;
   delete[] itemp;
   delete[] ftemp;
   return 0;
}

void MyFrame::PrintS38Fit(TGraphAsymmErrors *fitgraph, TF1 *fitfunc, float *fitpar, float *fitparErr, RootStyle *mystyle)
{
   string *stemp = new string[2];

   TLatex *ltext = new TLatex();

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   gStyle->SetEndErrorSize(2);
   
   mystyle->SetGraphColor(fitgraph, 0);
   mystyle->SetAxisTitles(fitgraph, "FD energy (EeV)", "S_{38} (VEM)");
   c1->SetLogx(kTRUE);
   c1->SetLogy(kTRUE);
   fitgraph->SetMarkerSize(0.8);
   fitgraph->GetXaxis()->SetMoreLogLabels(kTRUE);
   fitgraph->GetYaxis()->SetNoExponent(kTRUE);
   fitgraph->Draw("AP");
   mystyle->SetFuncColor(fitfunc, 2);
   fitfunc->Draw("L;SAME");

   ltext->SetTextAlign(11);
   stemp[0] = "A = (" + ToString(fitpar[0]*10.,2) + " #pm " + ToString(fitparErr[0]*10.,2) + ") #times 10^{17} eV";
   ltext->DrawLatex(3., 310., stemp[0].c_str());
   stemp[0] = "B = " + ToString(fitpar[1],3) + " #pm " + ToString(fitparErr[1],3);
   ltext->DrawLatex(3., 0.88*310., stemp[0].c_str());

   if(multipleEnergyBins)
      stemp[0] = "rm -fr " + string(*currentAnalysisDir) + "/../delta_conversion/s38_vs_energyFD_data.pdf";
   else
      stemp[0] = "rm -fr " + string(*currentAnalysisDir) + "/delta_conversion/s38_vs_energyFD_data.pdf";
   system(stemp[0].c_str());
   if(multipleEnergyBins)
      stemp[1] = string(*currentAnalysisDir) + "/../delta_conversion/s38_vs_energyFD_data.pdf";
   else
      stemp[1] = string(*currentAnalysisDir) + "/delta_conversion/s38_vs_energyFD_data.pdf";
   c1->SaveAs(stemp[1].c_str());

   delete ltext;
   delete[] stemp;
   delete c1;
}

void MyFrame::WriteoutS38Fits(int tree, int type, float minEn, float maxEn, int nrpar, float *fitpar, float *fitparErr)
{
   ofstream *fitFile = new ofstream;
   string *stemp = new string[2];

   if(type == 0)
   {
      if(multipleEnergyBins)
         stemp[0] = string(*currentAnalysisDir) + "/../delta_conversion/sectheta_s1000_fits_" + ToString(tree) + ".txt";
      else
         stemp[0] = string(*currentAnalysisDir) + "/delta_conversion/sectheta_s1000_fits_" + ToString(tree) + ".txt";
      fitFile->open(stemp[0].c_str(), ofstream::out | ofstream::app );

      *fitFile << minEn << "\t" << maxEn;
      for(int i = 0; i < nrpar; i++)
         *fitFile << "\t" << fitpar[i] << "\t" << fitparErr[i];
      *fitFile << endl;

      fitFile->close();
   }
   else if(type == 1)
   {
      if(multipleEnergyBins)
         stemp[0] = string(*currentAnalysisDir) + "/../delta_conversion/s38_data_fit.txt";
      else
         stemp[0] = string(*currentAnalysisDir) + "/delta_conversion/s38_data_fit.txt";
      fitFile->open(stemp[0].c_str(), ofstream::out | ofstream::app );

      *fitFile << minEn << "\t" << maxEn;
      for(int i = 0; i < nrpar; i++)
         *fitFile << "\t" << fitpar[i] << "\t" << fitparErr[i];
      *fitFile << endl;

      fitFile->close();
   }

   delete[] stemp;
   delete fitFile;
}

void MyFrame::CalculateS38(vector<double> *shwsize, vector<double> *zenith, float *fitpar, float *fitparErr, vector<double> *outVect)
{
   double *dtemp = new double[4];
   double *seczenith = new double[3];

   for(int i = 0; i < zenith->size()/3; i++)
   {
      // Calculate sec(theta) from zenith angle
      seczenith[0] = SecTheta(zenith->at(3*i),false);
      seczenith[1] = (TMath::Sin(zenith->at(3*i))*zenith->at(3*i+1))/TMath::Power(TMath::Cos(zenith->at(3*i)),2);
      seczenith[2] = (TMath::Sin(zenith->at(3*i))*zenith->at(3*i+2))/TMath::Power(TMath::Cos(zenith->at(3*i)),2);

      cout << i << ": sec(theta) = " << seczenith[0] << ", " << seczenith[1] << ", " << seczenith[2] << ", S1000 = " << shwsize->at(3*i) << ", " << shwsize->at(3*i+1) << ", " << shwsize->at(3*i+2);
      // x = (1/sec(theta))^2 - (cos(thetaref))^2
      dtemp[0] = TMath::Power(1./seczenith[0],2) - TMath::Power(TMath::Cos(DegToRad(38.)),2);
      cout << ", x = " << dtemp[0];
      // fCIC = 1 + a*x + b*x^2 + c*x^3
      dtemp[1] = 1.+fitpar[1]*dtemp[0]+fitpar[2]*TMath::Power(dtemp[0],2)+fitpar[3]*TMath::Power(dtemp[0],3);
      cout << ", fCIC = " << dtemp[1];
      // S38 = S1000/fCIC
      dtemp[2] = shwsize->at(3*i)/dtemp[1];
      cout << ", S38 = " << dtemp[2] << endl;
   
      outVect->push_back(dtemp[2]);
   
      // (dS1000/fCIC)^2
      dtemp[2] = shwsize->at(3*i+1)/dtemp[1];
      dtemp[3] = TMath::Power(dtemp[2],2);
      cout << "  S1000 err = " << dtemp[2];
      // (-S1000*da*x/fCIC^2)^2
      dtemp[2] = -(shwsize->at(3*i)*fitparErr[1]*dtemp[0])/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", a err = " << dtemp[2];
      // (-S1000*db*x^2/fCIC^2)^2
      dtemp[2] = -(shwsize->at(3*i)*fitparErr[2]*TMath::Power(dtemp[0],2))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", b err = " << dtemp[2];
      // (-S1000*dc*x^3/fCIC^2)^2
      dtemp[2] = -(shwsize->at(3*i)*fitparErr[3]*TMath::Power(dtemp[0],3))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", c err = " << dtemp[2];
      // (-S1000*dx*(a+2b*x+3c*x^2)/fCIC^2)^2
      // dx = 2*cos(theta)*sin(theta)*dtheta
      dtemp[2] = 2*seczenith[1]/TMath::Power(seczenith[0],3);
      cout << ", dx = " << dtemp[2];
      dtemp[2] = -(shwsize->at(3*i)*dtemp[2]*(fitpar[1]+2*fitpar[2]*dtemp[0]+3*fitpar[3]*TMath::Power(dtemp[0],2)))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", x err = " << dtemp[2];
   
      dtemp[3] = TMath::Sqrt(dtemp[3]);
      cout << ", dS38 neg = " << dtemp[3] << endl;
      outVect->push_back(dtemp[3]);
   
      // (dS1000/fCIC)^2
      dtemp[2] = shwsize->at(3*i+2)/dtemp[1];
      dtemp[3] = TMath::Power(dtemp[2],2);
      cout << "  S1000 err = " << dtemp[2];
      // (-S1000*da*x/fCIC^2)^2
      dtemp[2] = -(shwsize->at(3*i)*fitparErr[1]*dtemp[0])/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", a err = " << dtemp[2];
      // (-S1000*db*x^2/fCIC^2)^2
      dtemp[2] = -(shwsize->at(3*i)*fitparErr[2]*TMath::Power(dtemp[0],2))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", b err = " << dtemp[2];
      // (-S1000*dc*x^3/fCIC^2)^2
      dtemp[2] = -(shwsize->at(3*i)*fitparErr[3]*TMath::Power(dtemp[0],3))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", c err = " << dtemp[2];
      // (-S1000*dx*(a+2b*x+3c*x^2)/fCIC^2)^2
      // dx = 2*cos(theta)*sin(theta)*dtheta
      dtemp[2] = 2*seczenith[2]/TMath::Power(seczenith[0],3);
      cout << ", dx = " << dtemp[2];
      dtemp[2] = -(shwsize->at(3*i)*dtemp[2]*(fitpar[1]+2*fitpar[2]*dtemp[0]+3*fitpar[3]*TMath::Power(dtemp[0],2)))/TMath::Power(dtemp[1],2);
      dtemp[3] += TMath::Power(dtemp[2],2);
      cout << ", x err = " << dtemp[2];
   
      dtemp[3] = TMath::Sqrt(dtemp[3]);
      cout << ", dS38 pos = " << dtemp[3] << endl;
      outVect->push_back(dtemp[3]);
   }

   delete[] dtemp;
   delete[] seczenith;
}

//int MyFrame::IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, vector<int> *seleye, bool split, int splitbin)
int MyFrame::IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, bool split, int splitbin)
{
   bool *sepcut, *isinside;
   bool *btemp;
   float *ftemp;
   float *ftempaver;
   float *wQuantity, *quantitySum, *wQuantitySum;
   int *avercount;

   sepcut = new bool[3];
   isinside = new bool;
   btemp = new bool[2];
   ftemp = new float[4];
   ftempaver = new float[2];
   wQuantity = new float[2];
   quantitySum = new float[2];
   wQuantitySum = new float[2];
   avercount = new int[2];

   *isinside = false;
   sepcut[0] = false;
   sepcut[1] = false;
   sepcut[2] = false;

   int *selectedBin = new int[2];
   selectedBin[0] = (cutEnergyBins->widgetCB)->GetSelection();
   selectedBin[1] = (cutZenithBins->widgetCB)->GetSelection();

   // Determine the type of observables to cut on (SD or FD)
   if(!split)
      selcuttype = (cutObservables->widgetCB)->GetSelection();
   else
      selcuttype = (splitCutObservables->widgetCB)->GetSelection();

   for(int i = 0; i < 2; i++)
   {
      ftempaver[i] = 0;
      avercount[i] = 0;
      wQuantity[i] = 0;
      quantitySum[i] = 0;
      wQuantitySum[i] = 0;
   }

   if(!split)
   {
      btemp[0] = (cutEnergy->widgetChBox)->IsChecked();
      btemp[1] = (cutZenith->widgetChBox)->IsChecked();
      btemp[2] = (cutRisetime->widgetChBox)->IsChecked();
   }
   else
   {
      btemp[0] = (splitCutEnergy->widgetChBox)->IsChecked();
      btemp[1] = (splitCutZenith->widgetChBox)->IsChecked();
      btemp[2] = (splitCutRisetime->widgetChBox)->IsChecked();
   }

   if(DBGSIG > 1)
      cout << "# IsInsideCuts          #: ";

   // SD observables for cut
   if(selcuttype == 0)
   {
      if(btemp[0])
      {
         ftemp[0] = mean->GetValue("energySD");

         if(DBGSIG > 1)
            cout << "Energy = " << ftemp[0] << "\t";

         if( (ftemp[0] != -1) )
         {
            if(!split)
            {
               if((ftemp[0] > ecutBins[2*selectedBin[0]]) && (ftemp[0] <= ecutBins[2*selectedBin[0]+1]))
                  sepcut[0] = true;
               else
                  sepcut[0] = false;
            }
            else
            {
               if((ftemp[0] > esplitBins[2*splitbin]) && (ftemp[0] <= esplitBins[2*splitbin+1]))
                  sepcut[0] = true;
               else
                  sepcut[0] = false;
            }
         }
         else
            sepcut[0] = false;
      }
      else
         sepcut[0] = true;

      if(btemp[1]) 
      {
         ftemp[0] = mean->GetValue("zenithSD");

         if(DBGSIG > 1)
            cout << "Zenith = " << ftemp[0] << "\t";

         if( (ftemp[0] != -1) )
         {
            if(!split)
            {
//               if((ftemp[0] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
//               if((ftemp[0] > InvSecTheta(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= InvSecTheta(zcutBins[selectedBin[1]+1],false)))
               if((ftemp[0] > zcutBins[selectedBin[1]]) && (ftemp[0] <= zcutBins[selectedBin[1]+1]))
                  sepcut[1] = true;
               else
                  sepcut[1] = false;
            }
            else
            {
               ftemp[2] = DegToRad(splitCutZenith->GetNumber(splitCutZenith->widgetNE[0]));
               ftemp[3] = DegToRad(splitCutZenith->GetNumber(splitCutZenith->widgetNE[1]));
               if((ftemp[0] > ftemp[2]) && (ftemp[0] <= ftemp[3]))
                  sepcut[1] = true;
               else
                  sepcut[1] = false;
            }
         }
         else
            sepcut[1] = false;
      }
      else
         sepcut[1] = true;
   }
   // FD observables for cut
   else if(selcuttype == 1)
   {
      if(btemp[0]) 
      {
         ftemp[0] = mean->GetValue("energyFD");

         if(DBGSIG > 1)
            cout << ftemp[0] << " ";

         if( (ftemp[0] != -1) )
         {
            if(!split)
            {
               if((ftemp[0] > ecutBins[2*selectedBin[0]]) && (ftemp[0] <= ecutBins[2*selectedBin[0]+1]))
                  sepcut[0] = true;
               else
                  sepcut[0] = false;
            }
            else
            {
               if((ftemp[0] > esplitBins[2*splitbin]) && (ftemp[0] <= esplitBins[2*splitbin+1]))
                  sepcut[0] = true;
               else
                  sepcut[0] = false;
            }
         }
         else
            sepcut[0] = false;
      }
      else
         sepcut[0] = true;

      if(btemp[1]) 
      {
         ftemp[0] = mean->GetValue("zenithFD");

         if(DBGSIG > 1)
            cout << ftemp[0] << " ";

         if( (ftemp[0] != -1) )
         {
            if(!split)
            {
//               if((ftemp[0] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
//               if((ftemp[0] > InvSecTheta(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= InvSecTheta(zcutBins[selectedBin[1]+1],false)))
               if((ftemp[0] > zcutBins[selectedBin[1]]) && (ftemp[0] <= zcutBins[selectedBin[1]+1]))
                  sepcut[1] = true;
               else
                  sepcut[1] = false;
            }
            else
            {
               ftemp[2] = DegToRad(splitCutZenith->GetNumber(splitCutZenith->widgetNE[0]));
               ftemp[3] = DegToRad(splitCutZenith->GetNumber(splitCutZenith->widgetNE[1]));
               if((ftemp[0] > ftemp[2]) && (ftemp[0] <= ftemp[3]))
                  sepcut[1] = true;
               else
                  sepcut[1] = false;
            }
         }
         else
            sepcut[1] = false;
      }
      else
         sepcut[1] = true;
   }

   if(btemp[2]) 
   {
      ftemp[0] = mean->GetValue("risetimerecalc");
      ftemp[1] = neg->GetValue("risetimerecalc");

      if(DBGSIG > 1)
         cout << ftemp[1]/ftemp[0] << " ";

      if( ((ftemp[0] != -1) && (ftemp[1] != -1)) )
      {
         if(!split)
         {
            if(ftemp[1]/ftemp[0] <= (cutRisetime->widgetNE[0])->GetValue())
               sepcut[2] = true;
            else
               sepcut[2] = false;
         }
         else
         {
            if(ftemp[1]/ftemp[0] <= (splitCutRisetime->widgetNE[0])->GetValue())
               sepcut[2] = true;
            else
               sepcut[2] = false;
         }
      }
      else
         sepcut[2] = false;
   }
   else
      sepcut[2] = true;

   *isinside = sepcut[0] & sepcut[1] & sepcut[2];

   if(DBGSIG > 1)
      cout << "# IsInsideCuts          #: " << "isinside = " << (int)*isinside << ", sepcut[0] = " << (int)sepcut[0] << ", sepcut[1] = " << (int)sepcut[1] << ", sepcut[2] = " << (int)sepcut[2] << endl;

   delete[] sepcut;
   delete[] btemp;
   delete[] ftemp;
   delete[] ftempaver;
   delete[] avercount;
   delete[] wQuantity;
   delete[] quantitySum;
   delete[] wQuantitySum;
   delete[] selectedBin;

   if(*isinside)
   {
      if(DBGSIG > 1)
         cout << "# IsInsideCuts          #: " << "return code = 0" << endl;
      delete isinside;
      return 0;
   }
   else
   {
      if(DBGSIG > 1)
         cout << "# IsInsideCuts          #: " << "return code = -1" << endl;
      delete isinside;
      return -1;
   }
}

// Perform the MVA analysis on a collection of observables
int MyFrame::PerformMvaAnalysis(string *infilename, string *outfilename, int *curcount)
{
   int *nrTreeEvents;
   string *stemp;
   int *itemp;

   nrTreeEvents = new int[nrkeys];
   stemp = new string[4];
   itemp = new int;

   TFile *ifile = TFile::Open(infilename->c_str(), "READ");
   // Open the file to write out to
   TFile *ofile = TFile::Open(outfilename->c_str(), "RECREATE");
   // Prepare the MVA Factory
   // Factory usage:
   // - user-defined job name, reappearing in names of weight files for training results ("TMVAClassification")
   // - pointer to an output file (ofile)
   // - options
   // Factory has the following options:
   //        V = verbose
   //        Silent = batch mode
   //        Color = colored screen output
   //        DrawProgressBar = progress bar display during training and testing
   //        Transformations = the transformations to make (identity, decorrelation, PCA, uniform, gaussian, gaussian decorrelation)
   //        AnalysisType = setting the analysis type (Classification, Regression, Multiclass, Auto)
   // Default values = !V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Auto
   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",ofile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

   // Selecting the weights directory
   (TMVA::gConfig().GetIONames()).fWeightFileDir = ((*currentAnalysisDir) + "/weights").c_str();
   cout << "# PerformMvaAnalysis    #: " << "Weights directory after = " << (TMVA::gConfig().GetIONames()).fWeightFileDir << endl;
   (*curcount)++;
   progress->Update(*curcount);

   // Adding observables to the Factory
   nrTreeEvents[0] = 0;
   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
      {
         factory->AddVariable(observables[i].c_str(), 'F');
         cout << "# PerformMvaAnalysis    #: " << "Adding variable: " << observables[i] << " (" << (int)obssel[i] << ")" << endl;
         nrTreeEvents[0]++;
      }
   }
   nrselobs = MvaNoteObservables(nrTreeEvents[0]);
   (*curcount)++;
   progress->Update(*curcount);

   // Select signal and background trees (from the temporary input file)
   TTree *signalTree = new TTree;
   TTree *backgroundTree[mixednum];
   for(int j = 0; j < mixednum; j++)
      backgroundTree[j] = new TTree;
   *itemp = 0;

   nrTreeEvents[0] = -1;
   nrTreeEvents[1] = -1;
   TList *tempkeyslist = (TList*) ifile->GetListOfKeys();
   for(int j = 1; j <= ifile->GetNkeys(); j++)
   {
      stemp[0] = string((tempkeyslist->At(j-1))->GetName());
      stemp[1] = string(ifile->GetKey(stemp[0].c_str())->GetTitle());
      stemp[1] = RemovePath(&stemp[1]);

      // Signal tree setup
      if( string((signalSelect->widgetCB)->GetStringSelection()) == stemp[1] )
      {
         cout << "# PerformMvaAnalysis    #: " << "Using signal tree: " << stemp[1] << endl;
         signalTree = (TTree*)ifile->Get(stemp[0].c_str());
         nrTreeEvents[0] = signalTree->GetEntries();
      }

      // Background tree setup
      stemp[3] = string((backgroundSelect->widgetCB)->GetStringSelection());
      if( stemp[3].find(stemp[1]) != string::npos )
      {
         cout << "# PerformMvaAnalysis    #: " << "Using background tree " << *itemp << ": " << stemp[1] << endl;
         backgroundTree[*itemp] = (TTree*)ifile->Get(stemp[0].c_str());
         nrTreeEvents[1] = backgroundTree[*itemp]->GetEntries();
         (*itemp)++;
      }
   }

   (*curcount)++;
   progress->Update(*curcount);

   cout << "# PerformMvaAnalysis    #: " << "Number of entries in signal tree = " << nrTreeEvents[0] << endl;
   cout << "# PerformMvaAnalysis    #: " << "Number of entries in background tree = " << nrTreeEvents[1] << endl;

   // Add signal and background tree
   factory->AddSignalTree(signalTree, 1.0);
   for(int i = 0; i < *itemp; i++)
      factory->AddBackgroundTree(backgroundTree[i], 1.0);

   // Preparing and training from the trees:
   // - preselection cuts make cuts on input variables, before even starting the MVA
   // - options
   // These are the possible options:
   //        nTrain_Signal = number of training events of class Signal (0 takes all)
   //        nTrain_Background = number of training events of class Background (0 takes all)
   //        nTest_Signal = number of test events of class Signal (0 takes all)
   //        nTest_Background = number of test events of class Background (0 takes all)
   //        SplitMode = method of choosing training and testing events (Random, Alternate, Block)
   //        NormMode = renormalisation of event-by-event weights for training (NumEvents: average weight of 1 per event for signal and background, EqualNumEvents: average weight of 1 per event for signal and sum of weights for background equal to sum of weights for signal, None)
   //        V = verbose
   //        MixMode = method of mixing events of different classes into one dataset (SameAsSplitMode, Random, Alternate, Block)
   //        SplitSeed = seed for random event shuffling (default = 100)
   //        VerboseLevel = level of verbose (Debug, Verbose, Info)
   factory->PrepareTrainingAndTestTree("", "", "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
   (*curcount)++;
   progress->Update(*curcount);

   // Booking MVA methods:
   // - type of MVA method to be used (all defined in src/Types.h)
   // - the unique name for the MVA method suplied by the user
   // - options
   // The possible options for each method are defined here: http://tmva.sourceforge.net/optionRef.html
   // For example:
   //        H = print method-specific help message
   //        V = verbose
   //        NeuronType = neuron activation function type (default = sigmoid)
   //        VarTransform = list of variable transformations to do before training (D_Background,P_Signal,G,N_AllClasses -> N = Normalization for all classes)
   //        NCycles = number of training cycles
   //        HiddenLayers = hidden layer architecture (default = N,N-1)
   //        TestRate = test for overtraining at each #th epoch (default = 10)
   //        TrainingMethod = train with back propagation (BP), BFGS algorithm (BFGS) or generic algorithm (GA)
   //        UseRegulator = use regulator to avoid overtraining
   if(BookTheMethod(factory) == -1)
   {
      delete signalTree;
      for(int j = 0; j < mixednum; j++)
         delete backgroundTree[j];
      ifile->Close();
      delete factory;
      ofile->Close();

      AlertPopup("Invalid MVA method", "The selected MVA method is invalid. Please make sure it is correctly defined.");
      delete[] stemp;
      delete[] nrTreeEvents;
      delete itemp;
      return -1;
   }
   (*curcount)++;
   progress->Update(*curcount);

   // Train the selected methods and save them to the weights folder
   factory->TrainAllMethods();
   (*curcount)++;
   progress->Update(*curcount);
   // Test the selected methods by applying the trained data to the test data set -> outputs saved to TestTree output file and then to the output ROOT file
   factory->TestAllMethods();
   (*curcount)++;
   progress->Update(*curcount);
   // Evaluation of methods printed to stdout
   factory->EvaluateAllMethods();
   (*curcount)++;
   progress->Update(*curcount);

   // Close the open files
   delete signalTree;
   for(int j = 0; j < mixednum; j++)
      delete backgroundTree[j];
   ifile->Close();
   delete factory;
   ofile->Close();

   // Copy the training values from results/transformation_stats.dat to the current analysis directory
   stemp[3] = "cp -r " + string(rootdir) + "/results/transformation_stats.dat " + (*currentAnalysisDir) + "/";
   system(stemp[3].c_str());

   // Skip the GUI interface for best cut and automatically select signal/background
   if(((cutEnergyBins->widgetNE[0])->GetValue() > 1) && ((specialMva->widgetChBox[0])->IsChecked()))
   {
      TString *inname = new TString;
      *inname = (TString)(*outfilename);

      MvaEfficiency *effplot = new MvaEfficiency(nrTreeEvents[0], nrTreeEvents[1], currentAnalysisDir);
      effplot->RunMvaEfficiency(inname);
      // Selecting the signal/purity cut (not using signal/background or optimal)
      (cutMva->widgetNE[0])->SetValue(effplot->sigpurCut);

      delete effplot;
      delete inname;
   }
   // Open the MVA GUI to review the training and testing procedure
   else
   {
      if((specialMva->widgetChBox[1])->IsChecked())
      {
         TString *inname = new TString;
         *inname = (TString)(*outfilename);

         MvaEfficiency *effplot = new MvaEfficiency(nrTreeEvents[0], nrTreeEvents[1], currentAnalysisDir);
         effplot->RunMvaEfficiency(inname);

         stemp[2] = "Finished running MVA analysis. For Signal/Background values of: [" + ToString(nrTreeEvents[0]) + "/" + ToString(nrTreeEvents[1]) + "] best cuts are (sig/bgd/pur/SNR):\n";
         stemp[2] = stemp[2] + "- Optimal cut: " + ToString(effplot->optimalCut, 4) + " (" + ToString(100.*(effplot->GetHistValue(0, 0)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(0, 1)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(0, 2)), 2) + "%/" + ToString(effplot->GetHistValue(0, 3), 4) + ")\n";
         stemp[2] = stemp[2] + "- Cut with equal Signal/Background: " + ToString(effplot->sigbgdCut, 4) + " (" + ToString(100.*(effplot->GetHistValue(1, 0)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(1, 1)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(1, 2)), 2) + "%/" + ToString(effplot->GetHistValue(1, 3), 4) + ")\n";
         stemp[2] = stemp[2] + "- Cut with equal Signal/Purity: " + ToString(effplot->sigpurCut, 4) + " (" + ToString(100.*(effplot->GetHistValue(2, 0)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(2, 1)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(2, 2)), 2) + "%/" + ToString(effplot->GetHistValue(2, 3), 4) + ")\n";

         NEDialog *cutMvaDialog = new NEDialog(wxT("MVA cut"), wxSize(600,200), stemp[2], "Set MVA cut:", effplot->sigpurCut, &ID_MVACUTDIALOG);
         cutMvaDialog->SetNEntryFormat(cutMvaDialog->widgetNE, 4, 0.0001, 0, 0., 0.);
         if(cutMvaDialog->ShowModal() == wxID_OK)
            (cutMva->widgetNE[0])->SetValue(cutMvaDialog->GetNEValue());
         else
         {
            delete[] stemp;
            delete[] nrTreeEvents;
            delete itemp;
            delete effplot;
            delete inname;
            cutMvaDialog->Destroy();
            delete cutMvaDialog;
            return -1;
         }

         cutMvaDialog->Destroy();
         delete cutMvaDialog;
         delete effplot;
         delete inname;
      }
   }

   stemp[2] = (methodsSelect->widgetCB)->GetStringSelection();
   applymva = GetMethodName(stemp[2]);

   if( applymva == "All" )
   {
      AlertPopup("Multiple MVA methods", "Multiple MVA methods selected. To continue applying MVA cuts, please select only one method and rerun the analysis.");
      delete[] stemp;
      delete[] nrTreeEvents;
      delete itemp;
      return -1;
   }

   delete[] nrTreeEvents;
   delete[] stemp;
   delete itemp;

   return 0;
}

// Get the TMVA type
void MyFrame::SetTmvaType(TMVA::Factory *factory, int nr, string *formula)
{
   string *stemp;
   stemp = new string[2];

   // Cut optimization
   if(methods[nr].find("Cuts") != string::npos)
      factory->BookMethod(TMVA::Types::kCuts, methods[nr].c_str(), methodsOpt[nr].c_str());
   // 1D likelihood
   if(methods[nr].find("Likelihood") != string::npos)
      factory->BookMethod(TMVA::Types::kLikelihood, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Multidimensional likelihood
   if(methods[nr].find("PDERS") != string::npos)
      factory->BookMethod(TMVA::Types::kPDERS, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Likelihood estimator
   if(methods[nr].find("PDEFoam") != string::npos)
      factory->BookMethod(TMVA::Types::kPDEFoam, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Nearest neighbours
   if(methods[nr].find("KNN") != string::npos)
      factory->BookMethod(TMVA::Types::kKNN, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Linear discriminant
   if(methods[nr].find("LD") != string::npos)
      factory->BookMethod(TMVA::Types::kLD, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Fisher discriminants
   if(methods[nr].find("Fisher") != string::npos)
      factory->BookMethod(TMVA::Types::kFisher, methods[nr].c_str(), methodsOpt[nr].c_str());
   // H-Matrix discriminant
   if(methods[nr].find("HMatrix") != string::npos)
      factory->BookMethod(TMVA::Types::kHMatrix, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Functional discriminant
   if(methods[nr].find("FDA") != string::npos)
   {
      stemp[0] = methodsOpt[nr];
      cout << "Old method options: " << stemp[0] << endl;
      cout << "Formula: " << *formula << endl;
      stemp[0].replace(stemp[0].find("VARFORMULA"), 10, *formula);
      cout << "New method options: " << stemp[0] << endl;
      factory->BookMethod(TMVA::Types::kFDA, methods[nr].c_str(), stemp[0].c_str());
   }
   // Neural networks
   if(methods[nr].find("MLP") != string::npos)
      factory->BookMethod(TMVA::Types::kMLP, methods[nr].c_str(), methodsOpt[nr].c_str());
   if(methods[nr].find("CFMlpANN") != string::npos)
      factory->BookMethod(TMVA::Types::kCFMlpANN, methods[nr].c_str(), methodsOpt[nr].c_str());
   if(methods[nr].find("TMlpANN") != string::npos)
      factory->BookMethod(TMVA::Types::kTMlpANN, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Support vector machine
   if(methods[nr].find("SVM") != string::npos)
      factory->BookMethod(TMVA::Types::kSVM, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Boosted decision trees
   if(methods[nr].find("BDT") != string::npos)
      factory->BookMethod(TMVA::Types::kBDT, methods[nr].c_str(), methodsOpt[nr].c_str());
   // Friedman's rulefit
   if(methods[nr].find("RuleFit") != string::npos)
      factory->BookMethod(TMVA::Types::kRuleFit, methods[nr].c_str(), methodsOpt[nr].c_str());

   delete[] stemp;
}

// Book the method, depending on which was chosen
int MyFrame::BookTheMethod(TMVA::Factory *factory)
{
   string *stemp;
   vector<string>::iterator it;
   int *itemp;

   stemp = new string[3];
   itemp = new int[2];

   stemp[0] = (methodsSelect->widgetCB)->GetStringSelection();
   it = find(methods.begin(), methods.end(), GetMethodName(stemp[0]));
   itemp[0] = distance(methods.begin(), it);
   cout << "# BookTheMethod         #: " << stemp[0] << ", " << methods[itemp[0]] << endl;
   cout << "Method info: " << methods[itemp[0]] << ", " << methodsDesc[itemp[0]] << ", " << methodsOpt[itemp[0]] << endl;

   // Prepare formula for FDA methods
   for(int i = 0; i <= nrselobs; i++)
   {
      if(i == 0)
         stemp[2] = "(" + ToString(i) + ")";
      else
         stemp[2] = stemp[2] + "+(" + ToString(i) + ")*x" + ToString(i-1);
   }
   cout << "# BookTheMethod         #: Using formula " << stemp[2] << endl;
   
   if(GetMethodName(stemp[0]) == "All")
   {
      for(int i = 1; i <= nrmethods; i++)
         SetTmvaType(factory, i, &stemp[2]);
   }
   else
      SetTmvaType(factory, itemp[0], &stemp[2]);

   delete[] stemp;
   delete[] itemp;
   return 0;
}

int MyFrame::GetTrainingShift(string *mvafilename)
{
   string *stemp;
   char *ctemp;
   double *dtemp;

   stemp = new string[2];
   ctemp = new char[1024];
   dtemp = new double;

   TFile *mvafile = TFile::Open(mvafilename->c_str(), "READ");
   GetCorrelations(mvafile);
   mvafile->Close();

   // Check the stats used for the MVA training from file ./results/transformation_stats.dat and save them to vectors
   stemp[0] = (*currentAnalysisDir) + "/transformation_stats.dat";
   ifstream *fstats = new ifstream;
   fstats->open(stemp[0].c_str(), ifstream::in);

   if(!statsMin.empty()) statsMin.erase(statsMin.begin(), statsMin.end());
   if(!statsMax.empty()) statsMax.erase(statsMax.begin(), statsMax.end());

   if(fstats->is_open())
   {
      if(fstats->peek() == '#')
         fstats->getline(ctemp, 1024, '\n');

      for(int i = 0; i < nrobs; i++)
      {
         if(obssel[i])
         {
            if(DBGSIG > 0)
               cout << "# GetTrainingShift      #: " << "Stats for " << observables[i] << ": ";
            *fstats >> *dtemp >> *dtemp;	// Mean and RMS of the testing distribution
            *fstats >> *dtemp;		// Minimum of the testing distribution
            statsMin.push_back(*dtemp);
            if(DBGSIG > 0)
               cout << *dtemp << ", ";
            *fstats >> *dtemp;		// Maximum of the testing distribution
            statsMax.push_back(*dtemp);
            if(DBGSIG > 0)
               cout << *dtemp << endl;
            fstats->ignore(1, '\n');
         }
      }
   }

   fstats->close();
   
   delete fstats;
   delete[] stemp;
   delete[] ctemp;
   delete dtemp;

   return 0;
}

void MyFrame::GetCorrelations(TFile *corfile)
{
   TDirectoryFile *corrDir = new TDirectoryFile();
   corrDir = (TDirectoryFile*)corfile->Get("InputVariables_Id");
   corrDir = (TDirectoryFile*)corrDir->Get("CorrelationPlots");

   string *stemp;
   int *itemp;
   vector<string> *vtemp = new vector<string>;

   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
         vtemp->push_back(observables[i]);
   }

   stemp = new string;
   itemp = new int[2];
   itemp[0] = 0;

   for(int i = 0; i < nrselobs; i++)
   {
      for(int j = i+1; j < nrselobs; j++)
      {
         *stemp = "scat_" + vtemp->at(j) + "_vs_" + vtemp->at(i) + "_Signal_Id";
	 if(DBGSIG > 0)
	    cout << "# GetCorrelations       #: " << "Signal scatter plot: " << *stemp << " (" << ((TH2F*)corrDir->Get( (*stemp).c_str() ))->GetCorrelationFactor() << ")" << endl;
         vtemp->push_back(*stemp);

         *stemp = "scat_" + vtemp->at(j) + "_vs_" + vtemp->at(i) + "_Background_Id";
	 if(DBGSIG > 0)
   	    cout << "# GetCorrelations       #: " << "Background scatter plot: " << *stemp << " (" << ((TH2F*)corrDir->Get( (*stemp).c_str() ))->GetCorrelationFactor() << ")" << endl;
         vtemp->push_back(*stemp);
      }
   }

   vtemp->erase(vtemp->begin(), vtemp->begin()+nrselobs);


   for(int i = 0; i < nrselobs; i++)
   {
      for(int j = 0; j < nrselobs; j++)
      {
         if(i == j)
	 {
	    (*sigCorMat)(i,i) = 1;
	    (*backCorMat)(i,i) = 1;
	 }
	 else if(i < j)
	 {
	    (*sigCorMat)(i,j) = ((TH2F*)corrDir->Get( (vtemp->at(itemp[0])).c_str() ))->GetCorrelationFactor();
	    (*backCorMat)(i,j) = ((TH2F*)corrDir->Get( (vtemp->at(itemp[0]+1)).c_str() ))->GetCorrelationFactor();

	    (*sigCorMat)(j,i) = ((TH2F*)corrDir->Get( (vtemp->at(itemp[0])).c_str() ))->GetCorrelationFactor();
	    (*backCorMat)(j,i) = ((TH2F*)corrDir->Get( (vtemp->at(itemp[0]+1)).c_str() ))->GetCorrelationFactor();

	    itemp[0]+=2;
	 }

      }
   }

   if(DBGSIG > 0)
   {
      cout << "# GetCorrelations       #: " << "Signal correlation matrix: " << endl;

      for(int i = 0; i < nrselobs; i++)
      {
         for(int j = 0; j < nrselobs; j++)
         {
            cout << (*sigCorMat)(i,j) << " | ";
         }
         cout << endl;
      }
      cout << endl;

      cout << "# GetCorrelations       #: " << "Background correlation matrix: " << endl;

      for(int i = 0; i < nrselobs; i++)
      {
         for(int j = 0; j < nrselobs; j++)
         {
            cout << (*backCorMat)(i,j) << " | ";
         }
         cout << endl;
      }
   }

   delete vtemp;
   delete stemp;
   delete[] itemp;
   delete corrDir;
}

void MyFrame::GetApplyCorrelations(string *corname)
{
   string *stemp;
   stemp = new string[2];

   TCanvas *canvCor;
   TTree *TCor = new TTree;
   TGraph *gr1 = new TGraph;

   vector<string> *vtemp = new vector<string>;

   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
         vtemp->push_back(observables[i]);
   }

   TFile *fileCor = TFile::Open(corname->c_str(), "READ");

   canvCor = new TCanvas("canvCor", "", 800, 800);

   TList *tempkeyslist = (TList*)fileCor->GetListOfKeys();
   for(int k = 0; k < nrkeys; k++)
   {
      stemp[0] = string((tempkeyslist->At(k))->GetName());
      TCor = (TTree*)fileCor->Get(stemp[0].c_str());
      if(DBGSIG > 0)
         cout << "# GetApplyCorrelations  #: " << "Calculating correlation factor for tree " << stemp[0] << endl;
      
      for(int i = 0; i < nrselobs; i++)
      {
         for(int j = 0; j < nrselobs; j++)
	 {
            if(i == j)
	       (*otherCorMat[k])(i,i) = 1;
	    else if(i < j)
	    {
               stemp[1] = vtemp->at(i) + ":" + vtemp->at(j);
	       TCor->Draw(stemp[1].c_str());
               gr1 = (TGraph*)canvCor->GetPrimitive("Graph");
	       (*otherCorMat[k])(i,j) = gr1->GetCorrelationFactor();
	       (*otherCorMat[k])(j,i) = gr1->GetCorrelationFactor();
	    }
	 }
      }

      if(DBGSIG > 0)
      {
         for(int i = 0; i < nrselobs; i++)
         {
            for(int j = 0; j < nrselobs; j++)
            {
               cout << (*otherCorMat[k])(i,j) << " | ";
            }
            cout << endl;
         }
         cout << endl;
      }
   }

   fileCor->Close();

   // TEST

   delete vtemp;
//   delete TCor;
   delete gr1;
   delete[] stemp;
   delete canvCor;
   delete fileCor;
}

// Apply the MVA cut and save all events and MVA values into one output files (in order to plot anything from it)
void MyFrame::CreateOutput(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean)
{
   vector<string> *obs = new vector<string>;

   int *sigcount, *backcount;
   sigcount = new int[nrselobs+1];
   backcount = new int[nrselobs+1];

   for(int i = 0; i <= nrselobs; i++)
   {
      sigcount[i] = 0;
      backcount[i] = 0;
   }

   TBranch *mvabranch = new TBranch;
   if(mean == 0)
      mvabranch = app->Branch("MVA", &obsvars[nrobs], "MVA/F");

   if(DBGSIG > 0)
      cout << "# CreateOutput          #: " << "Selected observables:" << endl;
   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
      {
	 obs->push_back(observables[i]);
         cout << "# CreateOutput          #: " << "- " << observables[i] << endl;
      }
   }
   nrselobs = obs->size();

   if(DBGSIG > 1)
      cout << "# CreateOutput          #: " << "Number of events in the tree = " << app->GetEntries() << endl;

   if(!application)
      GetErrors(app, obsvars, obs, curtree);

   // Determine which values are signal and which are background
   for(int ievt = 0; ievt < app->GetEntries(); ievt++)
   {
      app->GetEntry(ievt);

      obsvars[nrobs] = reader->EvaluateMVA(mvamethod.c_str());
      if(mean == 0)
         mvabranch->Fill();

      for(int i = 0; i <= nrselobs; i++)
      {
	 if((i < nrselobs) && (DBGSIG > 1))
            cout << "# CreateOutput          #: " << "  Event = " << ievt << ": values = " << obsvars[3*i] << ", " << obsvars[3*i+1] << ", " << obsvars[3*i+2] << ", observable = " << obs->at(i) << ", mvavalue = " << obsvars[nrobs] << endl;

         if(obsvars[nrobs] >= (cutMva->widgetNE[0])->GetValue())
            sigcount[i]++;
	 else
	    backcount[i]++;
      }
   }

   if(DBGSIG > 0)
   {
      cout << "# CreateOutput          #: " << "Signal vs. background (" << app->GetTitle() << "):" << endl;
      for(int i = 0; i <= nrselobs; i++)
      {
         if(i < nrselobs)
            cout << " - " << obs->at(i) << " = " << sigcount[i] << " vs. " << backcount[i] << endl;
         else
            cout << " - MVA = " << sigcount[i] << " vs. " << backcount[i] << endl;
      }

      cout << endl;
   }

   cutResults.push_back(sigcount[0]);
   cutResults.push_back(backcount[0]);

   mvaprintout = mvaprintout + ToString(app->GetEntries()) + "\t" + ToString(sigcount[0]) + "\t" + ToString(backcount[0]) + "\t" + string(app->GetTitle()) + "\n";

   if(DBGSIG > 0)
      cout << "# CreateOutput          #: " << "Finished here..." << endl;

   delete obs;
   delete[] sigcount;
   delete[] backcount;
}

void MyFrame::GetErrors(TTree *app, float *obsvars, vector<string> *obs, int curtree)
{
   string *stemp;
   TFile *outFile;
   double *dtemp;

   stemp = new string[2];
   dtemp = new double[2];

   if(DBGSIG > 1)
      cout << "# GetErrors             #: " << "Current tree = " << curtree << ", Observables size = " << nrselobs << endl;

   if(curtree == 1)
      outFile = TFile::Open(((*currentAnalysisDir) + "/mva_error.root").c_str(), "RECREATE");
   else
      outFile = TFile::Open(((*currentAnalysisDir) + "/mva_error.root").c_str(), "UPDATE");
   
   vector<float> *errObs = new vector<float>;
   float *errVals = new float[2];

   TTree *outTree = new TTree(("mva_error" + ToString(curtree)).c_str(), "MVA errors");
   outTree->Branch("obs_errors", errObs);
   outTree->Branch("mva_errors", errVals, "mva_errors[2]/F");

   for(int ievt = 0; ievt < app->GetEntries(); ievt++)
   {
      errObs->clear();
     
      app->GetEntry(ievt);

      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "Entry printout " << ievt << ":" << endl;
      for(int j = 0; j < 3*nrselobs; j++)
      {
         if(DBGSIG > 1)
            cout << "# GetErrors             #: " << "j = " << j << ": observable value = " << obsvars[j] << endl;
         errObs->push_back(obsvars[j]);
      }

      double *norm1, *norm2;
      norm1 = new double[3];
      norm2 = new double[3];
      covMatNeg = new TMatrixD(nrselobs,nrselobs);
      covMatPos = new TMatrixD(nrselobs,nrselobs);
      eigenValMatNeg = new TMatrixD(nrselobs,nrselobs);
      eigenValMatPos = new TMatrixD(nrselobs,nrselobs);
     
      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "Covariance matrix:" << endl;

      for(int i = 0; i < nrselobs; i++)
      {
	 norm1[0] = ((obsvars[3*i] - statsMin[i])/(statsMax[i] - statsMin[i]))*2 - 1;
	 norm1[1] = ((obsvars[3*i] - obsvars[3*i+1] - statsMin[i])/(statsMax[i] - statsMin[i]))*2 - 1;
	 norm1[2] = ((obsvars[3*i] + obsvars[3*i+2] - statsMin[i])/(statsMax[i] - statsMin[i]))*2 - 1;
         if(DBGSIG > 1)
   	    cout << "# GetErrors             #: " << "Normalized values (" << obs->at(i) << "): " << norm1[0] << ", " << norm1[1] << ", " << norm1[2] << endl;

         for(int j = 0; j < nrselobs; j++)
         {
	    norm2[0] = ((obsvars[3*j] - statsMin[j])/(statsMax[j] - statsMin[j]))*2 - 1;
	    norm2[1] = ((obsvars[3*j] - obsvars[3*j+1] - statsMin[j])/(statsMax[j] - statsMin[j]))*2 - 1;
	    norm2[2] = ((obsvars[3*j] + obsvars[3*j+2] - statsMin[j])/(statsMax[j] - statsMin[j]))*2 - 1;
            if(DBGSIG > 1)
	    {
	       cout << "# GetErrors             #: " << "Xmin and Xmax (" << obs->at(j) << "): " << statsMin[j] << ", " << statsMax[j] << endl;
	       cout << "# GetErrors             #: " << "Original values (" << obs->at(j) << "): " << obsvars[3*j] << ", " << obsvars[3*j] - obsvars[3*j+1] << ", " << obsvars[3*j] + obsvars[3*j+2] << endl;
	       cout << "# GetErrors             #: " << "Normalized values (" << obs->at(j) << "): " << norm2[0] << ", " << norm2[1] << ", " << norm2[2] << endl;
	    }

	    if((signalSelect->widgetCB)->GetSelection() == curtree-1)
	    {
               (*covMatNeg)(i,j) = (*sigCorMat)(i,j)*(norm1[0]-norm1[1])*(norm2[0]-norm2[1]);
               (*covMatPos)(i,j) = (*sigCorMat)(i,j)*(norm1[2]-norm1[0])*(norm2[2]-norm2[0]);
               if(DBGSIG > 1)
	          cout << "# GetErrors             #: " << "Signal (neg: " << (*sigCorMat)(i,j) << ", " << (norm1[0]-norm1[1]) << ", " << (norm2[0]-norm2[1]) << " = " << (*covMatNeg)(i,j) << ")" << endl;
	    }
	    else if((backgroundSelect->widgetCB)->GetSelection() == curtree-1)
	    {
               (*covMatNeg)(i,j) = (*backCorMat)(i,j)*(norm1[0]-norm1[1])*(norm2[0]-norm2[1]);
               (*covMatPos)(i,j) = (*backCorMat)(i,j)*(norm1[2]-norm1[0])*(norm2[2]-norm2[0]);
               if(DBGSIG > 1)
	          cout << "# GetErrors             #: " << "Background (neg: " << (*backCorMat)(i,j) << ", " << (norm1[0]-norm1[1]) << ", " << (norm2[0]-norm2[1]) << " = " << (*covMatNeg)(i,j) << ")" << endl;
	    }
	    else
	    {
               (*covMatNeg)(i,j) = (*otherCorMat[curtree-1])(i,j)*(norm1[0]-norm1[1])*(norm2[0]-norm2[1]);
               (*covMatPos)(i,j) = (*otherCorMat[curtree-1])(i,j)*(norm1[2]-norm1[0])*(norm2[2]-norm2[0]);
               if(DBGSIG > 1)
	          cout << "# GetErrors             #: " << "Others (neg: " << (*otherCorMat[curtree-1])(i,j) << ", " << (norm1[0]-norm1[1]) << ", " << (norm2[0]-norm2[1]) << " = " << (*covMatNeg)(i,j) << ")" << endl;
	    }

            if(DBGSIG > 1)
               cout << "# GetErrors             #: " << " | " << (*covMatNeg)(i,j) << " | " << endl;

	    eigenCovMat = new TMatrixDEigen((const TMatrixD)*covMatNeg);
	    (*eigenValMatNeg) = eigenCovMat->GetEigenValues();
	    delete eigenCovMat;
	    eigenCovMat = new TMatrixDEigen((const TMatrixD)*covMatPos);
	    (*eigenValMatPos) = eigenCovMat->GetEigenValues();
	    delete eigenCovMat;
         }

         if(DBGSIG > 1)
	    cout << endl;
      }

      dtemp[0] = 0;

      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "Diagonalized matrix (negative error): ";
      for(int i = 0; i < nrselobs; i++)
      {
	 dtemp[0] += TMath::Power((*eigenValMatNeg)(i,i),2);
         if(DBGSIG > 1)
            cout << (*eigenValMatNeg)(i,i) << " ";
      }
      
      if(DBGSIG > 1)
         cout << endl;

      dtemp[1] = TMath::Sqrt(dtemp[0]);
      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "MVA variable error (negative error) = " << dtemp[1] << endl;

      errVals[0] = dtemp[1];

      dtemp[0] = 0;

      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "Diagonalized matrix (positive error): ";
      for(int i = 0; i < nrselobs; i++)
      {
	 dtemp[0] += TMath::Power((*eigenValMatPos)(i,i),2);
         if(DBGSIG > 1)
            cout << (*eigenValMatPos)(i,i) << " ";
      }
      
      if(DBGSIG > 1)
         cout << endl;

      dtemp[1] = TMath::Sqrt(dtemp[0]);
      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "MVA variable error (positive error) = " << dtemp[1] << endl << endl;

      errVals[1] = dtemp[1];
      outTree->Fill();

      delete covMatNeg;
      delete covMatPos;
      delete eigenValMatNeg;
      delete eigenValMatPos;
      delete[] norm2;
      delete[] norm1;
   }

   outFile->Write();
   delete outTree;
   outFile->Close();

   delete errObs;
   delete[] errVals;
   delete[] stemp;
   delete[] dtemp;
}

void MyFrame::GetMvaError(int selection, double *outvalue, string *inname)
{
   string *stemp;

   if(!(specialMva->widgetChBox[2])->IsChecked())
   {
      if(DBGSIG > 1)
      {
         cout << "# GetMvaError           #: " << "Selection: " << selection << endl;
         cout << "# GetMvaError           #: " << "Mean: " << outvalue[0] << endl;
      }

      stemp = new string;

      *stemp = (*currentAnalysisDir) + string("/mva_error.root");
      if(DBGSIG > 0)
         cout << "# GetMvaError           #: " << "File: " << *stemp << endl;

      int *vrstica = new int;
      *vrstica = 0;

      vector<double> *negarray = new vector<double>;
      vector<double> *posarray = new vector<double>;
      
      float *errVals = new float[2];

      TFile *inFile = TFile::Open(stemp->c_str(), "READ");
      TTree *inTree = new TTree;
      inTree = (TTree*)inFile->Get(("mva_error" + ToString(selection)).c_str());
      inTree->SetBranchAddress("mva_errors", errVals);

      for(int i = 0; i < inTree->GetEntries(); i++)
      {
         inTree->GetEntry(i);
         negarray->push_back(errVals[0]);
         posarray->push_back(errVals[1]);
         (*vrstica)++;
      }

      if(DBGSIG > 1)
         cout << "# GetMvaError           #: " << "Lines read = " << *vrstica << endl;

      double *dtemp;
      dtemp = new double[2];
      dtemp[0] = 0;
      dtemp[1] = 0;
      for(int i = 0; i < *vrstica; i++)
      {
         dtemp[0] += negarray->at(i);
         dtemp[1] += posarray->at(i);
      }
      dtemp[0] = dtemp[0]/(*vrstica);
      dtemp[1] = dtemp[1]/(*vrstica);

      if(DBGSIG > 1)
         cout << "# GetMvaError           #: " << "Mean values (negative error, positive error) = " << "(" << dtemp[0] << ", " << dtemp[1] << ")" << endl;

      for(int i = 0; i < *vrstica; i++)
      {
         outvalue[1] += (negarray->at(i) - dtemp[0])*(negarray->at(i) - dtemp[0]);
         outvalue[2] += (posarray->at(i) - dtemp[1])*(posarray->at(i) - dtemp[1]);
      }
      outvalue[1] = TMath::Sqrt(outvalue[1]/(*vrstica));
      outvalue[2] = TMath::Sqrt(outvalue[2]/(*vrstica));

      if(DBGSIG > 1)
         cout << "# GetMvaError           #: " << "Sigma values (negative error, positive error) = " << "(" << outvalue[1] << ", " << outvalue[2] << ")" << endl;

      outvalue[1] = outvalue[0] - outvalue[1];
      outvalue[2] = outvalue[0] + outvalue[2];

      inFile->Close();

      delete[] errVals;
//      delete inTree;
      delete negarray;
      delete posarray;
      delete vrstica;
      delete[] dtemp;
      delete stemp;
   }
   else
   {
      double *dtemp = new double[5];
      stemp = new string[3];

      TFile *ifile = TFile::Open(inname->c_str(), "READ");

      vector<string> *obser = new vector<string>;
      cout << "# GetMvaError           #: " << "Observables: ";
      for(int i = 0; i < nrobs; i++)
      {
         if(obssel[i])
         {
            obser->push_back(observables[i]);
            cout << observables[i] << ", ";
         }
      }
      cout << endl;
      nrselobs = obser->size();

      float *obsVals = new float[3*nrselobs];

      TTree *dataTree = new TTree;

      TList *tempkeyslist = (TList*)ifile->GetListOfKeys();
      for(int j = 1; j <= ifile->GetNkeys(); j++)
      {
         stemp[0] = string((tempkeyslist->At(j-1))->GetName());
         stemp[1] = string(ifile->GetKey(stemp[0].c_str())->GetTitle());
         stemp[1] = RemovePath(&stemp[1]);

         // Signal tree setup
         if( string((dataSelect->widgetCB)->GetStringSelection()) == stemp[1] )
         {
            cout << "# GetMvaError           #: " << "Using data tree: " << stemp[1] << endl;
            dataTree = (TTree*)ifile->Get(stemp[0].c_str());
         }
      }

      for(int i = 0; i < nrselobs; i++)
      {
         cout << "- " << obser->at(i) << endl;
         stemp[0] = obser->at(i);
         dataTree->SetBranchAddress(stemp[0].c_str(), &obsVals[3*i]);
         stemp[0] = obser->at(i) + "_neg";
         dataTree->SetBranchAddress(stemp[0].c_str(), &obsVals[3*i+1]);
         stemp[0] = obser->at(i) + "_pos";
         dataTree->SetBranchAddress(stemp[0].c_str(), &obsVals[3*i+2]);
      }

      TH1F *hist, *histNeg, *histPos;
      dtemp[3] = 0;
      dtemp[4] = 0;

      for(int i = 0; i < nrselobs; i++)
      {
         hist = new TH1F("h1","",100,-1.5,1.5);
         histNeg = new TH1F("h2","",100,0.,1.);
         histPos = new TH1F("h3","",100,0.,1.);

         for(int j = 0; j < dataTree->GetEntries(); j++)
         {
            dataTree->GetEntry(j);

            dtemp[0] = 2*(obsVals[3*i] - statsMin[i])/(statsMax[i] - statsMin[i]) - 1;
            dtemp[1] = dtemp[0] - (2*((obsVals[3*i]-obsVals[3*i+1]) - statsMin[i])/(statsMax[i] - statsMin[i]) - 1);
            dtemp[2] = 2*((obsVals[3*i]+obsVals[3*i+2]) - statsMin[i])/(statsMax[i] - statsMin[i]) - 1 - dtemp[0];

            hist->Fill(dtemp[0]);
            histNeg->Fill(dtemp[1]);
            histPos->Fill(dtemp[2]);
         }

         cout << "# GetMvaError           #: " << "0: Mean = " << hist->GetMean() << ", Sigma = " << hist->GetRMS() << endl;
         cout << "# GetMvaError           #: " << "-: Mean = " << histNeg->GetMean() << ", Sigma = " << histNeg->GetRMS() << endl;
         cout << "# GetMvaError           #: " << "+: Mean = " << histPos->GetMean() << ", Sigma = " << histPos->GetRMS() << endl;

         dtemp[3] += TMath::Power(histNeg->GetRMS(), 2);
         dtemp[4] += TMath::Power(histPos->GetRMS(), 2);

         delete hist;
         delete histNeg;
         delete histPos;
      }

      outvalue[1] = TMath::Sqrt(dtemp[3]);
      outvalue[2] = TMath::Sqrt(dtemp[4]);
      cout << "# GetMvaError           #: " << "Final combined negative error = " << outvalue[1] << endl;
      cout << "# GetMvaError           #: " << "Final combined positive error = " << outvalue[2] << endl;

      outvalue[1] = outvalue[0] - outvalue[1];
      outvalue[2] = outvalue[0] + outvalue[2];

      ifile->Close();

//      delete dataTree;
      delete obsVals;
      delete obser;
      delete[] stemp;
      delete[] dtemp;
   }
}
