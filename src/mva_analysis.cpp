#include "frame.h"
#include "separate_functions.h"
#include "adst_mva.h"
#include "popups.h"
#include "combine.h"

// Add observables to the MVA analysis
//int MyFrame::MvaAddObservables(TMVA::Factory *factory)
int MyFrame::MvaNoteObservables(int count)
{
//   int count = 0;
   string *stemp;

   // Write out the number of different selected observables
   stemp = new string;
   if(DBGSIG > 0)
      cout << "# MvaNoteObservables    #: " << "Number of different selected observables: " << count << endl;
   *stemp = string(rootdir) + "/results/observables_nr.dat";
   ofstream fobser;
   fobser.open(stemp->c_str(), ofstream::out | ofstream::trunc);
   if(fobser.is_open())
      fobser << count << endl;
   fobser.close();

   delete stemp;
   return count;
}

// Write out the MVA trees to contain only one of the eyes (best eye inside selection cuts or average of all eyes)
int MyFrame::MvaTreeFile(string *infilename, string *outfilename, int *nrEvents)
{
   string *stemp;
   stemp = new string[2];

   TFile *ifile = TFile::Open(infilename->c_str(), "READ");

   TTree *signalTempTree[nrkeys];

   int selectedBin[2];
   selectedBin[0] = (cutEnergyBins->widgetCB)->GetSelection();
   selectedBin[1] = (cutZenithBins->widgetCB)->GetSelection();

   if( (signalSelect->widgetCB)->GetSelection() == (backgroundSelect->widgetCB)->GetSelection() )
   {
      delete[] stemp;
      ifile->Close();
      return -1;
   }
   
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
   stemp = new string[3];
   vector<int> seleye;
   int *itemp;
   itemp = new int[4];
   float *ftemp;
   ftemp = new float[9];

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
//   outtree = new TTree(stemp[0].c_str(), treeName->c_str());
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
   TTree *tempTree = (TTree*)ifile->Get(stemp[0].c_str());

   Observables *invalues = new Observables(observables);
   Observables *invalues_neg = new Observables(observables);
   Observables *invalues_pos = new Observables(observables);

   for(int i = 0; i < nrobs; i++)
   {
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Setting mean variable: " << invalues->GetName(i) << endl;
      tempTree->SetBranchAddress((invalues->GetName(i)).c_str(), invalues->obsstruct[i].value);
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Setting neg variable: " << invalues_neg->GetName(i) << endl;
      tempTree->SetBranchAddress((invalues_neg->GetName(i) + "_neg").c_str(), invalues_neg->obsstruct[i].value);
      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << i << ": Setting pos variable: " << invalues_pos->GetName(i) << endl;
      tempTree->SetBranchAddress((invalues_pos->GetName(i) + "_pos").c_str(), invalues_pos->obsstruct[i].value);
   }

   itemp[0] = 0;
   itemp[1] = 0;
   itemp[2] = 0;
   itemp[3] = 0;

   for(int j = 0; j < tempTree->GetEntries(); j++)
   {
      tempTree->GetEntry(j);

      // TODO: Check if event is inside the selected cuts
      if(!seleye.empty()) seleye.erase(seleye.begin(), seleye.end());

      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << "Event = " << j << endl;
      ret = IsInsideCuts(invalues, invalues_neg, invalues_pos, &seleye);

      // Combine stereo FD events
      if(ret == 0)
      {
         itemp[0]++;

         for(int i = 0; i < nrobs; i++)
	 {
            ftemp[0] = 0;
            ftemp[1] = 0;
            ftemp[2] = 0;
            ftemp[3] = 0;
            ftemp[4] = 0;
            ftemp[5] = 0;
            ftemp[6] = 0;
	    outobs[i] = 0;

	    itemp[2] = 0;
	    itemp[3] = 0;

            for(int k = 0; k < seleye.size(); k++)
            {
	       itemp[2]++;
	       if(invalues->GetValue(i, k) != -1)
	       {
	          ftemp[6] = invalues->GetValue(i, k);
	          ftemp[7] = invalues_neg->GetValue(i, k);
	          ftemp[8] = invalues_pos->GetValue(i, k);

	          itemp[3]++;
		  if(invalues_neg->GetValue(i, k) == 0)
		  {
                     ftemp[0] = invalues->GetValue(i, k);
		     ftemp[2] = 1.;
		  }

		  if(invalues_pos->GetValue(i, k) == 0)
		  {
                     ftemp[1] = invalues->GetValue(i, k);
		     ftemp[3] = 1.;
		  }

		  if( (invalues_neg->GetValue(i, k) != 0) && (invalues_pos->GetValue(i, k) != 0) )
		  {
                     if(DBGSIG > 1)
	                cout << "# MvaSetTrees           #: " << "Input values (combined): " << ftemp[6] << ", " << ftemp[7] << ", " << ftemp[8] << endl;
                     ftemp[0] += ftemp[6]/TMath::Power(ftemp[7], 2);
                     if(DBGSIG > 1)
	                cout << "# MvaSetTrees           #: " << "sum of ftemp[0] = " << ftemp[0] << endl;
                     ftemp[1] += ftemp[6]/TMath::Power(ftemp[8], 2);
                     if(DBGSIG > 1)
	                cout << "# MvaSetTrees           #: " << "sum of ftemp[1] = " << ftemp[1] << endl;
                     ftemp[2] += 1./TMath::Power(ftemp[7], 2);
                     if(DBGSIG > 1)
	                cout << "# MvaSetTrees           #: " << "sum of ftemp[2] = " << ftemp[2] << endl;
                     ftemp[3] += 1./TMath::Power(ftemp[8], 2);
                     if(DBGSIG > 1)
	                cout << "# MvaSetTrees           #: " << "sum of ftemp[3] = " << ftemp[3] << endl;
		  }
	       }
	    }

            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Number of eyes = " << itemp[2] << ", Number of valid eyes = " << itemp[3] << endl;

	    ftemp[4] = ftemp[0]/ftemp[2];
            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "val of ftemp[4] = " << ftemp[4] << endl;
	    ftemp[5] = ftemp[1]/ftemp[3];
            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "val of ftemp[5] = " << ftemp[5] << endl;
	    outobs[i] = (ftemp[4] + ftemp[5])/2.;

	    if(ftemp[2] == 1.)	// negative error values were 0
	       outobs_neg[i] = 0;
	    else		// calculate the error value
	       outobs_neg[i] = TMath::Sqrt(1./ftemp[2]);

	    if(ftemp[3] == 1.)	// positive error values were 0
	       outobs_pos[i] = 0;
	    else		// calculate the error value
	       outobs_pos[i] = TMath::Sqrt(1./ftemp[3]);

	    // if all 4 values are the same, the errors should also be the same (for instance when going through 4 saved values for SD variables)
	    if( ((itemp[2] == ALLEYES) && (itemp[3] == ALLEYES)) || (itemp[3] == 1) )
            {
	       outobs[i] = ftemp[6];
	       outobs_neg[i] = ftemp[7];
	       outobs_pos[i] = ftemp[8];
	    }

	    // if all 4 values were negative (invalid variable)
	    if(itemp[3] == 0)
	    {
	       outobs[i] = -1;
	       outobs_neg[i] = -1;
	       outobs_pos[i] = -1;
	    }

	    if(DBGSIG > 1)
               cout << "# MvaSetTrees           #: " << "Writeout values (combined), " << observables[i] << ": " << outobs[i] << ", " << outobs_neg[i] << ", " << outobs_pos[i] << endl;
	 }
         outtree->Fill();

/*         if(DBGSIG > 1)
	    cout << "# MvaSetTrees           #: " << "Event " << j << " IsInsideCuts (" << seleye.size() << ")" << endl;
         itemp[2] = seleye[0];

         itemp[0]++;
         if(seleye.size() > 1)
         {
            itemp[1]++;

            if(DBGSIG > 1)
	    {
               cout << "# MvaSetTrees           #: " << "Event " << j << " has multiple eyes inside the cuts = ";
               for(int i = 0; i < seleye.size(); i++)
                  cout << seleye[i] << " ";
               cout << endl;
             }

            // Determine which FD eye has the smallest error on xmax
            ftemp[3] = 100.;

	    for(int k = 0; k < seleye.size(); k++)
	    {
               // Get error on xmax value (relative)
               ftemp[0] = invalues->GetValue(generalObservables->GetInt("xmax"), seleye[k]);
               ftemp[1] = invalues_neg->GetValue(generalObservables->GetInt("xmax"), seleye[k]);

               if(DBGSIG > 1)
	          cout << "# MvaSetTrees           #: " << "Values: " << ftemp[0] << ", " << ftemp[1] << endl;

	       ftemp[0] = ftemp[1]/ftemp[0];

               if(DBGSIG > 1)
	       {
   	          cout << "# MvaSetTrees           #: " << "Old relative error: " << ftemp[3] << endl;
   	          cout << "# MvaSetTrees           #: " << "New relative error: " << ftemp[0] << endl;
   	       }
	       
	       // If error is smaller, save them
	       if(ftemp[0] < ftemp[3])
	       {
                  ftemp[3] = ftemp[0];
		  itemp[2] = seleye[k];
	       }
	    }
	    
            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Best eye is: " << itemp[2] << endl;
	 }

         for(int i = 0; i < nrobs; i++)
	 {
            outobs[i] = invalues->GetValue(i, itemp[2]);
            outobs_neg[i] = invalues_neg->GetValue(i, itemp[2]);
            outobs_pos[i] = invalues_pos->GetValue(i, itemp[2]);
            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Saving data " << i << ": " << outobs[i] << " (" << outobs_neg[i] << ", " << outobs_pos[i] << ")" << endl;
	 }
	 outtree->Fill();*/
      }
      // Any of the eyes is inside the cut
      else if(ret == 1)
      {
         if(DBGSIG > 1)
	    cout << "# MvaSetTrees           #: " << "Event " << j << " IsInsideCuts (" << seleye.size() << ")" << endl;
         itemp[2] = seleye[0];

         itemp[0]++;
         if(seleye.size() > 1)
         {
            itemp[1]++;

            if(DBGSIG > 1)
	    {
               cout << "# MvaSetTrees           #: " << "Event " << j << " has multiple eyes inside the cuts = ";
               for(int i = 0; i < seleye.size(); i++)
                  cout << seleye[i] << " ";
               cout << endl;
            }

            // Determine which FD eye has the smallest error on xmax
            ftemp[3] = 100.;

	    for(int k = 0; k < seleye.size(); k++)
	    {
               // Get error on xmax value (relative)
               ftemp[0] = invalues->GetValue(generalObservables->GetInt("xmax"), seleye[k]);
               ftemp[1] = invalues_neg->GetValue(generalObservables->GetInt("xmax"), seleye[k]);

            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Values: " << ftemp[0] << ", " << ftemp[1] << endl;

	       ftemp[0] = ftemp[1]/ftemp[0];

            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Old relative error: " << ftemp[3] << endl;
	       cout << "# MvaSetTrees           #: " << "New relative error: " << ftemp[0] << endl;
	       
	       // If error is smaller, save them
	       if(ftemp[0] < ftemp[3])
	       {
                  ftemp[3] = ftemp[0];
		  itemp[2] = seleye[k];
	       }
	    }
	    
            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Best eye is: " << itemp[2] << endl;
	 }

         for(int i = 0; i < nrobs; i++)
	 {
            outobs[i] = invalues->GetValue(i, itemp[2]);
            outobs_neg[i] = invalues_neg->GetValue(i, itemp[2]);
            outobs_pos[i] = invalues_pos->GetValue(i, itemp[2]);
            if(DBGSIG > 1)
	       cout << "# MvaSetTrees           #: " << "Saving data " << i << ": " << outobs[i] << " (" << outobs_neg[i] << ", " << outobs_pos[i] << ")" << endl;
	 }
	 outtree->Fill();
      }
      // Average of eyes is inside the cut
      else if(ret == 2)
      {
         itemp[0]++;
         for(int i = 0; i < nrobs; i++)
	 {
	    outobs[i] = 0;
            for(int k = 0; k < seleye.size(); k++)
            {
	       if(invalues->GetValue(i, k) != -1)
	       {
                  outobs[i] += invalues->GetValue(i, k);
                  outobs_neg[i] += invalues_neg->GetValue(i, k);
                  outobs_pos[i] += invalues_pos->GetValue(i, k);
/*                  outobs[i] += invalues->obsstruct[i].value[k];
                  outobs_neg[i] += invalues->obsstruct[i].value[k];
                  outobs_pos[i] += invalues->obsstruct[i].value[k];*/
	       }
	    }

	    outobs[i] = outobs[i]/(seleye.size());
	    outobs_neg[i] = outobs_neg[i]/(seleye.size());
	    outobs_pos[i] = outobs_pos[i]/(seleye.size());

	    if(DBGSIG > 0)
               cout << "# MvaSetTrees           #: " << "Writeout values (average): " << outobs[i] << ", " << outobs_neg[i] << ", " << outobs_pos[i] << endl;
	 }
         outtree->Fill();
      }
   }

   cout << "# MvaSetTrees           #: " << "Number of events inside the cuts = " << itemp[0] << endl;
   cout << "# MvaSetTrees           #: " << "Number of events with multiple eyes inside the cuts = " << itemp[1] << endl;

   ret = itemp[0];

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

int MyFrame::IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, vector<int> *seleye)
{
   bool *sepcut, isinside;
   sepcut = new bool[3];
   float *ftemp;
   ftemp = new float[2];
   float *ftempaver;
   ftempaver = new float[2];
   float *wQuantity, *quantitySum, *wQuantitySum;
   wQuantity = new float[2];
   quantitySum = new float[2];
   wQuantitySum = new float[2];
   int *avercount;
   avercount = new int[2];
   isinside = false;
   sepcut[0] = false;
   sepcut[1] = false;
   sepcut[2] = false;

   int selectedBin[2];
   selectedBin[0] = (cutEnergyBins->widgetCB)->GetSelection();
   selectedBin[1] = (cutZenithBins->widgetCB)->GetSelection();

   // Determine the type of observables to cut on (SD or FD)
   selcuttype = (cutObservables->widgetCB)->GetSelection();
   // Determine how eye selection should be handled (any eye inside selection or average)
   seleyetype = (eyeSelection->widgetCB)->GetSelection();

   for(int i = 0; i < 2; i++)
   {
      ftempaver[i] = 0;
      avercount[i] = 0;
      wQuantity[i] = 0;
      quantitySum[i] = 0;
      wQuantitySum[i] = 0;
   }

   // Calculate averages for all eyes
   for(int i = 0; i < ALLEYES; i++)
   {
      if(selcuttype == 0)
      {
         if((cutEnergy->widgetChBox)->IsChecked())
         {
            if(mean->GetValue(mean->GetInt("energySD"), i) != -1)
            {
               ftempaver[0] += mean->GetValue(mean->GetInt("energySD"), i);
	       avercount[0]++;
            }
         }

         if((cutZenith->widgetChBox)->IsChecked())
         {
            if(mean->GetValue(mean->GetInt("zenithSD"), i) != -1)
            {
               ftempaver[1] += mean->GetValue(mean->GetInt("zenithSD"), i);
	       avercount[1]++;
            }
         }
      }
      else if(selcuttype == 1)
      {
         if((cutEnergy->widgetChBox)->IsChecked())
         {
            if(mean->GetValue(mean->GetInt("energyFD"), i) != -1)
            {
	       if(seleyetype == 0)
	       {
	          wQuantity[0] = 1/TMath::Power(neg->GetValue(mean->GetInt("energyFD"), i), 2);
	          quantitySum[0] += (mean->GetValue(mean->GetInt("energyFD"), i))*wQuantity[0];
	          wQuantitySum[0] += wQuantity[0];
	       }
	       else
                  ftempaver[0] += mean->GetValue(mean->GetInt("energyFD"), i);
	       avercount[0]++;
            }
         }

         if((cutZenith->widgetChBox)->IsChecked())
         {
            if(mean->GetValue(mean->GetInt("zenithFD"), i) != -1)
            {
	       if(seleyetype == 0)
	       {
	          wQuantity[1] = 1/TMath::Power(neg->GetValue(mean->GetInt("zenithFD"), i), 2);
	          quantitySum[1] += (mean->GetValue(mean->GetInt("zenithFD"), i))*wQuantity[1];
	          wQuantitySum[1] += wQuantity[1];
	       }
	       else
                  ftempaver[1] += mean->GetValue(mean->GetInt("zenithFD"), i);
	       avercount[1]++;
            }
         }
      }
   }

   if(avercount[0] > 0)
   {
      if(seleyetype == 0)
         quantitySum[0] = quantitySum[0]/wQuantitySum[0];
      else
         ftempaver[0] = ftempaver[0]/avercount[0];
   }
   else
   {
      if(seleyetype == 0)
         quantitySum[0] = -1;
      else
         ftempaver[0] = -1;
   }

   if(avercount[1] > 0)
   {
      if(seleyetype == 0)
         quantitySum[1] = quantitySum[1]/wQuantitySum[1];
      else
         ftempaver[1] = ftempaver[1]/avercount[1];
   }
   else
   {
      if(seleyetype == 0)
         quantitySum[1] = -1;
      else
         ftempaver[1] = -1;
   }

   if(DBGSIG > 1)
      cout << "# IsInsideCuts          #: " << "All values: Nr. eyes = " << avercount[0] << ", " << avercount[1] << ", Average = " << ftempaver[0] << ", " << ftempaver[1] << ", Combined = " << quantitySum[0] << ", " << quantitySum[1] << endl;

   if(DBGSIG > 1)
      cout << "# IsInsideCuts          #: ";

   for(int i = 0; i < ALLEYES; i++)
   {
      // SD observables for cut
      if(selcuttype == 0)
      {
         if((cutEnergy->widgetChBox)->IsChecked())
         {
	    if(seleyetype == 0)		// combined
               ftemp[0] = mean->GetValue(mean->GetInt("energySD"), i);
	    else if(seleyetype == 1)	// any eye
               ftemp[0] = mean->GetValue(mean->GetInt("energySD"), i);
	    else if(seleyetype == 2)	// average
               ftemp[0] = ftempaver[0];

	    if(DBGSIG > 1)
	       cout << "Energy = " << ftemp[0] << "\t";

            if( (ftemp[0] != -1) )
            {
               if((ftemp[0] > ecutBins[selectedBin[0]]) && (ftemp[0] <= ecutBins[selectedBin[0]+1]))
                  sepcut[0] = true;
	       else
                  sepcut[0] = false;
            }
	    else
               sepcut[0] = false;
         }
	 else
            sepcut[0] = true;

         if((cutZenith->widgetChBox)->IsChecked()) 
         {
	    if(seleyetype == 0)		// combined
               ftemp[0] = mean->GetValue(mean->GetInt("zenithSD"), i);
	    else if(seleyetype == 1)	// any eye
               ftemp[0] = mean->GetValue(mean->GetInt("zenithSD"), i);
	    else if(seleyetype == 2)	// average
               ftemp[0] = ftempaver[1];

	    if(DBGSIG > 1)
	       cout << "Zenith = " << ftemp[0] << "\t";

            if( (ftemp[0] != -1) )
            {
               if((ftemp[0] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
                  sepcut[1] = true;
	       else
                  sepcut[1] = false;
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
         if((cutEnergy->widgetChBox)->IsChecked()) 
         {
	    if(seleyetype == 0)		// combined
               ftemp[0] = quantitySum[0];
	    else if(seleyetype == 1)	// any eye
               ftemp[0] = mean->GetValue(mean->GetInt("energyFD"), i);
	    else if(seleyetype == 2)	// average
               ftemp[0] = ftempaver[0];

	    if(DBGSIG > 1)
	       cout << ftemp[0] << " ";

            if( (ftemp[0] != -1) )
            {
               if((ftemp[0] > ecutBins[selectedBin[0]]) && (ftemp[0] <= ecutBins[selectedBin[0]+1]))
                  sepcut[0] = true;
	       else
                  sepcut[0] = false;
            }
	    else
               sepcut[0] = false;
         }
	 else
            sepcut[0] = true;

         if((cutZenith->widgetChBox)->IsChecked()) 
         {
	    if(seleyetype == 0)		// combined
               ftemp[0] = quantitySum[1];
	    else if(seleyetype == 1)	// any eye
               ftemp[0] = mean->GetValue(mean->GetInt("zenithFD"), i);
	    else if(seleyetype == 2)	// average
               ftemp[0] = ftempaver[1];

	    if(DBGSIG > 1)
  	       cout << ftemp[0] << " ";

            if( (ftemp[0] != -1) )
            {
               if((ftemp[0] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
                  sepcut[1] = true;
	       else
                  sepcut[1] = false;
            }
	    else
               sepcut[1] = false;
         }
	 else
            sepcut[1] = true;
      }

      if((cutRisetime->widgetChBox)->IsChecked()) 
      {
         ftemp[0] = mean->GetValue(mean->GetInt("risetimerecalc"), i);
         ftemp[1] = neg->GetValue(neg->GetInt("risetimerecalc"), i);

	 if(DBGSIG > 1)
	    cout << ftemp[1]/ftemp[0] << " ";

         if( ((ftemp[0] != -1) && (ftemp[1] != -1)) )
         {
            if(ftemp[1]/ftemp[0] <= (cutRisetime->widgetNE[0])->GetValue())
               sepcut[2] = true;
	    else
               sepcut[2] = false;
         }
	 else
            sepcut[2] = false;
      }
      else
         sepcut[2] = true;

      isinside = sepcut[0] & sepcut[1] & sepcut[2];

      if(DBGSIG > 1)
         cout << endl << "# IsInsideCuts          #: " << "isinside = " << (int)isinside << ", sepcut[0] = " << (int)sepcut[0] << ", sepcut[1] = " << (int)sepcut[1] << ", sepcut[2] = " << (int)sepcut[2] << endl;

      if(isinside)
         seleye->push_back(i);
   }

   delete[] sepcut;
   delete[] ftemp;
   delete[] ftempaver;
   delete[] avercount;
   delete[] wQuantity;
   delete[] quantitySum;
   delete[] wQuantitySum;

   if(seleye->size() > 0)
   {
      if( (seleyetype >= 0) && (seleyetype < 3) )
         return seleyetype;
      else
         return -1;
   }
   else
      return -1;
}

// Book the method, depending on which was chosen
int MyFrame::BookTheMethod(TMVA::Factory *factory)
{
   string *stemp;
   stemp = new string[3];

   stemp[0] = (methodsSelect->widgetCB)->GetStringSelection();
   cout << "# BookTheMethod         #: " << stemp[0] << ", " << GetMethodName(stemp[0]) << endl;

   // Prepare formula for FDA methods
   for(int i = 0; i <= nrselobs; i++)
   {
      if(i == 0)
         stemp[2] = "(" + ToString(i) + ")";
      else
         stemp[2] = stemp[2] + "+(" + ToString(i) + ")*x" + ToString(i-1);
   }
   cout << "# BookTheMethod         #: Using formula " << stemp[2] << endl;

   if(GetMethodName(stemp[0]) == "Cuts")
      factory->BookMethod(TMVA::Types::kCuts, (GetMethodName(stemp[0])).c_str(), "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");
   else if(GetMethodName(stemp[0]) == "CutsD")
      factory->BookMethod(TMVA::Types::kCuts, (GetMethodName(stemp[0])).c_str(), "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");
   else if(GetMethodName(stemp[0]) == "CutsPCA")
      factory->BookMethod(TMVA::Types::kCuts, (GetMethodName(stemp[0])).c_str(), "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");
   else if(GetMethodName(stemp[0]) == "CutsGA")
      factory->BookMethod(TMVA::Types::kCuts, (GetMethodName(stemp[0])).c_str(), "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");
   else if(GetMethodName(stemp[0]) == "CutsSA")
      factory->BookMethod(TMVA::Types::kCuts, (GetMethodName(stemp[0])).c_str(), "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");
   else if(GetMethodName(stemp[0]) == "Likelihood")
      factory->BookMethod(TMVA::Types::kLikelihood, (GetMethodName(stemp[0])).c_str(), "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");
   else if(GetMethodName(stemp[0]) == "LikelihoodD")
      factory->BookMethod(TMVA::Types::kLikelihood, (GetMethodName(stemp[0])).c_str(), "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate");
   else if(GetMethodName(stemp[0]) == "LikelihoodPCA")
      factory->BookMethod(TMVA::Types::kLikelihood, (GetMethodName(stemp[0])).c_str(), "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"); 
   else if(GetMethodName(stemp[0]) == "LikelihoodKDE")
      factory->BookMethod(TMVA::Types::kLikelihood, (GetMethodName(stemp[0])).c_str(), "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50"); 
   else if(GetMethodName(stemp[0]) == "LikelihoodMIX")
      factory->BookMethod(TMVA::Types::kLikelihood, (GetMethodName(stemp[0])).c_str(), "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50"); 
   else if(GetMethodName(stemp[0]) == "PDERS")
      factory->BookMethod(TMVA::Types::kPDERS, (GetMethodName(stemp[0])).c_str(), "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");
   else if(GetMethodName(stemp[0]) == "PDERSD")
      factory->BookMethod(TMVA::Types::kPDERS, (GetMethodName(stemp[0])).c_str(), "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate");
   else if(GetMethodName(stemp[0]) == "PDERSPCA")
      factory->BookMethod(TMVA::Types::kPDERS, (GetMethodName(stemp[0])).c_str(), "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA");
   else if(GetMethodName(stemp[0]) == "PDEFoam")
      factory->BookMethod(TMVA::Types::kPDEFoam, (GetMethodName(stemp[0])).c_str(), "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");
   else if(GetMethodName(stemp[0]) == "PDEFoamBoost")
      factory->BookMethod(TMVA::Types::kPDEFoam, (GetMethodName(stemp[0])).c_str(), "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");
   else if(GetMethodName(stemp[0]) == "KNN")
      factory->BookMethod(TMVA::Types::kKNN, (GetMethodName(stemp[0])).c_str(), "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");
   else if(GetMethodName(stemp[0]) == "LD")
      factory->BookMethod(TMVA::Types::kLD, (GetMethodName(stemp[0])).c_str(), "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
   else if(GetMethodName(stemp[0]) == "Fisher")
      factory->BookMethod(TMVA::Types::kFisher, (GetMethodName(stemp[0])).c_str(), "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
   else if(GetMethodName(stemp[0]) == "FisherG")
      factory->BookMethod(TMVA::Types::kFisher, (GetMethodName(stemp[0])).c_str(), "H:!V:VarTransform=Gauss");
   else if(GetMethodName(stemp[0]) == "BoostedFisher")
      factory->BookMethod(TMVA::Types::kFisher, (GetMethodName(stemp[0])).c_str(), "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring");
   else if(GetMethodName(stemp[0]) == "HMatrix")
      factory->BookMethod(TMVA::Types::kHMatrix, (GetMethodName(stemp[0])).c_str(), "!H:!V:VarTransform=None");
   else if(GetMethodName(stemp[0]) == "FDA_GA")
   {
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1";
      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), stemp[1].c_str());
   }
   else if(GetMethodName(stemp[0]) == "FDA_SA")
   {
//      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale";
      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), stemp[1].c_str());
   }
   else if(GetMethodName(stemp[0]) == "FDA_MC")
   {
//      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1");
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1";
      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), stemp[1].c_str());
   }
   else if(GetMethodName(stemp[0]) == "FDA_MT")
   {
//      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch");
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch";
      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), stemp[1].c_str());
   }
   else if(GetMethodName(stemp[0]) == "FDA_GAMT")
   {
//      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim");
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim";
      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), stemp[1].c_str());
   }
   else if(GetMethodName(stemp[0]) == "FDA_MCMT")
   {
//      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20");
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20";
      factory->BookMethod(TMVA::Types::kFDA, (GetMethodName(stemp[0])).c_str(), stemp[1].c_str());
   }
   else if(GetMethodName(stemp[0]) == "MLP")
      factory->BookMethod(TMVA::Types::kMLP, (GetMethodName(stemp[0])).c_str(), "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");
   else if(GetMethodName(stemp[0]) == "MLPBFGS")
      factory->BookMethod(TMVA::Types::kMLP, (GetMethodName(stemp[0])).c_str(), "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator");
   else if(GetMethodName(stemp[0]) == "MLPBNN")
      factory->BookMethod(TMVA::Types::kMLP, (GetMethodName(stemp[0])).c_str(), "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator");
   else if(GetMethodName(stemp[0]) == "CFMlpANN")
      factory->BookMethod(TMVA::Types::kCFMlpANN, (GetMethodName(stemp[0])).c_str(), "!H:!V:NCycles=2000:HiddenLayers=N+1,N");
   else if(GetMethodName(stemp[0]) == "TMlpANN")
      factory->BookMethod(TMVA::Types::kTMlpANN, (GetMethodName(stemp[0])).c_str(), "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3");
   else if(GetMethodName(stemp[0]) == "SVM")
      factory->BookMethod(TMVA::Types::kSVM, (GetMethodName(stemp[0])).c_str(), "Gamma=0.25:Tol=0.001:VarTransform=Norm");
   else if(GetMethodName(stemp[0]) == "BDT")
      factory->BookMethod(TMVA::Types::kBDT, (GetMethodName(stemp[0])).c_str(), "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
   else if(GetMethodName(stemp[0]) == "BDTG")
      factory->BookMethod(TMVA::Types::kBDT, (GetMethodName(stemp[0])).c_str(), "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");
   else if(GetMethodName(stemp[0]) == "BDTB")
      factory->BookMethod(TMVA::Types::kBDT, (GetMethodName(stemp[0])).c_str(), "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");
   else if(GetMethodName(stemp[0]) == "BDTD")
      factory->BookMethod(TMVA::Types::kBDT, (GetMethodName(stemp[0])).c_str(), "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");
   else if(GetMethodName(stemp[0]) == "BDTF")
      factory->BookMethod(TMVA::Types::kBDT, "BDTMitFisher", "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");
   else if(GetMethodName(stemp[0]) == "RuleFit")
      factory->BookMethod(TMVA::Types::kRuleFit, (GetMethodName(stemp[0])).c_str(), "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");
   else if(GetMethodName(stemp[0]) == "All")
   {
      factory->BookMethod(TMVA::Types::kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");
      factory->BookMethod(TMVA::Types::kCuts, "CutsD", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");
      factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");
      factory->BookMethod(TMVA::Types::kLikelihood, "LikelihoodPCA", "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"); 
      factory->BookMethod(TMVA::Types::kPDERS, "PDERS", "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");
      factory->BookMethod(TMVA::Types::kKNN, "KNN", "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");
      factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
      factory->BookMethod(TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
      stemp[1] = "H:!V:Formula=" + stemp[2] + ":ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1";
      factory->BookMethod(TMVA::Types::kFDA, "FDA_GA", stemp[1].c_str());
//      factory->BookMethod(TMVA::Types::kFDA, "FDA_GA", "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1");
      factory->BookMethod(TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator");
      factory->BookMethod(TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm");
      factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
      factory->BookMethod(TMVA::Types::kRuleFit, "RuleFit", "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");
   }
   else
   {
      delete[] stemp;
      return -1;
   }

   delete[] stemp;
   return 0;
}

int MyFrame::GetTrainingShift(string *mvafilename)
{
   string *stemp;
   stemp = new string[2];
   char *ctemp;
   ctemp = new char[1024];
   double *dtemp;
   dtemp = new double;

   TFile *mvafile = TFile::Open(mvafilename->c_str(), "READ");
   GetCorrelations(mvafile);
   mvafile->Close();

/*   // Check the stats used for the MVA training from file ./results/transformation_stats.dat and save them to vectors
   stemp[0] = string(rootdir) + "/results/transformation_stats.dat";*/
   stemp[0] = (*currentAnalysisDir) + "/transformation_stats.dat";
   ifstream fstats;
   fstats.open(stemp[0].c_str(), ifstream::in);

   if(!statsMin.empty()) statsMin.erase(statsMin.begin(), statsMin.end());
   if(!statsMax.empty()) statsMax.erase(statsMax.begin(), statsMax.end());

   if(fstats.is_open())
   {
      if(fstats.peek() == '#')
         fstats.getline(ctemp, 1024, '\n');

      for(int i = 0; i < nrobs; i++)
      {
         if(obssel[i])
         {
            if(DBGSIG > 0)
               cout << "# GetTrainingShift      #: " << "Stats for " << observables[i] << ": ";
            fstats >> *dtemp >> *dtemp;	// Mean and RMS of the testing distribution
            fstats >> *dtemp;		// Minimum of the testing distribution
            statsMin.push_back(*dtemp);
            if(DBGSIG > 0)
               cout << *dtemp << ", ";
            fstats >> *dtemp;		// Maximum of the testing distribution
            statsMax.push_back(*dtemp);
            if(DBGSIG > 0)
               cout << *dtemp << endl;
            fstats.ignore(1, '\n');
         }
      }
   }

   fstats.close();
   
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
   stemp = new string;
   int *itemp;
   vector<string> vtemp;

   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
         vtemp.push_back(observables[i]);
   }

   itemp = new int[2];
   itemp[0] = 0;

   for(int i = 0; i < nrselobs; i++)
   {
      for(int j = i+1; j < nrselobs; j++)
      {
         *stemp = "scat_" + vtemp[j] + "_vs_" + vtemp[i] + "_Signal_Id";
	 if(DBGSIG > 0)
	    cout << "# GetCorrelations       #: " << "Signal scatter plot: " << *stemp << " (" << ((TH2F*)corrDir->Get( (*stemp).c_str() ))->GetCorrelationFactor() << ")" << endl;
         vtemp.push_back(*stemp);

         *stemp = "scat_" + vtemp[j] + "_vs_" + vtemp[i] + "_Background_Id";
	 if(DBGSIG > 0)
   	    cout << "# GetCorrelations       #: " << "Background scatter plot: " << *stemp << " (" << ((TH2F*)corrDir->Get( (*stemp).c_str() ))->GetCorrelationFactor() << ")" << endl;
         vtemp.push_back(*stemp);
      }
   }

   vtemp.erase(vtemp.begin(), vtemp.begin()+nrselobs);


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
	    (*sigCorMat)(i,j) = ((TH2F*)corrDir->Get( (vtemp[itemp[0]]).c_str() ))->GetCorrelationFactor();
	    (*backCorMat)(i,j) = ((TH2F*)corrDir->Get( (vtemp[itemp[0]+1]).c_str() ))->GetCorrelationFactor();

	    (*sigCorMat)(j,i) = ((TH2F*)corrDir->Get( (vtemp[itemp[0]]).c_str() ))->GetCorrelationFactor();
	    (*backCorMat)(j,i) = ((TH2F*)corrDir->Get( (vtemp[itemp[0]+1]).c_str() ))->GetCorrelationFactor();

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

      cout << "# GetCorrelations       #: " << "Bacnground correlation matrix: " << endl;

      for(int i = 0; i < nrselobs; i++)
      {
         for(int j = 0; j < nrselobs; j++)
         {
            cout << (*backCorMat)(i,j) << " | ";
         }
         cout << endl;
      }
   }

   delete stemp;
   delete[] itemp;
   delete corrDir;
}

void MyFrame::GetApplyCorrelations(string *corname)
{
   string *stemp;
   stemp = new string[2];

   TCanvas *canvCor;
   TTree *TCor;
   TGraph *gr1;

   vector<string> vtemp;

   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
         vtemp.push_back(observables[i]);
   }

   TFile *fileCor = new TFile(corname->c_str(), "READ");

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
               stemp[1] = vtemp[i] + ":" + vtemp[j];
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

   delete[] stemp;
   delete canvCor;
   delete fileCor;
}

// Apply the MVA cut and save all events and MVA values into one output files (in order to plot anything from it)
void MyFrame::CreateOutput(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean)
{
   vector<string> obs;

   int *sigcount, *backcount;
   sigcount = new int[nrselobs+1];
   backcount = new int[nrselobs+1];

   for(int i = 0; i <= nrselobs; i++)
   {
      sigcount[i] = 0;
      backcount[i] = 0;
   }

   int cnt;

   TBranch *mvabranch;
   if(mean == 0)
      mvabranch = app->Branch("MVA", &obsvars[nrobs], "MVA/F");

   if(DBGSIG > 0)
      cout << "# CreateOutput          #: " << "Selected observables:" << endl;
   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
      {
	 obs.push_back(observables[i]);
         cout << "# CreateOutput          #: " << "- " << observables[i] << endl;
      }
   }
   nrselobs = obs.size();

   if(DBGSIG > 1)
      cout << "# CreateOutput          #: " << "Number of events in the tree = " << app->GetEntries() << endl;

   if(!application)
      GetErrors(app, obsvars, obs, curtree);

   // Determine which values are signal and which are background
   for(int ievt = 0; ievt < app->GetEntries(); ievt++)
   {
      app->GetEntry(ievt);
      cnt = 0;

      obsvars[nrobs] = reader->EvaluateMVA(mvamethod.c_str());
      if(mean == 0)
         mvabranch->Fill();

      for(int i = 0; i <= nrselobs; i++)
      {
	 if((i < nrselobs) && (DBGSIG > 1))
            cout << "# CreateOutput          #: " << "  Event = " << ievt << ": values = " << obsvars[3*i] << ", " << obsvars[3*i+1] << ", " << obsvars[3*i+2] << ", observable = " << obs[i] << ", mvavalue = " << obsvars[nrobs] << endl;

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
            cout << " - " << obs[i] << " = " << sigcount[i] << " vs. " << backcount[i] << endl;
         else
            cout << " - MVA = " << sigcount[i] << " vs. " << backcount[i] << endl;
      }

      cout << endl;
   }

   cutResults.push_back(sigcount[0]);
   cutResults.push_back(backcount[0]);

   if(DBGSIG > 0)
      cout << "# CreateOutput          #: " << "Finished here..." << endl;

   delete[] sigcount;
   delete[] backcount;
}

/*// Apply the MVA cut and create the MVA plots
void MyFrame::CreateMVAPlots(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean)
{
   // Prepare colors for signal, background and MVA cut line
   static Int_t c_SignalLine     = TColor::GetColor("#0000ee");
   static Int_t c_SignalFill     = TColor::GetColor("#7d99d1");
   static Int_t c_AllLine        = TColor::GetColor("#ff0000");
   static Int_t c_AllFill        = TColor::GetColor("#ff0000");
   static Int_t c_MvaCut         = TColor::GetColor("#ffff66");

   vector<string> obs;

   cout << "Selected observables:" << endl;
   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
      {
	 obs.push_back(observables[i]);
         cout << "- " << observables[i] << endl;
      }
   }
   nrselobs = obs.size();

cout << "Number of events in the tree = " << app->GetEntries() << endl;
   
   // All additional things we need for plotting
   TLegend *legend;
   TLine *line;
   TH1F *basesig[nrselobs+1];
   TH1F *baseback[nrselobs+1];

   string *stemp;
   stemp = new string[2];

   int legendFill = 1001;
   float *yhistlimit;
   yhistlimit = new float[2];
   int cnt;
   double *dtemp;
   dtemp = new double[2];

   float *max;
   int *sigcount, *backcount;
   max = new float[nrselobs+1];
   sigcount = new int[nrselobs+1];
   backcount = new int[nrselobs+1];

   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);

   // Prepare the histograms
   for(int i = 0; i <= nrselobs; i++)
   {
      if(i < nrselobs)
      {
         SetupBinning(obs[i], yhistlimit);
         stemp[0] = "basesig" + ToString(i);
         basesig[i] = new TH1F(stemp[0].c_str(), obs[i].c_str(), 100, yhistlimit[0], yhistlimit[1]);
         basesig[i]->SetBit(TH1::kCanRebin);
         stemp[0] = "baseback" + ToString(i);
         baseback[i] = new TH1F(stemp[0].c_str(), obs[i].c_str(), 100, yhistlimit[0], yhistlimit[1]);
         baseback[i]->SetBit(TH1::kCanRebin);
      }
      else
      {
         SetupBinning("MVA", yhistlimit);
         stemp[0] = "basesig" + ToString(i);
         basesig[i] = new TH1F(stemp[0].c_str(), "mva", 100, yhistlimit[0], yhistlimit[1]);
         basesig[i]->SetBit(TH1::kCanRebin);
         stemp[0] = "baseback" + ToString(i);
         baseback[i] = new TH1F(stemp[0].c_str(), "mva", 100, yhistlimit[0], yhistlimit[1]);
         baseback[i]->SetBit(TH1::kCanRebin);
      }

      max[i] = 0.0;
      sigcount[i] = 0;
      backcount[i] = 0;
   }

   if(!application)
      GetErrors(app, obsvars, obs, curtree);

   // Parse values into histograms
   for(int ievt = 0; ievt < app->GetEntries(); ievt++)
   {
      app->GetEntry(ievt);
      cnt = 0;

      for(int i = 0; i <= nrselobs; i++)
      {
         if(i < nrselobs)
	 {
cout << "  Event = " << ievt << ": values = " << obsvars[3*i] << ", " << obsvars[3*i+1] << ", " << obsvars[3*i+2] << ", observable = " << obs[i] << endl;
            if(reader->EvaluateMVA(mvamethod.c_str()) >= (cutMva->widgetNE[0])->GetValue())
            {
               basesig[i]->Fill(obsvars[3*i]);
               sigcount[i]++;
	    }
	    else
	    {
               baseback[i]->Fill(obsvars[3*i]);
	       backcount[i]++;
	    }
	 }
	 else
	 {
cout << "  MVA Event = " << ievt << ": value = " << reader->EvaluateMVA(mvamethod.c_str()) << ", observable = " << "MVA" << endl;
            if(reader->EvaluateMVA(mvamethod.c_str()) >= (cutMva->widgetNE[0])->GetValue())
            {
               basesig[i]->Fill(reader->EvaluateMVA(mvamethod.c_str()));
               sigcount[i]++;
	    }
	    else
	    {
               baseback[i]->Fill(reader->EvaluateMVA(mvamethod.c_str()));
	       backcount[i]++;
	    }
	 }
      }
   }

   cout << "Signal vs. background (" << app->GetTitle() << "):" << endl;
   for(int i = 0; i <= nrselobs; i++)
   {
      if(i < nrselobs)
         cout << " - " << obs[i] << " = " << sigcount[i] << " vs. " << backcount[i] << endl;
      else
         cout << " - MVA = " << sigcount[i] << " vs. " << backcount[i] << endl;
   }

   cout << endl;
   cutResults.push_back(sigcount[0]);
   cutResults.push_back(backcount[0]);

   if(mean == 0)
   {
      stemp[1] = "mkdir -p " + (*currentAnalysisDir) + "/mean";
      system(stemp[1].c_str());
   }
   else if(mean == -1)
   {
      stemp[1] = "mkdir -p " + (*currentAnalysisDir) + "/negerror";
      system(stemp[1].c_str());
   }
   else if(mean == 1)
   {
      stemp[1] = "mkdir -p " + (*currentAnalysisDir) + "/poserror";
      system(stemp[1].c_str());
   }

   for(int i = 0; i <= nrselobs; i++)
   {
      if(basesig[i]->GetMaximum() > max[i])
         max[i] = basesig[i]->GetMaximum();
      if(baseback[i]->GetMaximum() > max[i])
         max[i] = baseback[i]->GetMaximum();

//      cout << "Maximum = " << max[i] << endl;

      basesig[i]->SetLineColor(c_SignalLine);
      basesig[i]->SetLineWidth(2);
      basesig[i]->SetFillColor(c_SignalFill);
      basesig[i]->SetFillStyle(1001);
      baseback[i]->SetLineColor(c_AllLine);
      baseback[i]->SetLineWidth(2);
      baseback[i]->SetFillColor(c_AllFill);
      baseback[i]->SetFillStyle(3554);

      basesig[i]->Draw();
      baseback[i]->Draw("same");

      basesig[i]->GetYaxis()->SetRangeUser(0.,max[i]*1.2);
      basesig[i]->SetMaximum(max[i]*1.2);

      if(i < nrselobs)
         SetupAxis(basesig[i], obs[i]);
      else
      {
         SetupAxis(basesig[i], "MVA");
         line = new TLine((cutMva->widgetNE[0])->GetValue(), 0., (cutMva->widgetNE[0])->GetValue(), basesig[i]->GetMaximum());
         line->SetLineWidth(2);
         line->SetLineStyle(7);
         line->SetLineColor(kOrange+2);
         line->Draw("same");
      }

      // Draw legend
      TLegend *legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-.12, gPad->GetLeftMargin()+.28, 1-gPad->GetTopMargin());
      legend->SetFillStyle(legendFill);
      legend->SetFillColor(c_MvaCut);
      stemp[0] = "MVA cut background events (" + ToString(backcount[i]) + ")";
      legend->AddEntry(baseback[i],stemp[0].c_str(),"f");
      stemp[0] = "MVA cut signal events (" + ToString(sigcount[i]) + ")";
      legend->AddEntry(basesig[i],stemp[0].c_str(),"f");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      if(i < nrselobs)
      {
         if(mean == 0)	// mean values
	 {
            stemp[1] = (*currentAnalysisDir) + "/mean/mva_analysis_" + signalName + "_" + obs[i] + ".pdf";
            c1->SaveAs(stemp[1].c_str());
	 }
	 else if(mean == -1)	// negative errors
	 {
            stemp[1] = (*currentAnalysisDir) + "/negerror/mva_analysis_" + signalName + "_" + obs[i] + ".pdf";
            c1->SaveAs(stemp[1].c_str());
	 }
	 else if(mean == 1)	// positive errors
	 {
            stemp[1] = (*currentAnalysisDir) + "/poserror/mva_analysis_" + signalName + "_" + obs[i] + ".pdf";
            c1->SaveAs(stemp[1].c_str());
	 }
      }
      else
      {
         if(mean == 0)	// mean values
	 {
            stemp[1] = (*currentAnalysisDir) + "/mean/mva_analysis_" + signalName + "_MVA.pdf";
            c1->SaveAs(stemp[1].c_str());
	 }
	 else if(mean == -1)	// negative errors
	 {
            stemp[1] = (*currentAnalysisDir) + "/negerror/mva_analysis_" + signalName + "_MVA.pdf";
            c1->SaveAs(stemp[1].c_str());
	 }
	 else if(mean == 1)	// positive errors
	 {
            stemp[1] = (*currentAnalysisDir) + "/poserror/mva_analysis_" + signalName + "_MVA.pdf";
            c1->SaveAs(stemp[1].c_str());
	 }
      }

      delete legend;
   }

   delete[] max;
   delete[] sigcount;
   delete[] backcount;
   for(int i = 0; i <= nrselobs; i++)
   {
      delete basesig[i];
      delete baseback[i];
   }
   delete c1;

   delete[] dtemp;
   delete[] stemp;
   delete[] yhistlimit;
}*/

void MyFrame::SetupBinning(string obs, float *limit)
{
   if(obs == "MVA")
   {
      limit[0] = mvalimit[0];
      limit[1] = mvalimit[1];
   }
   else
   {
      limit[0] = generalObservables->GetMin(obs);
      limit[1] = generalObservables->GetMax(obs);
   }
}

void MyFrame::SetupAxis(TH1F *hist, string obs)
{
   string *obsdesc;
   obsdesc = new string;

   if(obs == "MVA")
      *obsdesc = "MVA variable";
   else
      *obsdesc = generalObservables->GetLabel(obs);

   hist->GetYaxis()->SetTitle("Number of events");
   hist->GetXaxis()->SetTitle(obsdesc->c_str());

   hist->GetXaxis()->SetTitleOffset(1.2);
   hist->GetXaxis()->CenterTitle(kTRUE);
   hist->GetXaxis()->SetLabelSize(0.028);
   hist->GetXaxis()->SetLabelOffset(0.015);
   hist->GetYaxis()->SetTitleOffset(1.3);
   hist->GetYaxis()->CenterTitle(kTRUE);
   hist->GetYaxis()->SetLabelSize(0.028);
   hist->GetYaxis()->SetLabelOffset(0.015);

   hist->SetTitle("");

   delete obsdesc;
}

void MyFrame::GetErrors(TTree *app, float *obsvars, vector<string> obs, int curtree)
{
   string *stemp;
   stemp = new string[2];
   ofstream outFile;
   double *dtemp;
   dtemp = new double[2];

   int cnt = 0;

   // Delete the old MVA error files
   stemp[0] = "rm -fr " + (*currentAnalysisDir) + "/mva_error" + ToString(curtree) + ".dat";

   if(DBGSIG > 0)
      cout << "# GetErrors             #: " << "Deleting old mva_error file: " << stemp[0] << endl;
   system(stemp[0].c_str());

   if(DBGSIG > 1)
      cout << "# GetErrors             #: " << "Current tree = " << curtree << ", Observables size = " << nrselobs << endl;

   for(int ievt = 0; ievt < app->GetEntries(); ievt++)
   {
      outFile.open(((*currentAnalysisDir) + "/mva_error" + ToString(curtree) + ".dat").c_str(), ofstream::out | ofstream::app);
      outFile << nrselobs << "\t";
     
      app->GetEntry(ievt);

      if(DBGSIG > 1)
         cout << "# GetErrors             #: " << "Entry printout " << ievt << ":" << endl;
      for(int j = 0; j < 3*nrselobs; j++)
      {
         if(DBGSIG > 1)
            cout << "# GetErrors             #: " << "j = " << j << ": observable value = " << obsvars[j] << endl;
	 outFile << obsvars[j] << "\t";
      }
      cnt = 0;

      double *norm1, *norm2;
      norm1 = new double[3];
      norm2 = new double[3];
      covMatNeg = new TMatrixD(nrselobs,nrselobs);
      covMatPos = new TMatrixD(nrselobs,nrselobs);
      eigenValMatNeg = new TMatrixD(nrselobs,nrselobs);
      eigenValMatPos = new TMatrixD(nrselobs,nrselobs);
     
      if(DBGSIG > 0)
         cout << "# GetErrors             #: " << "Covariance matrix:" << endl;

      for(int i = 0; i < nrselobs; i++)
      {
	 norm1[0] = ((obsvars[3*i] - statsMin[i])/(statsMax[i] - statsMin[i]))*2 - 1;
	 norm1[1] = ((obsvars[3*i] - obsvars[3*i+1] - statsMin[i])/(statsMax[i] - statsMin[i]))*2 - 1;
	 norm1[2] = ((obsvars[3*i] + obsvars[3*i+2] - statsMin[i])/(statsMax[i] - statsMin[i]))*2 - 1;
         if(DBGSIG > 0)
   	    cout << "# GetErrors             #: " << "Normalized values (" << obs[i] << "): " << norm1[0] << ", " << norm1[1] << ", " << norm1[2] << endl;

         for(int j = 0; j < nrselobs; j++)
         {
	    norm2[0] = ((obsvars[3*j] - statsMin[j])/(statsMax[j] - statsMin[j]))*2 - 1;
	    norm2[1] = ((obsvars[3*j] - obsvars[3*j+1] - statsMin[j])/(statsMax[j] - statsMin[j]))*2 - 1;
	    norm2[2] = ((obsvars[3*j] + obsvars[3*j+2] - statsMin[j])/(statsMax[j] - statsMin[j]))*2 - 1;
            if(DBGSIG > 0)
	    {
	       cout << "# GetErrors             #: " << "Xmin and Xmax (" << obs[j] << "): " << statsMin[j] << ", " << statsMax[j] << endl;
	       cout << "# GetErrors             #: " << "Original values (" << obs[j] << "): " << obsvars[3*j] << ", " << obsvars[3*j] - obsvars[3*j+1] << ", " << obsvars[3*j] + obsvars[3*j+2] << endl;
	       cout << "# GetErrors             #: " << "Normalized values (" << obs[j] << "): " << norm2[0] << ", " << norm2[1] << ", " << norm2[2] << endl;
	    }

	    if((signalSelect->widgetCB)->GetSelection() == curtree-1)
	    {
               (*covMatNeg)(i,j) = (*sigCorMat)(i,j)*(norm1[0]-norm1[1])*(norm2[0]-norm2[1]);
               (*covMatPos)(i,j) = (*sigCorMat)(i,j)*(norm1[2]-norm1[0])*(norm2[2]-norm2[0]);
               if(DBGSIG > 0)
	          cout << "# GetErrors             #: " << "Signal (neg: " << (*sigCorMat)(i,j) << ", " << (norm1[0]-norm1[1]) << ", " << (norm2[0]-norm2[1]) << " = " << (*covMatNeg)(i,j) << ")" << endl;
	    }
	    else if((backgroundSelect->widgetCB)->GetSelection() == curtree-1)
	    {
               (*covMatNeg)(i,j) = (*backCorMat)(i,j)*(norm1[0]-norm1[1])*(norm2[0]-norm2[1]);
               (*covMatPos)(i,j) = (*backCorMat)(i,j)*(norm1[2]-norm1[0])*(norm2[2]-norm2[0]);
               if(DBGSIG > 0)
	          cout << "# GetErrors             #: " << "Background (neg: " << (*backCorMat)(i,j) << ", " << (norm1[0]-norm1[1]) << ", " << (norm2[0]-norm2[1]) << " = " << (*covMatNeg)(i,j) << ")" << endl;
	    }
	    else
	    {
               (*covMatNeg)(i,j) = (*otherCorMat[curtree-1])(i,j)*(norm1[0]-norm1[1])*(norm2[0]-norm2[1]);
               (*covMatPos)(i,j) = (*otherCorMat[curtree-1])(i,j)*(norm1[2]-norm1[0])*(norm2[2]-norm2[0]);
               if(DBGSIG > 0)
	          cout << "# GetErrors             #: " << "Others (neg: " << (*otherCorMat[curtree-1])(i,j) << ", " << (norm1[0]-norm1[1]) << ", " << (norm2[0]-norm2[1]) << " = " << (*covMatNeg)(i,j) << ")" << endl;
	    }

            if(DBGSIG > 0)
               cout << "# GetErrors             #: " << " | " << (*covMatNeg)(i,j) << " | " << endl;

	    eigenCovMat = new TMatrixDEigen((const TMatrixD)*covMatNeg);
	    (*eigenValMatNeg) = eigenCovMat->GetEigenValues();
	    delete eigenCovMat;
	    eigenCovMat = new TMatrixDEigen((const TMatrixD)*covMatPos);
	    (*eigenValMatPos) = eigenCovMat->GetEigenValues();
	    delete eigenCovMat;
         }

         if(DBGSIG > 0)
	    cout << endl;
      }

      dtemp[0] = 0;

      if(DBGSIG > 0)
         cout << "# GetErrors             #: " << "Diagonalized matrix (negative error): ";
      for(int i = 0; i < nrselobs; i++)
      {
	 dtemp[0] += TMath::Power((*eigenValMatNeg)(i,i),2);
         if(DBGSIG > 0)
            cout << (*eigenValMatNeg)(i,i) << " ";
      }
      
      if(DBGSIG > 0)
         cout << endl;

      dtemp[1] = TMath::Sqrt(dtemp[0]);
      if(DBGSIG > 0)
         cout << "# GetErrors             #: " << "MVA variable error (negative error) = " << dtemp[1] << endl;

      outFile << dtemp[1] << "\t";

      dtemp[0] = 0;

      if(DBGSIG > 0)
         cout << "# GetErrors             #: " << "Diagonalized matrix (positive error): ";
      for(int i = 0; i < nrselobs; i++)
      {
	 dtemp[0] += TMath::Power((*eigenValMatPos)(i,i),2);
         if(DBGSIG > 0)
            cout << (*eigenValMatPos)(i,i) << " ";
      }
      
      if(DBGSIG > 0)
         cout << endl;

      dtemp[1] = TMath::Sqrt(dtemp[0]);
      if(DBGSIG > 0)
         cout << "# GetErrors             #: " << "MVA variable error (positive error) = " << dtemp[1] << endl << endl;

      outFile << dtemp[1] << endl;
      outFile.close();

      delete covMatNeg;
      delete covMatPos;
      delete eigenValMatNeg;
      delete eigenValMatPos;
      delete[] norm2;
      delete[] norm1;
   }

   delete[] stemp;
   delete[] dtemp;
}

void MyFrame::GetMvaError(int selection, double *outvalue)
{
   if(DBGSIG > 1)
   {
      cout << "# GetMvaError           #: " << "Selection: " << selection << endl;
      cout << "# GetMvaError           #: " << "Mean: " << outvalue[0] << endl;
   }

   string *stemp;
   stemp = new string;

   *stemp = (*currentAnalysisDir) + string("/mva_error") + ToString(selection) + string(".dat");
   if(DBGSIG > 1)
      cout << "# GetMvaError           #: " << "File: " << *stemp << endl;

   ifstream infile;
   infile.open(stemp->c_str(), ifstream::in);
   int size, obssize;
   infile >> obssize;
   size = 3*obssize + 3;
   infile.seekg(0, infile.beg);
   double *readVal;
   readVal = new double[size];
   int vrstica[] = {0,0,0,0};

   vector<double> negarray;
   vector<double> posarray;
   
   if(infile.is_open())
   {
      while(1)
      {
         if(infile.eof()) break;

	 for(int i = 0; i < size; i++)
            infile >> readVal[i];

	 negarray.push_back(readVal[size-2]);
	 posarray.push_back(readVal[size-1]);

	 infile.ignore(1,' ');
	 vrstica[0]++;
      }
   }

   infile.close();

   if(DBGSIG > 1)
      cout << "# GetMvaError           #: " << "Lines read = " << vrstica[0] << endl;

   double *dtemp;
   dtemp = new double[2];
   dtemp[0] = 0;
   dtemp[1] = 0;
   for(int i = 0; i < vrstica[0]; i++)
   {
      dtemp[0] += negarray[i];
      dtemp[1] += posarray[i];
   }
   dtemp[0] = dtemp[0]/vrstica[0];
   dtemp[1] = dtemp[1]/vrstica[0];

   if(DBGSIG > 1)
      cout << "# GetMvaError           #: " << "Mean values (negative error, positive error) = " << "(" << dtemp[0] << ", " << dtemp[1] << ")" << endl;

   for(int i = 0; i < vrstica[0]; i++)
   {
      outvalue[1] += (negarray[i] - dtemp[0])*(negarray[i] - dtemp[0]);
      outvalue[2] += (posarray[i] - dtemp[1])*(posarray[i] - dtemp[1]);
   }
   outvalue[1] = TMath::Sqrt(outvalue[1]/vrstica[0]);
   outvalue[2] = TMath::Sqrt(outvalue[2]/vrstica[0]);

   if(DBGSIG > 1)
      cout << "# GetMvaError           #: " << "Sigma values (negative error, positive error) = " << "(" << outvalue[1] << ", " << outvalue[2] << ")" << endl;

   outvalue[1] = outvalue[0] - outvalue[1];
   outvalue[2] = outvalue[0] + outvalue[2];

   delete[] dtemp;
   delete stemp;
   delete[] readVal;
}
