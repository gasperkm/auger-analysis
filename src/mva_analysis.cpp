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
int MyFrame::MvaTreeFile(string *infilename, string *outfilename, int *nrEvents)
{
   string *stemp;
   stemp = new string[2];

   TFile *ifile = TFile::Open(infilename->c_str(), "READ");

   TTree *signalTempTree[nrkeys];

/*   int selectedBin[2];
   selectedBin[0] = (cutEnergyBins->widgetCB)->GetSelection();
   selectedBin[1] = (cutZenithBins->widgetCB)->GetSelection();*/

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

   vector<int> *seleye;
   seleye = new vector<int>;

   for(int j = 0; j < tempTree->GetEntries(); j++)
   {
      tempTree->GetEntry(j);

      // Check if event is inside the selected cuts
      if(!seleye->empty()) seleye->erase(seleye->begin(), seleye->end());

      if(DBGSIG > 1)
         cout << "# MvaSetTrees           #: " << "Event = " << j << endl;
      ret = IsInsideCuts(invalues, invalues_neg, invalues_pos, seleye, false);

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

            for(int k = 0; k < seleye->size(); k++)
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
      }
      // Any of the eyes is inside the cut
      else if(ret == 1)
      {
         if(DBGSIG > 1)
	    cout << "# MvaSetTrees           #: " << "Event " << j << " IsInsideCuts (" << seleye->size() << ")" << endl;
         itemp[2] = seleye->at(0);

         itemp[0]++;
         if(seleye->size() > 1)
         {
            itemp[1]++;

            if(DBGSIG > 1)
	    {
               cout << "# MvaSetTrees           #: " << "Event " << j << " has multiple eyes inside the cuts = ";
               for(int i = 0; i < seleye->size(); i++)
                  cout << seleye->at(i) << " ";
               cout << endl;
            }

            // Determine which FD eye has the smallest error on xmax
            ftemp[3] = 100.;

	    for(int k = 0; k < seleye->size(); k++)
	    {
               // Get error on xmax value (relative)
               ftemp[0] = invalues->GetValue(generalObservables->GetInt("xmax"), seleye->at(k));
               ftemp[1] = invalues_neg->GetValue(generalObservables->GetInt("xmax"), seleye->at(k));

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
		  itemp[2] = seleye->at(k);
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
            for(int k = 0; k < seleye->size(); k++)
            {
	       if(invalues->GetValue(i, k) != -1)
	       {
                  outobs[i] += invalues->GetValue(i, k);
                  outobs_neg[i] += invalues_neg->GetValue(i, k);
                  outobs_pos[i] += invalues_pos->GetValue(i, k);
	       }
	    }

	    outobs[i] = outobs[i]/(seleye->size());
	    outobs_neg[i] = outobs_neg[i]/(seleye->size());
	    outobs_pos[i] = outobs_pos[i]/(seleye->size());

	    if(DBGSIG > 0)
               cout << "# MvaSetTrees           #: " << "Writeout values (average): " << outobs[i] << ", " << outobs_neg[i] << ", " << outobs_pos[i] << endl;
	 }
         outtree->Fill();
      }
   }

   cout << "# MvaSetTrees           #: " << "Number of events inside the cuts = " << itemp[0] << endl;
   cout << "# MvaSetTrees           #: " << "Number of events with multiple eyes inside the cuts = " << itemp[1] << endl;

   ret = itemp[0];

   delete seleye;
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

int MyFrame::IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, vector<int> *seleye, bool split)
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
   // Determine how eye selection should be handled (any eye inside selection or average)
   if(!split)
      seleyetype = (eyeSelection->widgetCB)->GetSelection();
   else
      seleyetype = (splitEyeSelection->widgetCB)->GetSelection();

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

   // Calculate averages for all eyes
   for(int i = 0; i < ALLEYES; i++)
   {
      if(selcuttype == 0)
      {
         if(btemp[0])
         {
            if(mean->GetValue(mean->GetInt("energySD"), i) != -1)
            {
               ftempaver[0] += mean->GetValue(mean->GetInt("energySD"), i);
	       avercount[0]++;
            }
         }

         if(btemp[1])
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
         if(btemp[0])
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

         if(btemp[1])
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
         if(btemp[0])
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
               if(!split)
	       {
                  if((ftemp[0] > ecutBins[selectedBin[0]]) && (ftemp[0] <= ecutBins[selectedBin[0]+1]))
                     sepcut[0] = true;
	          else
                     sepcut[0] = false;
	       }
	       else
	       {
                  ftemp[2] = (float)pow(10, splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[0]));
                  ftemp[3] = (float)pow(10, splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[1]));
                  if((ftemp[0] > ftemp[2]) && (ftemp[0] <= ftemp[3]))
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
               if(!split)
	       {
                  if((ftemp[0] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
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
               if(!split)
	       {
                  if((ftemp[0] > ecutBins[selectedBin[0]]) && (ftemp[0] <= ecutBins[selectedBin[0]+1]))
                     sepcut[0] = true;
	          else
                     sepcut[0] = false;
	       }
	       else
	       {
                  ftemp[2] = (float)pow(10, splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[0]));
                  ftemp[3] = (float)pow(10, splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[1]));
                  if((ftemp[0] > ftemp[2]) && (ftemp[0] <= ftemp[3]))
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
               if(!split)
	       {
                  if((ftemp[0] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp[0] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
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
         ftemp[0] = mean->GetValue(mean->GetInt("risetimerecalc"), i);
         ftemp[1] = neg->GetValue(neg->GetInt("risetimerecalc"), i);

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
         cout << endl << "# IsInsideCuts          #: " << "isinside = " << (int)*isinside << ", sepcut[0] = " << (int)sepcut[0] << ", sepcut[1] = " << (int)sepcut[1] << ", sepcut[2] = " << (int)sepcut[2] << endl;

      if(*isinside)
         seleye->push_back(i);
   }

   delete[] sepcut;
   delete isinside;
   delete[] btemp;
   delete[] ftemp;
   delete[] ftempaver;
   delete[] avercount;
   delete[] wQuantity;
   delete[] quantitySum;
   delete[] wQuantitySum;
   delete[] selectedBin;

   if(DBGSIG > 1)
      cout << "# IsInsideCuts          #: " << "return code = " << seleyetype << endl;

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

// Perform the MVA analysis on a collection of observables
int MyFrame::PerformMvaAnalysis(string *infilename, string *outfilename, int type)
{
   int *nrTreeEvents;
   string *stemp;
   int *itemp;

   nrTreeEvents = new int[nrkeys];
   stemp = new string[4];
   itemp = new int[3];

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
   cout << "# StartMvaAnalysis      #: " << "Weights directory after = " << (TMVA::gConfig().GetIONames()).fWeightFileDir << endl;
//   itemp[1]++;
//   progress->Update(itemp[1]);

   // Adding observables to the Factory
   nrTreeEvents[0] = 0;
   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
      {
         factory->AddVariable(observables[i].c_str(), 'F');
         cout << "# StartMvaAnalysis      #: " << "Adding variable: " << observables[i] << " (" << (int)obssel[i] << ")" << endl;
         nrTreeEvents[0]++;
      }
   }
   nrselobs = MvaNoteObservables(nrTreeEvents[0]);
//   itemp[1]++;
//   progress->Update(itemp[1]);

   // Select signal and background trees (from the temporary input file)
   TTree *signalTree = new TTree;
   TTree *backgroundTree[mixednum];
   for(int j = 0; j < mixednum; j++)
      backgroundTree[j] = new TTree;
   itemp[2] = 0;

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
         cout << "# StartMvaAnalysis      #: " << "Using signal tree: " << stemp[1] << endl;
         signalTree = (TTree*)ifile->Get(stemp[0].c_str());
         nrTreeEvents[0] = signalTree->GetEntries();
/*         itemp[1]++;
         progress->Update(itemp[1]);*/
      }

      // Background tree setup
      stemp[3] = string((backgroundSelect->widgetCB)->GetStringSelection());
      if( stemp[3].find(stemp[1]) != string::npos )
      {
         cout << "# StartMvaAnalysis      #: " << "Using background tree " << itemp[2] << ": " << stemp[1] << endl;
         backgroundTree[itemp[2]] = (TTree*)ifile->Get(stemp[0].c_str());
         nrTreeEvents[1] = backgroundTree[itemp[2]]->GetEntries();
         itemp[2]++;
      }
   }

//   itemp[1]++;
//   progress->Update(itemp[1]);

   cout << "# StartMvaAnalysis      #: " << "Number of entries in signal tree = " << nrTreeEvents[0] << endl;
   cout << "# StartMvaAnalysis      #: " << "Number of entries in background tree = " << nrTreeEvents[1] << endl;

   // Add signal and background tree
   factory->AddSignalTree(signalTree, 1.0);
   for(int i = 0; i < itemp[2]; i++)
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
//   itemp[1]++;
//   progress->Update(itemp[1]);

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
      progress->Update(itemp[0]);
      delete signalTree;
      for(int j = 0; j < mixednum; j++)
         delete backgroundTree[j];
      ifile->Close();
      delete factory;
      ofile->Close();

      AlertPopup("Invalid MVA method", "The selected MVA method is invalid. Please make sure it is correctly defined.");
      delete[] stemp;
      delete[] nrTreeEvents;
      delete[] itemp;
      return -1;
   }
//   itemp[1]++;
//   progress->Update(itemp[1]);

   // Train the selected methods and save them to the weights folder
   factory->TrainAllMethods();
//   itemp[1]++;
//   progress->Update(itemp[1]);
   // Test the selected methods by applying the trained data to the test data set -> outputs saved to TestTree output file and then to the output ROOT file
   factory->TestAllMethods();
//   itemp[1]++;
//   progress->Update(itemp[1]);
   // Evaluation of methods printed to stdout
   factory->EvaluateAllMethods();
//   itemp[1]++;
//   progress->Update(itemp[1]);

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
            delete[] itemp;
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
      delete[] itemp;
      return -1;
   }

   delete[] nrTreeEvents;
   delete[] stemp;
   delete[] itemp;

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

void MyFrame::GetMvaError(int selection, double *outvalue)
{
   if(DBGSIG > 1)
   {
      cout << "# GetMvaError           #: " << "Selection: " << selection << endl;
      cout << "# GetMvaError           #: " << "Mean: " << outvalue[0] << endl;
   }

   string *stemp;
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
//   delete inTree;
   delete negarray;
   delete posarray;
   delete vrstica;
   delete[] dtemp;
   delete stemp;
}
