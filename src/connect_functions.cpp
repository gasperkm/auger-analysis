#include "frame.h"
#include "separate_functions.h"
#include "adst_mva.h"
#include "popups.h"
#include "combine.h"
#include "mvaefficiency.h"

#include <chrono>

// Double click on list to select a MVA file for input
void MyFrame::EnableMvaFile(wxCommandEvent& event)
{
   int *itemp;
   int cnt = 0;
   string *stemp;
   string signalName;

   if( (mvaList[2]->widgetLB)->GetSelection() != wxNOT_FOUND )
   {
      itemp = new int;
      stemp = new string[2];

      vector<string> backgroundType;

      // Place the file name into the text entry on the right side
      *itemp = (mvaList[2]->widgetLB)->GetSelection();
      stemp[0] = (mvaList[2]->widgetLB)->GetString(*itemp);
      (selectedMva->widgetTE)->SetValue(stemp[0]);

      // Open the file with root
      TFile *ifile = TFile::Open(stemp[0].c_str(), "READ");
      (signalSelect->widgetCB)->Clear();
      (backgroundSelect->widgetCB)->Clear();
      (dataSelect->widgetCB)->Clear();

      // Append all trees from root files
      nrkeys = (ifile->GetNkeys()-1);
      if(DBGSIG > 1)
         cout << "# EnableMvaFile         #: " << "Starting number of keys " << nrkeys << endl;
      keyslist = (TList*) ifile->GetListOfKeys();
      for(int i = 1; i <= nrkeys; i++)
      {
         if(DBGSIG > 1)
            cout << "# EnableMvaFile         #: " << "i = " << i << ", Searching for TreeA = " << strcmp((keyslist->At(i-1))->GetName(), "TreeA") << " (" << (keyslist->At(i-1))->GetName() << ")" << endl;
         if(strcmp((keyslist->At(i-1))->GetName(), "TreeA") != 0)
	 {
            signalName = "TreeS" + ToString(i);
	    stemp[1] = string(ifile->GetKey(signalName.c_str())->GetTitle());
/*	    stemp[0] = RemoveBeforeLast(&stemp[1], "/");
	    stemp[1] = RemoveFromLast(&stemp[0], "(");*/
	    (signalSelect->widgetCB)->Append(stemp[1]);
	    (backgroundSelect->widgetCB)->Append(stemp[1]);
	    (dataSelect->widgetCB)->Append(stemp[1]);
            cnt++;

	    if(stemp[1] != "Data")
 	       backgroundType.push_back(stemp[1]);
	 }
      }

      nrkeys = cnt;
      cout << "# EnableMvaFile         #: " << "Found " << nrkeys << " signal keys in the selected file." << endl;

      mixednum = backgroundType.size();
      // Add mixed combinations for background trees
      if(mixednum == 3)
      {
	 stemp[1] = backgroundType[0] + " + " + backgroundType[1];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[0] + " + " + backgroundType[2];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[1] + " + " + backgroundType[2];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
      }
      else if(mixednum == 4)
      {
	 stemp[1] = backgroundType[0] + " + " + backgroundType[1];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[0] + " + " + backgroundType[2];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[0] + " + " + backgroundType[3];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[1] + " + " + backgroundType[2];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[1] + " + " + backgroundType[3];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[2] + " + " + backgroundType[3];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);

	 stemp[1] = backgroundType[0] + " + " + backgroundType[1] + " + " + backgroundType[2];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[0] + " + " + backgroundType[1] + " + " + backgroundType[3];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[0] + " + " + backgroundType[2] + " + " + backgroundType[3];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 stemp[1] = backgroundType[1] + " + " + backgroundType[2] + " + " + backgroundType[3];
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
      }

      // Select the signal tree
      if(oldselect[0] >= cnt)
      {
         (signalSelect->widgetCB)->SetSelection(cnt-1);
	 oldselect[0] = cnt-1;
      }
      else
         (signalSelect->widgetCB)->SetSelection(oldselect[0]);

      // Select the background tree
      if(oldselect[1] >= cnt)
      {
         (backgroundSelect->widgetCB)->SetSelection(cnt-1);
	 oldselect[1] = cnt-1;
      }
      else
         (backgroundSelect->widgetCB)->SetSelection(oldselect[1]);

      // Select the data tree
      if(oldselect[2] >= cnt)
      {
         (dataSelect->widgetCB)->SetSelection(cnt-1);
	 oldselect[2] = cnt-1;
      }
      else
         (dataSelect->widgetCB)->SetSelection(oldselect[2]);

      ifile->Close();

      delete itemp;
      delete[] stemp;
   }
   else
      return;
}

// Update the dropbox options when we change number of bins or limits for energy
void MyFrame::UpdateEnergyBinSelect(wxSpinDoubleEvent& event)
{
   RunEnergyBinSelect();
}

void MyFrame::RunEnergyBinSelect()
{
   double *erange = new double;
   double *elimits = new double[2];
   elimits[0] = cutEnergy->GetNumber(cutEnergy->widgetNE[0]);
   elimits[1] = cutEnergy->GetNumber(cutEnergy->widgetNE[1]);

   string *stemp = new string;

   // Check the energy range (in logarithmic scale)
   *erange = (elimits[1] - elimits[0])/((double)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]));

   if(!ecutBins.empty()) ecutBins.erase(ecutBins.begin(), ecutBins.end());
   (cutEnergyBins->widgetCB)->Clear();

   // Save energy cut bins (in linear scale)
   for(int i = 0; i <= (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
   {
      ecutBins.push_back((double)pow(10, i*(*erange)+elimits[0]));
      if(DBGSIG > 0)
         cout << "# RunEnergyBinSelect    #: " << i << ", energy bin = " << ecutBins[i] << endl;
   }

   // Add all bins to the selection dropbox for energy
   for(int i = 0; i < (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
   {
      if(i < 9)
         *stemp = "0" + ToString(i+1) + " (" + ToSciString(ecutBins[i], 1) + " - " + ToSciString(ecutBins[i+1], 1) + ")";
      else
         *stemp = ToString(i+1) + " (" + ToSciString(ecutBins[i], 1) + " - " + ToSciString(ecutBins[i+1], 1) + ")";

      (cutEnergyBins->widgetCB)->Append(*stemp);
   }

   (cutEnergyBins->widgetCB)->SetSelection(0);

   delete stemp;
   delete erange;
   delete[] elimits;
}

// Check the amount of events inside the selected energy bin
void MyFrame::CheckEnergyBin(wxCommandEvent& event)
{
   // Determine the type of observables to cut on (SD or FD)
   selcuttype = (cutObservables->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckEnergyBin        #: " << "selcuttype = " << selcuttype << endl;
   // Determine how eye selection should be handled (any eye inside selection or average)
   seleyetype = (eyeSelection->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckEnergyBin        #: " << "seleyetype = " << seleyetype << endl;

   string *stemp = new string[2];
   stemp[0] = string((selectedMva->widgetTE)->GetLineText(0));

   // Open the file for MVA analysis input
   if(stemp[0] == "")
   {
      AlertPopup("No selected MVA file", "No MVA input file selected. Please select one of the files from the third listbox.");
      delete[] stemp;
      return;
   }

   float *ftemp;
   ftemp = new float[ALLEYES];
   float *ftempaver;
   int *avercount;
   int *evtcount;
   bool deleted = true;
   bool *uoflow;
  
   TFile *ifile = TFile::Open(stemp[0].c_str(),"READ");

   int selectedBin = (cutEnergyBins->widgetCB)->GetSelection();

   stemp[0] = "Number of entries inside energy bin:\n";
   for(int i = 1; i <= nrkeys; i++)
   {
      stemp[1] = "TreeS" + ToString(i);
      TTree *tempTree = (TTree*)ifile->Get(stemp[1].c_str());
      if(DBGSIG > 0)
         cout << "# CheckEnergyBin        #: " << stemp[1] << endl;

      // SD observables for cut
      if(selcuttype == 0)
      {
         ret = Find(observables, "energySD");
         if(ret != -1)
            tempTree->SetBranchAddress("energySD", ftemp);
         else
         {
            AlertPopup("No SD energy observable found", "No SD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energySD.");
            delete[] stemp;
            delete[] ftemp;
	    ifile->Close();
            return;
         }
      }
      // FD observables for cut
      else if(selcuttype == 1)
      {
         ret = Find(observables, "energyFD");
         if(ret != -1)
            tempTree->SetBranchAddress("energyFD", ftemp);
         else
         {
            AlertPopup("No FD energy observable found", "No FD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energyFD.");
            delete[] stemp;
            delete[] ftemp;
	    ifile->Close();
            return;
         }
      }

      evtcount = new int[4];
      evtcount[0] = 0;
      evtcount[1] = 0;
      evtcount[2] = 0;
      evtcount[3] = 0;

      uoflow = new bool[2];

      ftempaver = new float;
      avercount = new int;

      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

	 deleted = true;
	 uoflow[0] = false;
	 uoflow[1] = false;

	 // Calculate averages for all eyes
	 if(seleyetype == 2)
	 {
	    *ftempaver = 0;
	    *avercount = 0;

	    for(int i = 0; i < ALLEYES; i++)
	    {
               if(ftemp[i] != -1)
	       {
                  (*ftempaver) += ftemp[i];
	          (*avercount)++;
	       }
	    }

	    if(*avercount > 0)
               *ftempaver = (*ftempaver)/(*avercount);
	    else
               *ftempaver = -1;
	 }

	 if(DBGSIG > 1)
            cout << "# CheckEnergyBin        #: ";
	 
	 // Determine if event is inside energy bin
	 for(int i = 0; i < ALLEYES; i++)
	 {
            if(DBGSIG > 1)
               cout << ftemp[i] << " ";

            if(seleyetype == 2)
               ftemp[i] = *ftempaver;

	    if((ftemp[i] != -1) && deleted)
            {
               // Event is inside the energy cut
	       if((ftemp[i] > ecutBins[selectedBin]) && (ftemp[i] <= ecutBins[selectedBin+1]))
                  deleted = false;
               // Event is below the energy cut
               else if(ftemp[i] <= ecutBins[selectedBin])
                  uoflow[0] = true;
               // Event is above the energy cut
               else if(ftemp[i] > ecutBins[selectedBin+1])
                  uoflow[1] = true;
	    }
	 }

	 if(!deleted)
	 {
            evtcount[0]++;
	    uoflow[0] = false;
	    uoflow[1] = false;

            if(DBGSIG > 1)
	       cout << "\tgood";
	 }
	 else
	 {
            if(!uoflow[0] && !uoflow[1])
	    {
               evtcount[3]++;
               
	       if(DBGSIG > 1)
  	          cout << "\tnovalues";
	    }
	    else if(uoflow[0])
	    {
               evtcount[1]++;
               
	       if(DBGSIG > 1)
	          cout << "\tbelow";
	    }
	    else if(uoflow[1])
	    {
               evtcount[2]++;
               
	       if(DBGSIG > 1)
	          cout << "\tabove";
	    }
	 }

         if(DBGSIG > 1)
	    cout << endl;
      }

      stemp[0] = stemp[0] + stemp[1] + ":\n" + "- Number of total entries = " + ToString(tempTree->GetEntries()) + "\n- Number of selected entries = " + ToString(evtcount[0]) + "\n- Total dropped entries = " + ToString(tempTree->GetEntries() - evtcount[0]) + "\n- Number of entries outside of energy bin = " + ToString(evtcount[1]) + ", " + ToString(evtcount[2]) + "\n- Number of entries with no valid values = " + ToString(evtcount[3]) + "\n\n";

      delete[] evtcount;
      delete[] uoflow;
      delete ftempaver;
      delete avercount;
   }

   ifile->Close();

   InfoPopup("Check energy bins", stemp[0]);

   delete[] stemp;
   delete[] ftemp;
}

// Update the dropbox options when we change number of bins or limits for zenith angle
void MyFrame::UpdateZenithBinSelect(wxSpinDoubleEvent& event)
{
   RunZenithBinSelect();
}
void MyFrame::RunZenithBinSelect()
{
   double *arange = new double;
   double *asinlimits = new double[2];
   asinlimits[0] = SinSquare(cutZenith->GetNumber(cutZenith->widgetNE[0]), true);
   asinlimits[1] = SinSquare(cutZenith->GetNumber(cutZenith->widgetNE[1]), true);

   string *stemp = new string;

   // Check the zenith angle range (in sin square scale)
   *arange = (asinlimits[1] - asinlimits[0])/((double)cutZenithBins->GetNumber(cutZenithBins->widgetNE[0]));

   if(!zcutBins.empty()) zcutBins.erase(zcutBins.begin(), zcutBins.end());
   (cutZenithBins->widgetCB)->Clear();

   // Save zenith cut bins (in sin square scale)
   for(int i = 0; i <= (int)cutZenithBins->GetNumber(cutZenithBins->widgetNE[0]); i++)
   {
      zcutBins.push_back((double)i*(*arange)+asinlimits[0]);
      if(DBGSIG > 0)
         cout << "# RunZenithBinSelect    #: " << i << ", zenith bin = " << zcutBins[i] << endl;
   }

   // Add all bins to the selection dropbox for zenith angle
   for(int i = 0; i < (int)cutZenithBins->GetNumber(cutZenithBins->widgetNE[0]); i++)
   {
      if(i < 9)
         *stemp = "0" + ToString(i+1) + " (" + ToString(AsinSqrt(zcutBins[i], true), 1) + " - " + ToString(AsinSqrt(zcutBins[i+1], true), 1) + ")";
      else
         *stemp = ToString(i+1) + " (" + ToString(AsinSqrt(zcutBins[i], true), 1) + " - " + ToString(AsinSqrt(zcutBins[i+1], true), 1) + ")";

      (cutZenithBins->widgetCB)->Append(*stemp);
   }

   (cutZenithBins->widgetCB)->SetSelection(0);

   delete stemp;
   delete arange;
   delete[] asinlimits;
}

// Check the amount of events inside the selected energy bin
void MyFrame::CheckZenithBin(wxCommandEvent& event)
{
   // Determine the type of observables to cut on (SD or FD)
   selcuttype = (cutObservables->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckZenithBin        #: " << "selcuttype = " << selcuttype << endl;
   // Determine how eye selection should be handled (any eye inside selection or average)
   seleyetype = (eyeSelection->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckZenithBin        #: " << "seleyetype = " << seleyetype << endl;

   string *stemp = new string[2];
   stemp[0] = string((selectedMva->widgetTE)->GetLineText(0));

   // Open the file for MVA analysis input
   if(stemp[0] == "")
   {
      AlertPopup("No selected MVA file", "No MVA input file selected. Please select one of the files from the third listbox.");
      delete[] stemp;
      return;
   }

   float *ftemp;
   ftemp = new float[ALLEYES];
   float *ftempaver;
   int *avercount;
   int *evtcount;
   bool deleted = true;
   bool *uoflow;
  
   TFile *ifile = TFile::Open(stemp[0].c_str(),"READ");

   int selectedBin = (cutZenithBins->widgetCB)->GetSelection();

   stemp[0] = "Number of entries inside zenith angle bin:\n";
   for(int i = 1; i <= nrkeys; i++)
   {
      stemp[1] = "TreeS" + ToString(i);
      TTree *tempTree = (TTree*)ifile->Get(stemp[1].c_str());

      // SD observables for cut
      if(selcuttype == 0)
      {
         ret = Find(observables, "zenithSD");
         if(ret != -1)
            tempTree->SetBranchAddress("zenithSD", ftemp);
         else
         {
            AlertPopup("No SD zenith angle observable found", "No SD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithSD.");
            delete[] stemp;
            delete[] ftemp;
	    ifile->Close();
            return;
         }
      }
      // FD observables for cut
      else if(selcuttype == 1)
      {
         ret = Find(observables, "zenithFD");
         if(ret != -1)
            tempTree->SetBranchAddress("zenithFD", ftemp);
         else
         {
            AlertPopup("No FD zenith angle observable found", "No FD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithFD.");
            delete[] stemp;
            delete[] ftemp;
	    ifile->Close();
            return;
         }
      }

      evtcount = new int[4];
      evtcount[0] = 0;
      evtcount[1] = 0;
      evtcount[2] = 0;
      evtcount[3] = 0;

      uoflow = new bool[2];

      ftempaver = new float;
      avercount = new int;

      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

	 deleted = true;
	 uoflow[0] = false;
	 uoflow[1] = false;
	 
	 // Calculate averages for all eyes
	 if(seleyetype == 2)
	 {
	    *ftempaver = 0;
	    *avercount = 0;

	    for(int i = 0; i < ALLEYES; i++)
	    {
               if(ftemp[i] != -1)
	       {
                  (*ftempaver) += ftemp[i];
	          (*avercount)++;
	       }
	    }

	    if(*avercount > 0)
               *ftempaver = (*ftempaver)/(*avercount);
	    else
               *ftempaver = -1;
	 }

	 if(DBGSIG > 1) 
            cout << "# CheckZenithBin        #: ";

	 // Determine if event is inside zenith angle bin
	 for(int i = 0; i < ALLEYES; i++)
	 {
	    if(DBGSIG > 1)
               cout << ftemp[i] << " ";

            if(seleyetype == 2)
               ftemp[i] = *ftempaver;

	    if( (ftemp[i] != -1) && deleted )
            {
               // Event is inside the zenith angle cut
	       if((ftemp[i] > AsinSqrt(zcutBins[selectedBin],false)) && (ftemp[i] <= AsinSqrt(zcutBins[selectedBin+1],false)))
                  deleted = false;
               // Event is below the zenith angle cut
               else if(ftemp[i] <= AsinSqrt(zcutBins[selectedBin],false))
                  uoflow[0] = true;
               // Event is above the zenith angle cut
               else if(ftemp[i] > AsinSqrt(zcutBins[selectedBin+1],false))
                  uoflow[1] = true;
	    }
	 }

	 if(!deleted)
	 {
            evtcount[0]++;
	    uoflow[0] = false;
	    uoflow[1] = false;

	    if(DBGSIG > 1)
	       cout << "\tgood";
	 }
	 else
	 {
	    if(!uoflow[0] && !uoflow[1])
	    {
               evtcount[3]++;

	       if(DBGSIG > 1)
	          cout << "\tnovalues";
	    }
	    else if(uoflow[0])
	    {
               evtcount[1]++;

	       if(DBGSIG > 1)
	          cout << "\tbelow";
	    }
	    else if(uoflow[1])
	    {
               evtcount[2]++;

	       if(DBGSIG > 1)
	          cout << "\tabove";
	    }
	 }

	 if(DBGSIG > 1)
	    cout << endl;
      }

      stemp[0] = stemp[0] + stemp[1] + ":\n" + "- Number of total entries = " + ToString(tempTree->GetEntries()) + "\n- Number of selected entries = " + ToString(evtcount[0]) + "\n- Total dropped entries = " + ToString(tempTree->GetEntries() - evtcount[0]) + "\n- Number of entries outside of zenith bin = " + ToString(evtcount[1]) + ", " + ToString(evtcount[2]) + "\n- Number of entries with no valid values = " + ToString(evtcount[3]) + "\n\n";

      delete[] evtcount;
      delete[] uoflow;
      delete ftempaver;
      delete avercount;
   }

   ifile->Close();

   InfoPopup("Check zenith angle bins", stemp[0]);

   delete[] stemp;
   delete[] ftemp;
}

// Check the amount of events inside the selected energy, zenith angle and risetime bins
void MyFrame::CheckBothBins(wxCommandEvent& event)
{
   // Determine the type of observables to cut on (SD or FD)
   selcuttype = (cutObservables->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckBothBins         #: " << "selcuttype = " << selcuttype << endl;
   // Determine how eye selection should be handled (any eye inside selection or average)
   seleyetype = (eyeSelection->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckBothBins         #: " << "seleyetype = " << seleyetype << endl;

   string *stemp = new string[2];
   stemp[0] = string((selectedMva->widgetTE)->GetLineText(0));

   // Open the file for MVA analysis input
   if(stemp[0] == "")
   {
      AlertPopup("No selected MVA file", "No MVA input file selected. Please select one of the files from the third listbox.");
      delete[] stemp;
      return;
   }

   // Stop execution if all cuts and bins are disabled
   if( (!(cutEnergy->widgetChBox)->IsChecked()) && (!(cutZenith->widgetChBox)->IsChecked()) && (!(cutRisetime->widgetChBox)->IsChecked()) )
   {
      AlertPopup("No cut options enabled", "All cut options are disabled, meaning all events survive the binning cuts. If energy, zenith angle and/or risetime cuts are needed, please enable them by ticking the checkboxes in front of the cuts.");
      delete[] stemp;
      return;
   }

   float *ftemp1, *ftemp2, *ftemp3, *ftemp3_err;
   ftemp1 = new float[ALLEYES];
   ftemp2 = new float[ALLEYES];
   ftemp3 = new float[ALLEYES];
   ftemp3_err = new float[ALLEYES];
   float *ftempaver;
   int *avercount;
   int *evtcount;
   bool *sepcut, *uoflow;
   bool sepcutend = false;
  
   TFile *ifile = TFile::Open(stemp[0].c_str(),"READ");

   int selectedBin[2];
   selectedBin[0] = (cutEnergyBins->widgetCB)->GetSelection();
   selectedBin[1] = (cutZenithBins->widgetCB)->GetSelection();

   stemp[0] = "Number of entries inside cuts:\n";
   for(int i = 1; i <= nrkeys; i++)
   {
      stemp[1] = "TreeS" + ToString(i);
      TTree *tempTree = (TTree*)ifile->Get(stemp[1].c_str());

      // SD observables for cut
      if(selcuttype == 0)
      {
         ret = Find(observables, "energySD");
         if(ret != -1)
            tempTree->SetBranchAddress("energySD", ftemp1);
         else
         {
            AlertPopup("No SD energy observable found", "No SD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energySD.");
            delete[] stemp;
            delete[] ftemp1;
            delete[] ftemp2;
            delete[] ftemp3;
            delete[] ftemp3_err;
	    ifile->Close();
            return;
         }

         ret = Find(observables, "zenithSD");
         if(ret != -1)
            tempTree->SetBranchAddress("zenithSD", ftemp2);
         else
         {
            AlertPopup("No SD zenith angle observable found", "No SD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithSD.");
            delete[] stemp;
            delete[] ftemp1;
            delete[] ftemp2;
            delete[] ftemp3;
            delete[] ftemp3_err;
	    ifile->Close();
            return;
         }
      }
      // FD observables for cut
      else if(selcuttype == 1)
      {
         ret = Find(observables, "energyFD");
         if(ret != -1)
            tempTree->SetBranchAddress("energyFD", ftemp1);
         else
         {
            AlertPopup("No FD energy observable found", "No FD energy observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable energyFD.");
            delete[] stemp;
            delete[] ftemp1;
            delete[] ftemp2;
            delete[] ftemp3;
            delete[] ftemp3_err;
	    ifile->Close();
            return;
         }

         ret = Find(observables, "zenithFD");
         if(ret != -1)
            tempTree->SetBranchAddress("zenithFD", ftemp2);
         else
         {
            AlertPopup("No FD zenith angle observable found", "No FD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithFD.");
            delete[] stemp;
            delete[] ftemp1;
            delete[] ftemp2;
            delete[] ftemp3;
            delete[] ftemp3_err;
	    ifile->Close();
            return;
         }
      }

      ret = Find(observables, "risetimerecalc");
      if(ret != -1)
      {
         tempTree->SetBranchAddress("risetimerecalc", ftemp3);
         tempTree->SetBranchAddress("risetimerecalc_neg", ftemp3_err);
      }
      else
      {
         AlertPopup("No risetime observable found", "No risetime observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable risetimerecalc.");
         delete[] stemp;
         delete[] ftemp1;
         delete[] ftemp2;
         delete[] ftemp3;
         delete[] ftemp3_err;
	 ifile->Close();
	 return;
      }

      evtcount = new int[7];
      evtcount[0] = 0;	// All dropped
      evtcount[1] = 0;	// Below energy
      evtcount[2] = 0;	// Above energy
      evtcount[3] = 0;	// Below zenith
      evtcount[4] = 0;	// Above zenith
      evtcount[5] = 0;	// Above risetime
      evtcount[6] = 0;	// No VEM signal (no risetime)

      sepcut = new bool[3];
      uoflow = new bool[6];

      ftempaver = new float[2];
      avercount = new int[2];

      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

         sepcut[0] = false;	// Energy
         sepcut[1] = false;	// Zenith
         sepcut[2] = false;	// Risetime
         uoflow[0] = false;	// Below energy
         uoflow[1] = false;	// Above energy
         uoflow[2] = false;	// Below zenith
         uoflow[3] = false;	// Above zenith
	 uoflow[4] = false;	// Above risetime
	 uoflow[5] = false;	// No VEM signal (no risetime)
	 sepcutend = false;	// All

	 // Calculate averages for all eyes
	 if(seleyetype == 2)
	 {
	    ftempaver[0] = 0;
	    ftempaver[1] = 0;
	    avercount[0] = 0;
	    avercount[1] = 0;

	    for(int k = 0; k < ALLEYES; k++)
	    {
               if((cutEnergy->widgetChBox)->IsChecked()) 
	       {
                  if(ftemp1[k] != -1)
	          {
                     (ftempaver[0]) += ftemp1[k];
	             (avercount[0])++;
	          }
	       }

               if((cutZenith->widgetChBox)->IsChecked()) 
	       {
                  if(ftemp2[k] != -1)
	          {
                     (ftempaver[1]) += ftemp2[k];
	             (avercount[1])++;
	          }
	       }
	    }

	    if(avercount[0] > 0)
               ftempaver[0] = (ftempaver[0])/(avercount[0]);
	    else
               ftempaver[0] = -1;

	    if(avercount[1] > 0)
               ftempaver[1] = (ftempaver[1])/(avercount[1]);
	    else
               ftempaver[1] = -1;

	    if(DBGSIG > 1)
	       cout << "# CheckBothBins         #: " << "Number of eyes = " << avercount[0] << ", " << avercount[1] << endl;
	 }
	 
	 if(DBGSIG > 1) 
            cout << "# CheckBothBins         #: ";

	 // Determine if event is inside energy bin
	 for(int k = 0; k < ALLEYES; k++)
         {
	    if(DBGSIG > 1)
	    {
               cout << ftemp1[k] << " ";
               cout << ftemp2[k] << " ";
	       cout << TMath::Abs(ftemp3[k]-ftemp3_err[k])/ftemp3[k] << ", ";
	    }

            if(seleyetype == 2)
	    {
               ftemp1[k] = ftempaver[0];
               ftemp2[k] = ftempaver[1];
	    }

            if((cutEnergy->widgetChBox)->IsChecked()) 
	    {
	       if( (ftemp1[k] != -1) )
               {
                  // Event is inside the energy cut
                  if((ftemp1[k] > ecutBins[selectedBin[0]]) && (ftemp1[k] <= ecutBins[selectedBin[0]+1]))
                     sepcut[0] = true;
		  // Event is below the energy cut
	          else if(ftemp1[k] <= ecutBins[selectedBin[0]])
		  {
		     uoflow[0] = true;
		     sepcut[0] = false;
		  }
		  // Event is above the energy cut
	          else if(ftemp1[k] > ecutBins[selectedBin[0]+1])
		  {
		     uoflow[1] = true;
		     sepcut[0] = false;
		  }
               }
	       else
	          sepcut[0] = false;
	    }
	    else
               sepcut[0] = true;
			 
	    if((cutZenith->widgetChBox)->IsChecked())
	    {
	       if( (ftemp2[k] != -1) )
               {
                  // Event is inside the zenith angle cut
                  if((ftemp2[k] > AsinSqrt(zcutBins[selectedBin[1]],false)) && (ftemp2[k] <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
                     sepcut[1] = true;
		  // Event is below the zenith angle cut
	          else if(ftemp2[k] <= AsinSqrt(zcutBins[selectedBin[1]],false))
		  {
		     uoflow[2] = true;
		     sepcut[1] = false;
		  }
		  // Event is above the zenith angle cut
	          else if(ftemp2[k] > AsinSqrt(zcutBins[selectedBin[1]+1],false))
		  {
		     uoflow[3] = true;
		     sepcut[1] = false;
		  }
               }
	       else
		  sepcut[1] = false;
	    }
	    else
               sepcut[1] = true;
			 
	    if((cutRisetime->widgetChBox)->IsChecked())
	    {
	       if( ((ftemp3[k] != -1) && (ftemp3_err[k] != -1)) )
               {
                  // Event is below the risetime cut
                  if(TMath::Abs(ftemp3[k]-ftemp3_err[k])/ftemp3[k] <= (cutRisetime->widgetNE[0])->GetValue())
                     sepcut[2] = true;
		  // Event is above the risetime cut
	          else
		  {
		     uoflow[4] = true;
		     sepcut[2] = false;
		  }
               }
	       else
	       {
                  // Event has no valid VEM signal
                  uoflow[5] = true;
		  sepcut[2] = false;
	       }
	    }
	    else
               sepcut[2] = true;

	    // Only keep eye, if both energy, zenith angle and risetime are inside the cut
	    if(!sepcutend)
  	       sepcutend = sepcut[0] & sepcut[1] & sepcut[2];
	 }

	 if(sepcutend)
	 {
            evtcount[0]++;
            uoflow[0] = false;
	    uoflow[1] = false;
            uoflow[2] = false;
	    uoflow[3] = false;
	    uoflow[4] = false;
	    uoflow[5] = false;

	    if(DBGSIG > 1)
	       cout << "\tgood";
	 }
	 else
	 {
	    if(uoflow[0])
	    {
               evtcount[1]++;

 	       if(DBGSIG > 1)
	          cout << "\tbelow energy";
	    }
	    else if(uoflow[1])
	    {
               evtcount[2]++;

 	       if(DBGSIG > 1)
	          cout << "\tabove energy";
	    }

	    if(uoflow[2])
	    {
               evtcount[3]++;

 	       if(DBGSIG > 1)
	          cout << "\tbelow zenith";
	    }
	    else if(uoflow[3])
	    {
               evtcount[4]++;

 	       if(DBGSIG > 1)
	          cout << "\tabove zenith";
	    }

	    if(uoflow[4])
	    {
               evtcount[5]++;

 	       if(DBGSIG > 1)
	          cout << "\tbelow risetime";
	    }

	    if(uoflow[5])
	    {
               evtcount[6]++;

 	       if(DBGSIG > 1)
	          cout << "\tno vem";
	    }
	 }

 	 if(DBGSIG > 1)
 	    cout << endl;
      }

      stemp[0] = stemp[0] + stemp[1] + ":\n" + "- Number of total entries = " + ToString(tempTree->GetEntries()) + "\n- Number of selected entries = " + ToString(evtcount[0]) + "\n- Total dropped entries = " + ToString(tempTree->GetEntries() - evtcount[0]) + "\n";
      if((cutEnergy->widgetChBox)->IsChecked())
         stemp[0] = stemp[0] + "- Number of entries outside of energy bin = " + ToString(evtcount[1]) + ", " + ToString(evtcount[2]) + "\n";
      if((cutZenith->widgetChBox)->IsChecked())
         stemp[0] = stemp[0] + "- Number of entries outside of zenith angle bin = " + ToString(evtcount[3]) + ", " + ToString(evtcount[4]) + "\n";
      if((cutRisetime->widgetChBox)->IsChecked())
         stemp[0] = stemp[0] + "- Number of entries above the risetime cut = " + ToString(evtcount[5]) + "\n- Number of entries without VEM signal = " + ToString(evtcount[6]) + "\n";
      stemp[0] = stemp[0] + "\n";

      delete[] evtcount;
      delete[] sepcut;
      delete[] uoflow;
      delete[] ftempaver;
      delete[] avercount;
   }

   ifile->Close();

   InfoPopup("Check all bins", stemp[0]);

   delete[] stemp;
   delete[] ftemp1;
   delete[] ftemp2;
   delete[] ftemp3;
   delete[] ftemp3_err;
}

// Update the list of selected observables, when we select new ones
void MyFrame::UpdateObservableSelection(wxCommandEvent& event)
{
   for(int i = 0; i < observables.size(); i++)
      obssel[i] = (selectObservables->widgetLB)->IsSelected(i);

   freshAnalysis = false;
}

// Don't start analysis, but just rewrite the input file as a temporary event file
void MyFrame::CreateTempEventFile(wxCommandEvent& event)
{
   int *nrTreeEvents;

   string tempfile;
   tempfile = (selectedMva->widgetTE)->GetLineText(0);
   cout << "# StartMvaAnalysis      #: " << "Opening file " << tempfile << " to rewrite it into a temporary event file." << endl;

   tempAnalysisFile = RemoveFilename(&tempfile) + "/temporary_event_file.root";
   nrTreeEvents = new int[nrkeys];
   ret = MvaTreeFile(&tempfile, &tempAnalysisFile, nrTreeEvents);
   if(ret == -1)
   {
      AlertPopup("Invalid selection of signal and background trees", "The selected signal or background trees are invalid. Please make sure to correctly select them and that they are present inside the input file.");
      delete[] nrTreeEvents;
      return;
   }

   InfoPopup("Finished creating temporary event file", "Finished creating a temporary event file.");

   delete[] nrTreeEvents;
}

// Start the MVA analysis with the selected options
void MyFrame::StartMvaAnalysis(wxCommandEvent& event)
{
   string *mvafile;
   mvafile = new string[2];
   string *stemp;
   stemp = new string[4];
   int *nrTreeEvents;
   int *itemp;
   itemp = new int[3];
   double *dtemp;

   mvaresults = "Finished applying MVA cut.";

   mvafile[0] = SelectMva();

   if(!mvafile[0].empty())
   {
      itemp[0] = 9;
      itemp[1] = 0;

      stemp[2] = "Currently performing the MVA analysis, please wait for it to finish.";
      ShowProgress(wxT("Performing MVA analysis"), stemp[2], itemp[0]);

/*      // DELETE
      // Deleting all unneeded files
      stemp[2] = "rm -fr " + RemoveFilename(&mvafile[0]) + "/mva_error*.dat " + RemoveFilename(&mvafile[0]) + "/mean " + RemoveFilename(&mvafile[0]) + "/negerror " + RemoveFilename(&mvafile[0]) + "/poserror ";//+ RemoveFilename(&mvafile[0]) + "/weights " + RemoveFilename(&mvafile[0]) + "/temporary_mvatree_file.root " + mvafile[0];
      cout << stemp[2] << endl;
      system(stemp[2].c_str());
//      cin >> stemp[2];
      // DELETE*/

      cout << "# StartMvaAnalysis      #: " << "Opening file " << mvafile[0] << " for MVA analysis." << endl;

      mvafile[1] = (selectedMva->widgetTE)->GetLineText(0);
      tempAnalysisFile = RemoveFilename(&mvafile[0]) + "/temporary_mvatree_file.root";
      nrTreeEvents = new int[nrkeys];
      ret = MvaTreeFile(&mvafile[1], &tempAnalysisFile, nrTreeEvents);
      if(ret == -1)
      {
	 progress->Update(itemp[0]);
         AlertPopup("Invalid selection of signal and background trees", "The selected signal or background trees are invalid. Please make sure to correctly select them and that they are present inside the input file.");
         delete[] mvafile;
         delete[] stemp;
         delete[] nrTreeEvents;
         delete[] itemp;
	 return;
      }
      itemp[1]++;
      progress->Update(itemp[1]);

      TFile *ifile = TFile::Open(tempAnalysisFile.c_str(), "READ");

      // Open the file to write out to
      TFile *ofile = TFile::Open(mvafile[0].c_str(), "RECREATE");
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
      itemp[1]++;
      progress->Update(itemp[1]);

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
      itemp[1]++;
      progress->Update(itemp[1]);

      // Select signal and background trees (from the temporary input file)
      TTree *signalTree, *backgroundTree[mixednum];
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
/*            itemp[1]++;
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
/*         if( string((backgroundSelect->widgetCB)->GetStringSelection()) == stemp[1] )
	 {
            cout << "# StartMvaAnalysis      #: " << "Using background tree: " << stemp[1] << endl;
	    backgroundTree = (TTree*)ifile->Get(stemp[0].c_str());
            nrTreeEvents[1] = backgroundTree->GetEntries();
            itemp[1]++;
            progress->Update(itemp[1]);
	 }*/
      }

      itemp[1]++;
      progress->Update(itemp[1]);

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
      itemp[1]++;
      progress->Update(itemp[1]);

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
         ifile->Close();
         delete factory;
         ofile->Close();
         AlertPopup("Invalid MVA method", "The selected MVA method is invalid. Please make sure it is correctly defined.");
         delete[] mvafile;
         delete[] stemp;
         delete[] nrTreeEvents;
         delete[] itemp;
	 return;
      }
      itemp[1]++;
      progress->Update(itemp[1]);

      // Train the selected methods and save them to the weights folder
      factory->TrainAllMethods();
      itemp[1]++;
      progress->Update(itemp[1]);
      // Test the selected methods by applying the trained data to the test data set -> outputs saved to TestTree output file and then to the output ROOT file
      factory->TestAllMethods();
      itemp[1]++;
      progress->Update(itemp[1]);
      // Evaluation of methods printed to stdout
      factory->EvaluateAllMethods();
      itemp[1]++;
      progress->Update(itemp[1]);

      // Close the open files
      ifile->Close();
      delete factory;
      ofile->Close();

      // Copy the training values from results/transformation_stats.dat to the current analysis directory
      stemp[3] = "cp -r " + string(rootdir) + "/results/transformation_stats.dat " + (*currentAnalysisDir) + "/";
      system(stemp[3].c_str());

      // Get MVA training values and variable correlations
      sigCorMat = new TMatrixD(nrselobs, nrselobs);
      backCorMat = new TMatrixD(nrselobs, nrselobs);
      ret = GetTrainingShift(&mvafile[0]);
      for(int i = 0; i < nrkeys; i++)
         otherCorMat[i] = new TMatrixD(nrselobs, nrselobs);
      GetApplyCorrelations(&tempAnalysisFile);

      // Open the MVA GUI to review the training and testing procedure
      if((specialMva->widgetChBox[0])->IsChecked())
      {
         TString *inname = new TString;
	 *inname = (TString)mvafile[0];

	 MvaEfficiency *effplot = new MvaEfficiency(nrTreeEvents[0], nrTreeEvents[1], currentAnalysisDir);
	 effplot->RunMvaEfficiency(inname);

         stemp[2] = "Finished running MVA analysis. For Signal/Background values of: [" + ToString(nrTreeEvents[0]) + "/" + ToString(nrTreeEvents[1]) + "] best cuts are (sig/bgd/pur/SNR):\n";
	 stemp[2] = stemp[2] + "- Optimal cut: " + ToString(effplot->optimalCut, 4) + " (" + ToString(100.*(effplot->GetHistValue(0, 0)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(0, 1)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(0, 2)), 2) + "%/" + ToString(effplot->GetHistValue(0, 3), 4) + ")\n";
	 stemp[2] = stemp[2] + "- Cut with equal Signal/Background: " + ToString(effplot->sigbgdCut, 4) + " (" + ToString(100.*(effplot->GetHistValue(1, 0)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(1, 1)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(1, 2)), 2) + "%/" + ToString(effplot->GetHistValue(1, 3), 4) + ")\n";
	 stemp[2] = stemp[2] + "- Cut with equal Signal/Purity: " + ToString(effplot->sigpurCut, 4) + " (" + ToString(100.*(effplot->GetHistValue(2, 0)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(2, 1)), 2) + "%/" + ToString(100.*(effplot->GetHistValue(2, 2)), 2) + "%/" + ToString(effplot->GetHistValue(2, 3), 4) + ")\n";

         NEDialog cutMvaDialog(wxT("MVA cut"), wxSize(600,200), stemp[2], "Set MVA cut:", effplot->sigbgdCut, &ID_MVACUTDIALOG);
         cutMvaDialog.SetNEntryFormat(cutMvaDialog.widgetNE, 4, 0.0001, 0, 0., 0.);
	 if(cutMvaDialog.ShowModal() == wxID_OK)
            (cutMva->widgetNE[0])->SetValue(cutMvaDialog.GetNEValue());
	 else
	 {
            delete[] mvafile;
            delete[] stemp;
            delete[] nrTreeEvents;
            delete[] itemp;
	    delete effplot;
	    delete inname;
	    return;
	 }

	 delete effplot;
	 delete inname;
/*         stemp[2] = "Finished running MVA analysis. Currently the cut is set to: " + ToString((cutMva->widgetNE[0])->GetValue(), 4) + "\nMake sure that the MVA cut is correct. Check for correct cuts by plotting \"(5a) Classifier Cut Efficiency\" at S/B values of: [" + ToString(nrTreeEvents[0]) + "/" + ToString(nrTreeEvents[1]) + "].\n";
	 InfoPopup("MVA cut", stemp[2]);

	 stemp[2] = string(rootdir) + "/bin/tmvagui " + mvafile[0];
	 system(stemp[2].c_str());

         NEDialog cutMvaDialog(wxT("MVA cut"), wxSize(250,110), "Set MVA cut:", (cutMva->widgetNE[0])->GetValue(), &ID_MVACUTDIALOG);
         cutMvaDialog.SetNEntryFormat(cutMvaDialog.widgetNE, 4, 0.0001, 0, 0., 0.);
	 if(cutMvaDialog.ShowModal() == wxID_OK)
            (cutMva->widgetNE[0])->SetValue(cutMvaDialog.GetNEValue());*/
      }

      stemp[2] = (methodsSelect->widgetCB)->GetStringSelection();
      applymva = GetMethodName(stemp[2]);

      if( applymva == "All" )
      {
         AlertPopup("Multiple MVA methods", "Multiple MVA methods selected. To continue applying MVA cuts, please select only one method and rerun the analysis.");
         delete[] mvafile;
         delete[] stemp;
         delete[] nrTreeEvents;
         delete[] itemp;
	 return;
      }

      dtemp = new double[3];
      dtemp[0] = (cutMva->widgetNE[0])->GetValue();
      dtemp[1] = 0;
      dtemp[2] = 0;
      mvaresults = mvaresults + "\n\nResults for mean cut (" + ToString(dtemp[0], 4) + ")\n";
      MvaApplication(&tempAnalysisFile, false, 0);

      freshAnalysis = true;

      GetMvaError((dataSelect->widgetCB)->GetSelection()+1, dtemp);
      (cutMva->widgetNE[0])->SetValue(dtemp[1]);
      mvaresults = mvaresults + "\nResults for negative error cut (" + ToString(dtemp[1], 4) + ")\n";
      MvaApplication(&tempAnalysisFile, freshAnalysis, -1);

      (cutMva->widgetNE[0])->SetValue(dtemp[2]);
      mvaresults = mvaresults + "\nResults for positive error cut (" + ToString(dtemp[2], 4) + ")\n";
      MvaApplication(&tempAnalysisFile, freshAnalysis, 1);

      (cutMva->widgetNE[0])->SetValue(dtemp[0]);
      InfoPopup("Finished applying MVA cut", mvaresults);

      cout << "# StartMvaAnalysis      #: " << "MVA results:" << endl << mvaresults << endl;

      delete[] dtemp;
      delete[] nrTreeEvents;

      delete sigCorMat;
      delete backCorMat;
      for(int i = 0; i < nrkeys; i++)
         delete otherCorMat[i];
   }

   delete[] mvafile;
   delete[] stemp;
   delete[] itemp;
}

// Apply the MVA cut
void MyFrame::ApplyMvaCut(wxCommandEvent& event)
{
   mvaresults = "Finished applying MVA cut.";

   if(freshAnalysis)
   {
      double *dtemp;
      dtemp = new double[3];
      dtemp[0] = (cutMva->widgetNE[0])->GetValue();
      dtemp[1] = 0;
      dtemp[2] = 0;
      mvaresults = mvaresults + "\n\nResults for mean cut (" + ToString(dtemp[0], 4) + ")\n";
      MvaApplication(&tempAnalysisFile, freshAnalysis, 0);

      GetMvaError((dataSelect->widgetCB)->GetSelection()+1, dtemp);
      (cutMva->widgetNE[0])->SetValue(dtemp[1]);
      mvaresults = mvaresults + "\nResults for negative error cut (" + ToString(dtemp[1], 4) + ")\n";
      MvaApplication(&tempAnalysisFile, freshAnalysis, -1);

      (cutMva->widgetNE[0])->SetValue(dtemp[2]);
      mvaresults = mvaresults + "\nResults for positive error cut (" + ToString(dtemp[2], 4) + ")\n";
      MvaApplication(&tempAnalysisFile, freshAnalysis, 1);

      (cutMva->widgetNE[0])->SetValue(dtemp[0]);
      InfoPopup("Finished applying MVA cut", mvaresults);

      cout << "# ApplyMvaCut           #: " << "MVA results:" << endl << mvaresults << endl;

      delete[] dtemp;
   }
   else
      AlertPopup("Can't apply MVA cut", "Can not apply MVA cut. Please run a fresh analysis first (Start MVA analysis).");
}

// Application of MVA analysis
void MyFrame::MvaApplication(string *infilename, bool application, int mean)
{
   int *itemp;
   itemp = new int[4];
   string *stemp;
   stemp = new string[5];

   // Open a reader and add variables (must be the same as for the training)
   TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
   
   float *obsvars;
   obsvars = new float[(3*nrobs)+1];
   
   itemp[0] = 0;

   // Adding observables to the MVA reader
   for(int i = 0; i < nrobs; i++)
   {
      if(obssel[i])
      {
	 if(DBGSIG > 0)
            cout << "# MvaApplication        #: " << "Selected observable: " << observables[i] << endl;
         reader->AddVariable(observables[i].c_str(), &obsvars[3*itemp[0]]);
	 if(DBGSIG > 0)
            cout << "# MvaApplication        #: " << "Tracked variable: obsvars[3*" << itemp[0] << "] = obsvars[" << 3*itemp[0] << "]" << endl;
         itemp[0]++;
      }
   }

   // Book the MVA with the produced weights file
   stemp[0] = (*currentAnalysisDir) + "/weights/TMVAClassification_" + applymva + ".weights.xml";
   if(DBGSIG > 0)
      cout << "# MvaApplication        #: " << "Weights file: " << stemp[0] << endl;
   stemp[1] = applymva + " method";
   if(DBGSIG > 0)
      cout << "# MvaApplication        #: " << "MVA method: " << stemp[1] << endl;
   reader->BookMVA(stemp[1].c_str(), stemp[0].c_str());

   // Open the input file and prepare the TTree
   TFile *ifile;
   if(mean == 0)
   {
      stemp[3] = "cp -r " + *infilename + " " + RemoveFilename(infilename) + "/mvatree_file.root";
      if(DBGSIG > 0)
         cout << "# MvaApplication        #: " << "Copying temporary file: " << stemp[3] << endl;
      ret = system(stemp[3].c_str());
      stemp[3] = RemoveFilename(infilename) + "/mvatree_file.root";
      ifile = TFile::Open(stemp[3].c_str(), "UPDATE");
   }
   else
      ifile = TFile::Open(infilename->c_str(), "READ");


   if(!(cutResults.empty()))
      cutResults.erase(cutResults.begin(), cutResults.end());

   TList *tempkeyslist = (TList*) ifile->GetListOfKeys();
   for(int j = 1; j <= ifile->GetNkeys(); j++)
   {
      stemp[2] = string((tempkeyslist->At(j-1))->GetName());
      cout << "# MvaApplication        #: " << endl << stemp[2] << " has been selected for evaluation." << endl;

      TTree *signalApp;
      signalApp = (TTree*)ifile->Get(stemp[2].c_str());

      itemp[1] = 0;
      for(int i = 0; i < nrobs; i++)
      {
         if(obssel[i])
         {
            cout << "# MvaApplication        #: " << "Setting observable " << observables[i] << endl;
            signalApp->SetBranchAddress((observables[i]).c_str(), &obsvars[3*itemp[1]]); // mean
            if(DBGSIG > 0)
               cout << "# MvaApplication        #: " << "Mean variable: obsvars[3*" << itemp[1] << "] = obsvars[" << 3*itemp[1] << "]" << endl;
            signalApp->SetBranchAddress((observables[i] + "_neg").c_str(), &obsvars[3*itemp[1]+1]); // neg error
            if(DBGSIG > 0)
               cout << "# MvaApplication        #: " << "Negerror variable: obsvars[3*" << itemp[1] << "+1] = obsvars[" << 3*itemp[1]+1 << "]" << endl;
            signalApp->SetBranchAddress((observables[i] + "_pos").c_str(), &obsvars[3*itemp[1]+2]); // pos error
            if(DBGSIG > 0)
               cout << "# MvaApplication        #: " << "Poserror variable: obsvars[3*" << itemp[1] << "+1] = obsvars[" << 3*itemp[1]+2 << "]" << endl;
            itemp[1]++;
         }
      }

      CreateOutput(signalApp, reader, stemp[1], obsvars, stemp[2], j, application, mean);
   }

   if(mean == 0)
   {
      ifile->Write();
      itemp[3] = ifile->GetNkeys()/2;
   }
   else
      itemp[3] = ifile->GetNkeys();

   itemp[2] = 0;
   for(int j = 0; j < itemp[3]; j++)
   {
      if(mean == 0)
         stemp[4] = (tempkeyslist->At(2*j))->GetName();
      else
         stemp[4] = (tempkeyslist->At(j))->GetName();
      mvaresults = mvaresults + " - " + stemp[4] + " = " + ToString(cutResults[itemp[2]]) + " vs. " + ToString(cutResults[itemp[2]+1]) + " (" + ToString(100.*((double)cutResults[itemp[2]])/((double)cutResults[itemp[2]]+(double)cutResults[itemp[2]+1]), 2) + "%, " + ToString(100.*((double)cutResults[itemp[2]+1])/((double)cutResults[itemp[2]]+(double)cutResults[itemp[2]+1]), 2) + "%)\n";
      itemp[2] += 2;
   }

   ifile->Close();

   delete reader;
   delete[] itemp;
   delete[] obsvars;
   delete[] stemp;
}

// Set the default MVA options
void MyFrame::SetDefaultMva(wxCommandEvent& event)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsorigsel[i])
         (selectObservables->widgetLB)->SetSelection(i);
      else
         (selectObservables->widgetLB)->Deselect(i);

      obssel[i] = obsorigsel[i];
   }
   (selectObservables->widgetLB)->SetFirstItem(0);

   (methodsSelect->widgetCB)->SetSelection((methodsSelect->widgetCB)->FindString(wxT("Neural network (MLPBNN)")));
   (cutMva->widgetNE[0])->SetValue(0.);
   (cutObservables->widgetCB)->SetSelection(1);
   (cutEnergy->widgetChBox)->SetValue(true);
   (cutEnergy->widgetNE[0])->SetValue(17.);
   (cutEnergy->widgetNE[1])->SetValue(21.);
   (cutEnergyBins->widgetNE[0])->SetValue(1);
   RunEnergyBinSelect();
   (cutZenith->widgetChBox)->SetValue(true);
   (cutZenith->widgetNE[0])->SetValue(0.);
   (cutZenith->widgetNE[1])->SetValue(60.);
   (cutZenithBins->widgetNE[0])->SetValue(1);
   RunZenithBinSelect();
   (cutRisetime->widgetChBox)->SetValue(true);
   (cutRisetime->widgetNE[0])->SetValue(0.3);
   (eyeSelection->widgetCB)->SetSelection(0);

   (specialMva->widgetChBox[0])->SetValue(true);

   freshAnalysis = false;
}

// Rewrite the observables from an ADST file into a root file
int MyFrame::StartRewrite(string *outfile)
{
   // Check for all selections in the first listbox
   wxArrayInt selections;
   (mvaList[0]->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
   {
      // Prepare the AdstMva object, signal & sum observables, and a tree for all (must be outside the for loop)
      AdstMva *mvatool = new AdstMva;
      mvatool->outname = *outfile;

      Observables *obssig[3];
      Observables *obsall[3];

      for(int i = 0; i < 3; i++)
      {
         obssig[i] = new Observables(observables);
         obsall[i] = new Observables(observables);
      }

      mvatool->outfile = TFile::Open((mvatool->outname).c_str(), "RECREATE");
      mvatool->all_tree = new TTree("TreeA", "Background tree with all events, including signal events.");
      mvatool->rewritecode = 0;

      for(int i = 0; i < selections.GetCount(); i++)
      {
         mvatool->inname = (mvaList[0]->widgetLB)->GetString(selections[i]);
	 ShowProgress(wxT("Rewriting observables"), wxT("Currently rewriting observables, please wait for it to finish."), 100);
	 ret = mvatool->RewriteObservables(selections.GetCount(), i, obssig, obsall, progress);
	 if(ret == -1)
	 {
            for (int j = 0; j < 3; j++)
            {
               delete obssig[j];
               delete obsall[j];
            }

            delete mvatool;
            AlertPopup("No events in selected input file", "File selected inside first listbox has no events. Please change input file and select one or more files (holding Ctrl or Shift while clicking) to rewrite.");

            return -1;
	 }
	 mvatool->PrepareOtherTrees(selections.GetCount(), i, observables);
      }

      (mvatool->all_tree)->Write();
      (mvatool->outfile)->Close();

      for (int i = 0; i < 3; i++)
      {
         delete obssig[i];
         delete obsall[i];
      }

      delete mvatool;

      return 0;
   }
   else
   {
      AlertPopup("No selected files", "No files from the first listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to rewrite.");
      return -1;
   }
}

// Select a file to split depending on the fraction of events we need -> and write it back to rewritten ADST format
int MyFrame::StartFileSplit(string infile)
{
   float *ftemp;
   int *itemp;
   string *stemp;
   float *obsvars;
   int rewritecode;

   stemp = new string[4];
   stemp[0] = RemoveExtension(&infile) + "_split-1.root";
   stemp[1] = RemoveExtension(&infile) + "_split-2.root";

   ftemp = new float[2];
   itemp = new int[4];
   
   ftemp[0] = (startCombining->widgetNE[0])->GetValue();
   
   TFile *input = TFile::Open(infile.c_str(), "READ");
   TTree *tempTree = (TTree*)input->Get("TreeA");
   TList *tempkeyslist = (TList*)input->GetListOfKeys();
   TTree *readTree;
   TTree *writeTree;
   
   itemp[0] = tempTree->GetEntries();
   itemp[1] = (int)TMath::Ceil(itemp[0]*ftemp[0]);
   itemp[2] = (int)(itemp[0]-itemp[1]);

/*   Observables *obstemp = new Observables(observables);
   for(int j = 0; j < nrobs; j++)
      tempTree->SetBranchAddress((obstemp->GetName(j)).c_str(), obstemp->obsstruct[j].value);
   for(int i = 0; i < 40; i++)
   {
      cout << "# i = " << i << endl;
      tempTree->GetEntry(i);
      ftemp[1] = 0;
      for(int j = 0; j < ALLEYES; j++)
      {
         cout << "EYE " << j << ": " << obstemp->GetValue("energyFD", j) << endl;
         ftemp[1] += obstemp->GetValue("energyFD", j);
      }
      cout << "Testing value: " << ftemp[1] << endl;
   }

//   return 1;*/
   
   if( (itemp[0] == 0) || (itemp[1] == 0) || (itemp[2] == 0) )
   {
      AlertPopup("No events in split file", "One of the split files will have no events (" + ToString(itemp[1]) + "/" + ToString(itemp[2]) + "). Please adjust the fraction accordingly and restart.");

      delete[] ftemp;
      delete[] itemp;
      delete[] stemp;

      return 1;
   }
   else
   {
      itemp[3] = 0;
      stemp[2] = "Currently splitting file \"" + RemovePath(&infile) + "\" into two files: \"" + RemovePath(&stemp[0]) + "\" and \"" + RemovePath(&stemp[1])+ "\". Please wait for it to finish.";
      ShowProgress(wxT("Splitting rewritten ADST file"), stemp[2].c_str(), 2*(input->GetNkeys())*itemp[0]);

      cout << "# Will start splitting file " << RemovePath(&infile) << " into:" << endl;
      cout << "# - First file:  " << RemovePath(&stemp[0]) << " (" << itemp[1] << " of " << itemp[0] << " total events)" << endl;
      cout << "# - Second file: " << RemovePath(&stemp[1]) << " (" << itemp[2] << " of " << itemp[0] << " total events)" << endl;
   
      // Preparing shuffled list for sampling
      vector<int> shuflist;
      for(int i = 0; i < itemp[0]; i++)
         shuflist.push_back(i);
   
      shuffle(shuflist.begin(), shuflist.end(), default_random_engine((unsigned int)chrono::system_clock::now().time_since_epoch().count()));

      vector<int> split1list;
      vector<int> split2list;

      for(int i = 0; i < itemp[0]; i++)
      {
         if(i < itemp[1])
            split1list.push_back(shuflist[i]);
	 else
            split2list.push_back(shuflist[i]);
      }

      // Loop over the two split files
      TFile *output;
      for(int i = 0; i < 2; i++)
      {
         cout << "# Currently writing out to: " << stemp[i] << endl;
         output = TFile::Open(stemp[i].c_str(), "RECREATE");

         // Loop over all trees
         for(int k = 0; k < input->GetNkeys(); k++)
         {
            // Prepare observables for reading and writing
            Observables *obser = new Observables(observables);
            Observables *obser_neg = new Observables(observables);
            Observables *obser_pos = new Observables(observables);

            // Name and title of the current tree
            stemp[2] = string((tempkeyslist->At(k))->GetName());
            stemp[3] = string((tempkeyslist->At(k))->GetTitle());

            cout << "#   Currently selected tree: " << stemp[2] << "; " << stemp[3] << endl;

            // Tree for reading
            readTree = (TTree*)input->Get(stemp[2].c_str());
            readTree->SetBranchAddress("rewritecode", &rewritecode);
            for(int j = 0; j < nrobs; j++)
            {
               readTree->SetBranchAddress((obser->GetName(j)).c_str(), obser->obsstruct[j].value);
               readTree->SetBranchAddress((obser->GetName(j) + "_neg").c_str(), obser_neg->obsstruct[j].value);
               readTree->SetBranchAddress((obser->GetName(j) + "_pos").c_str(), obser_pos->obsstruct[j].value);
            }

            // Tree for writing
            writeTree = new TTree(stemp[2].c_str(), stemp[3].c_str());
            writeTree->Branch("rewritecode", &rewritecode, "rewritecode/I");
            for(int j = 0; j < nrobs; j++)
            {
               writeTree->Branch((obser->GetName(j)).c_str(), &(obser->obsstruct[j].value), (obser->GetName(j) + "[" + ToString(ALLEYES) + "]/F").c_str());
               writeTree->Branch((obser->GetName(j) + "_neg").c_str(), &(obser_neg->obsstruct[j].value), (obser->GetName(j) + "_neg[" + ToString(ALLEYES) + "]/F").c_str());
               writeTree->Branch((obser->GetName(j) + "_pos").c_str(), &(obser_pos->obsstruct[j].value), (obser->GetName(j) + "_pos[" + ToString(ALLEYES) + "]/F").c_str());
            }

	    itemp[2] = 0;
            // Loop over all events 
            for(int j = 0; j < itemp[0]; j++)
            {
               readTree->GetEntry(j);

	       // Check if the values in this tree are valid
	       ftemp[1] = 0;
               for(int j = 0; j < ALLEYES; j++)
                  ftemp[1] += obser->GetValue("energyFD", j);

	       if( (ftemp[1] == -4) || (ftemp[1] > 1.e+15) )
	       {
                  // Select the correct tree to write to
                  if( (find(split1list.begin(), split1list.end(), j) != split1list.end()) && (i == 0) )
                     writeTree->Fill();
                  if( (find(split2list.begin(), split2list.end(), j) != split2list.end()) && (i == 1) )
                     writeTree->Fill();
	       }
	       else
		  itemp[2]++;

	       // Update the progress bar
	       itemp[3]++;
	       if(itemp[3]%((int)(2*(input->GetNkeys())*itemp[0]*0.05)) == 0)
                  progress->Update(itemp[3]);
            }

	    cout << "#   There were " << itemp[2] << " invalid events." << endl;

	    writeTree->Write();

	    delete obser;
	    delete obser_neg;
	    delete obser_pos;
	    delete writeTree;
         }

	 output->Close();
      }

      progress->Update(2*(input->GetNkeys())*itemp[0]);
   }

   delete[] ftemp;
   delete[] itemp;
   delete[] stemp;

   return 0;
}

// Combine the observables from multiple ADST files into a single root file
int MyFrame::StartCombine(string *outfile)
{
   // Check for all selections in the second listbox
   wxArrayInt selections;
   (mvaList[1]->widgetLB)->GetSelections(selections);

   if(selections.GetCount() > 1)
   {
      string *flist = new string[selections.GetCount()];
      for(int i = 0; i < selections.GetCount(); i++)
         flist[i] = (mvaList[1]->widgetLB)->GetString(selections[i]);

      // Combine the files with the hadd function
      hadd(selections.GetCount(), flist, outfile);
      InfoPopup("Finished combining files", "The selected files have been combined and can be used as input files to the MVA analysis.");

      delete[] flist;
      return 0;
   }
   else
   {
      AlertPopup("No selected files", "At least two files must be selected from the second listbox. Please select at least two files (holding Ctrl or Shift while clicking) to combine them to a single rewritten file.");
      return -1;
   }
}

// Merge the observables from one or multiple ADST files into a single root file
int MyFrame::StartMerge(string *outfile)
{
   // Check for all selections in the second listbox
   wxArrayInt selections;
   (mvaList[1]->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
   {
      string *flist = new string[selections.GetCount()];
      for(int i = 0; i < selections.GetCount(); i++)
         flist[i] = (mvaList[1]->widgetLB)->GetString(selections[i]);

      // Combine the files with the hadd function
      hmerge(selections.GetCount(), flist, outfile);
      InfoPopup("Finished merging files", "The selected files have been merged together.");

      delete[] flist;
      return 0;
   }
   else
   {
      AlertPopup("No selected files", "At least two files must be selected from the second listbox. Please select at least two files (holding Ctrl or Shift while clicking) to combine them to a single rewritten file.");
      return -1;
   }
}
