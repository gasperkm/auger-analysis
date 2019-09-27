#include "frame.h"
#include "separate_functions.h"
#include "adst_mva.h"
#include "popups.h"
#include "combine.h"

#include <chrono>
#include <random>

// Double click on list to select a MVA file for input
void MyFrame::EnableMvaFile(wxCommandEvent& event)
{
   int *itemp;
   int *cnt;
   string *stemp;
   string *signalName;
   vector<string> *backgroundType;

   if( (mvaList[2]->widgetLB)->GetSelection() != wxNOT_FOUND )
   {
      itemp = new int;
      stemp = new string[2];
      signalName = new string;
      backgroundType = new vector<string>;
      cnt = new int;
      *cnt = 0;

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
      nrkeys = GetRootKeys(ifile, "TreeS");
      if(DBGSIG > 1)
         cout << "# EnableMvaFile         #: " << "Starting number of keys " << nrkeys << endl;
      for(int i = 1; i <= nrkeys; i++)
      {
         *signalName = "TreeS" + ToString(i);
	 stemp[1] = string(ifile->GetKey(signalName->c_str())->GetTitle());
	 (signalSelect->widgetCB)->Append(stemp[1]);
	 (backgroundSelect->widgetCB)->Append(stemp[1]);
	 (dataSelect->widgetCB)->Append(stemp[1]);
         (*cnt)++;

	 if(FindStringPart(stemp[1], "Data") == 0)
 	    backgroundType->push_back(stemp[1]);
      }
/*      keyslist = (TList*) ifile->GetListOfKeys();
      for(int i = 1; i <= nrkeys; i++)
      {
         if(strcmp((keyslist->At(i-1))->GetName(), "TreeA") != 0)
	 {
            *signalName = "TreeS" + ToString(i);
	    stemp[1] = string(ifile->GetKey(signalName->c_str())->GetTitle());
	    (signalSelect->widgetCB)->Append(stemp[1]);
	    (backgroundSelect->widgetCB)->Append(stemp[1]);
	    (dataSelect->widgetCB)->Append(stemp[1]);
            (*cnt)++;

	    if(FindStringPart(stemp[1], "Data") == 0)
 	       backgroundType->push_back(stemp[1]);
	 }
      }*/

      nrkeys = *cnt;
      cout << "# EnableMvaFile         #: " << "Found " << nrkeys << " signal keys in the selected file." << endl;

      mixednum = backgroundType->size();
      if(mixednum > 2)
      {
         for(int i = 0; i < mixednum; i++)
	 {
            for(int j = 0; j < mixednum; j++)
	    {
               if(i < j)
	       {
                  stemp[1] = backgroundType->at(i) + " + " + backgroundType->at(j);
                  (backgroundSelect->widgetCB)->Append(stemp[1]);
	       }
	    }
	 }
      }
      if(mixednum > 3)
      {
         for(int i = 0; i < mixednum; i++)
	 {
            for(int j = 0; j < mixednum; j++)
	    {
               for(int k = 0; k < mixednum; k++)
	       {
                  if( (i < j) && (j < k) )
	          {
                     stemp[1] = backgroundType->at(i) + " + " + backgroundType->at(j) + " + " + backgroundType->at(k);
                     (backgroundSelect->widgetCB)->Append(stemp[1]);
	          }
	       }
	    }
	 }
      }

      // Select the signal tree
      if(oldselect[0] >= *cnt)
      {
         (signalSelect->widgetCB)->SetSelection((*cnt)-1);
	 oldselect[0] = (*cnt)-1;
      }
      else
         (signalSelect->widgetCB)->SetSelection(oldselect[0]);

      // Select the background tree
      if(oldselect[1] >= *cnt)
      {
         (backgroundSelect->widgetCB)->SetSelection((*cnt)-1);
	 oldselect[1] = (*cnt)-1;
      }
      else
         (backgroundSelect->widgetCB)->SetSelection(oldselect[1]);

      // Select the data tree
      if(oldselect[2] >= *cnt)
      {
         (dataSelect->widgetCB)->SetSelection((*cnt)-1);
	 oldselect[2] = (*cnt)-1;
      }
      else
         (dataSelect->widgetCB)->SetSelection(oldselect[2]);

      ifile->Close();

      delete itemp;
      delete[] stemp;
      delete signalName;
      delete cnt;
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
   string *stemp = new string;

   if(!customBinning)
   {
      elimits[0] = cutEnergy->GetNumber(cutEnergy->widgetNE[0]);
      elimits[1] = cutEnergy->GetNumber(cutEnergy->widgetNE[1]);

      // Check the energy range (in logarithmic scale)
      *erange = (elimits[1] - elimits[0])/((double)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]));

      ecutBins.clear();
      (cutEnergyBins->widgetCB)->Clear();

      // Save energy cut bins (in linear scale)
      for(int i = 0; i < (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
      {
         ecutBins.push_back((double)pow(10, i*(*erange)+elimits[0]));
         ecutBins.push_back((double)pow(10, (i+1)*(*erange)+elimits[0]));
         if(DBGSIG > 0)
            cout << "# RunEnergyBinSelect    #: " << i << ", energy bin = " << ecutBins[2*i] << ", " << ecutBins[2*i+1] << endl;
      }

      // Add all bins to the selection dropbox for energy
      for(int i = 0; i < (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
      {
         if(i < 9)
            *stemp = "0" + ToString(i+1) + " (" + ToSciString(ecutBins[2*i], 2) + " - " + ToSciString(ecutBins[2*i+1], 2) + ")";
         else
            *stemp = ToString(i+1) + " (" + ToSciString(ecutBins[2*i], 2) + " - " + ToSciString(ecutBins[2*i+1], 2) + ")";

         (cutEnergyBins->widgetCB)->Append(*stemp);
      }

      (cutEnergyBins->widgetCB)->SetSelection(0);

      oldBinSettings[0] = elimits[0];
      oldBinSettings[1] = elimits[1];
      oldBinSettings[2] = (double)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]);
   }
   else
   {
      (cutEnergyBins->widgetCB)->Clear();

      for(int i = 0; i < (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
      {
         if(DBGSIG > 0)
            cout << "# RunEnergyBinSelect    #: " << i << ", energy bin = " << ecutBins[2*i] << ", " << ecutBins[2*i+1] << endl;
      }

      // Add all bins to the selection dropbox for energy
      for(int i = 0; i < (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
      {
         if(i < 9)
            *stemp = "0" + ToString(i+1) + " (" + ToSciString(ecutBins[2*i], 2) + " - " + ToSciString(ecutBins[2*i+1], 2) + ")";
         else
            *stemp = ToString(i+1) + " (" + ToSciString(ecutBins[2*i], 2) + " - " + ToSciString(ecutBins[2*i+1], 2) + ")";

         (cutEnergyBins->widgetCB)->Append(*stemp);
      }

      (cutEnergyBins->widgetCB)->SetSelection(0);
   }

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
   ftemp = new float;
   int *avercount;
   int *evtcount;
   bool *deleted;
   bool *uoflow;

   deleted = new bool;
   *deleted = true;
  
   TFile *ifile = TFile::Open(stemp[0].c_str(),"READ");

   int *selectedBin = new int;
   *selectedBin = (cutEnergyBins->widgetCB)->GetSelection();

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
            delete ftemp;
	    ifile->Close();
	    delete deleted;
	    delete selectedBin;
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
            delete ftemp;
	    ifile->Close();
	    delete deleted;
	    delete selectedBin;
            return;
         }
      }

      evtcount = new int[4];
      evtcount[0] = 0;
      evtcount[1] = 0;
      evtcount[2] = 0;
      evtcount[3] = 0;

      uoflow = new bool[2];

      avercount = new int;

      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

	 *deleted = true;
	 uoflow[0] = false;
	 uoflow[1] = false;

	 if(DBGSIG > 1)
            cout << "# CheckEnergyBin        #: ";
	 
	 // Determine if event is inside energy bin
	 for(int i = 0; i < ALLEYES; i++)
	 {
            if(DBGSIG > 1)
               cout << *ftemp << " ";

	    if((*ftemp != -1) && *deleted)
            {
               // Event is inside the energy cut
	       if((*ftemp > ecutBins[2*(*selectedBin)]) && (*ftemp <= ecutBins[2*(*selectedBin+1)]))
                  *deleted = false;
               // Event is below the energy cut
               else if(*ftemp <= ecutBins[2*(*selectedBin)])
                  uoflow[0] = true;
               // Event is above the energy cut
               else if(*ftemp > ecutBins[2*(*selectedBin+1)])
                  uoflow[1] = true;
	    }
	 }

	 if(!(*deleted))
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
      delete avercount;
   }

   ifile->Close();

   InfoPopup("Check energy bins", stemp[0]);

   delete deleted;
   delete selectedBin;
   delete[] stemp;
   delete ftemp;
}
// Use a file to define custom energy bins
void MyFrame::CustomEnergyBins(wxCommandEvent& event)
{
   string *tempstring;
   int *retval;
   double *dtemp;
   int *itemp;

   wxFileDialog *customBinsDlg = new wxFileDialog(this, wxT("Open custom energy bins file"), string(rootdir) + string("/input"), "", wxT("Txt files (*.txt)|*.txt"), wxFD_OPEN);

   if(customBinsDlg->ShowModal() == wxID_OK)
   {
      tempstring = new string;
      *tempstring = customBinsDlg->GetPath();

      if(!tempstring->empty())
      {
         (cutEnergyBinsCustom->widgetTE)->SetValue(*tempstring);
	 (cutEnergy->widgetNE[0])->Disable();
	 (cutEnergy->widgetNE[1])->Disable();
	 (cutEnergyBins->widgetNE[0])->Disable();
         (splitCutEnergy->widgetNE[0])->Disable();
         (splitCutEnergy->widgetNE[1])->Disable();
         (splitCutEnergy->widgetNE[2])->Disable();
         (splitCutEnergy->widgetChBox)->SetValue(1);
	 customBinning = true;
         ecutBins.clear();

         ifstream *ifs = new ifstream;
	 dtemp = new double[4];
	 itemp = new int[2];

         ifs->open(tempstring->c_str(), ifstream::in);

	 dtemp[2] = 1e+20;
	 dtemp[3] = -1;

	 if(ifs->is_open())
	 {
            itemp[0] = 0;
	    cout << "Writing out custom binning:" << endl;
            while(1)
	    {
               *ifs >> dtemp[0] >> dtemp[1];

	       if(dtemp[0] < dtemp[2])
                  dtemp[2] = dtemp[0];
	       if(dtemp[1] > dtemp[3])
                  dtemp[3] = dtemp[1];

               ecutBins.push_back((double)pow(10, dtemp[0]));
               ecutBins.push_back((double)pow(10, dtemp[1]));
               esplitBins.push_back((double)pow(10, dtemp[0]));
               esplitBins.push_back((double)pow(10, dtemp[1]));

               itemp[0]++;
	       ifs->ignore(1,' ');
               if(ifs->eof()) break;
	    }

	    cout << "Number of bins = " << itemp[0] << ", Minimum energy = " << dtemp[2] << ", Maximum energy = " << dtemp[3] << endl;
	    (cutEnergy->widgetNE[0])->SetValue(dtemp[2]);
	    (cutEnergy->widgetNE[1])->SetValue(dtemp[3]);
	    (cutEnergyBins->widgetNE[0])->SetValue(itemp[0]);
	    (splitCutEnergy->widgetNE[0])->SetValue(dtemp[2]);
	    (splitCutEnergy->widgetNE[1])->SetValue(dtemp[3]);
	    (splitCutEnergy->widgetNE[2])->SetValue(itemp[0]);
	    RunEnergyBinSelect();
	    RunEnergyBinSplitSelect();
	 }

	 ifs->close();

	 delete[] itemp;
	 delete[] dtemp;
	 delete ifs;
      }

      delete tempstring;
   }

   delete customBinsDlg;
}

// Disable custom energy binning
void MyFrame::CustomEnergyBinsDisable(wxCommandEvent& event)
{
   (cutEnergyBinsCustom->widgetTE)->SetValue("");
   (cutEnergy->widgetNE[0])->Enable();
   (cutEnergy->widgetNE[1])->Enable();
   (cutEnergyBins->widgetNE[0])->Enable();
   (cutEnergy->widgetNE[0])->SetValue(oldBinSettings[0]);
   (cutEnergy->widgetNE[1])->SetValue(oldBinSettings[1]);
   (cutEnergyBins->widgetNE[0])->SetValue(oldBinSettings[2]);
   (splitCutEnergy->widgetNE[0])->Enable();
   (splitCutEnergy->widgetNE[1])->Enable();
   (splitCutEnergy->widgetNE[2])->Enable();
   (splitCutEnergy->widgetNE[2])->SetValue(1);
   (splitCutEnergy->widgetChBox)->SetValue(0);
   customBinning = false;
   RunEnergyBinSelect();
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
//   asinlimits[0] = SinSquare(cutZenith->GetNumber(cutZenith->widgetNE[0]), true);
//   asinlimits[1] = SinSquare(cutZenith->GetNumber(cutZenith->widgetNE[1]), true);
   asinlimits[0] = SecTheta(cutZenith->GetNumber(cutZenith->widgetNE[0]), true);
   asinlimits[1] = SecTheta(cutZenith->GetNumber(cutZenith->widgetNE[1]), true);

   string *stemp = new string;

   // Check the zenith angle range (in sin square scale)
   *arange = (asinlimits[1] - asinlimits[0])/((double)cutZenithBins->GetNumber(cutZenithBins->widgetNE[0]));

   zcutBins.clear();
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
//      if(i < 9)
//         *stemp = "0" + ToString(i+1) + " (" + ToString(AsinSqrt(zcutBins[i], true), 1) + " - " + ToString(AsinSqrt(zcutBins[i+1], true), 1) + ")";
//      else
//         *stemp = ToString(i+1) + " (" + ToString(AsinSqrt(zcutBins[i], true), 1) + " - " + ToString(AsinSqrt(zcutBins[i+1], true), 1) + ")";
      if(i < 9)
         *stemp = "0" + ToString(i+1) + " (" + ToString(InvSecTheta(zcutBins[i], true), 1) + " - " + ToString(InvSecTheta(zcutBins[i+1], true), 1) + ")";
      else
         *stemp = ToString(i+1) + " (" + ToString(InvSecTheta(zcutBins[i], true), 1) + " - " + ToString(InvSecTheta(zcutBins[i+1], true), 1) + ")";

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
   ftemp = new float;
   int *avercount;
   int *evtcount;
   bool *deleted;
   bool *uoflow;

   deleted = new bool;
   *deleted = true;
  
   TFile *ifile = TFile::Open(stemp[0].c_str(),"READ");

   int *selectedBin = new int;
   *selectedBin = (cutZenithBins->widgetCB)->GetSelection();

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
            delete ftemp;
	    ifile->Close();
	    delete deleted;
	    delete selectedBin;
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
            delete ftemp;
	    ifile->Close();
	    delete deleted;
	    delete selectedBin;
            return;
         }
      }

      evtcount = new int[4];
      evtcount[0] = 0;
      evtcount[1] = 0;
      evtcount[2] = 0;
      evtcount[3] = 0;

      uoflow = new bool[2];

      avercount = new int;

      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);

	 *deleted = true;
	 uoflow[0] = false;
	 uoflow[1] = false;

	 if(DBGSIG > 1) 
            cout << "# CheckZenithBin        #: ";

	 // Determine if event is inside zenith angle bin
	 for(int i = 0; i < ALLEYES; i++)
	 {
	    if(DBGSIG > 1)
               cout << *ftemp << " ";

	    if( (*ftemp != -1) && *deleted )
            {
               // Event is inside the zenith angle cut
	       if((*ftemp > InvSecTheta(zcutBins[*selectedBin],false)) && (*ftemp <= InvSecTheta(zcutBins[*selectedBin+1],false)))
                  *deleted = false;
               // Event is below the zenith angle cut
               else if(*ftemp <= InvSecTheta(zcutBins[*selectedBin],false))
                  uoflow[0] = true;
               // Event is above the zenith angle cut
               else if(*ftemp > InvSecTheta(zcutBins[*selectedBin+1],false))
                  uoflow[1] = true;
//               // Event is inside the zenith angle cut
//               if((*ftemp > AsinSqrt(zcutBins[*selectedBin],false)) && (*ftemp <= AsinSqrt(zcutBins[*selectedBin+1],false)))
//                  *deleted = false;
//               // Event is below the zenith angle cut
//               else if(*ftemp <= AsinSqrt(zcutBins[*selectedBin],false))
//                  uoflow[0] = true;
//               // Event is above the zenith angle cut
//               else if(*ftemp > AsinSqrt(zcutBins[*selectedBin+1],false))
//                  uoflow[1] = true;
	    }
	 }

	 if(!(*deleted))
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
      delete avercount;
   }

   ifile->Close();

   InfoPopup("Check zenith angle bins", stemp[0]);

   delete deleted;
   delete selectedBin;
   delete[] stemp;
   delete ftemp;
}

// Check the amount of events inside the selected energy, zenith angle and risetime bins
void MyFrame::CheckBothBins(wxCommandEvent& event)
{
   // Determine the type of observables to cut on (SD or FD)
   selcuttype = (cutObservables->widgetCB)->GetSelection();
   if(DBGSIG > 0)
      cout << "# CheckBothBins         #: " << "selcuttype = " << selcuttype << endl;

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
   ftemp1 = new float;
   ftemp2 = new float;
   ftemp3 = new float;
   ftemp3_err = new float;
   int *avercount;
   int *evtcount;
   bool *sepcut, *uoflow;
   bool *sepcutend;
  
   TFile *ifile = TFile::Open(stemp[0].c_str(),"READ");

   int *selectedBin = new int[2];
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
            delete ftemp1;
            delete ftemp2;
            delete ftemp3;
            delete ftemp3_err;
	    ifile->Close();
	    delete[] selectedBin;
            return;
         }

         ret = Find(observables, "zenithSD");
         if(ret != -1)
            tempTree->SetBranchAddress("zenithSD", ftemp2);
         else
         {
            AlertPopup("No SD zenith angle observable found", "No SD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithSD.");
            delete[] stemp;
            delete ftemp1;
            delete ftemp2;
            delete ftemp3;
            delete ftemp3_err;
	    ifile->Close();
	    delete[] selectedBin;
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
            delete ftemp1;
            delete ftemp2;
            delete ftemp3;
            delete ftemp3_err;
	    ifile->Close();
	    delete[] selectedBin;
            return;
         }

         ret = Find(observables, "zenithFD");
         if(ret != -1)
            tempTree->SetBranchAddress("zenithFD", ftemp2);
         else
         {
            AlertPopup("No FD zenith angle observable found", "No FD zenith angle observable found in the list of observables (" + string(rootdir) + "/input/observables.txt). Please name the observable zenithFD.");
            delete[] stemp;
            delete ftemp1;
            delete ftemp2;
            delete ftemp3;
            delete ftemp3_err;
	    ifile->Close();
	    delete[] selectedBin;
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
         delete ftemp1;
         delete ftemp2;
         delete ftemp3;
         delete ftemp3_err;
	 ifile->Close();
	 delete[] selectedBin;
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

      sepcutend = new bool;
      *sepcutend = false;

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
	 *sepcutend = false;	// All
	 
	 if(DBGSIG > 1) 
            cout << "# CheckBothBins         #: ";

	 // Determine if event is inside energy bin
	 for(int k = 0; k < ALLEYES; k++)
         {
	    if(DBGSIG > 1)
	    {
               cout << *ftemp1 << " ";
               cout << *ftemp2 << " ";
	       cout << TMath::Abs((*ftemp3)-(*ftemp3_err))/(*ftemp3) << ", ";
	    }

            if((cutEnergy->widgetChBox)->IsChecked()) 
	    {
	       if( (*ftemp1 != -1) )
               {
                  // Event is inside the energy cut
                  if((*ftemp1 > ecutBins[2*(selectedBin[0])]) && (*ftemp1 <= ecutBins[2*(selectedBin[0])+1]))
                     sepcut[0] = true;
		  // Event is below the energy cut
	          else if(*ftemp1 <= ecutBins[2*(selectedBin[0])])
		  {
		     uoflow[0] = true;
		     sepcut[0] = false;
		  }
		  // Event is above the energy cut
	          else if(*ftemp1 > ecutBins[2*(selectedBin[0])+1])
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
	       if( (*ftemp2 != -1) )
               {
                  // Event is inside the zenith angle cut
                  if((*ftemp2 > InvSecTheta(zcutBins[selectedBin[1]],false)) && (*ftemp2 <= InvSecTheta(zcutBins[selectedBin[1]+1],false)))
                     sepcut[1] = true;
		  // Event is below the zenith angle cut
	          else if(*ftemp2 <= InvSecTheta(zcutBins[selectedBin[1]],false))
		  {
		     uoflow[2] = true;
		     sepcut[1] = false;
		  }
		  // Event is above the zenith angle cut
	          else if(*ftemp2 > InvSecTheta(zcutBins[selectedBin[1]+1],false))
		  {
		     uoflow[3] = true;
		     sepcut[1] = false;
		  }
//                  // Event is inside the zenith angle cut
//                  if((*ftemp2 > AsinSqrt(zcutBins[selectedBin[1]],false)) && (*ftemp2 <= AsinSqrt(zcutBins[selectedBin[1]+1],false)))
//                     sepcut[1] = true;
//        	  // Event is below the zenith angle cut
//                  else if(*ftemp2 <= AsinSqrt(zcutBins[selectedBin[1]],false))
//        	  {
//        	     uoflow[2] = true;
//        	     sepcut[1] = false;
//        	  }
//        	  // Event is above the zenith angle cut
//                  else if(*ftemp2 > AsinSqrt(zcutBins[selectedBin[1]+1],false))
//        	  {
//        	     uoflow[3] = true;
//        	     sepcut[1] = false;
//        	  }
               }
	       else
		  sepcut[1] = false;
	    }
	    else
               sepcut[1] = true;
			 
	    if((cutRisetime->widgetChBox)->IsChecked())
	    {
	       if( ((*ftemp3 != -1) && (*ftemp3_err != -1)) )
               {
                  // Event is below the risetime cut
                  if(TMath::Abs((*ftemp3)-(*ftemp3_err))/(*ftemp3) <= (cutRisetime->widgetNE[0])->GetValue())
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
	    if(!(*sepcutend))
  	       *sepcutend = sepcut[0] & sepcut[1] & sepcut[2];
	 }

	 if(*sepcutend)
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

      delete sepcutend;
      delete[] evtcount;
      delete[] sepcut;
      delete[] uoflow;
      delete[] avercount;
   }

   ifile->Close();

   InfoPopup("Check all bins", stemp[0]);

   delete[] stemp;
   delete[] selectedBin;
   delete ftemp1;
   delete ftemp2;
   delete ftemp3;
   delete ftemp3_err;
}

// Update the list of selected observables, when we select new ones
void MyFrame::UpdateObservableSelection(wxCommandEvent& event)
{
   for(int i = 0; i < observables.size(); i++)
      obssel[i] = (selectObservables->widgetLB)->IsSelected(i);

   freshAnalysis = false;
}

// Set the negative uncertainty shift
void MyFrame::SetNegativeUncertainty(wxCommandEvent& event)
{
   (specialMva->widgetChBox[5])->SetValue(true);
   (specialMva->widgetChBox[6])->SetValue(false);
}

// Set the positive uncertainty shift
void MyFrame::SetPositiveUncertainty(wxCommandEvent& event)
{
   (specialMva->widgetChBox[5])->SetValue(false);
   (specialMva->widgetChBox[6])->SetValue(true);
}

// Disable uncertainty shift
void MyFrame::DisableUncertainty(wxCommandEvent& event)
{
   (specialMva->widgetChBox[5])->SetValue(false);
   (specialMva->widgetChBox[6])->SetValue(false);
}

// Handle events concerning lists
void MyFrame::EditList(wxCommandEvent& event)
{
   int *itemp;
   string *stemp;

   itemp = new int[3];
   stemp = new string[2];

//   cout << "# EditList: " << event.GetId() << endl;
   for(int i = 0; i < nrlists; i++)
   {
      itemp[0] = (allLBE[i]->widgetLB)->GetCount();
//      cout << "List has " << itemp[0] << " events" << endl;
      (allLBE[i]->widgetLB)->GetSelections(selections);
      itemp[1] = selections.GetCount();

      if(event.GetId() == (1001 + 4*i))	// delete the selection
      {
//         cout << "Clicked the Delete command, " << i << endl;

	 if(!selections.IsEmpty())
	 {
	    for(int j = itemp[1]-1; j >= 0; j--)
	    {
//	       cout << "Deleting entry " << selections[j] << endl;
	       (allLBE[i]->widgetLB)->Delete(selections[j]);
	    }

	    if((allLBE[i]->widgetLB)->GetCount() > 0)
	    {
               if(selections[0] > 0)
	       {
                  (allLBE[i]->widgetLB)->SetSelection(selections[0]-1);
//	          cout << "Selecting entry " << selections[0]-1 << endl;
	       }
	       else
	       {
                  (allLBE[i]->widgetLB)->SetSelection(0);
//	          cout << "Selecting entry " << 0 << endl;
	       }
	    }
	 }
	 else
	 {
            AlertPopup("No selection", "There were no entries selected in the listbox.");
            delete[] itemp;
            delete[] stemp;
            return;
	 }

	 break;
      }
      else if(event.GetId() == (1002 + 4*i))	// move up the list
      {
//         cout << "Clicked the Move Up command, " << i << endl;

	 if(!selections.IsEmpty())
	 {
	    for(int j = 0; j < itemp[1]; j++)
	    {
	       if(selections[0] > 0)
	       {
//	          cout << "Moving entry " << selections[j] << endl;
                  stemp[0] = (allLBE[i]->widgetLB)->GetString(selections[j]-1);
                  stemp[1] = (allLBE[i]->widgetLB)->GetString(selections[j]);

	          (allLBE[i]->widgetLB)->SetString(selections[j]-1, stemp[1]);
	          (allLBE[i]->widgetLB)->SetString(selections[j], stemp[0]);
//	          cout << "Up to " << selections[j]-1 << endl;

                  (allLBE[i]->widgetLB)->SetSelection(selections[j]-1);
		  (allLBE[i]->widgetLB)->Deselect(selections[j]);
	       }
	       else
	       {
                  AlertPopup("Entry at top", "One of the selected entries is already at the top of the listbox. Please reselect.");
                  delete[] itemp;
                  delete[] stemp;
                  return;
	       }
	    }
	 }
	 else
	 {
            AlertPopup("No selection", "There were no entries selected in the listbox.");
            delete[] itemp;
            delete[] stemp;
            return;
	 }

	 break;
      }
      else if(event.GetId() == (1003 + 4*i))	// move down the list
      {
//         cout << "Clicked the Move Down command, " << i << endl;

	 if(!selections.IsEmpty())
	 {
	    for(int j = itemp[1]-1; j >= 0; j--)
	    {
	       if(selections[itemp[1]-1] < itemp[0]-1)
	       {
//	          cout << "Moving entry " << selections[j] << endl;
                  stemp[0] = (allLBE[i]->widgetLB)->GetString(selections[j]);
                  stemp[1] = (allLBE[i]->widgetLB)->GetString(selections[j]+1);

	          (allLBE[i]->widgetLB)->SetString(selections[j], stemp[1]);
	          (allLBE[i]->widgetLB)->SetString(selections[j]+1, stemp[0]);
//	          cout << "Down to " << selections[j]+1 << endl;

                  (allLBE[i]->widgetLB)->SetSelection(selections[j]+1);
		  (allLBE[i]->widgetLB)->Deselect(selections[j]);
	       }
	       else
	       {
                  AlertPopup("Entry at bottom", "One of the selected entries is already at the bottom of the listbox. Please reselect.");
                  delete[] itemp;
                  delete[] stemp;
                  return;
	       }
	    }
	 }
	 else
	 {
            AlertPopup("No selection", "There were no entries selected in the listbox.");
            delete[] itemp;
            delete[] stemp;
            return;
	 }

	 break;
      }
      else if(event.GetId() == (1004 + 4*i))	// clear the list
      {
//         cout << "Clicked the Clear command " << i << endl;
	 (allLBE[i]->widgetLB)->Clear();

	 break;
      }
   }

   delete[] itemp;
   delete[] stemp;
}

// Don't start analysis, but just rewrite the input file as a temporary event file
void MyFrame::CreateTempEventFile(wxCommandEvent& event)
{
   int *nrTreeEvents;

   string *tempfile = new string[2];
   tempfile[0] = (selectedMva->widgetTE)->GetLineText(0);
   cout << "# CreateTempEventFile   #: " << "Opening file " << tempfile[0] << " to rewrite it into a temporary event file." << endl;

   tempfile[1] = RemovePath(&tempfile[0]);
   tempAnalysisFile = RemoveFilename(&tempfile[0]) + "/temporary_event_file_" + tempfile[1];
   nrTreeEvents = new int[nrkeys];
   ret = MvaTreeFile(&tempfile[0], &tempAnalysisFile, nrTreeEvents, 0);
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
   string *stemp;
   int *nrTreeEvents;
   int *itemp;
   double *dtemp;

   if( ((specialMva->widgetChBox[0])->IsChecked()) && ((cutEnergyBins->widgetNE[0])->GetValue() == 1) )
   {
      AlertPopup("Automatic run failure", "If automatic run is selected, there should be more than 1 energy bins. Please adjust one of them accordingly.");
      return;
   }

   if( ((signalSelect->widgetCB)->GetSelection() == (backgroundSelect->widgetCB)->GetSelection()) || ((signalSelect->widgetCB)->GetSelection() == (dataSelect->widgetCB)->GetSelection()) || ((backgroundSelect->widgetCB)->GetSelection() == (dataSelect->widgetCB)->GetSelection()) )
   {
      AlertPopup("Wrong tree selections", "There is an error in tree selections. Any two from signal tree select, background tree select and data tree select combo boxes should not be equal. Always select different trees as signal, background and data.");
      return;
   }

   mvafile = new string[3];
   stemp = new string[4];
   itemp = new int[3];

   mvaresults = "Finished applying MVA cut.";

   mvafile[0] = SelectMva();

   if((cutEnergyBins->widgetNE[0])->GetValue() > 1)
   {
      if((specialMva->widgetChBox[0])->IsChecked())
         mvafile[2] = mvafile[0];

      multipleEnergyBins = true;
   }
   else
      multipleEnergyBins = false;

   // Printout the special settings for this run
   cout << "specialMva settings:" << endl
        << "- Automatically run over energy bins = " << (int)(specialMva->widgetChBox[0])->IsChecked() << endl
	<< "- Manually select MVA cut = " << (int)(specialMva->widgetChBox[1])->IsChecked() << endl
	<< "- Determine MVA cut error = " << (int)(specialMva->widgetChBox[2])->IsChecked() << endl
	<< "- Apply FD standard correction = " << (int)(specialMva->widgetChBox[3])->IsChecked() << endl
	<< "- Apply HECO correction = " << (int)(specialMva->widgetChBox[4])->IsChecked() << endl
	<< "- Apply negative uncertainty = " << (int)(specialMva->widgetChBox[5])->IsChecked() << endl
	<< "- Apply positive uncertainty = " << (int)(specialMva->widgetChBox[6])->IsChecked() << endl
	<< "- Use published S38 values = " << (int)(specialMva->widgetChBox[7])->IsChecked() << endl
	<< "- Apply smearing to MC trees = " << (int)(specialMva->widgetChBox[8])->IsChecked() << endl
	<< "- Treat selected data as MC = " << (int)(specialMva->widgetChBox[9])->IsChecked() << endl;

   // Run over all energy bins automatically
   for(int imva = 0; imva < (cutEnergyBins->widgetNE[0])->GetValue(); imva++)
   {
      // Only perform analysis if input rewritten ADST file exists
      if(!mvafile[0].empty())
      {
	 // Prepare directory structure and naming for outputs of different energy bins
         if(((cutEnergyBins->widgetNE[0])->GetValue() > 1) && ((specialMva->widgetChBox[0])->IsChecked()))
	 {
            mvafile[0] = mvafile[2];

            stemp[2] = RemoveExtension(&mvafile[0]);
            stemp[3] = ToString(TMath::Log10(ecutBins[2*imva]),2) + "-" + ToString(TMath::Log10(ecutBins[2*imva+1]),2);
	    stemp[2] = stemp[2] + "_" + stemp[3];
            *currentAnalysisDir = stemp[2];
            stemp[2] = "mkdir " + stemp[2];
            cout << "Creating directory: " << *currentAnalysisDir << endl;
            system(stemp[2].c_str());

            stemp[2] = RemoveExtension(&mvafile[0]);
            stemp[3] = ToString(TMath::Log10(ecutBins[2*imva]),2) + "-" + ToString(TMath::Log10(ecutBins[2*imva+1]),2);
	    stemp[2] = stemp[2] + "_" + stemp[3];
            stemp[2] = CheckExtension(&stemp[2], "root");
	    mvafile[0] = (*currentAnalysisDir) + "/" + RemovePath(&stemp[2]);
            cout << "Running analysis with output: " << mvafile[0] << endl;

	    (cutEnergyBins->widgetCB)->SetSelection(imva);
	 }

         itemp[0] = 9;
         itemp[1] = 0;

         stemp[2] = "Currently performing the MVA analysis, please wait for it to finish.";
         ShowProgress(wxT("Performing MVA analysis"), stemp[2], itemp[0]);

	 // Open input file and rewrite it into a temporary file temporary_mvatree_file.root
         cout << "# StartMvaAnalysis      #: " << "Opening file " << mvafile[0] << " for MVA analysis." << endl;

         mvafile[1] = (selectedMva->widgetTE)->GetLineText(0);
         tempAnalysisFile = RemoveFilename(&mvafile[0]) + "/temporary_mvatree_file.root";
         nrTreeEvents = new int[nrkeys];
         ret = MvaTreeFile(&mvafile[1], &tempAnalysisFile, nrTreeEvents, imva);
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

/*-------------TODO-------------*/
	 ret = PerformMvaAnalysis(&tempAnalysisFile, &mvafile[0], &itemp[1]);
         if(ret == -1)
	 {
            progress->Update(itemp[0]);
            delete[] mvafile;
            delete[] stemp;
            delete[] nrTreeEvents;
            delete[] itemp;
	    return;
	 }
/*-------------TODO-------------*/
/*---Errors: PerformMvaAnalysis(&tempAnalysisFile, &mvafile[0], 1)--------*/
/*-------------TODO-------------*/

/*         // Copy the training values from results/transformation_stats.dat to the current analysis directory
         stemp[3] = "cp -r " + string(rootdir) + "/results/transformation_stats.dat " + (*currentAnalysisDir) + "/";
         system(stemp[3].c_str());*/

	 /* NEWREMOVE - TODO
         // Get MVA training values and variable correlations
         sigCorMat = new TMatrixD(nrselobs, nrselobs);
         backCorMat = new TMatrixD(nrselobs, nrselobs);
         ret = GetTrainingShift(&mvafile[0]);
         for(int i = 0; i < nrkeys; i++)
            otherCorMat[i] = new TMatrixD(nrselobs, nrselobs);
         GetApplyCorrelations(&tempAnalysisFile);
	 */

/*	 // Skip the GUI interface for best cut and automatically select signal/background
         if(((cutEnergyBins->widgetNE[0])->GetValue() > 1) && ((specialMva->widgetChBox[0])->IsChecked()))
	 {
            TString *inname = new TString;
            *inname = (TString)mvafile[0];

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
               *inname = (TString)mvafile[0];

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
                  delete[] mvafile;
                  delete[] stemp;
                  delete[] nrTreeEvents;
                  delete[] itemp;
                  delete effplot;
                  delete inname;
                  delete sigCorMat;
                  delete backCorMat;
                  for(int i = 0; i < nrkeys; i++)
                     delete otherCorMat[i];
		  cutMvaDialog->Destroy();
		  delete cutMvaDialog;
                  return;
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
            delete[] mvafile;
            delete[] stemp;
            delete[] nrTreeEvents;
            delete[] itemp;
            delete sigCorMat;
            delete backCorMat;
            for(int i = 0; i < nrkeys; i++)
               delete otherCorMat[i];
            return;
         }*/

         dtemp = new double[3];
         dtemp[0] = (cutMva->widgetNE[0])->GetValue();
         dtemp[1] = 0;
         dtemp[2] = 0;
         mvaresults = mvaresults + "\n\nResults for MVA cut (" + ToString(dtemp[0], 4) + ")\n";
	 int *selectedBin = new int;
         *selectedBin = (cutEnergyBins->widgetCB)->GetSelection();
         mvaprintout = ToSciString(ecutBins[2*(*selectedBin)], 4) + "\t" + ToSciString(ecutBins[2*(*selectedBin)+1], 4) + "\n";
	 /* NEWREMOVE - TODO
         mvaprintout = mvaprintout + ToString(0) + "\t" + ToString(dtemp[0],4) + "\t";
	 */
         mvaprintout = mvaprintout + ToString(dtemp[0],4) + "\t";
         MvaApplication(&tempAnalysisFile, false, 0);

         freshAnalysis = true;

	 /* NEWREMOVE - TODO
         GetMvaError((dataSelect->widgetCB)->GetSelection()+1, dtemp, &tempAnalysisFile);
	 
         (cutMva->widgetNE[0])->SetValue(dtemp[1]);
         mvaresults = mvaresults + "\nResults for negative error cut (" + ToString(dtemp[1], 4) + ")\n";
         mvaprintout = mvaprintout + ToString(-1) + "\t" + ToString(dtemp[1],4) + "\t";
         MvaApplication(&tempAnalysisFile, freshAnalysis, -1);

         (cutMva->widgetNE[0])->SetValue(dtemp[2]);
         mvaresults = mvaresults + "\nResults for positive error cut (" + ToString(dtemp[2], 4) + ")\n";
         mvaprintout = mvaprintout + ToString(1) + "\t" + ToString(dtemp[2],4) + "\t";
         MvaApplication(&tempAnalysisFile, freshAnalysis, 1);
	 */

         (cutMva->widgetNE[0])->SetValue(dtemp[0]);

         // Show the results popup (only if running manually)
         if(((cutEnergyBins->widgetNE[0])->GetValue() == 1) || (!(specialMva->widgetChBox[0])->IsChecked()))
            InfoPopup("Finished applying MVA cut", mvaresults);

         cout << "# StartMvaAnalysis      #: " << "MVA results:" << endl << mvaresults << endl;
         cout << "# StartMvaAnalysis      #: " << "MVA printout:" << endl << mvaprintout << endl;

         ofstream *printoutFile = new ofstream;
         printoutFile->open(((*currentAnalysisDir) + "/application_results.txt").c_str(), ofstream::out | ofstream::trunc );
         *printoutFile << mvaprintout;
         printoutFile->close();

	 delete printoutFile;

         delete[] dtemp;
         delete[] nrTreeEvents;
	 delete selectedBin;

         delete sigCorMat;
         delete backCorMat;
         for(int i = 0; i < nrkeys; i++)
            delete otherCorMat[i];
      }
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
      int *selectedBin = new int;
      *selectedBin = (cutEnergyBins->widgetCB)->GetSelection();
      mvaprintout = ToSciString(ecutBins[2*(*selectedBin)], 4) + "\t" + ToSciString(ecutBins[2*(*selectedBin)+1], 4) + "\n";
      mvaprintout = mvaprintout + ToString(0) + "\t" + ToString(dtemp[0],4) + "\t";
      MvaApplication(&tempAnalysisFile, freshAnalysis, 0);

      /* NEWREMOVE - TODO
      GetMvaError((dataSelect->widgetCB)->GetSelection()+1, dtemp, &tempAnalysisFile);
      */
      (cutMva->widgetNE[0])->SetValue(dtemp[1]);
      mvaresults = mvaresults + "\nResults for negative error cut (" + ToString(dtemp[1], 4) + ")\n";
      mvaprintout = mvaprintout + ToString(-1) + "\t" + ToString(dtemp[1],4) + "\t";
      MvaApplication(&tempAnalysisFile, freshAnalysis, -1);

      (cutMva->widgetNE[0])->SetValue(dtemp[2]);
      mvaresults = mvaresults + "\nResults for positive error cut (" + ToString(dtemp[2], 4) + ")\n";
      mvaprintout = mvaprintout + ToString(1) + "\t" + ToString(dtemp[2],4) + "\t";
      MvaApplication(&tempAnalysisFile, freshAnalysis, 1);

      (cutMva->widgetNE[0])->SetValue(dtemp[0]);
      InfoPopup("Finished applying MVA cut", mvaresults);

      cout << "# ApplyMvaCut           #: " << "MVA results:" << endl << mvaresults << endl;
      cout << "# ApplyMvaCut           #: " << "MVA printout:" << endl << mvaprintout << endl;

      ofstream *printoutFile = new ofstream;
      printoutFile->open(((*currentAnalysisDir) + "/application_results.txt").c_str(), ofstream::out | ofstream::trunc );
      *printoutFile << mvaprintout;
      printoutFile->close();

      delete printoutFile;
      delete[] dtemp;
   }
   else
      AlertPopup("Can't apply MVA cut", "Can not apply MVA cut. Please run a fresh analysis first (Start MVA analysis).");
}

// Application of MVA analysis
void MyFrame::MvaApplication(string *infilename, bool application, int mean)
{
   int *itemp;
   string *stemp;

   itemp = new int[4];
   stemp = new string[6];

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

   // Open the output file and prepare the TTree
   TFile *ofile;
   if(mean == 0)
   {
      stemp[3] = "cp -r " + *infilename + " " + RemoveFilename(infilename) + "/mvatree_file.root";
      if(DBGSIG > 0)
         cout << "# MvaApplication        #: " << "Copying temporary file: " << stemp[3] << endl;
      ret = system(stemp[3].c_str());
      stemp[3] = RemoveFilename(infilename) + "/mvatree_file.root";
      ofile = TFile::Open(stemp[3].c_str(), "UPDATE");
   }
   else
      ofile = TFile::Open(infilename->c_str(), "READ");

   cutResults.clear();

   TList *tempkeyslist = (TList*) ofile->GetListOfKeys();
   mvaprintout = mvaprintout + ToString(ofile->GetNkeys()) + "\n";
   for(int j = 1; j <= ofile->GetNkeys(); j++)
   {
      stemp[2] = string((tempkeyslist->At(j-1))->GetName());
      cout << "# MvaApplication        #: " << stemp[2] << " has been selected for evaluation." << endl;

      stemp[4] = string((tempkeyslist->At(j-1))->GetTitle());

      cout << stemp[4] << "\t" << j-1 << endl;
      cout << string((signalSelect->widgetCB)->GetStringSelection()) << "\t";
      cout << (signalSelect->widgetCB)->GetSelection() << endl;
      cout << string((backgroundSelect->widgetCB)->GetStringSelection()) << "\t";
      cout << (backgroundSelect->widgetCB)->GetSelection() << endl;
      cout << string((dataSelect->widgetCB)->GetStringSelection()) << "\t";
      cout << (dataSelect->widgetCB)->GetSelection() << endl;

      if( (signalSelect->widgetCB)->GetSelection() == j-1 )
         mvaprintout = mvaprintout + ToString(1) + "\t";
      else if( (backgroundSelect->widgetCB)->GetSelection() == j-1 )
         mvaprintout = mvaprintout + ToString(2) + "\t";
      else if( (dataSelect->widgetCB)->GetSelection() == j-1 )
         mvaprintout = mvaprintout + ToString(3) + "\t";
      else
         mvaprintout = mvaprintout + ToString(0) + "\t";

      TTree *signalApp = (TTree*)ofile->Get(stemp[2].c_str());

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
//      ofile->cd();
      ofile->Write();
      itemp[3] = ofile->GetNkeys()/2;
   }
   else
      itemp[3] = ofile->GetNkeys();

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

   ofile->Close();

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

   (methodsSelect->widgetCB)->SetSelection((methodsSelect->widgetCB)->FindString(wxT("Boosted_decision_trees_(BDTG)")));
   (cutMva->widgetNE[0])->SetValue(0.);
   (cutObservables->widgetCB)->SetSelection(1);
   (cutEnergy->widgetChBox)->SetValue(true);
   (cutEnergy->widgetNE[0])->SetValue(17.);
   (cutEnergy->widgetNE[1])->SetValue(21.);
   (cutEnergyBins->widgetNE[0])->SetValue(1);
   (cutEnergyBinsCustom->widgetTE)->SetValue("");
   (cutEnergy->widgetNE[0])->Enable();
   (cutEnergy->widgetNE[1])->Enable();
   (cutEnergyBins->widgetNE[0])->Enable();
   (splitCutEnergy->widgetNE[0])->Enable();
   (splitCutEnergy->widgetNE[1])->Enable();
   (splitCutEnergy->widgetNE[2])->Enable();
   customBinning = false;
   RunEnergyBinSelect();
   (cutZenith->widgetChBox)->SetValue(true);
   (cutZenith->widgetNE[0])->SetValue(0.);
   (cutZenith->widgetNE[1])->SetValue(60.);
   (cutZenithBins->widgetNE[0])->SetValue(1);
   RunZenithBinSelect();
   (cutRisetime->widgetChBox)->SetValue(false);
   (cutRisetime->widgetNE[0])->SetValue(0.3);

   (specialMva->widgetChBox[0])->SetValue(true);
   (specialMva->widgetChBox[1])->SetValue(true);
   (specialMva->widgetChBox[2])->SetValue(true);
   (specialMva->widgetChBox[3])->SetValue(false);
   (specialMva->widgetChBox[4])->SetValue(false);
   (specialMva->widgetChBox[5])->SetValue(false);
   (specialMva->widgetChBox[6])->SetValue(false);
   (specialMva->widgetChBox[7])->SetValue(false);
   (specialMva->widgetChBox[8])->SetValue(true);
   (specialMva->widgetChBox[9])->SetValue(false);

   freshAnalysis = false;
}

// Rewrite the observables from an ADST file into a root file
int MyFrame::StartRewrite(string *outfile)
{
#if ROOTADST == 1
   // Check for all selections in the first listbox
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
#endif
}

void MyFrame::RunEnergyBinSplitSelect()
{
   double *erange = new double;
   double *elimits = new double[2];
   string *stemp = new string;

   if(!customBinning)
   {
      elimits[0] = splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[0]);
      elimits[1] = splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[1]);

      // Check the energy range (in logarithmic scale)
      *erange = (elimits[1] - elimits[0])/((double)splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[2]));

      esplitBins.clear();

      // Save energy cut bins (in linear scale)
      for(int i = 0; i < (int)splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[2]); i++)
      {
         esplitBins.push_back((double)pow(10, i*(*erange)+elimits[0]));
         esplitBins.push_back((double)pow(10, (i+1)*(*erange)+elimits[0]));
         if(DBGSIG > 0)
            cout << "# RunEnergyBinSplitSelect#: " << i << ", energy bin = " << esplitBins[2*i] << ", " << esplitBins[2*i+1] << endl;
      }
   }
   else
   {
      for(int i = 0; i < (int)cutEnergyBins->GetNumber(cutEnergyBins->widgetNE[0]); i++)
      {
         if(DBGSIG > 0)
            cout << "# RunEnergyBinSplitSelect#: " << i << ", energy bin = " << esplitBins[2*i] << ", " << esplitBins[2*i+1] << endl;
      }
   }

   delete stemp;
   delete erange;
   delete[] elimits;
}

// Select a file to split depending on the fraction of events we need -> and write it back to rewritten ADST format
int MyFrame::StartFileSplit(string infile, int cursel, int nrsel)
{
   RunEnergyBinSplitSelect();
   int *nrsplits = new int;
   *nrsplits = splitCutEnergy->GetNumber(splitCutEnergy->widgetNE[2]);

   for(int ebin = 0; ebin < *nrsplits; ebin++)
   {
      float *ftemp;
      int *itemp;
      string *stemp;
      int rewritecode;
      // Set the random seed
      unsigned int *randSeed = new unsigned int;
      *randSeed = (unsigned int)chrono::system_clock::now().time_since_epoch().count();
      bool *singlefile = new bool;
      *singlefile = false;

      ftemp = new float[2];
      itemp = new int[6];
      stemp = new string[5];
      
      // Get the value for splitting (either fraction of events or number of events)
      ftemp[0] = (startSplitting->widgetNE[0])->GetValue();
      
      // Open input file and save All Tree and list of keys
      TFile *input = TFile::Open(infile.c_str(), "READ");
      TTree *tempTree = (TTree*)input->Get("TreeA");
      TList *tempkeyslist = (TList*)input->GetListOfKeys();
      TTree *readTree;
      TTree *writeTree;

      // Set the number of all events (itemp[0])
      itemp[0] = tempTree->GetEntries();

      // Determine the type of observables to cut on (SD or FD)
      selcuttype = (cutObservables->widgetCB)->GetSelection();

      // Prepare observables for reading and writing
      Observables *obser = new Observables(observables);
      Observables *obser_neg = new Observables(observables);
      Observables *obser_pos = new Observables(observables);

      // Name and title of the current tree
      stemp[2] = string(tempTree->GetName());
      stemp[3] = string(tempTree->GetTitle());

//      cout << "#   Applying split filtering to the all tree: " << stemp[2] << "; " << stemp[3] << endl;

      for(int i = 0; i < 3; i++)
      {
         stationDistance[i] = new vector<float>;
         stationRisetime[i] = new vector<float>;
    
         stationDistance[i]->clear();
         stationRisetime[i]->clear();
      }
      stationHSat = new vector<bool>;
      stationHSat->clear();
    
      bdist = 0;
      bdistneg = 0;
      bdistpos = 0;
      brise = 0;
      briseneg = 0;
      brisepos = 0;
      bsat = 0;

      // Prepare test tree for reading
      tempTree = (TTree*)input->Get(stemp[2].c_str());
      tempTree->SetBranchAddress("rewritecode", &rewritecode);
      for(int j = 0; j < nrobs; j++)
      {
         tempTree->SetBranchAddress((obser->GetName(j)).c_str(), &(obser->obsstruct[j].value));
         tempTree->SetBranchAddress((obser->GetName(j) + "_neg").c_str(), &(obser_neg->obsstruct[j].value));
         tempTree->SetBranchAddress((obser->GetName(j) + "_pos").c_str(), &(obser_pos->obsstruct[j].value));
      }

      tempTree->SetBranchAddress("stationdistance", &stationDistance[0], &bdist);
      tempTree->SetBranchAddress("stationdistance_neg", &stationDistance[1], &bdistneg);
      tempTree->SetBranchAddress("stationdistance_pos", &stationDistance[2], &bdistpos);
      tempTree->SetBranchAddress("stationrisetime", &stationRisetime[0], &brise);
      tempTree->SetBranchAddress("stationrisetime_neg", &stationRisetime[1], &briseneg);
      tempTree->SetBranchAddress("stationrisetime_pos", &stationRisetime[2], &brisepos);
      tempTree->SetBranchAddress("stationhighsat", &stationHSat, &bsat);

      // Preparing shuffled list for sampling
      vector<int> *shuflist = new vector<int>;
      vector<bool> *shuflistSDrec = new vector<bool>;

      // Loop over all events 
      for(int j = 0; j < itemp[0]; j++)
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

         ret = IsInsideCuts(obser, obser_neg, obser_pos, true, ebin);

         // If event is inside selected cuts, enter it into the list for shuffling
         if(ret != -1)
	 {
            // If event has a valid SD reconstruction, save it for later selection
            if((splitEnforceSD->widgetChBox[0])->IsChecked())
            {
               ret = Find(observables, "risetimerecalc");
               if(ret != -1)
	       {
                  if((obser->GetValue("energySD") != -1.) && (obser->GetValue("risetimerecalc") != -1.))
                     shuflistSDrec->push_back(true);
	          else
                     shuflistSDrec->push_back(false);
	       }
	       else
	       {
                  if(obser->GetValue("energySD") != -1.)
                     shuflistSDrec->push_back(true);
	          else
                     shuflistSDrec->push_back(false);
	       }
	    }

            shuflist->push_back(j);
	 }
         else
         {
            if(DBGSIG > 1)
               cout << "# StartFileSplit        #: " << "Event = " << j << " is outside the selected cuts." << endl;
         }
      }

      // If the we are only filtering, then don't write out if input file will be equal to output file
      if( *singlefile && ((itemp[0] - shuflist->size()) == 0) )
      {
         AlertPopup("Input and output equal", "The filtering of input file removed no events (" + ToString(itemp[0]) + "/" + ToString(shuflist->size()) + "). Stopping filtering, since input and output files will be equal. Please adjust the settings accordingly and restart.");

         for(int i = 0; i < 3; i++)
         {
            delete stationDistance[i];
            delete stationRisetime[i];
         }
         delete stationHSat;
         delete[] ftemp;
         delete[] itemp;
         delete[] stemp;
         delete obser;
         delete obser_neg;
         delete obser_pos;
         delete randSeed;
         delete shuflist;
         delete shuflistSDrec;
         delete singlefile;
         input->Close();

         return 1;
      }

      // After filtering set the number of all events (itemp[0]), events in one output file (itemp[1]) and events in other output file (itemp[2])
      itemp[3] = shuflist->size();
      cout << "# StartFileSplit        #: " << "There are " << itemp[3] << " events inside and " << (itemp[0] - itemp[3]) << " events outside the selected cuts." << endl;

      if(ftemp[0] < 0)
      {
         itemp[1] = itemp[3];
         *singlefile = true;
      }
      else if(ftemp[0] < 1)
      {
         itemp[1] = TMath::Nint(itemp[3]*ftemp[0]);
         *singlefile = false;
      }
      else
      {
         itemp[1] = TMath::Nint(ftemp[0]);
         *singlefile = false;
      }
      itemp[2] = (int)(itemp[3]-itemp[1]);

      // If any of these are 0, cancel splitting
      if( (!(*singlefile) && ((itemp[0] == 0) || (itemp[1] == 0) || (itemp[2] <= 0))) || (itemp[0] == 0) || (itemp[3] == 0) )
      {
         AlertPopup("No events in split file", "One of the split files will have no events (" + ToString(itemp[1]) + "/" + ToString(itemp[2]) + "). Total number of filtered events is " + ToString(itemp[3]) + ". Please adjust the split setting accordingly and restart.");

         for(int i = 0; i < 3; i++)
         {
            delete stationDistance[i];
            delete stationRisetime[i];
         }
         delete stationHSat;
         delete[] ftemp;
         delete[] itemp;
         delete[] stemp;
         delete obser;
         delete obser_neg;
         delete obser_pos;
         delete randSeed;
         delete shuflist;
         delete shuflistSDrec;
         delete singlefile;
         input->Close();

         return 1;
      }
      else
      {
         if( (*nrsplits > 1) || (nrsel > 1) )
	 {
            cout << "# StartFileSplit        #: Selected random seed is: " << *randSeed << endl;
            ofstream *ofs = new ofstream;
	    stemp[2] = string(rootdir) + "/results/split_seeds.txt";
	    if( (ebin == 0) && (cursel == 0) )
               ofs->open(stemp[2].c_str(), ofstream::trunc);
	    else
               ofs->open(stemp[2].c_str(), ofstream::app);
	    *ofs << *randSeed << endl;
	    if( (nrsel > 1) && (ebin == (*nrsplits)-1) )
               *ofs << endl;
	    ofs->close();
	    delete ofs;
	 }
	 else
	 {
            // Open a dialog to select the random seed
            stemp[2] = "Selecting random seed for splitting the rewritten ADST file into two parts.\nThe current seed is selected randomly based on time.\n";
            cout << "# StartFileSplit        #: Random seed = " << *randSeed << endl;
            NEDialog *randomseedDialog = new NEDialog(wxT("Random seed"), wxSize(500,200), stemp[2], "Set random seed:", *randSeed, &ID_RANDSEEDDIALOG);
            randomseedDialog->SetNEntryFormat(randomseedDialog->widgetNE, 0, 1, 2, 0, 10000000000);
            if(randomseedDialog->ShowModal() == wxID_OK)
            {
               *randSeed = (unsigned int)randomseedDialog->GetNEValue();
               cout << "# StartFileSplit        #: Selected random seed is: " << *randSeed << endl;
            }
            else
            {
               for(int i = 0; i < 3; i++)
               {
                  delete stationDistance[i];
                  delete stationRisetime[i];
               }
               delete stationHSat;
               delete[] ftemp;
               delete[] itemp;
               delete[] stemp;
               delete obser;
               delete obser_neg;
               delete obser_pos;
               delete randSeed;
               delete shuflist;
               delete shuflistSDrec;
               delete singlefile;
               randomseedDialog->Destroy();
               delete randomseedDialog;
               input->Close();

               return 1;
            }

            randomseedDialog->Destroy();
            delete randomseedDialog;
	 }

         // Prepare the output filenames
	 if(*nrsplits > 1)
	 {
            stemp[0] = RemoveExtension(&infile) + "_" + ToString(TMath::Log10(esplitBins[2*ebin]), 2) + "-" + ToString(TMath::Log10(esplitBins[2*ebin+1]), 2) + "_split-1.root";
            stemp[1] = RemoveExtension(&infile) + "_" + ToString(TMath::Log10(esplitBins[2*ebin]), 2) + "-" + ToString(TMath::Log10(esplitBins[2*ebin+1]), 2) + "_split-2.root";
	 }
	 else
	 {
            stemp[0] = RemoveExtension(&infile) + "_split-1.root";
            stemp[1] = RemoveExtension(&infile) + "_split-2.root";
	 }

         // Printout some information about splitting the files & prepare the output filenames
         itemp[4] = 0;
         // Splitting into two files
         if(!(*singlefile))
         {
	    if(*nrsplits > 1)
	    {
               stemp[0] = RemoveExtension(&infile) + "_" + ToString(TMath::Log10(esplitBins[2*ebin]), 2) + "-" + ToString(TMath::Log10(esplitBins[2*ebin+1]), 2) + "_split-1.root";
               stemp[1] = RemoveExtension(&infile) + "_" + ToString(TMath::Log10(esplitBins[2*ebin]), 2) + "-" + ToString(TMath::Log10(esplitBins[2*ebin+1]), 2) + "_split-2.root";
	    }
	    else
	    {
               stemp[0] = RemoveExtension(&infile) + "_split-1.root";
               stemp[1] = RemoveExtension(&infile) + "_split-2.root";
	    }

            stemp[2] = "Currently splitting file\n   " + RemovePath(&infile) + "        \t" + ToString(itemp[0]) + " events\ninto two files:\n   " + RemovePath(&stemp[0]) + "\t" + ToString(itemp[1]) + " events\n   " + RemovePath(&stemp[1])+ "\t" + ToString(itemp[2]) + " events\n\nPlease wait for it to finish.";
            ShowProgress(wxT("Splitting rewritten ADST file"), stemp[2].c_str(), 2*(input->GetNkeys())*itemp[0]);

            cout << "# StartFileSplit        #: Will start splitting file " << RemovePath(&infile) << " into:" << endl;
            cout << "# StartFileSplit        #: - First file:  " << RemovePath(&stemp[0]) << " (" << itemp[1] << " of " << itemp[3] << " total filtered events)" << endl;
            cout << "# StartFileSplit        #: - Second file: " << RemovePath(&stemp[1]) << " (" << itemp[2] << " of " << itemp[3] << " total filtered events)" << endl;
         }
         // Filtering a file
         else
         {
	    if(*nrsplits > 1)
	    {
               if(bCustomSplit)
                  stemp[0] = RemoveExtension(&infile) + "_" + ToString(TMath::Log10(esplitBins[2*ebin]), 2) + "-" + ToString(TMath::Log10(esplitBins[2*ebin+1]), 2) + "_split-2.root";
               else
                  stemp[0] = RemoveExtension(&infile) + "_" + ToString(TMath::Log10(esplitBins[2*ebin]), 2) + "-" + ToString(TMath::Log10(esplitBins[2*ebin+1]), 2) + "_filtered.root";
	    }
	    else
	    {
               if(bCustomSplit)
                  stemp[0] = RemoveExtension(&infile) + "_split-2.root";
               else
                  stemp[0] = RemoveExtension(&infile) + "_filtered.root";
	    }

            stemp[2] = "Currently filtering file\n   " + RemovePath(&infile) + "        \t" + ToString(itemp[0]) + " events\ninto file:\n   " + RemovePath(&stemp[0]) + "\t" + ToString(itemp[1]) + " events\n\nPlease wait for it to finish.";
            ShowProgress(wxT("Filtering rewritten ADST file"), stemp[2].c_str(), (input->GetNkeys())*itemp[0]);

            cout << "# StartFileSplit        #: Will start filtering file " << RemovePath(&infile) << " into:" << endl;
            cout << "# StartFileSplit        #: - Filtered file:  " << RemovePath(&stemp[0]) << " (" << itemp[1] << " of " << itemp[3] << " total filtered events)" << endl;
         }

         // Shuffle the event list (saving to two vectors with event numbers)
         shuffle(shuflist->begin(), shuflist->end(), default_random_engine(*randSeed));
         shuffle(shuflistSDrec->begin(), shuflistSDrec->end(), default_random_engine(*randSeed));

         vector<int> *split1list = new vector<int>;
         vector<int> *split2list = new vector<int>;

	 itemp[5] = 0;

         for(int i = 0; i < itemp[3]; i++)
         {
            if(DBGSIG > 1)
               cout << i << "\tList event:\t" << shuflist->at(i) << "\t";

            if(i < itemp[1]+itemp[5])
            {
               if((splitEnforceSD->widgetChBox[0])->IsChecked())
               {
                  if(shuflistSDrec->at(i))
                  {
                     split1list->push_back(shuflist->at(i));
                     if(DBGSIG > 1)
                        cout << " (1st list)" << endl;
//		     cout << i << ": Event has SD reconstruction - add to split1list (itemp[1] = " << itemp[1] << ", itemp[5] = " << itemp[5] << ")." << endl;
		  }
		  else
                  {
                     split2list->push_back(shuflist->at(i));
                     if(DBGSIG > 1)
                        cout << " (2nd list)" << endl;
//		     cout << i << ": Event has no SD reconstruction - add to split2list (itemp[1] = " << itemp[1] << ", itemp[5] = " << itemp[5] << ")." << endl;
		     itemp[5]++;
                  }
	       }
	       else
	       {
                  split1list->push_back(shuflist->at(i));
                  if(DBGSIG > 1)
                     cout << " (1st list)" << endl;
	       }
            }
            else
            {
               split2list->push_back(shuflist->at(i));
               if(DBGSIG > 1)
                  cout << " (2nd list)" << endl;
            }
         }

         // Loop over the two split files (break loop, if only filtering)
         TFile *output;
         for(int i = 0; i < 2; i++)
         {
            cout << "# StartFileSplit        #: Currently writing out to: " << stemp[i] << endl;
            output = TFile::Open(stemp[i].c_str(), "RECREATE");

            // Loop over all trees
            for(int k = 0; k < input->GetNkeys(); k++)
            {
               // Name and title of the current tree
               stemp[2] = string((tempkeyslist->At(k))->GetName());
               stemp[3] = string((tempkeyslist->At(k))->GetTitle());

               cout << "# StartFileSplit        #:   Currently selected tree: " << stemp[2] << "; " << stemp[3] << endl;

               // Prepare tree for reading
               for(int i = 0; i < 3; i++)
               {
                  stationDistance[i]->clear();
                  stationRisetime[i]->clear();
               }
               stationHSat->clear();
    
               bdist = 0;
               bdistneg = 0;
               bdistpos = 0;
               brise = 0;
               briseneg = 0;
               brisepos = 0;
               bsat = 0;

               readTree = (TTree*)input->Get(stemp[2].c_str());
               readTree->SetBranchAddress("rewritecode", &rewritecode);
               for(int j = 0; j < nrobs; j++)
               {
                  readTree->SetBranchAddress((obser->GetName(j)).c_str(), &(obser->obsstruct[j].value));
                  readTree->SetBranchAddress((obser->GetName(j) + "_neg").c_str(), &(obser_neg->obsstruct[j].value));
                  readTree->SetBranchAddress((obser->GetName(j) + "_pos").c_str(), &(obser_pos->obsstruct[j].value));
               }

               readTree->SetBranchAddress("stationdistance", &stationDistance[0], &bdist);
               readTree->SetBranchAddress("stationdistance_neg", &stationDistance[1], &bdistneg);
               readTree->SetBranchAddress("stationdistance_pos", &stationDistance[2], &bdistpos);
               readTree->SetBranchAddress("stationrisetime", &stationRisetime[0], &brise);
               readTree->SetBranchAddress("stationrisetime_neg", &stationRisetime[1], &briseneg);
               readTree->SetBranchAddress("stationrisetime_pos", &stationRisetime[2], &brisepos);
               readTree->SetBranchAddress("stationhighsat", &stationHSat, &bsat);

               // Prepare tree for writing
               if( (stemp[3].compare("Signal tree from old file.") == 0) || (stemp[3].compare("Signal tree from new file.") == 0) || (stemp[3].compare("Background tree with all events, including signal events.") == 0) )
                  stemp[4] = stemp[3];
               else
               {
                  // Splitting by fraction or number of events
                  if(!(*singlefile))
                  {
                     stemp[4] = stemp[3];
                  }
                  // Filtering a file
                  else
                  {
                     stemp[4] = stemp[3];
                  }
               }
               writeTree = new TTree(stemp[2].c_str(), stemp[4].c_str());
               writeTree->Branch("rewritecode", &rewritecode, "rewritecode/I");
               for(int j = 0; j < nrobs; j++)
               {
                  writeTree->Branch((obser->GetName(j)).c_str(), &(obser->obsstruct[j].value), (obser->GetName(j) + "/F").c_str());
                  writeTree->Branch((obser->GetName(j) + "_neg").c_str(), &(obser_neg->obsstruct[j].value), (obser->GetName(j) + "_neg/F").c_str());
                  writeTree->Branch((obser->GetName(j) + "_pos").c_str(), &(obser_pos->obsstruct[j].value), (obser->GetName(j) + "_pos/F").c_str());
               }

               writeTree->Branch("stationdistance", &stationDistance[0]);
               writeTree->Branch("stationdistance_neg", &stationDistance[1]);
               writeTree->Branch("stationdistance_pos", &stationDistance[2]);
               writeTree->Branch("stationrisetime", &stationRisetime[0]);
               writeTree->Branch("stationrisetime_neg", &stationRisetime[1]);
               writeTree->Branch("stationrisetime_pos", &stationRisetime[2]);
               writeTree->Branch("stationhighsat", &stationHSat);

               // Check if the values in this tree are valid
               if( (stemp[3].compare("Signal tree from old file.") == 0) || (stemp[3].compare("Signal tree from new file.") == 0) )
               {
                  // Update the progress bar
                  itemp[4]+=itemp[0];
                  progress->Update(itemp[4]);
               }
               else
               {
                  // Loop over all events 
                  for(int j = 0; j < itemp[0]; j++)
                  {
                     // Read all observables
                     readTree->GetEntry(j);

                     // Read all vectors for SD stations
                     tentry = readTree->LoadTree(j);
                     bdist->GetEntry(tentry);
                     bdistneg->GetEntry(tentry);
                     bdistpos->GetEntry(tentry);
                     brise->GetEntry(tentry);
                     briseneg->GetEntry(tentry);
                     brisepos->GetEntry(tentry);
                     bsat->GetEntry(tentry);

                     // Select the correct tree to write to
                     if( (find(split1list->begin(), split1list->end(), j) != split1list->end()) && (i == 0) )
                        writeTree->Fill();
                     if( (find(split2list->begin(), split2list->end(), j) != split2list->end()) && (i == 1) )
                        writeTree->Fill();

                     // Update the progress bar
                     itemp[4]++;
                     if(!(*singlefile))
           	     {
                        if(itemp[4]%((int)(2*(input->GetNkeys())*itemp[0]*0.05)) == 0)
                           progress->Update(itemp[4]);
           	     }
           	     else
           	     {
                        if(itemp[4]%((int)((input->GetNkeys())*itemp[0]*0.05)) == 0)
                           progress->Update(itemp[4]);
           	     }
                  }
               }

               writeTree->Write();

               delete writeTree;
            }

            output->Close();

            if(*singlefile && (i == 0))
               break;
         }

         if(!(*singlefile))
            progress->Update(2*(input->GetNkeys())*itemp[0]);
         else
            progress->Update((input->GetNkeys())*itemp[0]);

         delete split1list;
         delete split2list;
      }

      for(int i = 0; i < 3; i++)
      {
         delete stationDistance[i];
         delete stationRisetime[i];
      }
      delete stationHSat;

      delete[] ftemp;
      delete[] itemp;
      delete[] stemp;
      delete obser;
      delete obser_neg;
      delete obser_pos;
      delete randSeed;
      delete shuflist;
      delete shuflistSDrec;
      delete singlefile;
      input->Close();
   }

   return 0;
}

// Combine the observables from multiple ADST files into a single root file
int MyFrame::StartCombine(string *outfile)
{
   // Check for all selections in the second listbox
   (mvaList[1]->widgetLB)->GetSelections(selections);

   if(selections.GetCount() > 1)
   {
      string *stemp = new string;

      string *flist = new string[selections.GetCount()];
      for(int i = 0; i < selections.GetCount(); i++)
         flist[i] = (mvaList[1]->widgetLB)->GetString(selections[i]);

      // Open a dialog to select the titles
      string *slabels = new string[selections.GetCount()];
      for(int i = 0; i < selections.GetCount(); i++)
         slabels[i] = "Set title for TreeS" + ToString(i+1) + ": ";

      string *svals = new string[selections.GetCount()];
      TFile *ifile;
      TTree *temptree;
      TList *tempkeyslist;
      for(int i = 0; i < selections.GetCount(); i++)
      {
         ifile = TFile::Open(flist[i].c_str(), "READ");
	 tempkeyslist = (TList*)ifile->GetListOfKeys();

	 for(int j = 0; j < ifile->GetNkeys(); j++)
	 {
            if(strcmp((tempkeyslist->At(j))->GetName(), "TreeA") != 0)
	    {
               if( (strcmp((tempkeyslist->At(j))->GetTitle(), "Signal tree from old file.") != 0) && (strcmp((tempkeyslist->At(j))->GetTitle(), "Signal tree from new file.") != 0) )
                  svals[i] = (tempkeyslist->At(j))->GetTitle();
	    }
	 }

	 ifile->Close();
      }

      *stemp = "Selecting titles for combining rewritten ADST files into an MVA input file.\nImportant! All titles must be different.\n";
      FSDialog *titleselectDialog = new FSDialog(wxT("Title selection"), wxSize(600,150+30*(selections.GetCount())), *stemp, slabels, selections.GetCount(), svals, &ID_TITLEDIALOG);
      if(titleselectDialog->ShowModal() == wxID_OK)
      {
         for(int i = 0; i < selections.GetCount(); i++)
            svals[i] = titleselectDialog->GetFSValue(i);

	 ret = 1;
      }
      else
         ret = 0;
      
      titleselectDialog->Destroy();
      delete titleselectDialog;

      // Combine the files with the hadd function
      hadd(selections.GetCount(), flist, outfile, svals, (bool)ret);
      InfoPopup("Finished combining files", "The selected files have been combined and can be used as input files to the MVA analysis.");

      delete[] slabels;
      delete[] svals;

      delete[] flist;
      delete stemp;
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
