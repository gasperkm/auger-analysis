#include "frame.h"
#include "separate_functions.h"
#include "root_style.h"
#include "primary_type.h"
#include "popups.h"
#include "parse_directory_files.h"
#include "mva_fit_histogram.h"
#include "histograming.h"
#include "mva_methods.h"

#include <algorithm>

// Open a file when selecting it for plotting (File selection: ...)
void MyFrame::SelectPlotFile(wxCommandEvent& event)
{
   string *tempstring;
   int *retval;

   wxFileDialog *openFileDlg = new wxFileDialog(this, wxT("Open root file for plotting"), *currentPlotDir, "", wxT("Root files (*.root)|*.root"), wxFD_OPEN | wxFD_MULTIPLE);

   if(openFileDlg->ShowModal() == wxID_OK)
   {
      wxArrayString *filenameArray = new wxArrayString;
      openFileDlg->GetPaths(*filenameArray);

      if(!filenameArray->IsEmpty())
      {
         tempstring = new string;
	 retval = new int;

         for(int i = 0; i < filenameArray->GetCount(); i++)
	 {
            *tempstring = filenameArray->Item(i);
	    *retval = CheckFormat(tempstring);

	    if(*retval == 3)
	    {
	       SetPlotTreeSelect(tempstring, (firstHistTree->widgetCB));
               (firstHistTree->widgetCB)->Append("Disable");
	       SetPlotTreeSelect(tempstring, (secondHistTree->widgetCB));
               (secondHistTree->widgetCB)->Append("Disable");
	       SetPlotTreeSelect(tempstring, (thirdHistTree->widgetCB));
               (thirdHistTree->widgetCB)->Append("Disable");

               (openFileList->widgetLB)->Append(filenameArray->Item(i));
	       *currentPlotDir = RemoveFromLast(tempstring, "/");
	    }
	 }

	 delete tempstring;
	 delete retval;
      }

      delete filenameArray;
   }

   delete openFileDlg;
}

// Open a directory when selecting it for plotting (Directory selection: ...)
void MyFrame::SelectPlotDir(wxCommandEvent& event)
{
   int *retval;

   wxDirDialog *openDirDlg = new wxDirDialog(this, wxT("Open directory with files for plotting"), *currentPlotDir, wxDD_DEFAULT_STYLE | wxDD_DIR_MUST_EXIST);

   if(openDirDlg->ShowModal() == wxID_OK)
   {
      wxString *directoryPath = new wxString;
      *directoryPath = openDirDlg->GetPath();

      retval = new int;

      vector<string> *tempvectstring = new vector<string>;
      *retval = ParseDirectory(string(*directoryPath), tempvectstring);

      if( (*retval == 0) && (tempvectstring->size() > 0) )
      {
         sort(tempvectstring->begin(), tempvectstring->end());
	 SetPlotTreeSelect(&(tempvectstring->at(0)), (dataPlotSelect->widgetCB));

         for(int i = 0; i < tempvectstring->size(); i++)
	 {
            *retval = CheckFormat(&(tempvectstring->at(i)));

	    if(*retval == 3)
	    {
               (openFileList->widgetLB)->Append(tempvectstring->at(i));
	       *currentPlotDir = string(*directoryPath);
	    }
	 }
      }

      delete tempvectstring;
      delete retval;

      delete directoryPath;
   }

   delete openDirDlg;
}

// Set tree names to the plotting dropbox
void MyFrame::SetPlotTreeSelect(string *plotFile, wxChoice *combo)
{
   string *stemp = new string[2];
   string *signalName = new string;
   int *cnt = new int;
   *cnt = 0;
   int *cntdat = new int;
   *cntdat = -1;

   cout << "Plot file = " << *plotFile << endl;

   // Open the file with root
   TFile *ifile = TFile::Open(plotFile->c_str(), "READ");
   combo->Clear();

   // Append all trees from root files
   nrkeys = GetRootKeys(ifile);
   keyslist = (TList*) ifile->GetListOfKeys();
   for(int i = 1; i <= nrkeys; i++)
   {
      if(strcmp((keyslist->At(i-1))->GetName(), "TreeA") != 0)
      {
         *signalName = "TreeS" + ToString(i);
	 if(ifile->GetListOfKeys()->Contains(signalName->c_str()))
	 {
            stemp[1] = string(ifile->GetKey(signalName->c_str())->GetTitle());
	    if(stemp[1].find("Data") != string::npos)
               *cntdat = i-1;
            combo->Append(stemp[1]);
            (*cnt)++;
	 }
      }
   }

   cout << "# SetPlotTreeSelect     #: " << "Found " << *cnt << " signal keys in the selected file." << endl;

   // Select the data tree
   if(*cntdat == -1)
      combo->SetSelection(0);
   else
      combo->SetSelection(*cntdat);

   ifile->Close();

   delete[] stemp;
   delete signalName;
   delete cnt;
   delete cntdat;
}

// Start plotting histograms
void MyFrame::StartHistogramPlot(wxCommandEvent& event)
{
   (openFileList->widgetLB)->GetSelections(selections);
   if(!selections.IsEmpty())
      StartHistogramScatterPlot(0);
   else
      AlertPopup("No selected files", "No files from the plotting listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to use for plotting.");
}

// Start plotting scatter plots
void MyFrame::StartScatterPlot(wxCommandEvent& event)
{
   (openFileList->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
      StartHistogramScatterPlot(1);
   else
      AlertPopup("No selected files", "No files from the plotting listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to use for plotting.");
}

void MyFrame::StartHistogramScatterPlot(int type)
{
   int *itemp = new int[3];
   string *stemp = new string[2];

   // Check for selections in the open plot file listbox
   vector<string> *fileSelections = new vector<string>;
   (openFileList->widgetLB)->GetSelections(selections);
   for(int i = 0; i < selections.GetCount(); i++)
      fileSelections->push_back((string)(openFileList->widgetLB)->GetString(selections[i]));

   // Check for selections of the observables
   vector<string> *obsSelections = new vector<string>;
   (plottingList[1]->widgetLB)->GetSelections(plotSelections);
   for(int i = 0; i < plotSelections.GetCount(); i++)
      obsSelections->push_back((string)(plottingList[1]->widgetLB)->GetString(plotSelections[i]));

   // Check for selections of trees (when merging)
   vector<string> *treeSelections = new vector<string>;
   treeSelections->push_back((string)(firstHistTree->widgetCB)->GetStringSelection());
   treeSelections->push_back((string)(secondHistTree->widgetCB)->GetStringSelection());
   treeSelections->push_back((string)(thirdHistTree->widgetCB)->GetStringSelection());

   // Check for the number of bins
   int *nrbins = new int;
   *nrbins = (plotNrBins->widgetNE[0])->GetValue();

   // Check for X axis limits
   double *xlimit = new double[2];
   xlimit[0] = (plotXaxisRange->widgetNE[0])->GetValue();
   xlimit[1] = (plotXaxisRange->widgetNE[1])->GetValue();

   // Check for Y axis limits
   double *ylimit = new double[2];
   ylimit[0] = (plotYaxisRange->widgetNE[0])->GetValue();
   ylimit[1] = (plotYaxisRange->widgetNE[1])->GetValue();

   // Check other plotting settings
   bool *otherSet = new bool[4];
   otherSet[0] = (bool)(specialHistPlot->widgetChBox[0])->IsChecked();
   otherSet[1] = (bool)(specialHistPlot->widgetChBox[1])->IsChecked();

   if(fileSelections->size() > 0)
   {
      ScatHist *scathist = new ScatHist(&type);
      
      scathist->SetFilenames(fileSelections);
      scathist->SetObservables(obsSelections);
      scathist->SetTrees(treeSelections);
      scathist->SetNrBins(nrbins);
      scathist->SetOtherSettings(otherSet);
      scathist->SetXaxisLimits((plotXaxisRange->widgetChBox)->IsChecked(), xlimit);
      scathist->SetYaxisLimits((plotYaxisRange->widgetChBox)->IsChecked(), ylimit);

      scathist->StartPlotting();

      delete scathist;
   }
   else
   {
      AlertPopup("No selected files", "No files from the plotting listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to use for plotting.");
      return;
   }

   delete ylimit;
   delete xlimit;
   delete nrbins;
   delete[] otherSet;
   delete fileSelections;
   delete obsSelections;
   delete treeSelections;
   delete[] itemp;
   delete[] stemp;
}

// Start performing the MVA histogram fit
void MyFrame::StartMvaHistFit(wxCommandEvent& event)
{
   int *itemp = new int[3];
   string *stemp = new string[2];

   // Check for selections in the open plot file listbox
   vector<string> *fileSelections = new vector<string>;
   (openFileList->widgetLB)->GetSelections(selections);
   for(int i = 0; i < selections.GetCount(); i++)
      fileSelections->push_back((string)(openFileList->widgetLB)->GetString(selections[i]));

   stemp[0] = fileSelections->at(0);

   // Check for MVA method
   string *plotMethod = new string;
   *plotMethod = (methodsPlotSelect->widgetCB)->GetStringSelection();

   // Check for selections of the particle composition
   vector<int> *compSelections = new vector<int>;
   (plottingList[0]->widgetLB)->GetSelections(selections);
   for(int i = 0; i < selections.GetCount(); i++)
      compSelections->push_back((int)selections[i]);

   // Check for HI model
   int *plotHImodel = new int;
   *plotHImodel = (hiSelect->widgetCB)->GetSelection();

   // Check for simulation production
   int *plotProd = new int;
   *plotProd = (prodSelect->widgetCB)->GetSelection();

   // Check for data tree
   int *plotDataTree = new int;
   *plotDataTree = (dataPlotSelect->widgetCB)->GetSelection();

   // Check for the number of bins
   int *nrbins = new int;
   *nrbins = (plotNrBins->widgetNE[0])->GetValue();

   // Check for X axis limits
   double *xlimit = new double[2];
   xlimit[0] = (plotXaxisRange->widgetNE[0])->GetValue();
   xlimit[1] = (plotXaxisRange->widgetNE[1])->GetValue();

   // Check for Y axis limits
   double *ylimit = new double[2];
   ylimit[0] = (plotYaxisRange->widgetNE[0])->GetValue();
   ylimit[1] = (plotYaxisRange->widgetNE[1])->GetValue();

   // Check other plotting settings
   bool *otherSet = new bool[4];
   otherSet[0] = (bool)(specialPlot->widgetChBox[0])->IsChecked();
   otherSet[1] = (bool)(specialPlot->widgetChBox[1])->IsChecked();
   otherSet[2] = (bool)(specialPlot->widgetChBox[2])->IsChecked();
   otherSet[3] = (bool)(specialPlot->widgetChBox[3])->IsChecked();

   if(fileSelections->size() > 1)
   {
      if(compSelections->size() > 1)
      {
         vector<double> *finalLna = new vector<double>;
         vector<double> *finalChi2 = new vector<double>;
         vector<double> *finalPvalue = new vector<double>;
         vector<double> *finalNdf = new vector<double>;
         vector<double> *finalComp = new vector<double>;
         vector<double> *finalEnergy = new vector<double>;
         vector<double> *finalStep = new vector<double>;

         MvaFitHist *fithist = new MvaFitHist();

         fithist->SetPrimaries(compSelections);
	 fithist->SetMethod(plotMethod);
	 fithist->SetHImodel(plotHImodel);
	 fithist->SetSimProduction(plotProd);
	 fithist->SetData(plotDataTree);
         fithist->SetNrBins(nrbins);
         fithist->SetOtherSettings(otherSet);
         fithist->SetXaxisLimits((plotXaxisRange->widgetChBox)->IsChecked(), xlimit);
         fithist->SetYaxisLimits((plotYaxisRange->widgetChBox)->IsChecked(), ylimit);

	 itemp[1] = 0;

         for(int i = 0; i < fileSelections->size(); i++)
         {
            fithist->SetInputFile(&(fileSelections->at(i)));
	    fithist->PrepareHistograms(i);
	    itemp[0] = fithist->StartFitting();
	    fithist->ResetHistograms();

	    if(itemp[0] != -1)
	    {
	       for(int j = 0; j < compSelections->size(); j++)
	       {
	          finalComp->push_back(fithist->GetFinalComposition(j));
	          finalComp->push_back(fithist->GetFinalCompositionErr(j,0));		// 1 is low error
	          finalComp->push_back(fithist->GetFinalCompositionErr(j,1));		// 0 is high error
	       }

	       finalLna->push_back(fithist->GetFinalLna(0));
	       finalLna->push_back(fithist->GetFinalLna(1));
	       finalLna->push_back(fithist->GetFinalLna(2));

	       finalEnergy->push_back(fithist->GetFinalEnergy());

	       finalChi2->push_back(fithist->GetChi2());
	       finalPvalue->push_back(fithist->GetPvalue());
	       finalNdf->push_back((double)fithist->GetNdf());
	       finalStep->push_back(fithist->GetStep());

	       itemp[1]++;
	    }
         }

	 fithist->PrintoutFinalResults(&itemp[1], finalLna, finalChi2, finalPvalue, finalNdf, finalComp, finalEnergy, finalStep);
	 fithist->PlotLnaComposition();

         delete fithist;

	 delete finalLna;
         delete finalChi2;
         delete finalPvalue;
	 delete finalNdf;
	 delete finalComp;
	 delete finalEnergy;
	 delete finalStep;

	 stemp[1] = RemoveExtension(&stemp[0]);
	 stemp[0] = RemoveExtension(&stemp[1]) + "/plots";
         InfoPopup("MVA histogram fit finished", "MVA histogram fit has finished. Created plots have been saved to " + stemp[0] + ".");
      }
      else
      {
         AlertPopup("Incorrect composition", "At least two elements need to be selected from the composition listbox. Please select at least two elements (holding Ctrl or Shift while clicking) to use them for setting the composition.");
         return;
      }
   }
   else
   {
      AlertPopup("No selected files", "At least two files must be selected from the plotting listbox. Please select at least two files (holding Ctrl or Shift while clicking) to use them for MVA histogram fitting.");
      return;
   }

   delete plotDataTree;
   delete plotProd;
   delete plotHImodel;
   delete[] otherSet;
   delete ylimit;
   delete xlimit;
   delete nrbins;
   delete fileSelections;
   delete plotMethod;
   delete compSelections;
   delete[] itemp;
   delete[] stemp;
}

// Set default values for plotting
void MyFrame::SetDefaultPlot(wxCommandEvent& event)
{
}

// Set default values for plotting
void MyFrame::SetDefaultMvaPlot(wxCommandEvent& event)
{
}

// Set default values for histogram and scatter plot plotting
void MyFrame::SetDefaultHistPlot(wxCommandEvent& event)
{
}
