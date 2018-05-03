#include "frame.h"
#include "separate_functions.h"
#include "popups.h"

// Open a file when selecting it for MVA analysis (File selection: ...)
void MyFrame::SelectMvaFile(wxCommandEvent& event)
{
   string *tempstring;
   int *retval;

   wxFileDialog *openFileDlg = new wxFileDialog(this, wxT("Open Offline root file"), *currentMvaDir, "", wxT("MVA files (*.root)|*.root"), wxFD_OPEN | wxFD_MULTIPLE);

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

	    if(*retval == 1)
	    {
               (mvaList[0]->widgetLB)->Append(filenameArray->Item(i));
	       *currentMvaDir = RemoveFromLast(tempstring, "/");
	    }
	    else if(*retval == 2)
	    {
               (mvaList[1]->widgetLB)->Append(filenameArray->Item(i));
	       *currentMvaDir = RemoveFromLast(tempstring, "/");
	    }
	    else if(*retval == 3)
	    {
               (mvaList[2]->widgetLB)->Append(filenameArray->Item(i));
	       *currentMvaDir = RemoveFromLast(tempstring, "/");
	    }
	 }

	 delete tempstring;
	 delete retval;
      }

      delete filenameArray;
   }

   delete openFileDlg;
}

// Check to decide if the file format is correct
int MyFrame::CheckFormat(string *infile)
{
   TFile *tempfile = TFile::Open(infile->c_str(), "READ");
   if(tempfile->IsOpen())
   {
      if( (tempfile->GetListOfKeys()->Contains("recEvent")) || (tempfile->GetListOfKeys()->Contains("eventInfo")) )
      {
         cout << "# CheckFormat           #: " << "Input file " << *infile << " is in ADST format." << endl;
	 tempfile->Close();
	 return 1;
      }
      else if( (tempfile->GetListOfKeys()->Contains("TreeOldS1")) || (tempfile->GetListOfKeys()->Contains("TreeNewS1")) )
      {
         cout << "# CheckFormat           #: " << "Input file " << *infile << " is in MVA rewritten ADST format." << endl;
	 tempfile->Close();
	 return 2;
      }
      else if( (tempfile->GetListOfKeys()->Contains("TreeS1")) || (tempfile->GetListOfKeys()->Contains("TreeS1")) )
      {
         cout << "# CheckFormat           #: " << "Input file " << *infile << " is in MVA rewritten ADST format and will be used as input for MVA analysis." << endl;
	 tempfile->Close();
	 return 3;
      }
      else
      {
	 AlertPopup("Input file error", "The selected input file " + *infile + " is not in a readable format for this program. Please select ADST files, rewritten ADST files or combined rewritten ADST files.");
         tempfile->Close();
	 return -1;
      }
   }
   else
   {
      AlertPopup("Input file error", "The selected input file " + *infile + " can't be opened. Please make sure it is readable and that it exists.");
      tempfile->Close();
      return -1;
   }
}

// Select a file to split depending on the fraction of events we need -> and write it back to rewritten ADST format
void MyFrame::PrepareFileSplit(wxCommandEvent& event)
{
   (mvaList[1]->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
   {
      if(selections.GetCount() == 1)
      {
         ret = StartFileSplit(string((mvaList[1]->widgetLB)->GetString(selections[0])));
         if(ret == 0)
            InfoPopup("Finished spliting files", "The selected files have been split according to the set fraction.");
      }
      else
         AlertPopup("Multiple files selected", "Multiple files from the second listbox were selected. Please select only one file to split.");
   }
   else
      AlertPopup("No selected files", "No files from the second listbox were selected. Please select only one file to split.");
}

// Select a file to use for rewritten file
void MyFrame::SelectRewrite(wxCommandEvent& event)
{
   (mvaList[0]->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
   {
      string *stemp;
      stemp = new string[2];
      ret = 0;

      // Select the file to save the rewritten format
      wxFileDialog *saveFileDlg = new wxFileDialog(this, wxT("Save the rewritten ADST file"), *currentRewriteDir, "", wxT("Rewritten ADST (*.root)|*.root"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

      if(saveFileDlg->ShowModal() == wxID_CANCEL)
      {
         delete[] stemp;
	 delete saveFileDlg;
         return;
      }
      
      stemp[0] = saveFileDlg->GetPath();
      stemp[0] = CheckExtension(&stemp[0], "root");
      if(!stemp[0].empty())
      {
         *currentRewriteDir = RemoveFromLast(&stemp[0], "/");
         ret = StartRewrite(&stemp[0]);
      }

      delete[] stemp;
      delete saveFileDlg;
   }
   else
      AlertPopup("No selected files", "No files from the first listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to rewrite.");
}

// Select a file to use for combined file
void MyFrame::SelectCombine(wxCommandEvent& event)
{
   (mvaList[1]->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
   {
      string *stemp;
      stemp = new string[2];
      ret = 0;

      // Select the file to save the combined file
      wxFileDialog *saveFileDlg = new wxFileDialog(this, wxT("Save the combined ADST file"), *currentRewriteDir, "", wxT("Combined rewritten ADST (*.root)|*.root"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

      if(saveFileDlg->ShowModal() == wxID_CANCEL)
      {
         delete[] stemp;
         delete saveFileDlg;
         return;
      }
      
      stemp[0] = saveFileDlg->GetPath();
      stemp[0] = CheckExtension(&stemp[0], "root");
      if(!stemp[0].empty())
      {
         *currentRewriteDir = RemoveFromLast(&stemp[0], "/");
         ret = StartCombine(&stemp[0]);
      }

      delete[] stemp;
      delete saveFileDlg;
   }
   else
      AlertPopup("No selected files", "No files from the second listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to rewrite.");
}

// Select a file to use for combined file
void MyFrame::SelectMerge(wxCommandEvent& event)
{
   (mvaList[1]->widgetLB)->GetSelections(selections);

   if(!selections.IsEmpty())
   {
      string *stemp;
      stemp = new string[2];
      ret = 0;

      // Select the file to save the combined file
      wxFileDialog *saveFileDlg = new wxFileDialog(this, wxT("Save the merged ADST file"), *currentRewriteDir, "", wxT("Merged rewritten ADST (*.root)|*.root"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

      if(saveFileDlg->ShowModal() == wxID_CANCEL)
      {
         delete[] stemp;
         delete saveFileDlg;
         return;
      }
      
      stemp[0] = saveFileDlg->GetPath();
      stemp[0] = CheckExtension(&stemp[0], "root");
      if(!stemp[0].empty())
      {
         *currentRewriteDir = RemoveFromLast(&stemp[0], "/");
         ret = StartMerge(&stemp[0]);
      }

      delete[] stemp;
      delete saveFileDlg;
   }
   else
      AlertPopup("No selected files", "No files from the second listbox were selected. Please select one or more files (holding Ctrl or Shift while clicking) to rewrite.");
}

// Select a file to use for MVA analysis (Start analysis)
string MyFrame::SelectMva()
{
   if(!((selectedMva->widgetTE)->GetLineText(0)).IsEmpty())
   {
      string stemp;
      ret = 0;

      // Select the file to save the MVA file
      wxFileDialog *saveFileDlg = new wxFileDialog(this, wxT("Save the MVA analysis file"), *currentAnalysisDir, "", wxT("MVA save file (*.root)|*.root"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

      if(saveFileDlg->ShowModal() == wxID_CANCEL)
      {
         return "";
      }
      
      stemp = saveFileDlg->GetPath();
      stemp = CheckExtension(&stemp, "root");
      if(!stemp.empty())
      {
         *currentAnalysisDir = RemoveFromLast(&stemp, "/");
         delete saveFileDlg;
         return stemp;
      }
      else
      {
         delete saveFileDlg;
         return "";
      }
   }
   else
   {
      AlertPopup("No selected file", "No file was selected for MVA analysis. Please double click a file from the third listbox.");
      return "";
   }
}
