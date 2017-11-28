#include "frame.h"
#include "separate_functions.h"

// Function for reading the layout
void MyFrame::ReadLayout(int *width, int *height)
{
   ifstream ilayout;
   string *stemp;
   stemp = new string[2];
   int *itemp;
   itemp = new int[2];
   char readTemp[1024];

   // Check for the last used layout file
   stemp[0] = string(rootdir) + "/layout/selected_layout.txt";
   ilayout.open(stemp[0].c_str(), ifstream::in);
   if(ilayout.is_open())
   {
      ilayout >> stemp[1];
   }
   ilayout.close();
   if(DBGSIG > 1)
      cout << "# ReadLayout            #: " << "Loaded layout file is: " << stemp[1] << endl;

   // Read the current layout
   stemp[0] = string(rootdir) + "/layout/" + stemp[1];
   ilayout.open(stemp[0].c_str(), ifstream::in);
   if(ilayout.is_open())
   {
      while(1)
      {
         if(ilayout.peek() == '#')
	    ilayout.getline(readTemp, 1024, '\n');
	 else if(ilayout.peek() == '\n')
	    ilayout.ignore(1, '\n');
	 else
	 {
	    ilayout >> itemp[0] >> itemp[1] >> readTemp;
	    ilayout.ignore(1, '\n');
	    break;
	 }
      }
   }
   ilayout.close();

   // Set the window width and height
   *width = itemp[0];
   *height = itemp[1];

   delete[] stemp;
   delete[] itemp;
}

// Function for setting a saved layout
void MyFrame::SetLayout(wxCommandEvent& event)
{
   ofstream olayout;
   string *stemp;
   stemp = new string[3];

   // Select the layout file and save this to selected_layout.txt file for future use
   stemp[0] = string(rootdir) + "/layout";
   wxFileDialog *openFileDlg = new wxFileDialog(this, wxT("Open the layout file"), stemp[0], "", wxT("Layout files (*.layout)|*.layout"), wxFD_OPEN);

   if(openFileDlg->ShowModal() == wxID_OK)
   {
      stemp[1] = openFileDlg->GetPath();
      if(!stemp[1].empty())
      {
         stemp[2] = string(rootdir) + "/layout/selected_layout.txt";
         olayout.open(stemp[2].c_str(), ofstream::out | ofstream::trunc);
         stemp[0] = RemovePath(&stemp[1]);
	 olayout << stemp[0];
	 olayout.close();
      }

//      wxMessageDialog *dial = new wxMessageDialog(NULL, wxT("Please restart the program to enable the selected layout for future use."), wxT("Setting new layout"), wxOK);
//      dial->ShowModal();
   }

   delete[] stemp;
}

// Function for saving the current layout
void MyFrame::SaveLayout(wxCommandEvent& event)
{
   ofstream olayout;
   string *stemp;
   stemp = new string[2];
   int *itemp;
   itemp = new int[2];

   // Select the file to save the new layout
   stemp[0] = string(rootdir) + "/layout";
   wxFileDialog *saveFileDlg = new wxFileDialog(this, wxT("Save the layout file"), stemp[0], "", wxT("Layout files (*.layout)|*.layout"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

   if(saveFileDlg->ShowModal() == wxID_CANCEL)
   {
      delete[] stemp;
      delete[] itemp;
      return;
   }
   
   stemp[1] = saveFileDlg->GetPath();
   stemp[1] = CheckExtension(&stemp[1], "layout");
   if(!stemp[1].empty())
   {
      this->GetSize(&itemp[0], &itemp[1]);
      olayout.open(stemp[1].c_str(), ofstream::out);
      if(olayout.is_open())
      {
         olayout << "# Whole window width and height" << std::endl;
	 olayout << itemp[0] << "\t" << itemp[1] << endl;
      }
      else
         cout << "# SaveLayout            #: " << "Error! Save file can not be opened (please do not use default.layout since it is write protected)." << endl;
      olayout.close();
   }

   delete[] stemp;
   delete[] itemp;
}

// Function for creating a title
void MyFrame::MakeTitle(wxPanel *parent, wxBoxSizer *box, wxString title)
{
   wxStaticText *lab = new wxStaticText(parent, -1, title, wxDefaultPosition, wxSize(-1,-1), wxALIGN_CENTRE_HORIZONTAL);
   lab->SetForegroundColour(wxColour(252,252,252));
   lab->SetBackgroundColour(wxColour(46,90,134));
   lab->SetFont(wxFont(12, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD, false));
   box->Add(lab, 0, wxEXPAND, 0);
}

// Function for opening the About information
void MyFrame::OnAbout(wxCommandEvent& event)
{
   wxMessageBox("Software for analysis of CORSIKA simulations and reconstructions with Offline. Enables multivariate analysis of mass composition sensitive observables.\n\nCreated by Gasper Kukec Mezek (gasper.kukec@ung.si). Updated on July 17th 2015", "About auger-analysis-gui", wxOK | wxICON_INFORMATION);
}

// Function for quitting the application
void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
   Close(true);
}

void MyFrame::ShowProgress(wxString title, wxString message, int maxval)
{
   progress = new wxProgressDialog(title, message, maxval, NULL, wxPD_APP_MODAL | wxPD_AUTO_HIDE | wxPD_ELAPSED_TIME | wxPD_SMOOTH);
   progress->Update(0);
}
