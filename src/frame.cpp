#include "frame.h"
#include "separate_functions.h"
#include "observables.h"
#include "popups.h"
#include "primary_type.h"

// wxCollapsiblePane

MyFrame::MyFrame(const wxString& title) : wxFrame(NULL, -1, title, wxDefaultPosition, wxSize(1024,671))
{
   InitObservables();
   InitMethods();
   InitVariables();

   ReadLayout(&WW, &WH);
   if(DBGSIG > 1)
      cout << "# MyFrame               #: " << "Window size is: " << WW << "x" << WH << endl;
   this->SetSize(0,0,WW,WH);

   // Window size (default)
   int lwidth = WW/2;
   int rwidth = WW/2;

   wxStaticText *temptext;
   vector<string> *vstemp;
   vector<string> *vstemp2;
   vector<int> *vitemp;
   vector<int> *vitemp2;
   vector<double> *vdtemp;

   vstemp = new vector<string>;
   vstemp2 = new vector<string>;
   vitemp = new vector<int>;
   vitemp2 = new vector<int>;
   vdtemp = new vector<double>;

   // Prepare the menu bar
   menubar = new wxMenuBar;
   file = new wxMenu;
   file->Append(ID_SETLAYOUT, wxT("Set &user layout"));
   file->Append(ID_SAVELAYOUT, wxT("Save &current layout"));
   file->Append(wxID_EXIT, wxT("&Quit\tCtrl+W"));
   menubar->Append(file, wxT("&File"));
   help = new wxMenu;
   help->Append(ID_OPENHELP, wxT("Open &help in browser"));
   help->Append(wxID_ABOUT, wxT("&About"));
   menubar->Append(help, wxT("&Help"));
   SetMenuBar(menubar);

   // Connection functions for menu bar items
   Connect(ID_SAVELAYOUT, wxEVT_COMMAND_TOOL_CLICKED, wxCommandEventHandler(MyFrame::SaveLayout));
   Connect(ID_SETLAYOUT, wxEVT_COMMAND_TOOL_CLICKED, wxCommandEventHandler(MyFrame::SetLayout));
   Connect(wxID_EXIT, wxEVT_COMMAND_TOOL_CLICKED, wxCommandEventHandler(MyFrame::OnQuit));
   Connect(wxID_ABOUT, wxEVT_COMMAND_TOOL_CLICKED, wxCommandEventHandler(MyFrame::OnAbout));

   // Create a tab for MVA analysis
   wxNotebook *tabs = new wxNotebook(this, -1, wxDefaultPosition, wxDefaultSize, wxNB_BOTTOM);

// MVA analysis -----------------------------------------------------

   // Set the scrollable window inside the frame
   wxScrolledWindow *sw = new wxScrolledWindow(tabs);

   // Horizontally align panels inside the scrolled window
   wxBoxSizer *hboxhead = new wxBoxSizer(wxHORIZONTAL);

   // Set the left panel inside the scrollable window (MVA analysis) ------------------------------------
   wxPanel *leftmvapanel = new wxPanel(sw, -1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL | wxBORDER_RAISED);

   // Vertically align elements inside the panel
   wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

   // Make a title
   MakeTitle(leftmvapanel, vbox, wxT("Filelist"));

   // Label and button for opening root files
   selectMvaFile = new LabelButton(leftmvapanel, wxT("File selection:"), wxT("..."), ID_OPENFILE, lwidth);
   vbox->Add(selectMvaFile->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_OPENFILE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SelectMvaFile));

   // Label and listbox for displaying opened root files
   vitemp->clear();
   vitemp->push_back(ID_DELETELIST + 4*nrlists);
   vitemp->push_back(ID_UPLIST + 4*nrlists);
   vitemp->push_back(ID_DOWNLIST + 4*nrlists);
   vitemp->push_back(ID_CLEARLIST + 4*nrlists);

   mvaList[0] = new LabelListEdit(leftmvapanel, wxT("ADST files:"), lwidth, 70, 1, -1, vitemp);
   vbox->Add(mvaList[0]->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(vitemp->at(0), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(1), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(2), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(3), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   allLBE[nrlists] = mvaList[0];
   nrlists++;

   // Label and button to start rewriting ADST into simple form
   startRewriting = new LabelButton(leftmvapanel, wxT("Select ADST files above to rewrite into a single rewritten ADST file (first list):"), wxT("Rewrite"), ID_REWRITE, lwidth);
   vbox->Add(startRewriting->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_REWRITE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SelectRewrite));

   // Label and listbox for displaying rewriten root files
   vitemp->clear();
   vitemp->push_back(ID_DELETELIST + 4*nrlists);
   vitemp->push_back(ID_UPLIST + 4*nrlists);
   vitemp->push_back(ID_DOWNLIST + 4*nrlists);
   vitemp->push_back(ID_CLEARLIST + 4*nrlists);

   mvaList[1] = new LabelListEdit(leftmvapanel, wxT("Rewritten ADST files:"), lwidth, 70, 1, -1, vitemp);
   vbox->Add(mvaList[1]->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(vitemp->at(0), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(1), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(2), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(3), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   allLBE[nrlists] = mvaList[1];
   nrlists++;

   // Label + NEntry + three buttons
   vitemp->clear();
   vstemp->clear();
   vstemp->push_back("Combine (to MVA input)");
   vitemp->push_back(ID_COMBINE);
   vstemp->push_back("Merge (to rewritten ADST)");
   vitemp->push_back(ID_MERGE);

   startCombining = new LabelButton(leftmvapanel, wxT("Select rewritten ADST files above to combine into one rewritten ADST file (second list):"), vstemp, vitemp, lwidth);
   vbox->Add(startCombining->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_COMBINE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SelectCombine));
   Connect(ID_MERGE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SelectMerge));

   // Label and listbox for MVA input root files
   vitemp->clear();
   vitemp->push_back(ID_DELETELIST + 4*nrlists);
   vitemp->push_back(ID_UPLIST + 4*nrlists);
   vitemp->push_back(ID_DOWNLIST + 4*nrlists);
   vitemp->push_back(ID_CLEARLIST + 4*nrlists);

   mvaList[2] = new LabelListEdit(leftmvapanel, wxT("MVA input files:"), lwidth, 70, 0, ID_SELECTMVAFILE, vitemp);
   vbox->Add(mvaList[2]->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_SELECTMVAFILE, wxEVT_LISTBOX_DCLICK, wxCommandEventHandler(MyFrame::EnableMvaFile));
   Connect(vitemp->at(0), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(1), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(2), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(3), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   allLBE[nrlists] = mvaList[2];
   nrlists++;

   // Label for some extra information
   temptext = new wxStaticText(leftmvapanel, -1, wxT("Double click a file from the above list (third list) to use it in the MVA analysis."));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT | wxBOTTOM, 10);

   // Make a title
   MakeTitle(leftmvapanel, vbox, wxT("Rewritten ADST preparation"));

   // Label for splitting
   vbox->Add(-1, 10);
   temptext = new wxStaticText(leftmvapanel, -1, wxT("Settings for splitting rewritten ADST files into two separate files:"));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 10);

   // Label and combo box for selecting the type of observables for cuts SD or FD (splitting)
   vstemp->clear();
   vstemp->push_back("SD observables (energySD and zenithSD)");
   vstemp->push_back("FD observables (energyFD and zenithFD)");
   splitCutObservables = new LabelDrop(leftmvapanel, wxT("Select cut observables type:"), vstemp, vstemp->at(1), -1, rwidth);
   vbox->Add(splitCutObservables->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (splitCutObservables->widgetCB)->SetSelection(1);

   // Check + label + NEntry for choosing energy cut (splitting)
   vitemp->clear();
   vdtemp->clear();
   vdtemp->push_back(17.);
   vitemp->push_back(-1);
   vdtemp->push_back(21.);
   vitemp->push_back(-1);
   vdtemp->push_back(1.);
   vitemp->push_back(-1);
   splitCutEnergy = new CheckNEntry(leftmvapanel, false, wxT("Energy limits and bins:"), -1, vdtemp, vitemp, rwidth);
   splitCutEnergy->SetNEntryFormat(splitCutEnergy->widgetNE[0], 2, 0.01, 2, 0., 25.);
   splitCutEnergy->SetNEntryFormat(splitCutEnergy->widgetNE[1], 2, 0.01, 2, 0., 25.);
   splitCutEnergy->SetNEntryFormat(splitCutEnergy->widgetNE[2], 0, 1, 2, 1., 30.);
   vbox->Add(splitCutEnergy->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Check + label + NEntry for choosing zenith angle cut (splitting)
   vitemp->clear();
   vdtemp->clear();
   vdtemp->push_back(0.);
   vitemp->push_back(-1);
   vdtemp->push_back(60.);
   vitemp->push_back(-1);
   splitCutZenith = new CheckNEntry(leftmvapanel, false, wxT("Zenith angle limits:"), -1, vdtemp, vitemp, rwidth);
   splitCutZenith->SetNEntryFormat(splitCutZenith->widgetNE[0], 2, 0.01, 2, 0., 180.);
   splitCutZenith->SetNEntryFormat(splitCutZenith->widgetNE[1], 2, 0.01, 2, 0., 180.);
   vbox->Add(splitCutZenith->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Check + label + NEntry for choosing maximum risetime error cut (splitting)
   splitCutRisetime = new CheckNEntry(leftmvapanel, false, wxT("Maximum relative risetime limit:"), -1, 0.3, -1, rwidth);
   splitCutRisetime->SetNEntryFormat(splitCutRisetime->widgetNE[0], 3, 0.001, 2, 0., 10.);
   vbox->Add(splitCutRisetime->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Check + label for choosing if we want to enforce that the event also has a valid SD reconstruction
   vitemp->clear();
   vitemp2->clear();
   vstemp->clear();
   vstemp->push_back("Enforce that split events (in first file) have a valid SD reconstruction");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   splitEnforceSD = new CheckList(leftmvapanel, vitemp, vstemp, vitemp2, rwidth, "vertical");
   vbox->Add(splitEnforceSD->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

/*   // Label and combo box for selecting the eye selection method
     vstemp->clear();
   vstemp->push_back("Combine stereo FD events");
   vstemp->push_back("Any FD eye inside cut");
   vstemp->push_back("Average of active eyes");
   splitEyeSelection = new LabelDrop(leftmvapanel, wxT("Eye selection method, if more than one FD eye:"), vstemp, vstemp->at(0), -1, rwidth);
   vbox->Add(splitEyeSelection->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (splitEyeSelection->widgetCB)->SetSelection(0);*/

   // Label + NEntry + button for splitting rewritten ADST files
   startSplitting = new LabelNEntryButton(leftmvapanel, wxT("Split by a fraction (decimal) or number of events (integer):"), 0.5, -1, "Split", ID_SPLIT, lwidth);
   startSplitting->SetNEntryFormat(startSplitting->widgetNE[0], 3, 0.001, -1, -1.0, -1.);
   vbox->Add(startSplitting->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_SPLIT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::PrepareFileSplit));

   // Label + button for splitting rewritten ADST files according to a text file with split settings (numbers)
   startSplittingFile = new LabelButton(leftmvapanel, wxT("Split by file (decimals are fractions, integers are number of events):"), wxT("Split (File)"), ID_SPLITFILE, lwidth);
   vbox->Add(startSplittingFile->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_SPLITFILE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::PrepareFileSplitCustom));

   // Label for additional splitting instructions
   vbox->Add(-1, 5);
   temptext = new wxStaticText(leftmvapanel, -1, wxT("For splitting, use the following in the above number field:"));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 10);
   temptext = new wxStaticText(leftmvapanel, -1, wxT(" - decimal between 0.00 and 1.00 to split by fraction"));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT, 10);
   temptext = new wxStaticText(leftmvapanel, -1, wxT(" - integer to split by number of events"));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT, 10);
   temptext = new wxStaticText(leftmvapanel, -1, wxT(" - negative number to export to a single file"));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT, 10);
   vbox->Add(-1, 5);

   leftmvapanel->SetSizer(vbox);
   hboxhead->Add(leftmvapanel, 1, wxEXPAND);
   sw->SetSizer(hboxhead);
   // Set the left panel inside the scrollable window (MVA analysis) ------------------------------------
   
   // Set the right panel inside the scrollable window (MVA analysis) -----------------------------------
   wxPanel *rightmvapanel = new wxPanel(sw, -1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL | wxBORDER_RAISED);

   // Vertically align elements inside the panel
   vbox = new wxBoxSizer(wxVERTICAL);

   // Make a title
   MakeTitle(rightmvapanel, vbox, wxT("MVA analysis"));

   // Label and text entry for MVA selected input MVA file
   selectedMva = new LabelTEntry(rightmvapanel, wxT("Selected file for MVA analysis:"), wxT(""), -1, rwidth);
   (selectedMva->widgetTE)->SetEditable(false);
   vbox->Add(selectedMva->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Label and listbox for MVA input observables
   selectObservables = new LabelList(rightmvapanel, wxT("Observables to use in the MVA (Ctrl or Shift + Click to select multiple):"), 90, 1, ID_CHANGEOBSSELECT);
   vbox->Add(selectObservables->subsizer, 0, wxLEFT | wxRIGHT | wxTOP, 5);
   for(int i = 0; i < nrobs; i++)
   {
      (selectObservables->widgetLB)->Append(observables[i]);
      if(obsorigsel[i])
         (selectObservables->widgetLB)->SetSelection(i);
   }
   (selectObservables->widgetLB)->SetFirstItem(0);
   Connect(ID_CHANGEOBSSELECT, wxEVT_LISTBOX, wxCommandEventHandler(MyFrame::UpdateObservableSelection));

   // Label and combo box for selecting signal tree
   vitemp->clear();
   vstemp->clear();
   vstemp2->clear();
   for(int i = 0; i < nrobs; i++)
      vstemp->push_back(observables[i]);
   vstemp2->push_back("Apply negative");
   vitemp->push_back(ID_NUNCERTSELECT);
   vstemp2->push_back("Apply positive");
   vitemp->push_back(ID_PUNCERTSELECT);
   vstemp2->push_back("Disable");
   vitemp->push_back(ID_DISABLEUNCERT);
   uncertSelect = new LabelDropButton(rightmvapanel, wxT("Uncertainty shift:"), vstemp, vstemp->at(0), -1, vstemp2, vitemp, rwidth);
   vbox->Add(uncertSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (uncertSelect->widgetCB)->SetSelection(0);
   Connect(ID_NUNCERTSELECT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SetNegativeUncertainty));
   Connect(ID_PUNCERTSELECT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SetPositiveUncertainty));
   Connect(ID_DISABLEUNCERT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::DisableUncertainty));

   // Label and combo box for selecting signal tree
   vstemp->clear();
   vstemp->push_back("Select signal...");
   signalSelect = new LabelDrop(rightmvapanel, wxT("Select 'signal' tree:"), vstemp, vstemp->at(0), -1, rwidth);
   vbox->Add(signalSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (signalSelect->widgetCB)->SetSelection(0);

   // Label and combo box for selecting background tree
   vstemp->clear();
   vstemp->push_back("Select background...");
   backgroundSelect = new LabelDrop(rightmvapanel, wxT("Select 'background' tree:"), vstemp, vstemp->at(0), -1, rwidth);
   vbox->Add(backgroundSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (backgroundSelect->widgetCB)->SetSelection(0);

   // Label and combo box for selecting MVA method type
   vstemp->clear();
   for(int i = 0; i <= nrmethods; i++)
      vstemp->push_back(GetMethodName(methods[i]));
   methodsSelect = new LabelDrop(rightmvapanel, wxT("Choose MVA analysis method:"), vstemp, GetMethodName("Fisher"), -1, rwidth);
   vbox->Add(methodsSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (methodsSelect->widgetCB)->SetSelection((methodsSelect->widgetCB)->FindString(GetMethodName("Fisher")));

   // Label + NEntry for choosing the MVA cut
   cutMva = new LabelNEntry(rightmvapanel, wxT("Choose MVA cut value:"), 0., -1, rwidth);
   cutMva->SetNEntryFormat(cutMva->widgetNE[0], 4, 0.0001, 0, 0., 0.);
   vbox->Add(cutMva->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Label for cuts
   vbox->Add(-1, 10);
   temptext = new wxStaticText(rightmvapanel, -1, wxT("Cut and binning on energy, zenith angle and/or risetime:"));
   vbox->Add(temptext, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 10);

   // Label and combo box for selecting the type of observables for cuts SD or FD
   vstemp->clear();
   vstemp->push_back("SD observables (energySD and zenithSD)");
   vstemp->push_back("FD observables (energyFD and zenithFD)");
   cutObservables = new LabelDrop(rightmvapanel, wxT("Select cut observables type:"), vstemp, vstemp->at(1), -1, rwidth);
   vbox->Add(cutObservables->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (cutObservables->widgetCB)->SetSelection(1);

   // Check + label + NEntry for choosing energy cut
   vitemp->clear();
   vdtemp->clear();
   vdtemp->push_back(17.);
   vitemp->push_back(ID_ENERGYLIMITMIN);
   vdtemp->push_back(21.);
   vitemp->push_back(ID_ENERGYLIMITMAX);
   cutEnergy = new CheckNEntry(rightmvapanel, true, wxT("Energy limits:"), -1, vdtemp, vitemp, rwidth);
   cutEnergy->SetNEntryFormat(cutEnergy->widgetNE[0], 2, 0.01, 2, 0., 25.);
   cutEnergy->SetNEntryFormat(cutEnergy->widgetNE[1], 2, 0.01, 2, 0., 25.);
   vbox->Add(cutEnergy->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_ENERGYLIMITMIN, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(MyFrame::UpdateEnergyBinSelect));
   Connect(ID_ENERGYLIMITMAX, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(MyFrame::UpdateEnergyBinSelect));

   // Label + NEntry + dropbox + button for choosing energy binning options
   cutEnergyBins = new LabelNEntryDropButton(rightmvapanel, wxT("Energy binning:"), 1., ID_ENERGYBIN, vstemp, "", -1, wxT("Check bins"), ID_CHECKENERGYBINS, rwidth);
   cutEnergyBins->SetNEntryFormat(cutEnergyBins->widgetNE[0], 0, 1, 2, 1., 30.);
   vbox->Add(cutEnergyBins->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_ENERGYBIN, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(MyFrame::UpdateEnergyBinSelect));
   RunEnergyBinSelect();
   Connect(ID_CHECKENERGYBINS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::CheckEnergyBin));

   // Label + TEntry + button for choosing custom energy binning (from a file)
   vstemp->clear();
   vitemp->clear();
   vstemp->push_back("Select file...");
   vitemp->push_back(ID_ENERGYBINCUSTOM);
   vstemp->push_back("Disable");
   vitemp->push_back(ID_CUSTOMDISABLE);

   cutEnergyBinsCustom = new LabelTEntryButton(rightmvapanel, wxT("Custom energy binning:"), wxT(""), -1, vstemp, vitemp, rwidth);
   (cutEnergyBinsCustom->widgetTE)->SetEditable(false);
   vbox->Add(cutEnergyBinsCustom->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_ENERGYBINCUSTOM, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::CustomEnergyBins));
   Connect(ID_CUSTOMDISABLE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::CustomEnergyBinsDisable));
   customBinning = false;

   // Check + label + NEntry for choosing zenith angle cut
   vitemp->clear();
   vdtemp->clear();
   vdtemp->push_back(0.);
   vitemp->push_back(ID_ZENITHLIMITMIN);
   vdtemp->push_back(60.);
   vitemp->push_back(ID_ZENITHLIMITMAX);
   cutZenith = new CheckNEntry(rightmvapanel, true, wxT("Zenith angle limits:"), -1, vdtemp, vitemp, rwidth);
   cutZenith->SetNEntryFormat(cutZenith->widgetNE[0], 2, 0.01, 2, 0., 180.);
   cutZenith->SetNEntryFormat(cutZenith->widgetNE[1], 2, 0.01, 2, 0., 180.);
   vbox->Add(cutZenith->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_ZENITHLIMITMIN, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(MyFrame::UpdateZenithBinSelect));
   Connect(ID_ZENITHLIMITMAX, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(MyFrame::UpdateZenithBinSelect));

   // Label + NEntry + dropbox + button for choosing zenith angle binning options
   cutZenithBins = new LabelNEntryDropButton(rightmvapanel, wxT("Zenith binning:"), 1., ID_ZENITHBIN, vstemp, "", -1, wxT("Check bins"), ID_CHECKZENITHBINS, rwidth);
   cutZenithBins->SetNEntryFormat(cutZenithBins->widgetNE[0], 0, 1, 2, 1., 30.);
   vbox->Add(cutZenithBins->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_ZENITHBIN, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(MyFrame::UpdateZenithBinSelect));
   RunZenithBinSelect();
   Connect(ID_CHECKZENITHBINS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::CheckZenithBin));

   // Check + label + NEntry for choosing maximum risetime error cut
   cutRisetime = new CheckNEntry(rightmvapanel, false, wxT("Maximum relative risetime limit:"), -1, 0.3, -1, rwidth);
   cutRisetime->SetNEntryFormat(cutRisetime->widgetNE[0], 3, 0.001, 2, 0., 10.);
   vbox->Add(cutRisetime->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Label and combo box for selecting data tree
   vstemp->clear();
   vstemp->push_back("Select data...");
   dataSelect = new LabelDrop(rightmvapanel, wxT("Select 'data' tree:"), vstemp, vstemp->at(0), -1, rwidth);
   vbox->Add(dataSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (dataSelect->widgetCB)->SetSelection(0);

   vbox->Add(-1, 10);

   // Checkboxes for additional MVA settings
   vitemp->clear();
   vitemp2->clear();
   vstemp->clear();
   vstemp->push_back("Automatically run analysis over all energy bins");
   vitemp->push_back(1);
   vitemp2->push_back(-1);
   vstemp->push_back("Manually select an MVA cut after method training");
   vitemp->push_back(1);
   vitemp2->push_back(-1);
   vstemp->push_back("Determine MVA cut error through calculation of data tree standard deviation");
   vitemp->push_back(1);
   vitemp2->push_back(-1);
   vstemp->push_back("Apply bias correction for Auger FD standard data files (AugerWiki/XmaxHeatIcrc2017)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Apply bias correction for Auger HECO data files (AugerWiki/XmaxHeatIcrc2017)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Apply negative uncertainty to mean values (systematics estimation)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Apply positive uncertainty to mean values (systematics estimation)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Use published S38 fit (arXiv:1502.01323, page 65)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Apply atmospheric and alignment resolution smearing on Monte Carlo simulations");
   vitemp->push_back(1);
   vitemp2->push_back(-1);
   vstemp->push_back("Treat selected data tree as Monte Carlo simulations (for AugerMix)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   specialMva = new CheckList(rightmvapanel, vitemp, vstemp, vitemp2, rwidth, "vertical");
   vbox->Add(specialMva->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (specialMva->widgetChBox[5])->Disable();
   (specialMva->widgetChBox[6])->Disable();

   // Multiple buttons to start MVA analysis
   vitemp->clear();
   vstemp->clear();
   vstemp->push_back("Start MVA analysis");
   vitemp->push_back(ID_STARTMVA);
   vstemp->push_back("Create temp event file");
   vitemp->push_back(ID_TEMPFILE);
   vstemp->push_back("Apply MVA cut");
   vitemp->push_back(ID_MVACUT);
   vstemp->push_back("Check all bins");
   vitemp->push_back(ID_CHECKBINS);
   vstemp->push_back("Default options");
   vitemp->push_back(ID_DEFOPTIONS);
   startMva = new LabelButton(rightmvapanel, "", vstemp, vitemp, rwidth);
   vbox->Add(startMva->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_STARTMVA, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::StartMvaAnalysis));
   Connect(ID_TEMPFILE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::CreateTempEventFile));
   Connect(ID_MVACUT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::ApplyMvaCut));
   Connect(ID_CHECKBINS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::CheckBothBins));
   Connect(ID_DEFOPTIONS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SetDefaultMva));

   rightmvapanel->SetSizer(vbox);
   hboxhead->Add(rightmvapanel, 1, wxEXPAND);
   sw->SetSizer(hboxhead);
   // Set the right panel inside the scrollable window (MVA analysis) -----------------------------------

   // Event for when we close the window
//   Connect(wxEVT_CLOSE_WINDOW, wxCloseEventHandler(MyFrame::OnClose));

   // Set scroll bars
   sw->SetScrollbars(10, 10, 500/10, 350/10);
   sw->Scroll(0,0);

   // Add MVA analysis subframe to a tab
   tabs->AddPage(sw, wxT("MVA analysis"));

// MVA analysis -----------------------------------------------------

// Plotting page ----------------------------------------------------

   lwidth = 5*WW/6;

   // Set the scrollable window inside the frame
   sw = new wxScrolledWindow(tabs);

   // Horizontally align panels inside the scrolled window
   hboxhead = new wxBoxSizer(wxHORIZONTAL);

   // Set the plotting panel inside the scrollable window (plotting) ------------------------------------
   wxPanel *plotpanel = new wxPanel(sw, -1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL | wxBORDER_RAISED);

   // Vertically align elements inside the panel
   vbox = new wxBoxSizer(wxVERTICAL);

   // Make a title
   MakeTitle(plotpanel, vbox, wxT("General plotting settings"));

   // Label and button for opening files for plotting
   selectPlotFile = new LabelButton(plotpanel, wxT("File selection:"), wxT("..."), ID_OPENPLOTFILE, lwidth);
   vbox->Add(selectPlotFile->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_OPENPLOTFILE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SelectPlotFile));

   // Label and button for opening files for plotting
   selectPlotDir = new LabelButton(plotpanel, wxT("Directory selection:"), wxT("..."), ID_OPENPLOTDIR, lwidth);
   vbox->Add(selectPlotDir->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_OPENPLOTDIR, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SelectPlotDir));

   // Label and listbox for displaying opened plot files
   vitemp->clear();
   vitemp->push_back(ID_DELETELIST + 4*nrlists);
   vitemp->push_back(ID_UPLIST + 4*nrlists);
   vitemp->push_back(ID_DOWNLIST + 4*nrlists);
   vitemp->push_back(ID_CLEARLIST + 4*nrlists);

   openFileList = new LabelListEdit(plotpanel, wxT("Files for plotting:"), lwidth, 70, 1, -1, vitemp);
   vbox->Add(openFileList->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(vitemp->at(0), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(1), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(2), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   Connect(vitemp->at(3), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::EditList));
   allLBE[nrlists] = openFileList;
   nrlists++;

   // Label + NEntry for choosing the number of bins
   plotNrBins = new LabelNEntry(plotpanel, wxT("Choose number of bins (for histogram plots):"), 75., -1, lwidth);
   plotNrBins->SetNEntryFormat(plotNrBins->widgetNE[0], 0, 1, -1, 5., -1.);
   vbox->Add(plotNrBins->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Check + Label + NEntry for choosing x-axis range
   vitemp->clear();
   vdtemp->clear();
   vdtemp->push_back(-100.);
   vitemp->push_back(-1);
   vdtemp->push_back(100.);
   vitemp->push_back(-1);
   plotXaxisRange = new CheckNEntry(plotpanel, false, wxT("X-axis range:"), -1, vdtemp, vitemp, lwidth);
   plotXaxisRange->SetNEntryFormat(plotXaxisRange->widgetNE[0], 2, 0.01, 0, -1., -1.);
   plotXaxisRange->SetNEntryFormat(plotXaxisRange->widgetNE[1], 2, 0.01, 0, -1., -1.);
   (plotXaxisRange->widgetNE[0])->SetValue(-100.);
   vbox->Add(plotXaxisRange->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Check + Label + NEntry for choosing y-axis range
   vitemp->clear();
   vdtemp->clear();
   vdtemp->push_back(-100.);
   vitemp->push_back(-1);
   vdtemp->push_back(100.);
   vitemp->push_back(-1);
   plotYaxisRange = new CheckNEntry(plotpanel, false, wxT("Y-axis range:"), -1, vdtemp, vitemp, lwidth);
   plotYaxisRange->SetNEntryFormat(plotYaxisRange->widgetNE[0], 2, 0.01, 0, -1., -1.);
   plotYaxisRange->SetNEntryFormat(plotYaxisRange->widgetNE[1], 2, 0.01, 0, -1., -1.);
   (plotYaxisRange->widgetNE[0])->SetValue(-100.);
   vbox->Add(plotYaxisRange->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Label and combo box for selecting MVA method type to plot
   vstemp->clear();
   // All is set as the first method, so we skip it (but doesn't count in the nrmethods variable)
   for(int i = 0; i < nrmethods; i++)
      vstemp->push_back(GetMethodName(methods[i+1]));
   methodsPlotSelect = new LabelDrop(plotpanel, wxT("Chosen MVA analysis method:"), vstemp, GetMethodName("Fisher"), -1, lwidth);
   vbox->Add(methodsPlotSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (methodsPlotSelect->widgetCB)->SetSelection((methodsPlotSelect->widgetCB)->FindString(GetMethodName("Fisher")));

   vbox->Add(-1, 10);

   // Make a title
   MakeTitle(plotpanel, vbox, wxT("MVA histogram fit for lnA and composition plots"));

//   wxBoxSizer *hboxMid = new wxBoxSizer(wxHORIZONTAL);

   // Label and listbox for for selecting composition
   plottingList[0] = new LabelList(plotpanel, wxT("Composition to use during plotting (Ctrl or Shift + Click to select multiple):"), 90, 1, -1);
   vbox->Add(plottingList[0]->subsizer, 0, wxLEFT | wxRIGHT | wxTOP, 5);

   PrimPart *ptype = new PrimPart();

   vstemp->clear();
   for(int i = 0; i < ptype->Nr(); i++)
   {
      if(ptype->GetName(i) != "Other")
      {
         if(ptype->GetName(i) == "Proton")
            vstemp->push_back(ptype->GetName(i));
         if(ptype->GetName(i) == "Helium")
            vstemp->push_back(ptype->GetName(i));
         if(ptype->GetName(i) == "Oxygen")
            vstemp->push_back(ptype->GetName(i));
         if(ptype->GetName(i) == "Iron")
            vstemp->push_back(ptype->GetName(i));
         
         (plottingList[0]->widgetLB)->Append(ptype->GetName(i));
      }
   }

   for(int i = 0; i < vstemp->size(); i++)
      (plottingList[0]->widgetLB)->SetSelection(ptype->GetZ(&(vstemp->at(i))));
   (plottingList[0]->widgetLB)->SetFirstItem(0);

   delete ptype;

//   vbox->Add(hboxMid, 0, wxLEFT | wxRIGHT | wxTOP, 5);

   // Label and combo box for selecting the hadronic interaction model
   vstemp->clear();
   vstemp->push_back("EPOS_LHC");
   vstemp->push_back("QGSJET-II.04");
   vstemp->push_back("Sibyll-2.3");
   vstemp->push_back("None");
   hiSelect = new LabelDrop(plotpanel, wxT("Select hadronic interaction model:"), vstemp, vstemp->at(0), -1, lwidth);
   vbox->Add(hiSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (hiSelect->widgetCB)->SetSelection(0);

   // Label and combo box for selecting the hadronic interaction model
   vstemp->clear();
   vstemp->push_back("v2r9p5");
   vstemp->push_back("v3r3p4");
   prodSelect = new LabelDrop(plotpanel, wxT("Select simulation production:"), vstemp, vstemp->at(1), -1, lwidth);
   vbox->Add(prodSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (prodSelect->widgetCB)->SetSelection(1);

   // Label and combo box for selecting the data tree file
   vstemp->clear();
   vstemp->push_back("Select data...");
   dataPlotSelect = new LabelDrop(plotpanel, wxT("Select 'data' tree for histogram fit:"), vstemp, vstemp->at(0), -1, lwidth);
   vbox->Add(dataPlotSelect->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (dataPlotSelect->widgetCB)->SetSelection(0);

   // Checkboxes for additional plotting settings
   vitemp->clear();
   vitemp2->clear();
   vstemp->clear();
   vstemp->push_back("Use one of the fractions as constraint");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Fit Xmax instead of MVA");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Remove zero bins");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Use TFractionFitter");
   vitemp->push_back(1);
   vitemp2->push_back(-1);
   vstemp->push_back("Plot in a two column grid");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   specialPlot = new CheckList(plotpanel, vitemp, vstemp, vitemp2, lwidth, "vertical");
   vbox->Add(specialPlot->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Multiple buttons for plotting
   vitemp->clear();
   vstemp->clear();
/*   vstemp->push_back("Create histograms");
   vitemp->push_back(ID_PLOTHIST);
   vstemp->push_back("Create scatter plots");
   vitemp->push_back(ID_PLOTSCAT);*/
   vstemp->push_back("Start MVA histogram fit");
   vitemp->push_back(ID_PLOTMVAFIT);
   vstemp->push_back("Default options");
   vitemp->push_back(ID_PLOTDEFOPTIONS);
   startPlot = new LabelButton(plotpanel, "", vstemp, vitemp, lwidth);
   vbox->Add(startPlot->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
/*   Connect(ID_PLOTHIST, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::StartHistogramPlot));
   Connect(ID_PLOTSCAT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::StartScatterPlot));*/
   Connect(ID_PLOTMVAFIT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::StartMvaHistFit));
   Connect(ID_PLOTMVADEFOPTIONS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SetDefaultMvaPlot));

   vbox->Add(-1, 10);

   // Make a title
   MakeTitle(plotpanel, vbox, wxT("Distribution and scatter plot settings"));

   // Label and listbox for observables to plot
   plottingList[1] = new LabelList(plotpanel, wxT("Observables to plot (Ctrl or Shift + Click to select multiple):"), 90, 1, -1);
   vbox->Add(plottingList[1]->subsizer, 0, wxLEFT | wxRIGHT | wxTOP, 5);
   for(int i = 0; i < nrobs; i++)
   {
      (plottingList[1]->widgetLB)->Append(observables[i]);
      if(obsorigsel[i])
         (plottingList[1]->widgetLB)->SetSelection(i);
   }
   (plottingList[1]->widgetLB)->Append("MVA");
   (plottingList[1]->widgetLB)->SetFirstItem(0);

   // Label and combo box for selecting the first tree file
   vstemp->clear();
   vstemp->push_back("Disabled");
   dropHistTree[0] = new LabelDrop(plotpanel, wxT("Select first tree for plotting (signal colors):"), vstemp, vstemp->at(0), -1, lwidth);
   vbox->Add(dropHistTree[0]->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (dropHistTree[0]->widgetCB)->SetSelection(0);

   // Label and combo box for selecting the second tree file
   vstemp->clear();
   vstemp->push_back("Disabled");
   dropHistTree[1] = new LabelDrop(plotpanel, wxT("Select second tree for plotting (background colors):"), vstemp, vstemp->at(0), -1, lwidth);
   vbox->Add(dropHistTree[1]->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (dropHistTree[1]->widgetCB)->SetSelection(0);

   // Label and combo box for selecting the third tree file
   vstemp->clear();
   vstemp->push_back("Disabled");
   dropHistTree[2] = new LabelDrop(plotpanel, wxT("Select third tree for plotting (data colors):"), vstemp, vstemp->at(0), -1, lwidth);
   vbox->Add(dropHistTree[2]->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   (dropHistTree[2]->widgetCB)->SetSelection(0);

   // Checkboxes for additional plotting settings
   vitemp->clear();
   vitemp2->clear();
   vstemp->clear();
/*   vstemp->push_back("Merge multiple trees to one plot (select in dropdown menus)");
   vitemp->push_back(1);
   vitemp2->push_back(-1);*/
   vstemp->push_back("Perform a comparison of different observables (select in listbox)");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Use negative uncertainties instead of mean values");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Use positive uncertainties instead of mean values");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
/*   vstemp->push_back("Fit Xmax instead of MVA");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Remove zero bins");
   vitemp->push_back(0);
   vitemp2->push_back(-1);
   vstemp->push_back("Use FractionFitter");
   vitemp->push_back(1);
   vitemp2->push_back(-1);*/
   specialHistPlot = new CheckList(plotpanel, vitemp, vstemp, vitemp2, lwidth, "vertical");
   vbox->Add(specialHistPlot->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);

   // Multiple buttons for plotting
   vitemp->clear();
   vstemp->clear();
   vstemp->push_back("Create histograms");
   vitemp->push_back(ID_PLOTHIST);
   vstemp->push_back("Create scatter plots");
   vitemp->push_back(ID_PLOTSCAT);
/*   vstemp->push_back("Start MVA histogram fit");
   vitemp->push_back(ID_PLOTMVAFIT);*/
   vstemp->push_back("Default options");
   vitemp->push_back(ID_PLOTHISTDEFOPTIONS);
   startHistPlot = new LabelButton(plotpanel, "", vstemp, vitemp, lwidth);
   vbox->Add(startHistPlot->subsizer, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   Connect(ID_PLOTHIST, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::StartHistogramPlot));
   Connect(ID_PLOTSCAT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::StartScatterPlot));
   Connect(ID_PLOTHISTDEFOPTIONS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MyFrame::SetDefaultHistPlot));

   plotpanel->SetSizer(vbox);
   hboxhead->Add(plotpanel, 1, wxEXPAND);
   sw->SetSizer(hboxhead);

   // Set scroll bars
   sw->SetScrollbars(10, 10, 500/10, 350/10);
   sw->Scroll(0,0);

   // Add ADST cuts subframe to a tab
   tabs->AddPage(sw, wxT("Plotting"));

// Plotting page ----------------------------------------------------

//   CreateStatusBar();
   Centre();

   cout << "# MyFrame               #: " << " Deleting temporary variables" << endl;
   delete vstemp;
   delete vstemp2;
   delete vitemp;
   delete vitemp2;
   delete vdtemp;
}
