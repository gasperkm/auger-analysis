#include "substructure/subheader.h"
#include "substructure/CheckNEntry.h"
#include "workstation.h"
#include <iostream>

// Single NEntry
CheckNEntry::CheckNEntry(wxPanel *parent, bool checked, wxString label, const int checkID, double numval, const int nentryID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a checkbox
    widgetChBox = new wxCheckBox(parent, checkID, label);
    if(DBGSIG > 1)
       cout << "# CheckNEntry           #: " << "Check box width = " << GetCheckBoxWidth(widgetChBox) << ", Check box height = " << GetCheckBoxHeight(widgetChBox) << endl;
    widgetChBox->SetValue(checked);
    // Create a number entry
    widgetNE[0] = new wxSpinCtrlDouble(parent, nentryID);
    widgetNE[0]->SetValue(numval);
    if(DBGSIG > 1)
       cout << "# CheckNEntry           #: " << "Number entry width = " << GetNumEntryWidth(widgetNE[0]) << ", Number entry height = " << GetNumEntryHeight(widgetNE[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetCheckBoxWidth(widgetChBox) + GetNumEntryWidth(widgetNE[0]) + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetChBox, 0, wxLEFT | wxTOP, padding);

       if(GetCheckBoxWidth(widgetChBox) + 2*padding > maxsize)
          vbox->Add(-1, GetCheckBoxHeight(widgetChBox));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       hbox->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetChBox, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple NEntry
CheckNEntry::CheckNEntry(wxPanel *parent, bool checked, wxString label, const int checkID, vector<double> *numval, vector<int> *nentryID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a check box
    widgetChBox = new wxCheckBox(parent, checkID, label);
    if(DBGSIG > 1)
       cout << "# CheckNEntry           #: " << "Check box width = " << GetCheckBoxWidth(widgetChBox) << ", Check box height = " << GetCheckBoxHeight(widgetChBox) << endl;
    widgetChBox->SetValue(checked);
    // Create all number entries
    int nesize = 0;
    for(int i = 0; i < numval->size(); i++)
    {
       widgetNE[i] = new wxSpinCtrlDouble(parent, (const int)nentryID->at(i));
       widgetNE[i]->SetValue(numval->at(i));
       if(DBGSIG > 1)
          cout << "# CheckNEntry           #: " << "Number entry " << i+1 << " width = " << GetNumEntryWidth(widgetNE[i]) << ", Number entry " << i+1 << " height = " << GetNumEntryHeight(widgetNE[i]) << endl;
       nesize += GetNumEntryWidth(widgetNE[i]) + 2*padding;
    }

    // Determine their sizes and arrange them accordingly
    if(GetCheckBoxWidth(widgetChBox) + nesize + 2*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetChBox, 0, wxLEFT | wxTOP, padding);

       if(GetCheckBoxWidth(widgetChBox) + 2*padding > maxsize)
          vbox->Add(-1, GetCheckBoxHeight(widgetChBox));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       for(int i = 0; i < numval->size(); i++)
          hbox->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetChBox, 0, wxLEFT | wxTOP, padding);
       for(int i = 0; i < numval->size(); i++)
          subsizer->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
