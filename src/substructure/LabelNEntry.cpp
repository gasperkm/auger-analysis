#include "substructure/subheader.h"
#include "substructure/LabelNEntry.h"
#include "workstation.h"
//#include "separate_functions.h"
#include <iostream>

// Single NEntry
LabelNEntry::LabelNEntry(wxPanel *parent, wxString label, double numval, const int nentryID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntry           #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a number entry
    widgetNE[0] = new wxSpinCtrlDouble(parent, nentryID);
    widgetNE[0]->SetValue(numval);
    if(DBGSIG > 1)
       cout << "# LabelNEntry           #: " << "Number entry width = " << GetNumEntryWidth(widgetNE[0]) << ", Number entry height = " << GetNumEntryHeight(widgetNE[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetNumEntryWidth(widgetNE[0]) + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       hbox->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple NEntry
LabelNEntry::LabelNEntry(wxPanel *parent, wxString label, vector<double> *numval, vector<int> *nentryID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntry           #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create all number entries
    int nesize = 0;
    for(int i = 0; i < numval->size(); i++)
    {
       widgetNE[i] = new wxSpinCtrlDouble(parent, (const int)nentryID->at(i));
       widgetNE[i]->SetValue(numval->at(i));
       if(DBGSIG > 1)
          cout << "# LabelNEntry           #: " << "Number entry " << i+1 << " width = " << GetNumEntryWidth(widgetNE[i]) << ", Number entry " << i+1 << " height = " << GetNumEntryHeight(widgetNE[i]) << endl;
       nesize += GetNumEntryWidth(widgetNE[i]) + 2*padding;
    }

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + nesize + 2*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       for(int i = 0; i < numval->size(); i++)
          hbox->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       for(int i = 0; i < numval->size(); i++)
          subsizer->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
