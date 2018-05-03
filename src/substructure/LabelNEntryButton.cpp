#include "substructure/subheader.h"
#include "substructure/LabelNEntryButton.h"
#include "workstation.h"
#include <iostream>

// Single NEntry + single button
LabelNEntryButton::LabelNEntryButton(wxPanel *parent, wxString label, double numval, const int nentryID, wxString buttext, const int butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a number entry
    widgetNE[0] = new wxSpinCtrlDouble(parent, nentryID);
    widgetNE[0]->SetValue(numval);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "Number entry width = " << GetNumEntryWidth(widgetNE[0]) << ", Number entry height = " << GetNumEntryHeight(widgetNE[0]) << endl;
    // Create a button
    widgetTB[0] = new wxButton(parent, butID, buttext);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << ", TextButton text height = " << GetButtonHeight(widgetTB[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetButtonWidth(widgetTB[0]) + GetNumEntryWidth(widgetNE[0]) + 6*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       hbox->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
       hbox->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
       subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple NEntry + single button
LabelNEntryButton::LabelNEntryButton(wxPanel *parent, wxString label, vector<double> *numval, vector<int> *nentryID, wxString buttext, const int butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a number entry
    int nesize = 0;
    for(int i = 0; i < numval->size(); i++)
    {
       widgetNE[i] = new wxSpinCtrlDouble(parent, (const int)nentryID->at(i));
       widgetNE[i]->SetValue(numval->at(i));
       if(DBGSIG > 1)
          cout << "# LabelNEntryButton     #: " << "Number entry " << i+1 << " width = " << GetNumEntryWidth(widgetNE[i]) << ", Number entry " << i+1 << " height = " << GetNumEntryHeight(widgetNE[i]) << endl;
       nesize += GetNumEntryWidth(widgetNE[i]) + 2*padding;
    }
    // Create a button
    widgetTB[0] = new wxButton(parent, butID, buttext);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << ", TextButton text height = " << GetButtonHeight(widgetTB[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetButtonWidth(widgetTB[0]) + nesize + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       for(int i = 0; i < numval->size(); i++)
          hbox->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);
       hbox->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       for(int i = 0; i < numval->size(); i++)
          subsizer->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);
       subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Single NEntry + multiple buttons
LabelNEntryButton::LabelNEntryButton(wxPanel *parent, wxString label, double numval, const int nentryID, vector<string> *buttext, vector<int> *butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a number entry
    widgetNE[0] = new wxSpinCtrlDouble(parent, nentryID);
    widgetNE[0]->SetValue(numval);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "Number entry width = " << GetNumEntryWidth(widgetNE[0]) << ", Number entry height = " << GetNumEntryHeight(widgetNE[0]) << endl;
    // Create all button
    int butsize = 0;
    for(int i = 0; i < buttext->size(); i++)
    {
       widgetTB[i] = new wxButton(parent, (const int)butID->at(i), buttext->at(i));
       if(DBGSIG > 1)
          cout << "# LabelNEntryButton     #: " << "TextButton " << i+1 << " text width = " << GetButtonWidth(widgetTB[i]) << ", TextButton " << i+1 << " text height = " << GetButtonHeight(widgetTB[i]) << endl;
       butsize += GetButtonWidth(widgetTB[i]) + 2*padding;
    }

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetNumEntryWidth(widgetNE[0]) + butsize + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
       hbox->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext->size(); i++)
          hbox->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext->size(); i++)
          subsizer->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple NEntry + multiple buttons
LabelNEntryButton::LabelNEntryButton(wxPanel *parent, wxString label, vector<double> *numval, vector<int> *nentryID, vector<string> *buttext, vector<int> *butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntryButton     #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a number entry
    int nesize = 0;
    for(int i = 0; i < numval->size(); i++)
    {
       widgetNE[i] = new wxSpinCtrlDouble(parent, (const int)nentryID->at(i));
       widgetNE[i]->SetValue(numval->at(i));
       if(DBGSIG > 1)
          cout << "# LabelNEntryButton     #: " << "Number entry " << i+1 << " width = " << GetNumEntryWidth(widgetNE[i]) << ", Number entry " << i+1 << " height = " << GetNumEntryHeight(widgetNE[i]) << endl;
       nesize += GetNumEntryWidth(widgetNE[i]) + 2*padding;
    }
    // Create all button
    int butsize = 0;
    for(int i = 0; i < buttext->size(); i++)
    {
       widgetTB[i] = new wxButton(parent, (const int)butID->at(i), buttext->at(i));
       if(DBGSIG > 1)
          cout << "# LabelNEntryButton     #: " << "TextButton " << i+1 << " text width = " << GetButtonWidth(widgetTB[i]) << ", TextButton " << i+1 << " text height = " << GetButtonHeight(widgetTB[i]) << endl;
       butsize += GetButtonWidth(widgetTB[i]) + 2*padding;
    }

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + nesize + butsize + 2*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
       for(int i = 0; i < numval->size(); i++)
          hbox->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext->size(); i++)
          hbox->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       for(int i = 0; i < numval->size(); i++)
          subsizer->Add(widgetNE[i], 0, wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext->size(); i++)
          subsizer->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
