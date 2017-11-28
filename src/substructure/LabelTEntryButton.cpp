#include "substructure/subheader.h"
#include "substructure/LabelTEntryButton.h"
#include "workstation.h"
#include <iostream>

// Single button
LabelTEntryButton::LabelTEntryButton(wxPanel *parent, wxString label, wxString deftext, const int textID, wxString buttext, const int butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelTEntryButton     #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a text entry
    widgetTE = new wxTextCtrl(parent, textID, deftext);
    if(DBGSIG > 1)
       cout << "# LabelTEntryButton     #: " << "TextEntry text width = " << GetTextEntryWidth(widgetTE) << endl;
    // Create a button
    widgetTB[0] = new wxButton(parent, butID, buttext);
    if(DBGSIG > 1)
       cout << "# LabelTEntryButton     #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetTextEntryWidth(widgetTE) + GetButtonWidth(widgetTB[0]) + 6*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       vbox->Add(widgetTE, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);
       vbox->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetTE, 1, wxLEFT | wxRIGHT, padding);
       subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple buttons
LabelTEntryButton::LabelTEntryButton(wxPanel *parent, wxString label, wxString deftext, const int textID, vector<string> buttext, vector<int> butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelTEntryButton     #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a text entry
    widgetTE = new wxTextCtrl(parent, textID, deftext);
    if(DBGSIG > 1)
       cout << "# LabelTEntryButton     #: " << "TextEntry text width = " << GetTextEntryWidth(widgetTE) << endl;
    // Create a button
    int butsize = 0;
    for(int i = 0; i < buttext.size(); i++)
    {
       widgetTB[i] = new wxButton(parent, (const int)butID[i], buttext[i]);
       if(DBGSIG > 1)
          cout << "# LabelTEntryButton     #: " << "TextButton " << i+1 << " text width = " << GetButtonWidth(widgetTB[i]) << ", TextButton " << i+1 << " text height = " << GetButtonHeight(widgetTB[i]) << endl;
       butsize += GetButtonWidth(widgetTB[i]) + 2*padding;
    }

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetTextEntryWidth(widgetTE) + butsize + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       vbox->Add(widgetTE, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
       for(int i = 0; i < buttext.size(); i++)
          hbox->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetTE, 1, wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext.size(); i++)
          subsizer->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
