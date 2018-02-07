#include "substructure/subheader.h"
#include "substructure/LabelNEntryDropButton.h"
#include "workstation.h"
#include <iostream>

// Single NEntry and single button
LabelNEntryDropButton::LabelNEntryDropButton(wxPanel *parent, wxString label, double numval, const int nentryID, vector<string> entrytext, string selecttext, const int dropID, wxString buttext, const int butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelNEntryDropButton #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a number entry
    widgetNE[0] = new wxSpinCtrlDouble(parent, nentryID);
    widgetNE[0]->SetValue(numval);
    if(DBGSIG > 1)
       cout << "# LabelNEntryDropButton #: " << "Number entry width = " << GetNumEntryWidth(widgetNE[0]) << ", Number entry height = " << GetNumEntryHeight(widgetNE[0]) << endl;
    // Create a dropbox
    widgetCB = new wxChoice(parent, dropID);
    if(DBGSIG > 1)
       cout << "# LabelNEntryDropButton #: " << "DropBox width = " << GetDropBoxWidth(widgetCB) << ", DropBox height = " << GetDropBoxHeight(widgetCB) << endl;
    for(int i = 0; i < entrytext.size(); i++)
       widgetCB->Append(entrytext[i]);
    // Create a button
    widgetTB[0] = new wxButton(parent, butID, buttext);
    if(DBGSIG > 1)
       cout << "# LabelNEntryDropButton #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << ", TextButton text height = " << GetButtonHeight(widgetTB[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetNumEntryWidth(widgetNE[0]) + GetDropBoxWidth(widgetCB) + GetButtonWidth(widgetTB[0]) + 8*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       hbox->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
       hbox->Add(widgetCB, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);
       hbox->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetNE[0], 0, wxLEFT | wxRIGHT, padding);
       subsizer->Add(widgetCB, 1, wxLEFT | wxRIGHT, padding);
       subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
