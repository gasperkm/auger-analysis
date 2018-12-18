#include "substructure/subheader.h"
#include "substructure/LabelDropButton.h"
#include "workstation.h"

// Single button
LabelDropButton::LabelDropButton(wxPanel *parent, wxString label, vector<string> *entrytext, string selecttext, const int dropID, wxString buttext, const int butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);
    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelDrop             #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a dropbox
    widgetCB = new wxChoice(parent, dropID);
    if(DBGSIG > 1)
       cout << "# LabelDrop             #: " << "DropBox width = " << GetDropBoxWidth(widgetCB) << ", DropBox height = " << GetDropBoxHeight(widgetCB) << endl;
    for(int i = 0; i < entrytext->size(); i++)
       widgetCB->Append(entrytext->at(i));
    // Create a button
    widgetTB[0] = new wxButton(parent, butID, buttext);
    if(DBGSIG > 1)
       cout << "# LabelButton           #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << ", TextButton text height = " << GetButtonHeight(widgetTB[0]) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetDropBoxWidth(widgetCB) + GetButtonWidth(widgetTB[0]) + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       hbox->Add(widgetCB, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);
       hbox->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetCB, 1, wxLEFT | wxRIGHT, padding);
       subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple buttons
LabelDropButton::LabelDropButton(wxPanel *parent, wxString label, vector<string> *entrytext, string selecttext, const int dropID, vector<string> *buttext, vector<int> *butID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);
    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelDrop             #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a dropbox
    widgetCB = new wxChoice(parent, dropID);
    if(DBGSIG > 1)
       cout << "# LabelDrop             #: " << "DropBox width = " << GetDropBoxWidth(widgetCB) << ", DropBox height = " << GetDropBoxHeight(widgetCB) << endl;
    for(int i = 0; i < entrytext->size(); i++)
       widgetCB->Append(entrytext->at(i));
    // Create all buttons
    int butsize = 0;
    for(int i = 0; i < buttext->size(); i++)
    {
       widgetTB[i] = new wxButton(parent, (const int)butID->at(i), buttext->at(i));
       if(DBGSIG > 1)
          cout << "# LabelButton           #: " << "TextButton " << i+1 << " text width = " << GetButtonWidth(widgetTB[i]) << ", TextButton " << i+1 << " text height = " << GetButtonHeight(widgetTB[i]) << endl;
       butsize += GetButtonWidth(widgetTB[i]) + 2*padding;
    }

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetDropBoxWidth(widgetCB) + butsize + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

       hbox->Add(widgetCB, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext->size(); i++)
          hbox->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);

       vbox->Add(hbox, 1);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetCB, 1, wxLEFT | wxRIGHT, padding);
       for(int i = 0; i < buttext->size(); i++)
          subsizer->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
