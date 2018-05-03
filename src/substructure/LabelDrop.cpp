#include "substructure/subheader.h"
#include "substructure/LabelDrop.h"
#include "workstation.h"

LabelDrop::LabelDrop(wxPanel *parent, wxString label, vector<string> *entrytext, string selecttext, const int dropID, int maxsize)
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

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetDropBoxWidth(widgetCB) + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       vbox->Add(widgetCB, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetCB, 1, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
