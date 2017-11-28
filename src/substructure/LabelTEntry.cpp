#include "substructure/subheader.h"
#include "substructure/LabelTEntry.h"
#include "workstation.h"

LabelTEntry::LabelTEntry(wxPanel *parent, wxString label, wxString deftext, const int textID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    if(DBGSIG > 1)
       cout << "# LabelTEntry           #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
    // Create a text entry
    widgetTE = new wxTextCtrl(parent, textID, deftext);
    if(DBGSIG > 1)
       cout << "# LabelTEntry           #: " << "TextEntry text width = " << GetTextEntryWidth(widgetTE) << ", TextEntry text height = " << GetTextEntryHeight(widgetTE) << endl;

    // Determine their sizes and arrange them accordingly
    if(GetLabelWidth(widgetLabel) + GetTextEntryWidth(widgetTE) + 4*padding > maxsize)
    {
       *isbigger = true;
       vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

       if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
          vbox->Add(-1, GetLabelHeight(widgetLabel));

       vbox->Add(widgetTE, 1, wxEXPAND | wxLEFT | wxRIGHT, padding);
       subsizer->Add(vbox, 1);
    }
    else
    {
       subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
       subsizer->Add(widgetTE, 1, wxLEFT | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
