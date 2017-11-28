#include "substructure/subheader.h"
#include "substructure/LabelList.h"
#include <iostream>

LabelList::LabelList(wxPanel *parent, wxString label, int height, int multi, const int listID)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    // Create a listbox
    if(multi == 1)
       widgetLB = new wxListBox(parent, listID, wxDefaultPosition, wxSize(-1, height), 0, NULL, wxLB_EXTENDED);
    else
       widgetLB = new wxListBox(parent, listID, wxDefaultPosition, wxSize(-1, height), 0, NULL, wxLB_SINGLE);

    vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
    vbox->Add(widgetLB, 0, wxEXPAND | wxALL, padding);
    subsizer->Add(vbox, 1);
}
