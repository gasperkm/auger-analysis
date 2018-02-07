#include "substructure/subheader.h"
#include "substructure/LabelListEdit.h"
#include "workstation.h"
#include <iostream>

LabelListEdit::LabelListEdit(wxPanel *parent, wxString label, int width, int height, int multi, const int listID, vector<int> butID)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *hbox1 = new wxBoxSizer(wxHORIZONTAL);

    // Create a label
    widgetLabel = new wxStaticText(parent, -1, label);
    hbox1->Add(widgetLabel, 1, wxLEFT | wxTOP, padding);
    // Create a button
    wxBoxSizer *hbox2 = new wxBoxSizer(wxHORIZONTAL);
    widgetBM = new wxBitmap();
    for(int i = 0; i < 4; i++)
    {
       widgetTB[i] = new wxButton(parent, butID[i], "", wxDefaultPosition, wxSize(28,28));

       if(i == 0)
          widgetBM->LoadFile(string(rootdir) + "/input/delete.png");
       else if(i == 1)
          widgetBM->LoadFile(string(rootdir) + "/input/up-arrow.png");
       else if(i == 2)
          widgetBM->LoadFile(string(rootdir) + "/input/down-arrow.png");
       else if(i == 3)
          widgetBM->LoadFile(string(rootdir) + "/input/clear.png");

       widgetTB[i]->SetBitmap(*widgetBM, wxLEFT);
       hbox2->Add(widgetTB[i], 1, wxLEFT | wxRIGHT, 1);
    }
    hbox1->Add(hbox2, 1, wxALIGN_RIGHT, 0);
    // Create a listbox
    if(multi == 1)
       widgetLB = new wxListBox(parent, listID, wxDefaultPosition, wxSize(-1, height), 0, NULL, wxLB_EXTENDED);
    else
       widgetLB = new wxListBox(parent, listID, wxDefaultPosition, wxSize(-1, height), 0, NULL, wxLB_SINGLE);

    vbox->Add(hbox1, 0, wxEXPAND | wxALL, padding);
    vbox->Add(widgetLB, 0, wxEXPAND | wxLEFT | wxBOTTOM | wxRIGHT, padding);
    subsizer->Add(vbox, 1);
}
