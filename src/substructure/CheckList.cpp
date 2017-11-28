#include "substructure/subheader.h"
#include "substructure/CheckList.h"
#include "workstation.h"
#include <iostream>

// Single checkbox
CheckList::CheckList(wxPanel *parent, bool checked, wxString label, const int checkID)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);

    // Create a checkbox
    widgetChBox[0] = new wxCheckBox(parent, checkID, label);
    widgetChBox[0]->SetValue(checked);
    if(DBGSIG > 1)
       cout << "# CheckList             #: " << "Check box width = " << GetCheckBoxWidth(widgetChBox[0]) << ", Check box height = " << GetCheckBoxHeight(widgetChBox[0]) << endl;

    subsizer->Add(widgetChBox[0], 0, wxLEFT | wxTOP | wxRIGHT, padding);
}

// Multiple checkboxes
CheckList::CheckList(wxPanel *parent, vector<int> checked, vector<string> label, vector<int> checkID, int maxsize)
{
    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    // Create a checkbox
    int checksize = 0;
    for(int i = 0; i < checked.size(); i++)
    {
       widgetChBox[i] = new wxCheckBox(parent, (const int)checkID[i], label[i]);
       widgetChBox[i]->SetValue(checked[i]);
       if(DBGSIG > 1)
          cout << "# CheckList             #: " << "Check box width = " << GetCheckBoxWidth(widgetChBox[i]) << ", Check box height = " << GetCheckBoxHeight(widgetChBox[i]) << endl;
       checksize += GetCheckBoxWidth(widgetChBox[i]) + 2*padding;
    }

    if(checksize > maxsize)
    {
       *isbigger = true;

       if(checksize/2 > maxsize)
       {
          wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
          for(int i = 0; i < checked.size()/3; i++)
             hbox->Add(widgetChBox[i], 0, wxLEFT | wxRIGHT, padding);
          vbox->Add(hbox, 1);

          hbox = new wxBoxSizer(wxHORIZONTAL);
          for(int i = checked.size()/3; i < 2*checked.size()/3; i++)
             hbox->Add(widgetChBox[i], 0, wxLEFT | wxRIGHT, padding);
          vbox->Add(hbox, 1);

          hbox = new wxBoxSizer(wxHORIZONTAL);
          for(int i = 2*checked.size()/3; i < checked.size(); i++)
             hbox->Add(widgetChBox[i], 0, wxLEFT | wxRIGHT, padding);
          vbox->Add(hbox, 1);
       }
       else
       {
          wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
          for(int i = 0; i < checked.size()/2; i++)
             hbox->Add(widgetChBox[i], 0, wxLEFT | wxRIGHT, padding);
          vbox->Add(hbox, 1);

          hbox = new wxBoxSizer(wxHORIZONTAL);
          for(int i = checked.size()/2; i < checked.size(); i++)
             hbox->Add(widgetChBox[i], 0, wxLEFT | wxRIGHT, padding);
          vbox->Add(hbox, 1);
       }

       subsizer->Add(vbox, 1);
    }
    else
    {
       for(int i = 0; i < checked.size(); i++)
          subsizer->Add(widgetChBox[i], 0, wxLEFT | wxTOP | wxRIGHT, padding);
    }

    if(!(*isbigger))
       delete vbox;
}
