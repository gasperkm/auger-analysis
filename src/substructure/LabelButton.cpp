#include "substructure/subheader.h"
#include "substructure/LabelButton.h"
#include "workstation.h"
#include <iostream>

// Single button
LabelButton::LabelButton(wxPanel *parent, wxString label, wxString buttext, const int butID, int maxsize)
{
    bool nolabel = false;
    if(label.IsEmpty())
       nolabel = true;

    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    if(nolabel)
    {
       // Create a button
       widgetTB[0] = new wxButton(parent, butID, buttext);
       if(DBGSIG > 1)
          cout << "# LabelButton           #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << ", TextButton text height = " << GetButtonHeight(widgetTB[0]) << endl;

       subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
    }
    else
    {
       // Create a label
       widgetLabel = new wxStaticText(parent, -1, label);
       if(DBGSIG > 1)
          cout << "# LabelButton           #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
       // Create a button
       widgetTB[0] = new wxButton(parent, butID, buttext);
       if(DBGSIG > 1)
          cout << "# LabelButton           #: " << "TextButton text width = " << GetButtonWidth(widgetTB[0]) << ", TextButton text height = " << GetButtonHeight(widgetTB[0]) << endl;

       // Determine their sizes and arrange them accordingly
       if(GetLabelWidth(widgetLabel) + GetButtonWidth(widgetTB[0]) + 4*padding > maxsize)
       {
          *isbigger = true;
          vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
    
          if(GetLabelWidth(widgetLabel) + 4*padding > maxsize)
             vbox->Add(-1, GetLabelHeight(widgetLabel));
    
          vbox->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
          subsizer->Add(vbox, 1);
       }
       else
       {
          subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
          subsizer->Add(widgetTB[0], 0, wxLEFT | wxRIGHT, padding);
       }
    }

    if(!(*isbigger))
       delete vbox;
}

// Multiple buttons
LabelButton::LabelButton(wxPanel *parent, wxString label, vector<string> *buttext, vector<int> *butID, int maxsize)
{
    bool nolabel = false;
    if(label.IsEmpty())
       nolabel = true;

    subsizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

    if(nolabel)
    {
       // Create all buttons
       int butsize = 0;
       for(int i = 0; i < buttext->size(); i++)
       {
          widgetTB[i] = new wxButton(parent, (const int)butID->at(i), buttext->at(i));
          if(DBGSIG > 1)
             cout << "# LabelButton           #: " << "TextButton " << i+1 << " text width = " << GetButtonWidth(widgetTB[i]) << ", TextButton " << i+1 << " text height = " << GetButtonHeight(widgetTB[i]) << endl;
          butsize += GetButtonWidth(widgetTB[i]) + 2*padding;

          subsizer->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
       }
    }
    else
    {
       // Create a label
       widgetLabel = new wxStaticText(parent, -1, label);
       if(DBGSIG > 1)
          cout << "# LabelButton           #: " << "Label text width = " << GetLabelWidth(widgetLabel) << ", Label text height = " << GetLabelHeight(widgetLabel) << endl;
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
       if(GetLabelWidth(widgetLabel) + butsize + 2*padding > maxsize)
       {
          *isbigger = true;
          vbox->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);

          if(GetLabelWidth(widgetLabel) + 2*padding > maxsize)
             vbox->Add(-1, GetLabelHeight(widgetLabel));

          wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);

          for(int i = 0; i < buttext->size(); i++)
             hbox->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);

          vbox->Add(hbox, 1);
          subsizer->Add(vbox, 1);
       }
       else
       {
          subsizer->Add(widgetLabel, 0, wxLEFT | wxTOP, padding);
          for(int i = 0; i < buttext->size(); i++)
             subsizer->Add(widgetTB[i], 0, wxLEFT | wxRIGHT, padding);
       }
    }

    if(!(*isbigger))
       delete vbox;
}
