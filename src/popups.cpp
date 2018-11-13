#include "popups.h"

using namespace std;

void ErrorPopup(string title, string descr)
{
/*   wxMessageDialog *dial = new wxMessageDialog(NULL, descr, title, wxOK | wxICON_ERROR);
   dial->ShowModal();*/
   wxMessageDialog dial(NULL, descr, title, wxOK | wxICON_ERROR);
   dial.ShowModal();
}

void InfoPopup(string title, string descr)
{
/*   wxMessageDialog *dial = new wxMessageDialog(NULL, descr, title, wxOK);
   dial->ShowModal();*/
   wxMessageDialog dial(NULL, descr, title, wxOK);
   dial.ShowModal();
}

void AlertPopup(string title, string descr)
{
/*   wxMessageDialog *dial = new wxMessageDialog(NULL, descr, title, wxOK | wxICON_EXCLAMATION);
   dial->ShowModal();*/
   wxMessageDialog dial(NULL, descr, title, wxOK | wxICON_EXCLAMATION);
   dial.ShowModal();
}

void YesNoPopup(string title, string descr)
{
   wxMessageDialog dial(NULL, descr, title, wxYES_NO);
   dial.ShowModal();
}

// Popup for number settings
NEDialog::NEDialog(const wxString &title, const wxSize &size, string text, string label, double value, const int *numID) : wxDialog(NULL, wxID_ANY, title, wxDefaultPosition, size)
{
   nevalue = new double;
   *nevalue = value;

   panel = new wxPanel(this, -1);

   vbox = new wxBoxSizer(wxVERTICAL);

   widgetText = new wxStaticText(panel, -1, text);
   vbox->Add(widgetText, 1, wxALL, 5);

   hbox = new wxBoxSizer(wxHORIZONTAL);
   widgetLabel = new wxStaticText(panel, -1, label);
   widgetNE = new wxSpinCtrlDouble(panel, *numID);
   widgetNE->SetRange(-9999999999., 9999999999.);
   widgetNE->SetValue(*nevalue);

   hbox->Add(widgetLabel, 0, wxLEFT | wxTOP, 5);
   hbox->Add(widgetNE, 0, wxRIGHT | wxTOP, 5);

   Connect(*numID, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(NEDialog::UpdateNEValue));

   vbox->Add(hbox, 1, wxEXPAND);

   hbox2 = new wxBoxSizer(wxHORIZONTAL);
   okButton = new wxButton(panel, wxID_OK, wxT("OK"));
   cancelButton = new wxButton(panel, wxID_CANCEL, wxT("Cancel"));
   hbox2->Add(okButton, 0, wxLEFT | wxRIGHT | wxBOTTOM, 5);
   hbox2->Add(cancelButton, 0, wxLEFT | wxRIGHT | wxBOTTOM, 5);

   vbox->Add(hbox2, 1, wxEXPAND);

   panel->SetSizer(vbox);

   Centre();
}

NEDialog::~NEDialog()
{
   delete nevalue;
}

void NEDialog::SetNEntryFormat(wxSpinCtrlDouble *sc, int nrdig, double incr, int limit, double minlim, double maxlim)
{
   sc->SetDigits((unsigned int)nrdig);
   sc->SetIncrement(incr);

   if(limit == 0)
      sc->SetRange(-9999999999., 9999999999.);
   else if(limit == 1)
      sc->SetRange(-9999999999., maxlim);
   else if(limit == 2)
      sc->SetRange(minlim, maxlim);
   else if(limit == -1)
      sc->SetRange(minlim, 99999999999.);
}

void NEDialog::UpdateNEValue(wxSpinDoubleEvent &event)
{
   *nevalue = widgetNE->GetValue();
   if(DBGSIG > 1)
      cout << "# UpdateNEValue         #: " << "Value updated to " << *nevalue << endl;
}

double NEDialog::GetNEValue()
{
   return *nevalue;
}

// Popup for text/name/title settings
FSDialog::FSDialog(const wxString &title, const wxSize &size, std::string text, std::string *label, int nrval, std::string *values, const int *textID) : wxDialog(NULL, wxID_ANY, title, wxDefaultPosition, size)
{
   svalue = new std::string[nrval];
   nrs = new int;

   *nrs = nrval;

   for(int i = 0; i < nrval; i++)
      svalue[i] = values[i];

   panel = new wxPanel(this, -1);

   vbox = new wxBoxSizer(wxVERTICAL);

   widgetText = new wxStaticText(panel, -1, text);
   vbox->Add(widgetText, 0, wxLEFT | wxRIGHT | wxTOP, 5);

   for(int i = 0; i < nrval; i++)
   {
      hbox = new wxBoxSizer(wxHORIZONTAL);
      widgetLabel = new wxStaticText(panel, -1, label[i]);
      widgetTE[i] = new wxTextCtrl(panel, (*textID)+i, svalue[i]);

      hbox->Add(widgetLabel, 0, wxLEFT, 5);
      hbox->Add(widgetTE[i], 1, wxLEFT | wxRIGHT, 5);

      Connect((*textID)+i, wxEVT_TEXT, wxCommandEventHandler(FSDialog::UpdateFSValue));

      vbox->Add(hbox, 0, wxEXPAND | wxLEFT | wxRIGHT | wxTOP, 5);
   }

   hbox2 = new wxBoxSizer(wxHORIZONTAL);
   okButton = new wxButton(panel, wxID_OK, wxT("OK"));

   bool *btest = new bool;
   *btest = CheckFSValues();
   if(*btest)
      okButton->Enable();
   else
      okButton->Disable();
   delete btest;

   hbox2->Add(okButton, 0, wxALL, 5);

   vbox->Add(hbox2, 0, wxEXPAND);

   panel->SetSizer(vbox);

   Centre();
}

FSDialog::~FSDialog()
{
   delete[] svalue;
   delete nrs;
}

std::string FSDialog::GetFSValue(int type)
{
   return svalue[type];
}

void FSDialog::UpdateFSValue(wxCommandEvent &event)
{
   for(int i = 0; i < *nrs; i++)
      svalue[i] = widgetTE[i]->GetLineText(0);

   bool *btest = new bool;
   *btest = CheckFSValues();

   if(*btest)
      okButton->Enable();
   else
      okButton->Disable();

   if(DBGSIG > 1)
      cout << "# UpdateFSValue         #: " << "Updating the FS values" << endl;

   delete btest;
}

bool FSDialog::CheckFSValues()
{
   for(int i = 0; i < *nrs; i++)
   {
      for(int j = 0; j < *nrs; j++)
      {
         if( (i != j) && (GetFSValue(i) == GetFSValue(j)) )
            return false;
      }
   }

   return true;
}
