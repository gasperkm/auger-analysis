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

NEDialog::NEDialog(const wxString &title, const wxSize &size, string text, string label, double value, const int *numID) : wxDialog(NULL, wxID_ANY, title, wxDefaultPosition, size)
{
   nevalue = value;

   wxPanel *panel = new wxPanel(this, -1);

   wxBoxSizer *vbox = new wxBoxSizer(wxVERTICAL);

   wxStaticText *widgetText = new wxStaticText(panel, -1, text);
   vbox->Add(widgetText, 1, wxALL, 5);

   wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
   wxStaticText *widgetLabel = new wxStaticText(panel, -1, label);
   widgetNE = new wxSpinCtrlDouble(panel, *numID);
   widgetNE->SetValue(nevalue);

   hbox->Add(widgetLabel, 0, wxLEFT | wxTOP, 5);
   hbox->Add(widgetNE, 0, wxRIGHT | wxTOP, 5);

   Connect(*numID, wxEVT_SPINCTRLDOUBLE, wxSpinDoubleEventHandler(NEDialog::UpdateNEValue));

   vbox->Add(hbox, 1, wxEXPAND);

   wxBoxSizer *hbox2 = new wxBoxSizer(wxHORIZONTAL);
   wxButton *okButton = new wxButton(panel, wxID_OK, wxT("OK"));
   wxButton *cancelButton = new wxButton(panel, wxID_CANCEL, wxT("Cancel"));
   hbox2->Add(okButton, 0, wxLEFT | wxRIGHT | wxBOTTOM, 5);
   hbox2->Add(cancelButton, 0, wxLEFT | wxRIGHT | wxBOTTOM, 5);

   vbox->Add(hbox2, 1, wxEXPAND);

   panel->SetSizer(vbox);

   Centre();
}

NEDialog::~NEDialog()
{
//   delete widgetNE;
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
   nevalue = widgetNE->GetValue();
   if(DBGSIG > 1)
      cout << "# UpdateNEValue         #: " << "Value updated to " << nevalue << endl;
}

double NEDialog::GetNEValue()
{
   return nevalue;
}
