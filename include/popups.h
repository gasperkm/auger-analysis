#ifndef _POPUPS_H_
#define _POPUPS_H_

#include <wx/wx.h>
#include <wx/spinctrl.h>
#include "workstation.h"
#include <string>

void ErrorPopup(std::string title, std::string descr);
void InfoPopup(std::string title, std::string descr);
void AlertPopup(std::string title, std::string descr);
void YesNoPopup(std::string title, std::string descr);

class NEDialog : public wxDialog
{
private:
   double nevalue;
public:
   NEDialog(const wxString &title, const wxSize &size, std::string text, std::string label, double value, const int *numID);
   virtual ~NEDialog();

   wxSpinCtrlDouble *widgetNE;

   void SetNEntryFormat(wxSpinCtrlDouble *sc, int nrdig, double incr, int limit, double minlim, double maxlim);
   void UpdateNEValue(wxSpinDoubleEvent &event);
   double GetNEValue();
};	

#endif
