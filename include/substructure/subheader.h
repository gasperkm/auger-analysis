#ifndef _SUBHEADER_H_
#define _SUBHEADER_H_

#include <wx/wx.h>
#include <wx/spinctrl.h>
#include <wx/editlbox.h>
#include <iostream>
#include <string>
#include <vector>
#include <limits>

using namespace std;

class Subheader
{
public:
   Subheader();
   virtual ~Subheader();

   int GetLabelWidth(wxStaticText *st);
   int GetLabelHeight(wxStaticText *st);
   int GetTextEntryWidth(wxTextCtrl *tc);
   int GetTextEntryHeight(wxTextCtrl *tc);
   int GetNumEntryWidth(wxSpinCtrlDouble *sc);
   int GetNumEntryHeight(wxSpinCtrlDouble *sc);
   int GetButtonWidth(wxButton *bt);
   int GetButtonHeight(wxButton *bt);
   int GetDropBoxWidth(wxChoice *db);
   int GetDropBoxHeight(wxChoice *db);
   int GetCheckBoxWidth(wxCheckBox *cb);
   int GetCheckBoxHeight(wxCheckBox *cb);

   void SetLabelSize(wxStaticText *st, int width, int height);

   void SetNEntryFormat(wxSpinCtrlDouble *sc, int nrdig, double incr, int limit, double minlim, double maxlim);

   double GetNumber(wxSpinCtrlDouble *sc);
   void SetNumber(wxSpinCtrlDouble *sc, double value);

   int padding;
   bool *isbigger;
};

#endif
