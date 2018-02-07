#include "substructure/subheader.h"

Subheader::Subheader()
{
   padding = 5;
   isbigger = new bool;
   *isbigger = false;
}

Subheader::~Subheader()
{
   delete isbigger;
}

int Subheader::GetLabelWidth(wxStaticText *st)
{
//   return (st->GetTextExtent(st->GetLabelText())).GetWidth();
   return (st->GetSize()).GetWidth();
}

int Subheader::GetLabelHeight(wxStaticText *st)
{
//   return (st->GetTextExtent(st->GetLabelText())).GetHeight();
   return (st->GetSize()).GetHeight();
}

int Subheader::GetTextEntryWidth(wxTextCtrl *tc)
{
//   return (tc->GetTextExtent(tc->GetLineText(0))).GetWidth();
   return (tc->GetSize()).GetWidth();
}

int Subheader::GetTextEntryHeight(wxTextCtrl *tc)
{
//   return (tc->GetTextExtent(tc->GetLineText(0))).GetHeight();
   return (tc->GetSize()).GetHeight();
}

int Subheader::GetNumEntryWidth(wxSpinCtrlDouble *sc)
{
   return (sc->GetSize()).GetWidth();
}

int Subheader::GetNumEntryHeight(wxSpinCtrlDouble *sc)
{
   return (sc->GetSize()).GetHeight();
}

int Subheader::GetButtonWidth(wxButton *bt)
{
//   return (bt->GetTextExtent(bt->GetLabelText())).GetWidth();
   return (bt->GetSize()).GetWidth();
}

int Subheader::GetButtonHeight(wxButton *bt)
{
//   return (bt->GetTextExtent(bt->GetLabelText())).GetHeight();
   return (bt->GetSize()).GetHeight();
}

int Subheader::GetDropBoxWidth(wxChoice *db)
{
//   return (bt->GetTextExtent(bt->GetLabelText())).GetWidth();
   return (db->GetSize()).GetWidth();
}

int Subheader::GetDropBoxHeight(wxChoice *db)
{
//   return (bt->GetTextExtent(bt->GetLabelText())).GetHeight();
   return (db->GetSize()).GetHeight();
}

int Subheader::GetCheckBoxWidth(wxCheckBox *cb)
{
//   return (bt->GetTextExtent(bt->GetLabelText())).GetWidth();
   return (cb->GetSize()).GetWidth();
}

int Subheader::GetCheckBoxHeight(wxCheckBox *cb)
{
//   return (bt->GetTextExtent(bt->GetLabelText())).GetHeight();
   return (cb->GetSize()).GetHeight();
}

void Subheader::SetLabelSize(wxStaticText *st, int width, int height)
{
   if(width == -1)
      st->SetSize(GetLabelWidth(st), height);
   else if(height == -1)
      st->SetSize(width, GetLabelHeight(st));
   else if( (width != -1) && (height != -1) )
      st->SetSize(width, height);
}

void Subheader::SetNEntryFormat(wxSpinCtrlDouble *sc, int nrdig, double incr, int limit, double minlim, double maxlim)
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

double Subheader::GetNumber(wxSpinCtrlDouble *sc)
{
   return sc->GetValue();
}

void Subheader::SetNumber(wxSpinCtrlDouble *sc, double value)
{
   sc->SetValue(value);
}
