#include "subheader.h"

/*
 *  Format for number entry fields:
 *  [nr of digits] [increment] [limits] [min limit] [max limit]
 *
 *  [nr of digits]: number of digits -> SetDigits()
 *  [increment]: 0 = integer, 1 = real
 *  [limits]: 0 = no limits, 1 = maximum limit, 2 = min+max limits, -1 = minimum limit
 *  [min limit], [max limit]: minimum and maximum limit values
 */

class LabelNEntryDropButton : public Subheader
{
public:
    // Single NEntry and single button constructor
    LabelNEntryDropButton(wxPanel *parent, wxString label, double numval, const int nentryID, vector<string> entrytext, string selecttext, const int dropID, wxString buttext, const int butID, int maxsize);
    // TODO: add other constructors
/*    // Multiple NEntry constructors
    CheckNEntry(wxPanel *parent, bool checked, wxString label, int nrnentries, vector<double> numval, int maxsize);*/

    wxStaticText *widgetLabel;
    wxSpinCtrlDouble *widgetNE[2];
    wxComboBox *widgetCB;
    wxButton *widgetTB[6];

    wxBoxSizer *subsizer;
};

