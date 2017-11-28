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

class CheckNEntry : public Subheader
{
public:
    // Single NEntry constructor
    CheckNEntry(wxPanel *parent, bool checked, wxString label, const int checkID, double numval, const int nentryID, int maxsize);
    // Multiple NEntry constructors
    CheckNEntry(wxPanel *parent, bool checked, wxString label, const int checkID, vector<double> numval, vector<int> nentryID, int maxsize);

    wxCheckBox *widgetChBox;
    wxSpinCtrlDouble *widgetNE[2];

    wxBoxSizer *subsizer;
};
