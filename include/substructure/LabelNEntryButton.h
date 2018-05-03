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

class LabelNEntryButton : public Subheader
{
public:
    // Single NEntry + single button constructor
    LabelNEntryButton(wxPanel *parent, wxString label, double numval, const int nentryID, wxString buttext, const int butID, int maxsize);
    // Multiple NEntry + single button constructor
    LabelNEntryButton(wxPanel *parent, wxString label, vector<double> *numval, vector<int> *nentryID, wxString buttext, const int butID, int maxsize);
    // Single NEntry + multiple button constructor
    LabelNEntryButton(wxPanel *parent, wxString label, double numval, const int nentryID, vector<string> *buttext, vector<int> *butID, int maxsize);
    // Multiple NEntry + multiple button constructor
    LabelNEntryButton(wxPanel *parent, wxString label, vector<double> *numval, vector<int> *nentryID, vector<string> *buttext, vector<int> *butID, int maxsize);

    wxStaticText *widgetLabel;
    wxSpinCtrlDouble *widgetNE[2];
    wxButton *widgetTB[6];

    wxBoxSizer *subsizer;
};
