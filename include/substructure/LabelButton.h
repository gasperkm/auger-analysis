#include "subheader.h"

class LabelButton : public Subheader
{
public:
    // Single button constructor
    LabelButton(wxPanel *parent, wxString label, wxString buttext, const int butID, int maxsize);
    // Multiple button constructor
    LabelButton(wxPanel *parent, wxString label, vector<string> *buttext, vector<int> *butID, int maxsize);

    wxStaticText *widgetLabel;
    wxButton *widgetTB[6];

    wxBoxSizer *subsizer;
};
