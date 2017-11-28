#include "subheader.h"

class LabelTEntryButton : public Subheader
{
public:
    // Single button constructor
    LabelTEntryButton(wxPanel *parent, wxString label, wxString deftext, const int textID, wxString buttext, const int butID, int maxsize);
    // Multiple button constructor
    LabelTEntryButton(wxPanel *parent, wxString label, wxString deftext, const int textID, vector<string> buttext, vector<int> butID, int maxsize);

    wxStaticText *widgetLabel;
    wxTextCtrl *widgetTE;
    wxButton *widgetTB[6];

    wxBoxSizer *subsizer;
};
