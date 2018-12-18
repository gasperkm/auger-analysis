#include "subheader.h"

class LabelDropButton : public Subheader
{
public:
    // Single button constructor
    LabelDropButton(wxPanel *parent, wxString label, vector<string> *entrytext, string selecttext, const int dropID, wxString buttext, const int butID, int maxsize);
    // Multiple button constructor
    LabelDropButton(wxPanel *parent, wxString label, vector<string> *entrytext, string selecttext, const int dropID, vector<string> *buttext, vector<int> *butID, int maxsize);

    wxStaticText *widgetLabel;
    wxChoice *widgetCB;
    wxButton *widgetTB[6];

    wxBoxSizer *subsizer;
};
