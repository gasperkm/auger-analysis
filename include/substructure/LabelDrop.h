#include "subheader.h"

class LabelDrop : public Subheader
{
public:
    LabelDrop(wxPanel *parent, wxString label, vector<string> entrytext, string selecttext, const int dropID, int maxsize);

    wxStaticText *widgetLabel;
    wxChoice *widgetCB;

    wxBoxSizer *subsizer;
};
