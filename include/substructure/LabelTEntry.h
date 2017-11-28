#include "subheader.h"

class LabelTEntry : public Subheader
{
public:
    LabelTEntry(wxPanel *parent, wxString label, wxString deftext, const int textID, int maxsize);

    wxStaticText *widgetLabel;
    wxTextCtrl *widgetTE;

    wxBoxSizer *subsizer;
};
