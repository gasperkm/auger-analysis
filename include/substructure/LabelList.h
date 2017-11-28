#include "subheader.h"

class LabelList : public Subheader
{
public:
    LabelList(wxPanel *parent, wxString label, int height, int multi, const int listID);

    wxStaticText *widgetLabel;
    wxListBox *widgetLB;

    wxBoxSizer *subsizer;
};
