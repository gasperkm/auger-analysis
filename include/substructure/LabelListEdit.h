#include "subheader.h"

class LabelListEdit : public Subheader
{
public:
    LabelListEdit(wxPanel *parent, wxString label, int width, int height, int multi, const int listID, vector<int> *butID);

    wxStaticText *widgetLabel;
//    wxEditableListBox *widgetLB;
    wxListBox *widgetLB;
    wxButton *widgetTB[4];
    wxBitmap *widgetBM;

    wxBoxSizer *subsizer;
};
