#include "subheader.h"

class CheckList : public Subheader
{
public:
    // Single Checkbox
    CheckList(wxPanel *parent, bool checked, wxString label, const int checkID);
    // Multiple Checkboxes
    CheckList(wxPanel *parent, vector<int> checked, vector<string> label, vector<int> checkID, int maxsize);

    wxCheckBox *widgetChBox[18];

    wxBoxSizer *subsizer;
};
