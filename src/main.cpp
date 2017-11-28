#include "main.h"
#include "frame.h"

IMPLEMENT_APP(MyApp);

bool MyApp::OnInit()
{
    MyFrame *frame = new MyFrame(wxT("Auger analysis software"));
    frame->Show(true);
    frame->Centre();

    return true;
}
