#ifndef _FRAME_H_
#define _FRAME_H_

#define _STANDALONE_ 0

#include <wx/wx.h>
#include <wx/notebook.h>
#include <wx/progdlg.h>
#include "workstation.h"
#include "substr.h"
#include "root_include.h"
#include "adst_mva.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class MyFrame : public wxFrame
{
private:
    int ret;
public:
    MyFrame(const wxString& title);

    void SetLayout(wxCommandEvent& event);
    void SaveLayout(wxCommandEvent& event);
    void ReadLayout(int *width, int *height);
    void OnAbout(wxCommandEvent& event);
//    void OnClose(wxCloseEvent& event);
    void OnQuit(wxCommandEvent& event);
    void ShowProgress(wxString title, wxString message, int maxval);

    void InitObservables();
    void InitMethods();
    void InitVariables();
    void MakeTitle(wxPanel *parent, wxBoxSizer *box, wxString title);
    void SetTooltips();
    string GetMethodName(string name);

    wxMenuBar *menubar;
    wxMenu *file;
    wxMenu *help;
    wxProgressDialog *progress;

    // Window variables
    int WW;
    int WH;

    // Variables
    vector<string> observables;
    int nrobs;
    vector<bool> obssel;
    vector<bool> obsorigsel;
    int nrselobs;
    vector<string> methods;
    vector<string> methodsOpt;
    vector<string> methodsDesc;
    int nrmethods;
    int oldselect[3];
    vector<double> ecutBins;
    vector<double> zcutBins;
    int selcuttype;
    int seleyetype;
    vector<int> sigSeleye;
    vector<int> backSeleye;

    int nrlists;
    LabelListEdit *allLBE[20];

    int nrkeys;
    TList *keyslist;

    TMatrixD *sigCorMat;
    TMatrixD *backCorMat;
    TMatrixD *otherCorMat[10];
    TMatrixD *covMatNeg;
    TMatrixD *covMatPos;
    TMatrixDEigen *eigenCovMat;
    TMatrixD *eigenValMatNeg;
    TMatrixD *eigenValMatPos;

    string applymva;
    Observables *generalObservables;
    float mvalimit[2];

    vector<double> statsMin;
    vector<double> statsMax;

    vector<int> cutResults;
    string tempAnalysisFile;
    bool freshAnalysis;

    string mvaresults;
    string mvaprintout;

    int mixednum;

    wxArrayInt selections;

    vector<double> esplitBins;
    bool customBinning;
    double oldBinSettings[3];

    // File open variables
    string *currentMvaDir;
    string *currentRewriteDir;
    string *currentAnalysisDir;
    string *currentCutsDir;
    string *currentCutsInputDir;

    // Control substructures
    LabelButton *selectMvaFile;
    LabelListEdit *mvaList[3];
//    LabelListEdit *krnekiList;
    LabelButton *startRewriting;
    LabelButton *startCombining;

    LabelNEntryButton *startSplitting;
    LabelDrop *splitCutObservables;
    CheckNEntry *splitCutEnergy;
    CheckNEntry *splitCutZenith;
    CheckNEntry *splitCutRisetime;
    LabelDrop *splitEyeSelection;

    LabelTEntry *selectedMva;
    LabelList *selectObservables;
    LabelDrop *signalSelect;
    LabelDrop *backgroundSelect;
    LabelDrop *methodsSelect;
    LabelNEntry *cutMva;
    LabelDrop *cutObservables;
    CheckNEntry *cutEnergy;
    LabelNEntryDropButton *cutEnergyBins;
    LabelTEntryButton *cutEnergyBinsCustom;
    CheckNEntry *cutZenith;
    CheckNEntry *cutRisetime;
    LabelNEntryDropButton *cutZenithBins;
//    LabelDrop *eyeSelection;
    LabelDrop *dataSelect;
    CheckList *specialMva;
    LabelButton *startMva;

    // File operation functions
    void SelectMvaFile(wxCommandEvent& event);
    int CheckFormat(string *infile);
    void SelectRewrite(wxCommandEvent& event);
    void PrepareFileSplit(wxCommandEvent& event);
    void SelectCombine(wxCommandEvent& event);
    void SelectMerge(wxCommandEvent& event);

    // Connect functions
    void EnableMvaFile(wxCommandEvent& event);
    void UpdateEnergyBinSelect(wxSpinDoubleEvent& event);
    void RunEnergyBinSelect();
    void RunEnergyBinSplitSelect();
    void CheckEnergyBin(wxCommandEvent& event);
    void CustomEnergyBins(wxCommandEvent& event);
    void CustomEnergyBinsDisable(wxCommandEvent& event);
    void UpdateZenithBinSelect(wxSpinDoubleEvent& event);
    void RunZenithBinSelect();
    void CheckZenithBin(wxCommandEvent& event);
    void CheckBothBins(wxCommandEvent& event);
    void UpdateObservableSelection(wxCommandEvent& event);
    void StartMvaAnalysis(wxCommandEvent& event);
    void CreateTempEventFile(wxCommandEvent& event);
    string SelectMva();
    void ApplyMvaCut(wxCommandEvent& event);
    void SetDefaultMva(wxCommandEvent& event);
    int StartRewrite(string *outfile);
    int StartFileSplit(string infile, int cursel, int nrsel);
    int StartCombine(string *outfile);
    int StartMerge(string *outfile);
    void EditList(wxCommandEvent& event);

    // Functions for the MVA analysis
//    int MvaNoteObservables(TMVA::Factory *factory);
    int MvaNoteObservables(int count);
    int MvaTreeFile(string *infilename, string *outfilename, int *nrEvents);
    int MvaSetTrees(int type, TFile *ifile, TTree *outtree);
//    int IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, vector<int> *seleye, bool split, int splitbin);
    int IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, bool split, int splitbin);
    int PerformMvaAnalysis(string *infilename, string *outfilename, int *curcount);
    void SetTmvaType(TMVA::Factory *factory, int nr, string *formula);
    int BookTheMethod(TMVA::Factory *factory);
    int GetTrainingShift(string *mvafilename);
    void GetCorrelations(TFile *corfile);
    void GetApplyCorrelations(string *corname);

    // Functions for the application of MVA analysis
    void MvaApplication(string *infilename, bool application, int mean);
    void CreateMVAPlots(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean);
    void CreateOutput(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean);
    void GetErrors(TTree *app, float *obsvars, vector<string> *obs, int curtree);
    void GetMvaError(int selection, double *outvalue, string *inname);
};

// Menu IDs
const int ID_SETLAYOUT 		= 101;
const int ID_SAVELAYOUT 	= 102;
const int ID_OPENHELP 		= 103;

// MVA analysis left panel IDs
const int ID_OPENFILE 		= 201;
const int ID_REWRITE 		= 202;
const int ID_SPLIT 		= 203;
const int ID_COMBINE 		= 204;
const int ID_MERGE 		= 205;
const int ID_SELECTMVAFILE 	= 206;

// MVA analysis right panel IDs
const int ID_CHANGEOBSSELECT 	= 300;
const int ID_ENERGYLIMITMIN 	= 301;
const int ID_ENERGYLIMITMAX 	= 302;
const int ID_ENERGYBIN 		= 303;
const int ID_ENERGYBINCUSTOM	= 304;
const int ID_CUSTOMDISABLE	= 305;
const int ID_CHECKENERGYBINS 	= 306;
const int ID_ZENITHLIMITMIN 	= 307;
const int ID_ZENITHLIMITMAX 	= 308;
const int ID_ZENITHBIN 		= 309;
const int ID_CHECKZENITHBINS 	= 310;
const int ID_STARTMVA 		= 311;
const int ID_TEMPFILE 		= 312;
const int ID_MVACUT 		= 313;
const int ID_CHECKBINS 		= 314;
const int ID_DEFOPTIONS 	= 315;

// Additional IDs for any custom dialogs
const int ID_MVACUTDIALOG	= 401;
const int ID_RANDSEEDDIALOG	= 402;

// Additional IDs for custom options regarding listboxes
const int ID_DELETELIST		= 1001;
const int ID_UPLIST		= 1002;
const int ID_DOWNLIST		= 1003;
const int ID_CLEARLIST		= 1004;

#endif
