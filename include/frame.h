#ifndef _FRAME_H_
#define _FRAME_H_

#define _STANDALONE_ 0

#include <wx/wx.h>
#include <wx/notebook.h>
#include <wx/progdlg.h>
#include "workstation.h"
#include "substr.h"
#include "root_include.h"
#include "root_style.h"
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
    wxArrayInt plotSelections;

    vector<double> esplitBins;
    bool customBinning;
    double oldBinSettings[3];
    bool multipleEnergyBins;

    bool setdeltas[2];

    bool bCustomSplit;

    Long64_t tentry;
    vector<float> *stationDistance[3];
    vector<float> *stationRisetime[3];
    vector<bool> *stationHSat;
    TBranch *bdist;
    TBranch *bdistneg;
    TBranch *bdistpos;
    TBranch *brise;
    TBranch *briseneg;
    TBranch *brisepos;
    TBranch *bsat;

    // File open variables
    string *currentMvaDir;
    string *currentRewriteDir;
    string *currentAnalysisDir;
    string *currentCutsDir;
    string *currentCutsInputDir;
    string *currentPlotDir;

    // Control substructures
    LabelButton *selectMvaFile;
    LabelListEdit *mvaList[3];
//    LabelListEdit *krnekiList;
    LabelButton *startRewriting;
    LabelButton *startCombining;

    LabelNEntryButton *startSplitting;
    LabelButton *startSplittingFile;
    LabelDrop *splitCutObservables;
    CheckNEntry *splitCutEnergy;
    CheckNEntry *splitCutZenith;
    CheckNEntry *splitCutRisetime;
    CheckList *splitEnforceSD;
    LabelDrop *splitEyeSelection;

    LabelTEntry *selectedMva;
    LabelList *selectObservables;
    LabelDropButton *uncertSelect;
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

    // Structures for plotting
    LabelButton *selectPlotFile;
    LabelButton *selectPlotDir;
    LabelListEdit *openFileList;
    LabelDrop *methodsPlotSelect;
    LabelList *plottingList[2];
    LabelNEntry *plotNrBins;
    CheckNEntry *plotXaxisRange;
    CheckNEntry *plotYaxisRange;
    LabelDrop *hiSelect;
    LabelDrop *prodSelect;
    LabelDrop *dataPlotSelect;
    CheckList *specialPlot;
    LabelButton *startPlot;
    LabelDrop *dropHistTree[3];
/*    LabelDrop *secondHistTree;
    LabelDrop *thirdHistTree;*/
    CheckList *specialHistPlot;
    LabelButton *startHistPlot;

    // File operation functions
    void SelectMvaFile(wxCommandEvent& event);
    int CheckFormat(string *infile);
    void SelectRewrite(wxCommandEvent& event);
    void PrepareFileSplit(wxCommandEvent& event);
    void PrepareFileSplitCustom(wxCommandEvent& event);
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
    void SetNegativeUncertainty(wxCommandEvent& event);
    void SetPositiveUncertainty(wxCommandEvent& event);
    void DisableUncertainty(wxCommandEvent& event);
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
    int MvaTreeFile(string *infilename, string *outfilename, int *nrEvents, int curenbin);
    int MvaSetTrees(int type, TFile *ifile, TTree *outtree);
    int SetDeltas(int s38rise, int type, TFile *ifile, bool isdata);
    void CalculateS38(vector<double> *shwsize, vector<double> *zenith, float *fitpar, float *fitparErr, vector<double> *outVect);
    void WriteoutS38Fits(int tree, int type, float minEn, float maxEn, int nrpar, float *fitpar, float *fitparErr);
    void WriteoutDeltaFits(int tree, float minEn, float maxEn, float minZen, float maxZen, int nrpar, float *fitpar, float *fitparErr);
    void PrintS1000Fit(int *ebinS1000, TGraphAsymmErrors *fitgraph, TF1 *fitfunc, float *fitpar, float *fitparErr, RootStyle *mystyle, int seltype);
    void PrintS38Fit(TGraphAsymmErrors *fitgraph, TF1 *fitfunc, float *fitpar, float *fitparErr, RootStyle *mystyle, int seltype);
    void PrintRisetimeFit(int *zbinRise, TGraphAsymmErrors *fitgraphHG, TGraphAsymmErrors *fitgraph, TF1 *fitfuncHG, TF1 *fitfunc, float *fitpar, float *fitparErr, RootStyle *mystyle);
//    int IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, vector<int> *seleye, bool split, int splitbin);
    int IsInsideCuts(Observables *mean, Observables *neg, Observables *pos, bool split, int splitbin);
    int PerformMvaAnalysis(string *infilename, string *outfilename, int *curcount);
#if ROOTVER == 5
    void SetTmvaType(TMVA::Factory *factory, int nr, string *formula);
    int BookTheMethod(TMVA::Factory *factory);
#elif ROOTVER == 6
    void SetTmvaType(TMVA::Factory *factory, TMVA::DataLoader *dataloader, int nr, string *formula);
    int BookTheMethod(TMVA::Factory *factory, TMVA::DataLoader *dataloader);
#endif
    int GetTrainingShift(string *mvafilename);
    void GetCorrelations(TFile *corfile);
    void GetApplyCorrelations(string *corname);

    // Functions for the application of MVA analysis
    void MvaApplication(string *infilename, bool application, int mean);
    void CreateMVAPlots(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean);
    void CreateOutput(TTree *app, TMVA::Reader *reader, string mvamethod, float *obsvars, string signalName, int curtree, bool application, int mean);
    void GetErrors(TTree *app, float *obsvars, vector<string> *obs, int curtree);
    void GetMvaError(int selection, double *outvalue, string *inname);

    // Functions for plotting
    void SelectPlotFile(wxCommandEvent& event);
    void SelectPlotDir(wxCommandEvent& event);
    bool CheckPlotTreeSelect(string *plotFile, wxChoice *combo);
    void SetPlotTreeSelect(string *plotFile, wxChoice *combo, bool wdisable);
    void StartHistogramPlot(wxCommandEvent& event);
    void StartScatterPlot(wxCommandEvent& event);
    void StartROCPlot(wxCommandEvent& event);
    void StartHistogramScatterPlot(int type);
    void StartMvaHistFit(wxCommandEvent& event);
    void SetDefaultPlot(wxCommandEvent& event);
    void SetDefaultMvaPlot(wxCommandEvent& event);
    void SetDefaultHistPlot(wxCommandEvent& event);
};

// Menu IDs
const int ID_SETLAYOUT 		= 101;
const int ID_SAVELAYOUT 	= 102;
const int ID_OPENHELP 		= 103;

// MVA analysis left panel IDs
const int ID_OPENFILE 		= 201;
const int ID_REWRITE 		= 202;
const int ID_SPLIT 		= 203;
const int ID_SPLITFILE 		= 204;
const int ID_COMBINE 		= 205;
const int ID_MERGE 		= 206;
const int ID_SELECTMVAFILE 	= 207;

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
const int ID_NUNCERTSELECT 	= 316;
const int ID_PUNCERTSELECT 	= 317;
const int ID_DISABLEUNCERT 	= 318;

// Plotting panel IDs
const int ID_OPENPLOTFILE	= 401;
const int ID_OPENPLOTDIR	= 402;
const int ID_PLOTHIST		= 403;
const int ID_PLOTSCAT		= 404;
const int ID_PLOTROC		= 405;
const int ID_PLOTMVAFIT		= 406;
const int ID_PLOTDEFOPTIONS	= 407;
const int ID_PLOTMVADEFOPTIONS	= 408;
const int ID_PLOTHISTDEFOPTIONS	= 409;

// Additional IDs for any custom dialogs
const int ID_MVACUTDIALOG	= 901;
const int ID_RANDSEEDDIALOG	= 902;
const int ID_TITLEDIALOG	= 910;
const int ID_LEGENDDIALOG	= 920;

// Additional IDs for custom options regarding listboxes
const int ID_DELETELIST		= 1001;
const int ID_UPLIST		= 1002;
const int ID_DOWNLIST		= 1003;
const int ID_CLEARLIST		= 1004;

#endif
