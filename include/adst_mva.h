#ifndef _ADST_MVA_H_
#define _ADST_MVA_H_

#include "observables.h"
#include "workstation.h"
#if _STANDALONE_ == 0
   #include <wx/progdlg.h>
#endif
#if OFFVER == 0
   #include "OfflineIncludeOld.h"
#elif OFFVER == 1
   #include "OfflineIncludeNew.h"
#endif

class AdstMva
{
public:
   AdstMva();
   virtual ~AdstMva();

   // Variables for files and their structure
   string outname;
   string inname;
   TFile *outfile;
   TTree *all_tree, *sig_tree;

   // Objects for ADST
   RecEvent *fRecEvent;
   RecEventFile *fFile;
   DetectorGeometry *fDetGeo;
   SdRecShower *sdrecshw;
   GenShower *genshw;

   // Variables for SD stations and FD eyes
   vector<SdRecStation> actstations;
   int nreyes;
   vector<FDEvent> acteyes;
   bool isheco;
//   int besteye;
   bool goodrec, goodrecsim, goodrecsd, goodrecfd;
   double risemean, risemin, risemax, riseerr;
   double shfootmean, shfootmin, shfootmax, shfootlimit;
   double aopmean, aopmin, aopmax;
   /* Rewrite code gives information, if any of the observables have not been correctly rewritten (missing SD or FD data):
     - Written as a binary number
     - First five bits (from right) determine if Sim, SD and FD reconstructions are done, if event is hybrid and if VEM signal exists
     - The following bits determine if calculation for each observable is correct
   */
   int rewritecode;

   // In addition to observables, station specific vectors are also given
   int nractstations;
   vector<float> stationDistance[3];
   vector<float> stationRisetime[3];
   vector<bool> stationHSat;

   // Settings and variables for risetime calculations
   TFormula *fRTWeights;
   double limitTankDistance[2];
   double minSignal;
   bool includeSaturated;
   int minPoints;
   double evalDistance;

   // Settings and variables for combining stereo events
   double *wQuantity;
   double *quantitySum;
   double *wQuantitySum;
   double *quantitySumErr;

   // Functions to get all active SD stations and FD eyes in an event
/*   int GetActiveStations();
   int GetActiveEyes();*/

   // Functions for setting write out trees for rewriting ADST files for MVA input
   void PrepareOtherTrees(unsigned int nrfiles, int innr, vector<string> obs);
#if _STANDALONE_ == 0
   int RewriteObservables(int nrfiles, int innr, Observables **sig, Observables **all, wxProgressDialog *progress);
#elif _STANDALONE_ == 1
   int RewriteObservables(int nrfiles, int innr, Observables **sig, Observables **all);
#endif

   // Functions for rewriting and calculating observables
   int SetSimObservables(Observables **cursig);
   void ZeroSimObservables(Observables **cursig);
   int SetSdObservables(Observables **cursig);
   void ZeroSdObservables(Observables **cursig);
   int SetFdObservables(Observables **cursig);
   void ZeroFdObservables(Observables **cursig);
   int SetStationValues();
   float TotalSignalFromPMT(SdRecStation *station, bool silent);
//   int GetBestEye();
   float GetXmax(int eye, int type);
   float GetX0(int eye, int type);
   float GetLambda(int eye, int type);
   float GetFdEnergy(int eye, int type);
   float GetFdZenith(int eye, int type);
   float GetFdAzimuth(int eye, int type);
   float GetFdLatitude(int eye, int type);
   float GetFdLongitude(int eye, int type);
   void CalculateShowerFoot(int eye);
   float GetShowerFoot(int eye, int type);

   float GetNrStations(int type);
   float GetShowerSize(int type);
   float GetSdEnergy(int type);
   float GetBeta(int type);
   float GetCurvature(int type);
   void InitRisetimeVariables();
   void RisetimeFunction(double zenith, double energy, TF1 *risetimeFit);
   void CalculateRisetime();
   float GetRisetime(int type, bool recalc);
   float GetSdZenith(int type);
   float GetSdAzimuth(int type);
   float GetSdLatitude(int type);
   float GetSdLongitude(int type);
   void CalculateAoP();
   float GetAoP(int type);

   float GetNrMuons(int type);
};

#endif
