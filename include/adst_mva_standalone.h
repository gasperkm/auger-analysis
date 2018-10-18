#ifndef _ADST_MVA_H_
#define _ADST_MVA_H_

#define _STANDALONE_ 1
#include "observables.h"
#include "workstation.h"
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
//   int besteye;
   bool goodrec;
   double risemean, risemin, risemax;
   double shfootmean, shfootmin, shfootmax, shfootlimit;
   double aopmean, aopmin, aopmax;
   /* Rewrite code gives information, if any of the observables have not been correctly rewritten (missing SD or FD data):
     - Written as a binary number
     - First five bits (from right) determine if Sim, SD and FD reconstructions are done, if event is hybrid and if VEM signal exists
     - The following bits determine if calculation for each observable is correct
   */
   int rewritecode;

   // Settings and variables for risetime calculations
   TFormula *fRTWeights;
   double limitTankDistance[2];
   double minSignal;
   bool includeSaturated;
   int minPoints;
   double evalDistance;

   // Functions to get all active SD stations and FD eyes in an event
/*   int GetActiveStations();
   int GetActiveEyes();*/

   // Functions for setting write out trees for rewriting ADST files for MVA input
   void PrepareOtherTrees(unsigned int nrfiles, int innr, vector<string> obs);
   int RewriteObservables(int nrfiles, int innr, Observables **sig, Observables **all);

   // Functions for rewriting and calculating observables
   int SetSimObservables(Observables **cursig);
   int SetSdObservables(Observables **cursig);
   int SetFdObservables(Observables **cursig);
//   int GetBestEye();
   double GetXmax(int eye, int type);
   double GetX0(int eye, int type);
   double GetLambda(int eye, int type);
   double GetFdEnergy(int eye, int type);
   double GetFdZenith(int eye, int type);
   double GetFdAzimuth(int eye, int type);
   double GetFdLatitude(int eye, int type);
   double GetFdLongitude(int eye, int type);
   void CalculateShowerFoot(int eye);
   double GetShowerFoot(int eye, int type);

   double GetShowerSize(int eye, int type);
   double GetSdEnergy(int eye, int type);
   double GetBeta(int eye, int type);
   double GetCurvature(int eye, int type);
   void InitRisetimeVariables();
   void RisetimeFunction(double zenith, double energy, TF1 *risetimeFit);
   void CalculateRisetime();
   double GetRisetime(int eye, int type, bool recalc);
   double GetSdZenith(int eye, int type);
   double GetSdAzimuth(int eye, int type);
   double GetSdLatitude(int eye, int type);
   double GetSdLongitude(int eye, int type);
   void CalculateAoP();
   double GetAoP(int eye, int type);

   double GetNrMuons(int eye, int type);
};

#endif
