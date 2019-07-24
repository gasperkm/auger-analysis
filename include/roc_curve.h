#ifndef _ROCCURVE_H_
#define _ROCCURVE_H_

#include "workstation.h"
#include "root_include.h"
#include "root_style.h"

class RocCurve
{
private:
   double ret;
   double *dtemp;
   TH1 *fsigH, *fbgH;
   TH1D *fcumsigH, *fcumbgH;
   TSpline *fsplcumsigH, *fsplcumbgH;
   TH1D *rocH;
   bool *invert;
   double *fNbins;
   double *fXrange;
   double *p1, *p2;
   int *sbin;
//   double *rootVal;

   double Root(double refVal, int fMaxIter, double fAbsTol);
   void GetCumulative(TH1 *inhist, TH1 *outhist);
   double GetEffForRoot(double cut);

/*   void Root(double refVal);
   void GetCumulative(TH1 *inhist, TH1 *outhist);
   double InterpolationX(double yref, double *point1, double *point2);
   double InterpolationY(double xref, double *point1, double *point2);*/
public:
   RocCurve(TH1 *sigH, TH1 *bgH);
   virtual ~RocCurve();

   TH1D* GetROC();
   int GetOutBins();
   bool IsGood();
};

#endif
