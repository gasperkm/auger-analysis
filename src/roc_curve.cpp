#include <cstdlib>
#include <iomanip>
#include "roc_curve.h"

using namespace std;

RocCurve::RocCurve(TH1 *sigH, TH1 *bgH)
{
   dtemp = new double[4];
   invert = new bool;
   fNbins = new double[2];
   fXrange = new double[2];

   fsigH = sigH;
   fbgH = bgH;

   if(fsigH->GetMean() >= fbgH->GetMean())
      *invert = false;
   else
      *invert = true;

   cout << "Invert = " << (int)*invert << endl;

   fNbins[0] = fsigH->GetNbinsX();
   fNbins[1] = fNbins[0]*10;

   fXrange[0] = fsigH->GetXaxis()->GetXmin();
   fXrange[1] = fsigH->GetXaxis()->GetXmax();
   cout << "X range = (" << fXrange[0] << ", " << fXrange[1] << ")" << endl;

   // Cumulative distributions
   fcumsigH = new TH1D("fcumsigH", "", fNbins[0], fXrange[0], fXrange[1]);
   fcumbgH = new TH1D("fcumbgH", "", fNbins[0], fXrange[0], fXrange[1]);

   // Return ROC histogram
   rocH = new TH1D("rocH", "", fNbins[1], 0, 1);
}

RocCurve::~RocCurve()
{
   delete[] dtemp;
   delete invert;
   delete[] fNbins;
   delete[] fXrange;
   delete fcumsigH;
   delete fcumbgH;
   delete fsplcumsigH;
   delete fsplcumbgH;
   delete rocH;
}

bool RocCurve::IsGood()
{
   if( (TMath::Abs(fXrange[1] - fbgH->GetXaxis()->GetXmax()) > 0.000001) || (TMath::Abs(fXrange[0] - fbgH->GetXaxis()->GetXmin()) > 0.000001) || (fbgH->GetNbinsX() != fsigH->GetNbinsX()) )
   {
      cout << "Error! Histograms don't have the same number of bins or range." << endl;
      return false;
   }
   else
      return true;
}

void RocCurve::GetCumulative(TH1 *inhist, TH1 *outhist)
{
   cout << "Entering GetCumulative function: Returns the cumulative version of a distribution" << endl;

   dtemp[0] = 0;
   dtemp[1] = 0;
   dtemp[2] = 0;
   dtemp[3] = 0;

   dtemp[3] = inhist->GetEntries();
   cout << "  Total number of entries = " << dtemp[3] << endl;

   for(int i = 1; i <= fNbins[0]; i++)
   {
      dtemp[0] = inhist->GetBinCenter(i);
      dtemp[1] = inhist->GetBinContent(i)/dtemp[3];
      dtemp[2] += dtemp[1];
      cout << "  " << i << "\t" << dtemp[0] << "\t" << dtemp[1] << "\t" << 1. - dtemp[2] << endl;

      if(!(*invert))
         outhist->SetBinContent(i, 1. - dtemp[2]);
      else
         outhist->SetBinContent(i, dtemp[2]);
   }
}

TH1D* RocCurve::GetROC()
{
   cout << "Entering GetROC function: Returns the ROC curve from signal and background." << endl;

   // Return cumulative distributions
   GetCumulative(fsigH, fcumsigH);
   GetCumulative(fbgH, fcumbgH);

   // Transform cumulative distributions to splines
   fsplcumsigH = new TMVA::TSpline1("spline_signal", new TGraph(fcumsigH));
   fsplcumbgH = new TMVA::TSpline1("spline_background", new TGraph(fcumbgH));

   // Get the ROC curve
   for(int i = 1; i <= fNbins[1]; i++)
   {
      cout << endl << "i = " << i << ":" << endl;
      dtemp[0] = rocH->GetBinCenter(i);
      dtemp[1] = Root(dtemp[0], 100, 0.0);
      dtemp[2] = fsplcumbgH->Eval(dtemp[1]);
      cout << "  dtemp[0] = " << dtemp[0] << ", dtemp[1] = " << dtemp[1] << ", dtemp[2] = " << dtemp[2] << endl;

      rocH->SetBinContent(i, dtemp[0], 1. - dtemp[2]);
   }

   return rocH;
}

double RocCurve::GetEffForRoot(double cut)
{
   return fsplcumsigH->Eval(cut);
}

double RocCurve::Root(double refVal, int fMaxIter, double fAbsTol)
{
   // Root finding using Brents algorithm; taken from CERNLIB function RZERO
   double a  = fXrange[0], b = fXrange[1];
   double fa = GetEffForRoot(a) - refVal;
   double fb = GetEffForRoot(b) - refVal;
   if(fb*fa > 0)
   {
      cout << "<ROCCalc::Root> initial interval w/o root: "
            << "(a=" << a << ", b=" << b << "),"
            << " (Eff_a=" << GetEffForRoot(a) 
            << ", Eff_b=" << GetEffForRoot(b) << "), "
            << "(fa=" << fa << ", fb=" << fb << "), "
            << "refVal = " << refVal << endl;
      return 1;
   }

   bool ac_equal(kFALSE);
   double fc = fb;
   double c  = 0, d = 0, e = 0;
   for(int iter= 0; iter <= fMaxIter; iter++)
   {
      if((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
      {
         // Rename a,b,c and adjust bounding interval d
         ac_equal = kTRUE;
         c = a;
	 fc = fa;
         d = b - a;
	 e = b - a;
      }
  
      if(TMath::Abs(fc) < TMath::Abs(fb))
      {
         ac_equal = kTRUE;
         a = b;
	 b = c;
	 c = a;
         fa = fb;
	 fb = fc;
	 fc = fa;
      }

      double tol = 0.5*2.2204460492503131e-16*TMath::Abs(b);
      double m = 0.5*(c - b);
      if(fb == 0 || TMath::Abs(m) <= tol || TMath::Abs(fb) < fAbsTol)
         return b;
  
      // Bounds decreasing too slowly: use bisection
      if(TMath::Abs(e) < tol || TMath::Abs(fa) <= TMath::Abs(fb))
      {
         d = m;
	 e = m;
      }      
      else
      {
         // Attempt inverse cubic interpolation
         double p, q, r;
         double s = fb/fa;
      
         if(ac_equal)
	 {
            p = 2*m*s;
	    q = 1 - s;
	 }
         else
	 {
            q = fa/fc;
	    r = fb/fc;
            p = s*(2*m*q*(q - r) - (b - a)*(r - 1));
            q = (q - 1)*(r - 1)*(s - 1);
         }
         // Check whether we are in bounds
         if(p > 0)
            q = -q;
         else
            p = -p;
      
         double min1 = 3*m*q - TMath::Abs(tol*q);
         double min2 = TMath::Abs(e*q);
         if(2*p < (min1 < min2 ? min1 : min2))
	 {
            // Accept the interpolation
            e = d;
	    d = p/q;
         }
         else
	 {
            d = m;
	    e = m;
	 } // Interpolation failed: use bisection.
      }
      // Move last best guess to a
      a = b;
      fa = fb;
      // Evaluate new trial root
      if(TMath::Abs(d) > tol)
         b += d;
      else
         b += (m > 0 ? +tol : -tol);

      fb = GetEffForRoot(b) - refVal;
   }

   // Return our best guess if we run out of iterations
   cout << "<ROCCalc::Root> maximum iterations (" << fMaxIter << ") reached before convergence" << endl;

   return b;
}

int RocCurve::GetOutBins()
{
   return fNbins[1];
}

/*RocCurve::RocCurve(TH1 *sigH, TH1 *bgH)
{
   dtemp = new double[4];
   invert = new bool;
   fNbins = new double[2];
   fXrange = new double[2];

   fsigH = sigH;
   fbgH = bgH;

   if(fsigH->GetMean() >= fbgH->GetMean())
      *invert = false;
   else
      *invert = true;

   cout << "Invert = " << (int)*invert << endl;

   fNbins[0] = fsigH->GetNbinsX();
   fNbins[1] = fNbins[0]*10;

   fXrange[0] = fsigH->GetXaxis()->GetXmin();
   fXrange[1] = fsigH->GetXaxis()->GetXmax();
   cout << "X range = (" << fXrange[0] << ", " << fXrange[1] << ")" << endl;

   // Cumulative distributions
   fcumsigH = new TH1D("fcumsigH", "", fNbins[0], fXrange[0], fXrange[1]);
   fcumbgH = new TH1D("fcumbgH", "", fNbins[0], fXrange[0], fXrange[1]);

   // Return ROC histogram
   rocH = new TH1D("rocH", "", fNbins[1], 0, 1);

   rootVal = new double[2];
}

RocCurve::~RocCurve()
{
   delete[] dtemp;
   delete invert;
   delete[] fNbins;
   delete[] fXrange;
   delete fcumsigH;
   delete fcumbgH;
   delete[] rootVal;
   delete rocH;
}

bool RocCurve::IsGood()
{
   if( (TMath::Abs(fXrange[1] - fbgH->GetXaxis()->GetXmax()) > 0.000001) || (TMath::Abs(fXrange[0] - fbgH->GetXaxis()->GetXmin()) > 0.000001) || (fbgH->GetNbinsX() != fsigH->GetNbinsX()) )
   {
      cout << "Error! Histograms don't have the same number of bins or range." << endl;
      return false;
   }
   else
      return true;
}

void RocCurve::GetCumulative(TH1 *inhist, TH1 *outhist)
{
   cout << "Entering GetCumulative function: Returns the cumulative version of a distribution" << endl;

   dtemp[0] = 0;
   dtemp[1] = 0;
   dtemp[2] = 0;
   dtemp[3] = 0;

   dtemp[3] = inhist->GetEntries();
   cout << "  Total number of entries = " << dtemp[3] << endl;

   for(int i = 1; i <= fNbins[0]; i++)
   {
      dtemp[0] = inhist->GetBinCenter(i);
      dtemp[1] = inhist->GetBinContent(i)/dtemp[3];
      dtemp[2] += dtemp[1];
      cout << "  " << i << "\t" << dtemp[0] << "\t" << dtemp[1] << "\t" << 1. - dtemp[2] << endl;

      if(!(*invert))
         outhist->SetBinContent(i, 1. - dtemp[2]);
      else
         outhist->SetBinContent(i, dtemp[2]);
   }
}

TH1D* RocCurve::GetROC()
{
   cout << "Entering GetROC function: Returns the ROC curve from signal and background." << endl;

   // Return cumulative distributions
   GetCumulative(fsigH, fcumsigH);
   GetCumulative(fbgH, fcumbgH);

   // Get the ROC curve
   for(int i = 1; i <= fNbins[1]; i++)
   {
      cout << endl << "i = " << i << ":" << endl;
      dtemp[0] = rocH->GetBinCenter(i);
      Root(dtemp[0]);
      dtemp[1] = rootVal[0];
      dtemp[2] = rootVal[1];
      cout << "  dtemp[0] = " << dtemp[0] << ", dtemp[1] = " << dtemp[1] << ", dtemp[2] = " << dtemp[2] << endl;

      rocH->SetBinContent(i, dtemp[0], 1. - dtemp[2]);
   }

   return rocH;
}

void RocCurve::Root(double refVal)
{
   sbin = new int;
   *sbin = -1;
   p1 = new double[2];
   p2 = new double[2];

   if(!(*invert))
   {
      // Search, where signal cumulative distribution is lower than reference signal efficiency value
      for(int i = 1; i <= fNbins[0]; i++)
      {
         if(fcumsigH->GetBinContent(i) < refVal)
         {
            *sbin = i;
            break;
         }
      }
   }
   else
   {
      // Search, where signal cumulative distribution is higher than reference signal efficiency value
      for(int i = 1; i <= fNbins[0]; i++)
      {
         if(fcumsigH->GetBinContent(i) > refVal)
         {
            *sbin = i;
            break;
         }
      }
   }

   cout << "  Using bin " << *sbin << " (value = " << fcumsigH->GetBinContent(*sbin) << ", reference value = " << refVal << ")" << endl;

   p1[0] = fcumsigH->GetBinCenter(*sbin-1);
   p1[1] = fcumsigH->GetBinContent(*sbin-1);
   p2[0] = fcumsigH->GetBinCenter(*sbin);
   p2[1] = fcumsigH->GetBinContent(*sbin);
   cout << "  Signal: Point 1 = (" << p1[0] << ", " << p1[1] << "), Point 2 = (" << p2[0] << ", " << p2[1] << ")" << endl;

   rootVal[0] = InterpolationX(refVal, p1, p2);

   p1[0] = fcumbgH->GetBinCenter(*sbin-1);
   p1[1] = fcumbgH->GetBinContent(*sbin-1);
   p2[0] = fcumbgH->GetBinCenter(*sbin);
   p2[1] = fcumbgH->GetBinContent(*sbin);
   cout << "  Background: Point 1 = (" << p1[0] << ", " << p1[1] << "), Point 2 = (" << p2[0] << ", " << p2[1] << ")" << endl;

   rootVal[1] = InterpolationY(rootVal[0], p1, p2);

   cout << "  Return value corresponding to " << refVal << " signal = (" << rootVal[0] << ", " << rootVal[1] << ")" << endl;

   delete sbin;
   delete[] p1;
   delete[] p2;
}

double RocCurve::InterpolationX(double yref, double *point1, double *point2)
{
   double *linPar;
   linPar = new double[2];

   // slope of the linear curve
   linPar[0] = (point2[1] - point1[1])/(point2[0] - point1[0]);
   // starting value of the linear curve
   linPar[1] = point1[1] - point1[0]*linPar[0];

   // x value corresponding to reference y value
   ret = (yref - linPar[1])/linPar[0];

   cout << "  a = " << linPar[0] << ", b = " << linPar[1] << ", xret = " << ret << endl;

   delete[] linPar;

   return ret;
}

double RocCurve::InterpolationY(double xref, double *point1, double *point2)
{
   double *linPar;
   linPar = new double[2];

   // slope of the linear curve
   linPar[0] = (point2[1] - point1[1])/(point2[0] - point1[0]);
   // starting value of the linear curve
   linPar[1] = point1[1] - point1[0]*linPar[0];

   // y value corresponding to reference x value
   ret = linPar[0]*xref + linPar[1];

   cout << "  a = " << linPar[0] << ", b = " << linPar[1] << ", yret = " << ret << endl;

   delete[] linPar;

   return ret;
}

int RocCurve::GetOutBins()
{
   return fNbins[1];
}*/
