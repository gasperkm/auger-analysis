#define _STANDALONE_ 1
#include "workstation.h"
#include "separate_functions.h"
#include "root_style.h"
#include "mva_result_read.h"

using namespace std;

void SetColor(TGraphAsymmErrors *gr, int cur, int nrbins)
{
   int ci = 1738+cur;
   double colormix = (double)cur/(double)(nrbins-1.);
   TColor *color = new TColor(ci, colormix, 0, 1.-colormix, "", 1);

   gr->SetMarkerColor(ci);
   gr->SetMarkerSize(0.9);
   gr->SetMarkerStyle(20+cur);
   gr->SetLineColor(ci);
   gr->SetLineWidth(2);
}

int main(int argc, char **argv)
{
   float yrange[2];
   float fraction;
   int onlyplot = -1;
   float *ftemp;
   char *ctemp;
   string *stemp;
   int *itemp;
   double *dtemp;
   bool runnorm;
   bool autoFrac;
   bool xErrors;
   string plotInstr[2];

   int nrbins = -1;

   ResultRead *analRes = new ResultRead();

   string filename;
   if(argc > 1)
   {
      ftemp = new float[3];
      ctemp = new char;
      stemp = new string[3];
      itemp = new int[2];
      dtemp = new double[2];

      bool plottest[2] = {false, false};

      for(int i = 0; i < argc-1; i++)
      {
         stemp[0] = string(argv[i+1]);
	 itemp[0] = stemp[0].find("individual_results");
	 if(RemoveBeforeLast(&stemp[0], ".") == "root")
	    plottest[0] = plottest[0] || true;

	 itemp[0] = stemp[0].find("application_results.txt");
	 if(itemp[0] != string::npos)
	    plottest[1] = plottest[1] || true;

	 itemp[0] = stemp[0].find("individual_results");
	 if(itemp[0] != string::npos)
	    plottest[1] = plottest[1] || true;
      }

      if(plottest[0] && !plottest[1])
         onlyplot = 1;
      else if(!plottest[0] && plottest[1])
	 onlyplot = 0;
      else
	 onlyplot = -1;

      if(onlyplot == -1)
      {
         cout << "Run the fraction plots script (0) or combine plots from created fraction plots (1)? ";
         cin >> onlyplot;
      }

      if( (onlyplot == 0) || (onlyplot == 1) )
      {
         cout << "Represent x-axis bins with error bars (0) or not (1)? ";
	 cin >> itemp[0];
	 while( (itemp[0] != 0) && (itemp[0] != 1) )
	 {
            cout << "Represent x-axis bins with error bars (0) or not (1)? ";
	    cin >> itemp[0];
	 }

	 if(itemp[0] == 0)
	 {
            xErrors = true;
	    plotInstr[0] = "AP";
	    plotInstr[1] = "P;SAME";
	    cout << "Using errors." << endl;
	 }
	 else if(itemp[0] == 1)
	 {
            xErrors = false;
	    plotInstr[0] = "APL";
	    plotInstr[1] = "PL;SAME";
	    cout << "Not using errors." << endl;
	 }
      }

      if(onlyplot == 0)
      {
         cout << "For linear normalization, please enter the fraction of data events taken as signal (if unknown, enter -1): ";
         cin >> fraction;
	 if(fraction == -1)
            autoFrac = true;
	 else
            autoFrac = false;

         nrbins = argc-1;
         cout << "Number of bins on the plot = " << nrbins << endl;

         // Arrays for plotting
         float *xbin[3];
         float *ybinSig[3];
         float *ybinSigNorm[3];
         float *ybinBack[3];
         float *ybinBackNorm[3];
         float *ybinData[3];
         float *ybinDataNorm[3];
         float *ybinDataLna[3];
         float *ybinDataLnaNorm[3];
	 float *ylowfractest[3], *yhighfractest[3];
	 float *ylowlimitfrac, *yhighlimitfrac;

	 vector<string> treeName;

         for(int i = 0; i < 3; i++)
         {
            xbin[i] = new float[nrbins];
            ybinSig[i] = new float[nrbins];
            ybinSigNorm[i] = new float[nrbins];
            ybinBack[i] = new float[nrbins];
            ybinBackNorm[i] = new float[nrbins];
            ybinData[i] = new float[nrbins];
            ybinDataNorm[i] = new float[nrbins];
            ybinDataLna[i] = new float[nrbins];
            ybinDataLnaNorm[i] = new float[nrbins];

	    if(autoFrac)
	    {
	       ylowfractest[i] = new float[nrbins];
	       yhighfractest[i] = new float[nrbins];
	    }
         }

	 if(autoFrac)
         {
            ylowlimitfrac = new float[nrbins];
            yhighlimitfrac = new float[nrbins];
	 }

         // Going through all input files
         for(int i = 0; i < nrbins; i++)
         {
            filename = string(argv[i+1]);
            analRes->ReadFile(filename);

	    // Get tree name
	    treeName.push_back(analRes->GetTreeName(i));

            // Getting xbin values
            xbin[0][i] = analRes->GetEnergy();
            analRes->GetEnergyError(ftemp);
	    if(xErrors)
	    {
               xbin[1][i] = ftemp[0];
               xbin[2][i] = ftemp[1];
	    }
	    else
	    {
               xbin[1][i] = 0;
               xbin[2][i] = 0;
	    }

            // Getting signal ybin values
            ybinSig[0][i] = analRes->GetFraction(1, -1);
            analRes->GetFractionError(ftemp);
            ybinSig[1][i] = ftemp[0];
            ybinSig[2][i] = ftemp[1];

            // Getting background ybin values
            ybinBack[0][i] = analRes->GetFraction(0, -1);
            analRes->GetFractionError(ftemp);
            ybinBack[1][i] = ftemp[0];
            ybinBack[2][i] = ftemp[1];

            // Getting data ybin values
            ybinData[0][i] = analRes->GetFraction(2, -1);
            analRes->GetFractionError(ftemp);
            ybinData[1][i] = ftemp[0];
            ybinData[2][i] = ftemp[1];

            if(!autoFrac)
            {
               ybinSigNorm[0][i] = analRes->GetFraction(1, fraction);
               analRes->GetFractionError(ftemp);
               ybinSigNorm[1][i] = ftemp[0];
               ybinSigNorm[2][i] = ftemp[1];

               ybinBackNorm[0][i] = analRes->GetFraction(0, fraction);
               analRes->GetFractionError(ftemp);
               ybinBackNorm[1][i] = ftemp[0];
               ybinBackNorm[2][i] = ftemp[1];

               // Getting data ybin values (normalized)
               ybinDataNorm[0][i] = analRes->GetFraction(2, fraction);
               analRes->GetFractionError(ftemp);
               ybinDataNorm[1][i] = ftemp[0];
               ybinDataNorm[2][i] = ftemp[1];

	       runnorm = true;
            }
	    else
	    {
	       // Getting the extreme value of normalization for data -> fraction is 0
               ylowfractest[0][i] = analRes->GetFraction(2, 0.);
               analRes->GetFractionError(ftemp);
               ylowfractest[1][i] = ftemp[0];
               ylowfractest[2][i] = ftemp[1];

	       // Getting the extreme value of normalization for data -> fraction is 1
               yhighfractest[0][i] = analRes->GetFraction(2, 1.);
               analRes->GetFractionError(ftemp);
               yhighfractest[1][i] = ftemp[0];
               yhighfractest[2][i] = ftemp[1];

	       ybinSigNorm[0][i] = 0.;
	       ybinSigNorm[1][i] = 0.;
	       ybinSigNorm[2][i] = 0.;

	       ybinBackNorm[0][i] = 0.;
	       ybinBackNorm[1][i] = 0.;
	       ybinBackNorm[2][i] = 0.;

	       ybinDataNorm[0][i] = 0.;
	       ybinDataNorm[1][i] = 0.;
	       ybinDataNorm[2][i] = 0.;

	       runnorm = false;
	    }
         }

	 // Set normalization for data
	 if(autoFrac)
	 {
            // Printout original signal fraction for data and the two extreme normalizations (0 and 1)
	    cout << "Low and high fraction tests: " << endl;
            for(int i = 0; i < nrbins; i++)
	       cout << i << ": orig = " << ybinData[0][i] << " (" << ybinData[1][i] << "," << ybinData[2][i] << ")\t norm0 = " << ylowfractest[0][i] << " (" << ylowfractest[1][i] << "," << ylowfractest[2][i] << ")\t norm1 = " << yhighfractest[0][i] << " (" << yhighfractest[1][i] << "," << yhighfractest[2][i] << ")" << endl;

	    // The two normalizations will give a possible signal fraction (including errors) for each point
	    cout << "Lowest and highest possible fractions for each point: " << endl;
            for(int i = 0; i < nrbins; i++)
	    {
               ftemp[0] = ylowfractest[0][i];
               ftemp[1] = yhighfractest[0][i];

	       if(TMath::MinElement(2,ftemp) == ylowfractest[0][i])
	       {
                  ylowlimitfrac[i] = (ylowfractest[0][i]-ylowfractest[1][i]);
                  yhighlimitfrac[i] = (yhighfractest[0][i]+yhighfractest[2][i]);
	       }
	       else
	       {
                  ylowlimitfrac[i] = (yhighfractest[0][i]-yhighfractest[1][i]);
                  yhighlimitfrac[i] = (ylowfractest[0][i]+ylowfractest[2][i]);
	       }

	       cout << i << ": lowest (ylow) = " << ylowlimitfrac[i] << "\t highest (yhigh) = " << yhighlimitfrac[i] << endl;
	    }

	    // Rerun normalization, but using the new range of signal fractions for each point
            for(int i = 0; i < nrbins; i++)
	    {
               filename = string(argv[i+1]);
               analRes->ReadFile(filename);

	       // Limit the low limit fraction to 0 and 1
/*	       if(ylowlimitfrac[i] < 0.)
                  ylowlimitfrac[i] = 0.;
	       if(ylowlimitfrac[i] > 1.)
                  ylowlimitfrac[i] = 1.;*/

	       cout << i << ": Normalization using lower fraction " << ylowlimitfrac[i] << ": " << endl;
               ylowfractest[0][i] = analRes->GetFraction(2, ylowlimitfrac[i]);
               analRes->GetFractionError(ftemp);
               ylowfractest[1][i] = ftemp[0];
               ylowfractest[2][i] = ftemp[1];

	       cout << "   orig = " << ybinData[0][i] << " (" << ybinData[1][i] << "," << ybinData[2][i] << ")\t norm = " << ylowfractest[0][i] << " (" << ylowfractest[1][i] << "," << ylowfractest[2][i] << ")" << endl;

/*	       // Limit the high limit fraction to 0 and 1
	       if(yhighlimitfrac[i] < 0.)
                  yhighlimitfrac[i] = 0.;
	       if(yhighlimitfrac[i] > 1.)
                  yhighlimitfrac[i] = 1.;*/

	       cout << i << ": Normalization using higher fraction " << yhighlimitfrac[i] << ": " << endl;
               yhighfractest[0][i] = analRes->GetFraction(2, yhighlimitfrac[i]);
               analRes->GetFractionError(ftemp);
               yhighfractest[1][i] = ftemp[0];
               yhighfractest[2][i] = ftemp[1];

	       cout << "   orig = " << ybinData[0][i] << " (" << ybinData[1][i] << "," << ybinData[2][i] << ")\t norm = " << yhighfractest[0][i] << " (" << yhighfractest[1][i] << "," << yhighfractest[2][i] << ")" << endl;

	       // Determine the final mean value and errors
	       cout << i << ": Final normalization: " << endl;
               ybinDataNorm[0][i] = (ylowfractest[0][i] + yhighfractest[0][i])/2.;

	       ftemp[0] = ylowfractest[0][i];
	       ftemp[1] = yhighfractest[0][i];
	       if(TMath::MinElement(2,ftemp) == ylowfractest[0][i])
	       {
                  ybinDataNorm[1][i] = TMath::Abs((ylowfractest[0][i] - ylowfractest[1][i]) - ybinDataNorm[0][i]);
                  ybinDataNorm[2][i] = TMath::Abs((yhighfractest[0][i] + yhighfractest[2][i]) - ybinDataNorm[0][i]);
	       }
	       else
	       {
                  ybinDataNorm[1][i] = TMath::Abs((yhighfractest[0][i] - yhighfractest[1][i]) - ybinDataNorm[0][i]);
                  ybinDataNorm[2][i] = TMath::Abs((ylowfractest[0][i] + ylowfractest[2][i]) - ybinDataNorm[0][i]);
	       }

/*               ybinDataNorm[1][i] = (ylowfractest[1][i] + yhighfractest[1][i])/2.;
               ybinDataNorm[2][i] = (ylowfractest[2][i] + yhighfractest[2][i])/2.;*/
/*	       dtemp[0] = TMath::Power((ylowfractest[1][i] - ylowfractest[0][i]), 2) + TMath::Power((yhighfractest[1][i] - yhighfractest[0][i]), 2);
	       ybinDataNorm[1][i] = TMath::Sqrt(dtemp[0]/2.);
	       dtemp[1] = TMath::Power((ylowfractest[2][i] - ylowfractest[0][i]), 2) + TMath::Power((yhighfractest[2][i] - yhighfractest[0][i]), 2);
	       ybinDataNorm[2][i] = TMath::Sqrt(dtemp[1]/2.);
*/
	       cout << "   orig = " << ybinData[0][i] << " (" << ybinData[1][i] << "," << ybinData[2][i] << ")\t norm = " << ybinDataNorm[0][i] << " (" << ybinDataNorm[1][i] << "," << ybinDataNorm[2][i] << ")\t diff = " << (ybinDataNorm[0][i] - ybinData[0][i]) << endl;

	       fraction = ybinDataNorm[0][i];

               // Getting signal ybin values (normalized)
               ybinSigNorm[0][i] = analRes->GetFraction(1, fraction);
               analRes->GetFractionError(ftemp);
               ybinSigNorm[1][i] = ftemp[0];
               ybinSigNorm[2][i] = ftemp[1];

               // Getting signal ybin values (normalized)
               ybinBackNorm[0][i] = analRes->GetFraction(0, fraction);
               analRes->GetFractionError(ftemp);
               ybinBackNorm[1][i] = ftemp[0];
               ybinBackNorm[2][i] = ftemp[1];
	    }

	    // Decide if we want to create plots or not
	    runnorm = true;
	 }

	 // Calculate LnA values
	 for(int i = 0; i < nrbins; i++)
	 {
	    // Old calculation
            // Raw data (mean + propagation of errors)
            ybinDataLna[0][i] = TMath::Log(56)*(ybinSig[0][i] - ybinData[0][i])/(ybinSig[0][i] - (1.-ybinBack[0][i]));
	    ftemp[0] = ybinSig[1][i] + ybinData[1][i];
	    ftemp[1] = ybinSig[1][i] + ybinBack[1][i];
	    ybinDataLna[1][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);
	    ftemp[0] = ybinSig[2][i] + ybinData[2][i];
	    ftemp[1] = ybinSig[2][i] + ybinBack[2][i];
	    ybinDataLna[2][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);

//	    cerr << i << ": Raw = " << ybinDataLna[0][i] << "\t"  << ybinDataLna[1][i] << "\t" << ybinDataLna[2][i] << endl;

	    // Normalized data (mean + propagation of errors)
            ybinDataLnaNorm[0][i] = TMath::Log(56)*(1. - ybinDataNorm[0][i]);
            ybinDataLnaNorm[1][i] = TMath::Log(56)*(1. - (ybinDataNorm[0][i]-ybinDataNorm[1][i]));
            ybinDataLnaNorm[1][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[1][i]);
            ybinDataLnaNorm[2][i] = TMath::Log(56)*(1. - (ybinDataNorm[0][i]+ybinDataNorm[2][i]));
            ybinDataLnaNorm[2][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[2][i]);

/*            // New calculation
            // Raw data (mean + propagation of errors)
            ybinDataLna[0][i] = TMath::Log(56)*(ybinSig[0][i] - ybinData[0][i])/(ybinSig[0][i] - (1.-ybinBack[0][i]));
	    ftemp[0] = ybinSig[1][i] + ybinData[1][i];
	    ftemp[1] = ybinSig[1][i] + ybinBack[1][i];
	    ybinDataLna[1][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);
	    ftemp[0] = ybinSig[2][i] + ybinData[2][i];
	    ftemp[1] = ybinSig[2][i] + ybinBack[2][i];
	    ybinDataLna[2][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);

//	    cerr << i << ": Raw = " << ybinDataLna[0][i] << "\t"  << ybinDataLna[1][i] << "\t" << ybinDataLna[2][i] << endl;

	    // Normalized data (mean + propagation of errors)
	    if(ybinDataNorm[0][i] <= 1.)
               ybinDataLnaNorm[0][i] = TMath::Log(56)*TMath::Sqrt(1. - ybinDataNorm[0][i]);
	    else
               ybinDataLnaNorm[0][i] = TMath::Log(56)*(1. - ybinDataNorm[0][i]);

	    if(ybinDataNorm[1][i] <= 1.)
               ybinDataLnaNorm[1][i] = TMath::Log(56)*TMath::Sqrt(1. - (ybinDataNorm[0][i]-ybinDataNorm[1][i]));
	    else
               ybinDataLnaNorm[1][i] = TMath::Log(56)*(1. - (ybinDataNorm[0][i]-ybinDataNorm[1][i]));
            ybinDataLnaNorm[1][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[1][i]);

	    if(ybinDataNorm[2][i] <= 1.)
               ybinDataLnaNorm[2][i] = TMath::Log(56)*TMath::Sqrt(1. - (ybinDataNorm[0][i]+ybinDataNorm[2][i]));
	    else
               ybinDataLnaNorm[2][i] = TMath::Log(56)*(1. - (ybinDataNorm[0][i]+ybinDataNorm[2][i]));
            ybinDataLnaNorm[2][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[2][i]);*/

/*            ybinDataLnaNorm[1][i] = TMath::Log(56)*(ybinSigNorm[0][i] - (ybinDataNorm[0][i]-ybinDataNorm[1][i])/(ybinSigNorm[0][i] - (1.-ybinBackNorm[0][i]));
            ybinDataLnaNorm[2][i] = TMath::Log(56)*(ybinSigNorm[0][i] - ybinDataNorm[0][i])/(ybinSigNorm[0][i] - (1.-ybinBackNorm[0][i]));*/
/*	    ftemp[0] = ybinSigNorm[1][i] + ybinDataNorm[1][i];
	    ftemp[1] = ybinSigNorm[1][i] + ybinBackNorm[1][i];
	    ybinDataLnaNorm[1][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);
	    ftemp[0] = ybinSigNorm[2][i] + ybinDataNorm[2][i];
	    ftemp[1] = ybinSigNorm[2][i] + ybinBackNorm[2][i];
	    ybinDataLnaNorm[2][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);*/

//	    cerr << i << ": Normalized = " << ybinDataLnaNorm[0][i] << "\t"  << ybinDataLnaNorm[1][i] << "\t" << ybinDataLnaNorm[2][i] << endl;
	 }

	 cout << endl << "Point printout:" << endl;
	 ftemp[0] = 0;
	 ftemp[1] = 0;
	 ftemp[2] = 0;
	 // Skip the first point, just because
	 for(int i = 1; i < nrbins; i++)
	 {
            cerr << i << ": Signal     = " << ybinSig[0][i] << "\t" << ybinSig[1][i] << "\t" << ybinSig[2][i] << endl;
            cerr << i << ": Background = " << ybinBack[0][i] << "\t" << ybinBack[1][i] << "\t" << ybinBack[2][i] << endl;
            cerr << i << ": Raw data   = " << ybinData[0][i] << "\t" << ybinData[1][i] << "\t" << ybinData[2][i] << endl;
            cerr << i << ": Norm. data = " << ybinDataNorm[0][i] << "\t" << ybinDataNorm[1][i] << "\t" << ybinDataNorm[2][i] << endl;
	    ftemp[0] += ybinDataNorm[0][i];
	    ftemp[1] += ybinDataNorm[1][i];
	    ftemp[2] += ybinDataNorm[2][i];
	 }
	 cout << endl << "Mean norm. data value = " << ftemp[0]/(nrbins-1.) << "\t" << ftemp[1]/(nrbins-1.) << "\t" << ftemp[2]/(nrbins-1.) << endl;
	 cout << endl;

         // Minimum and maximum values for all bins and types
         cerr << "Minimum and maximum values for all bins (without error bars):" << endl;
	 ftemp[0] = TMath::MinElement(nrbins, ybinSig[0]);
	 ftemp[1] = TMath::MaxElement(nrbins, ybinSig[0]);
         cerr << "- Signal:           \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	 ftemp[0] = TMath::MinElement(nrbins, ybinBack[0]);
	 ftemp[1] = TMath::MaxElement(nrbins, ybinBack[0]);
         cerr << "- Background:       \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	 ftemp[0] = TMath::MinElement(nrbins, ybinData[0]);
	 ftemp[1] = TMath::MaxElement(nrbins, ybinData[0]);
         cerr << "- Data:             \t" << ftemp[0] << "\t" << ftemp[1] << endl;
         if(runnorm)
	 {
	    ftemp[0] = TMath::MinElement(nrbins, ybinSigNorm[0]);
	    ftemp[1] = TMath::MaxElement(nrbins, ybinSigNorm[0]);
            cerr << "- Norm. Signal:     \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	    ftemp[0] = TMath::MinElement(nrbins, ybinBackNorm[0]);
	    ftemp[1] = TMath::MaxElement(nrbins, ybinBackNorm[0]);
            cerr << "- Norm. Background: \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	    ftemp[0] = TMath::MinElement(nrbins, ybinDataNorm[0]);
	    ftemp[1] = TMath::MaxElement(nrbins, ybinDataNorm[0]);
            cerr << "- Norm. Data:       \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	 }

         cerr << "Largest possible minimum and maximum values for all bins (with error bars):" << endl;
	 ftemp[0] = TMath::MinElement(nrbins, ybinSig[0]) - TMath::MaxElement(nrbins, ybinSig[1]);
	 ftemp[1] = TMath::MaxElement(nrbins, ybinSig[0]) + TMath::MaxElement(nrbins, ybinSig[2]);
         cerr << "- Signal:           \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	 ftemp[0] = TMath::MinElement(nrbins, ybinBack[0]) - TMath::MaxElement(nrbins, ybinBack[1]);
	 ftemp[1] = TMath::MaxElement(nrbins, ybinBack[0]) + TMath::MaxElement(nrbins, ybinBack[2]);
         cerr << "- Background:       \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	 ftemp[0] = TMath::MinElement(nrbins, ybinData[0]) - TMath::MaxElement(nrbins, ybinData[1]);
	 ftemp[1] = TMath::MaxElement(nrbins, ybinData[0]) + TMath::MaxElement(nrbins, ybinData[2]);
         cerr << "- Data:             \t" << ftemp[0] << "\t" << ftemp[1] << endl;
         if(runnorm)
	 {
	    ftemp[0] = TMath::MinElement(nrbins, ybinSigNorm[0]) - TMath::MaxElement(nrbins, ybinSigNorm[1]);
	    ftemp[1] = TMath::MaxElement(nrbins, ybinSigNorm[0]) + TMath::MaxElement(nrbins, ybinSigNorm[2]);
            cerr << "- Norm. Signal:     \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	    ftemp[0] = TMath::MinElement(nrbins, ybinBackNorm[0]) - TMath::MaxElement(nrbins, ybinBackNorm[1]);
	    ftemp[1] = TMath::MaxElement(nrbins, ybinBackNorm[0]) + TMath::MaxElement(nrbins, ybinBackNorm[2]);
            cerr << "- Norm. Background: \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	    ftemp[0] = TMath::MinElement(nrbins, ybinDataNorm[0]) - TMath::MaxElement(nrbins, ybinDataNorm[1]);
	    ftemp[1] = TMath::MaxElement(nrbins, ybinDataNorm[0]) + TMath::MaxElement(nrbins, ybinDataNorm[2]);
            cerr << "- Norm. Data:       \t" << ftemp[0] << "\t" << ftemp[1] << endl;
	 }

         // y-range settings
         cerr << endl << "Please set the y-axis range for all plots (comma separated): ";
         cin >> yrange[0] >> *ctemp >> yrange[1];

         // Create directory structure for plots and delete old plots
         stemp[2] = RemoveFilename(&filename);
         itemp[0] = analRes->GetFileType();
         if(itemp[0] == 0)
            stemp[0] = RemoveFilename(&stemp[2]);
         else if(itemp[0] == 1)
            stemp[0] = stemp[2];

         stemp[1] = "mkdir -p " + RemoveFilename(&stemp[0]) + "/plots";
         system(stemp[1].c_str());

         if(itemp[0] == 0)
         {
            stemp[2] = analRes->GetObservableType();
            stemp[1] = "rm -fr " + RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[2] + "*";
            system(stemp[1].c_str());
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_all.root";
         }
         else if(itemp[0] == 1)
         {
            stemp[1] = "rm -fr " + RemoveFilename(&stemp[0]) + "/plots/fraction_plot_*";
            system(stemp[1].c_str());
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_" + stemp[1] + "_all.root";
         }

         TFile *printRes = TFile::Open(stemp[2].c_str(), "RECREATE");;

         // Prepare canvases and clear any ROOT default statistics
         RootStyle *mystyle = new RootStyle();
         mystyle->SetBaseStyle();

         TCanvas *c1 = new TCanvas("c1","",1200,900);

         // Preparing all graphs
	 yrange[1] += 0.25*(yrange[1]-yrange[0]);
         // Signal
         TGraphAsymmErrors *grSig = new TGraphAsymmErrors(nrbins, xbin[0], ybinSig[0], xbin[1], xbin[2], ybinSig[1], ybinSig[2]);
         mystyle->SetGraphColor(grSig, 1);
         grSig->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

         // Background
         TGraphAsymmErrors *grBack = new TGraphAsymmErrors(nrbins, xbin[0], ybinBack[0], xbin[1], xbin[2], ybinBack[1], ybinBack[2]);
         mystyle->SetGraphColor(grBack, 0);
         grBack->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

         // Data
         TGraphAsymmErrors *grData = new TGraphAsymmErrors(nrbins, xbin[0], ybinData[0], xbin[1], xbin[2], ybinData[1], ybinData[2]);
         mystyle->SetGraphColor(grData, 2);
         grData->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

         TGraphAsymmErrors *grSigNorm, *grBackNorm, *grDataNorm;
         // Normalized values
         if(runnorm)
         {
            grSigNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinSigNorm[0], xbin[1], xbin[2], ybinSigNorm[1], ybinSigNorm[2]);
            mystyle->SetGraphColor(grSigNorm, 1);
            grSigNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

            grBackNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinBackNorm[0], xbin[1], xbin[2], ybinBackNorm[1], ybinBackNorm[2]);
            mystyle->SetGraphColor(grBackNorm, 0);
            grBackNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

            grDataNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataNorm[0], xbin[1], xbin[2], ybinDataNorm[1], ybinDataNorm[2]);
            mystyle->SetGraphColor(grDataNorm, 3);
            grDataNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);
         }

	 // LnA graphs for data and normalized data
	 TGraphAsymmErrors *grDataLna, *grDataLnaNorm;
         grDataLna = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataLna[0], xbin[1], xbin[2], ybinDataLna[1], ybinDataLna[2]);
         mystyle->SetGraphColor(grDataLna, 2);
         grDataLna->GetYaxis()->SetRangeUser(-0.7, 5.);
         if(runnorm)
         {
            grDataLnaNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataLnaNorm[0], xbin[1], xbin[2], ybinDataLnaNorm[1], ybinDataLnaNorm[2]);
            mystyle->SetGraphColor(grDataLnaNorm, 3);
            grDataLnaNorm->GetYaxis()->SetRangeUser(-0.7, 5.);
         }

         // Plotting each graph separately
	 // Signal only
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of signal simulation events");
         grSig->Draw(plotInstr[0].c_str());
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_sig-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            if(autoFrac)
	       stemp[1] = "frac-auto";
	    else
               stemp[1] = "frac-" + ToString(fraction, 3);
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_sig-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         grSig->Write(("sig_" + stemp[1]).c_str());

	 // Background only
         mystyle->SetAxisTitles(grBack, "FD energy [log(E/eV)]", "Purity of background simulation events");
         grBack->Draw(plotInstr[0].c_str());
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_back-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            if(autoFrac)
	       stemp[1] = "frac-auto";
	    else
               stemp[1] = "frac-" + ToString(fraction, 3);
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_back-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         grBack->Write(("back_" + stemp[1]).c_str());

	 // Raw data only
         mystyle->SetAxisTitles(grData, "FD energy [log(E/eV)]", "Signal fraction of data events");
         grData->Draw(plotInstr[0].c_str());
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_data-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            if(autoFrac)
	       stemp[1] = "frac-auto";
	    else
               stemp[1] = "frac-" + ToString(fraction, 3);
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_data-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         grData->Write(("data_" + stemp[1]).c_str());

         if(runnorm)
         {
            // Normalized signal only
            mystyle->SetAxisTitles(grSigNorm, "FD energy [log(E/eV)]", "Purity of signal simulation events");
            grSigNorm->Draw(plotInstr[0].c_str());
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_sig-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               if(autoFrac)
	          stemp[1] = "frac-auto";
	       else
                  stemp[1] = "frac-" + ToString(fraction, 3);
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_sig-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            grSigNorm->Write(("signorm_" + stemp[1]).c_str());

	    // Normalized background only
            mystyle->SetAxisTitles(grBackNorm, "FD energy [log(E/eV)]", "Purity of background simulation events");
            grBackNorm->Draw(plotInstr[0].c_str());
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_back-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               if(autoFrac)
	          stemp[1] = "frac-auto";
	       else
                  stemp[1] = "frac-" + ToString(fraction, 3);
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_back-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            grBackNorm->Write(("backnorm_" + stemp[1]).c_str());

	    // Normalized data only
            mystyle->SetAxisTitles(grDataNorm, "FD energy [log(E/eV)]", "Signal fraction of data events");
            grDataNorm->Draw(plotInstr[0].c_str());
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_data-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               if(autoFrac)
	          stemp[1] = "frac-auto";
	       else
                  stemp[1] = "frac-" + ToString(fraction, 3);
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_data-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            grDataNorm->Write(("datanorm_" + stemp[1]).c_str());
         }

         // Plotting combinations (signal/background, signal_norm/background_norm, data/data-norm)
	 TLegend *legend;
	 int c_Legend = TColor::GetColor("#ffff66");
	 int nrfigs = 2;
	 // Signal/background
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events");
         grSig->Draw(plotInstr[0].c_str());
         grBack->Draw(plotInstr[1].c_str());
	 legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	 legend->SetFillStyle(1001);
	 legend->SetFillColor(c_Legend);
	 legend->AddEntry(grSig, "Signal", "lp");
	 legend->AddEntry(grBack, "Background", "lp");
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_sig-back.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_sig-back.pdf";
            c1->SaveAs(stemp[2].c_str());
         }

         if(runnorm)
         {
	    // Normalized signal/background
            mystyle->SetAxisTitles(grSigNorm, "FD energy [log(E/eV)]", "Purity of simulation events");
            grSigNorm->Draw(plotInstr[0].c_str());
            grBackNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grSigNorm, "Signal (normalized)", "lp");
	    legend->AddEntry(grBackNorm, "Background (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_sig-back_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_sig-back_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }

	    // Data/data-norm
            mystyle->SetAxisTitles(grData, "FD energy [log(E/eV)]", "Signal fraction of data events");
//	    grData->SetLineStyle(9);
            grData->Draw(plotInstr[0].c_str());
            grDataNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grData, "Data (raw)", "lp");
	    legend->AddEntry(grDataNorm, "Data (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
//	    grData->SetLineStyle(1);

	    // LnA data/data-norm
            mystyle->SetAxisTitles(grDataLna, "FD energy [log(E/eV)]", "<lnA> of data events");
//	    grDataLna->SetLineStyle(9);
            grDataLna->Draw(plotInstr[0].c_str());
            grDataLnaNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grDataLna, "Data (raw)", "lp");
	    legend->AddEntry(grDataLnaNorm, "Data (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_lnA_" + stemp[1] + "_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_lnA_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
//	    grDataLna->SetLineStyle(1);
	 }

         // Collections of three (signal/background/data, signal/background/data-norm, signal-norm/background-norm/data, signal-norm/background-norm/data-norm, signal/background/data/data-norm)
	 nrfigs = 3;
	 // Signal/background/data
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
         grSig->Draw(plotInstr[0].c_str());
         grBack->Draw(plotInstr[1].c_str());
         grData->Draw(plotInstr[1].c_str());
	 legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	 legend->SetFillStyle(1001);
	 legend->SetFillColor(c_Legend);
	 legend->AddEntry(grSig, "Signal", "lp");
	 legend->AddEntry(grBack, "Background", "lp");
	 legend->AddEntry(grData, "Data (raw)", "lp");
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_all.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_all.pdf";
            c1->SaveAs(stemp[2].c_str());
         }

         if(runnorm)
         {
            // Signal/background/data-norm
            mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
            grSig->Draw(plotInstr[0].c_str());
            grBack->Draw(plotInstr[1].c_str());
            grDataNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grSig, "Signal", "lp");
	    legend->AddEntry(grBack, "Background", "lp");
	    legend->AddEntry(grDataNorm, "Data (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_all-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_all-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }

            // Signal-norm/background-norm/data
            mystyle->SetAxisTitles(grSigNorm, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
            grSigNorm->Draw(plotInstr[0].c_str());
            grBackNorm->Draw(plotInstr[1].c_str());
            grData->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grSigNorm, "Signal (normalized)", "lp");
	    legend->AddEntry(grBackNorm, "Background (normalized)", "lp");
	    legend->AddEntry(grDataNorm, "Data (raw)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_all-sim_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_all-sim_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }

	    nrfigs = 4;
            // Signal-norm/background-norm/data-norm
            mystyle->SetAxisTitles(grSigNorm, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
            grSigNorm->Draw(plotInstr[0].c_str());
            grBackNorm->Draw(plotInstr[1].c_str());
            grDataNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grSigNorm, "Signal (normalized)", "lp");
	    legend->AddEntry(grBackNorm, "Background (normalized)", "lp");
	    legend->AddEntry(grDataNorm, "Data (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_all_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_all_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }

            // Signal/background/data/data-norm
            mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
            grSig->Draw(plotInstr[0].c_str());
            grBack->Draw(plotInstr[1].c_str());
//	    grData->SetLineStyle(9);
            grData->Draw(plotInstr[1].c_str());
            grDataNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(4*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grSig, "Signal", "lp");
	    legend->AddEntry(grBack, "Background", "lp");
	    legend->AddEntry(grData, "Data (raw)", "lp");
	    legend->AddEntry(grDataNorm, "Data (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_all-data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_all-data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
//	    grData->SetLineStyle(1);
         }

	 TText *t = new TText();
         t->SetTextAlign(12);
         t->SetTextColor(28);
         t->SetTextSize(24);
	 TLine *l = new TLine();
         l->SetLineWidth(2);
         l->SetLineStyle(7);
	 l->SetLineColor(28);

	 nrfigs = 2;

	 // LnA raw data only
         mystyle->SetAxisTitles(grDataLna, "FD energy [log(E/eV)]", "<lnA> of data events");
         grDataLna->Draw(plotInstr[0].c_str());

	 // Right side markings
         c1->Update();
         t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
         t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
         t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
         t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(16), "O");
         t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");
         // Draw lines marking the particle types
         l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
         l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
         l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
         l->DrawLine(gPad->GetUxmin(), TMath::Log(16), gPad->GetUxmax(), TMath::Log(16));
         l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_lnA_" + stemp[1] + "_data-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
/*            if(autoFrac)
	       stemp[1] = "frac-auto";
	    else
               stemp[1] = "frac-" + ToString(fraction, 3);*/
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_lnA_data-only.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
//         grDataLna->Write(("datalna_" + stemp[1]).c_str());

         if(runnorm)
         {
	    // LnA normalized data only
            mystyle->SetAxisTitles(grDataLnaNorm, "FD energy [log(E/eV)]", "<lnA> of data events");
            grDataLnaNorm->Draw(plotInstr[0].c_str());

	    // Right side markings
            c1->Update();
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(16), "O");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");
            // Draw lines marking the particle types
            l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(16), gPad->GetUxmax(), TMath::Log(16));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_lnA_" + stemp[1] + "_data-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
/*               if(autoFrac)
	          stemp[1] = "frac-auto";
	       else
                  stemp[1] = "frac-" + ToString(fraction, 3);*/
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_lnA_data-only_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
//            grDataLnaNorm->Write(("datalnanorm_" + stemp[1]).c_str());

	    // LnA data/data-norm
            mystyle->SetAxisTitles(grDataLna, "FD energy [log(E/eV)]", "<lnA> of data events");
//	    grDataLna->SetLineStyle(9);
            grDataLna->Draw(plotInstr[0].c_str());
            grDataLnaNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(nrfigs*.03), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grDataLna, "Data (raw)", "lp");
	    legend->AddEntry(grDataLnaNorm, "Data (normalized)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");

	    // Right side markings
            c1->Update();
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(1),  "p");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(4),  "He");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(14), "N");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(16), "O");
            t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(56), "Fe");
            // Draw lines marking the particle types
            l->DrawLine(gPad->GetUxmin(), TMath::Log(1), gPad->GetUxmax(), TMath::Log(1));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(4), gPad->GetUxmax(), TMath::Log(4));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(14), gPad->GetUxmax(), TMath::Log(14));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(16), gPad->GetUxmax(), TMath::Log(16));
            l->DrawLine(gPad->GetUxmin(), TMath::Log(56), gPad->GetUxmax(), TMath::Log(56));

            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_lnA_" + stemp[1] + "_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_lnA_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
//	    grDataLna->SetLineStyle(1);
	 }

         printRes->Close();

         for(int i = 0; i < 3; i++)
         {
            delete[] xbin[i];
            delete[] ybinSig[i];
            delete[] ybinSigNorm[i];
            delete[] ybinBack[i];
            delete[] ybinBackNorm[i];
            delete[] ybinData[i];
            delete[] ybinDataNorm[i];
            delete[] ybinDataLna[i];
            delete[] ybinDataLnaNorm[i];

            if(autoFrac)
            {
               delete[] ylowfractest[i];
               delete[] yhighfractest[i];
            }
         }      

         if(autoFrac)
         {
            delete[] ylowlimitfrac;
            delete[] yhighlimitfrac;
         }
      }         
      else
      {
         nrbins = argc-1; // gives the number of graphs

         RootStyle *mystyle = new RootStyle();
         mystyle->SetBaseStyle();

         TCanvas *c1 = new TCanvas("c1","",1200,900);
         gStyle->SetEndErrorSize(6);

	 TGraphAsymmErrors *ingr;
	 TGraphAsymmErrors *outgr[40];
         TList *tempkeyslist;
	 TLegend *legend;

	 vector<string> treeType;
	 vector<string> treeObs;
	 vector<string> treeNames;
	 vector<int> treePoints;

	 double yrange[2];

	 string names[2];

	 double x[3], y[3];
	 int nrkeys;

	 // Loop over all files to get number and names of observables
	 for(int i = 0; i < nrbins; i++)
	 {
            filename = string(argv[i+1]);
	    cout << "Opening file: " << filename << endl;
	    TFile *infile = TFile::Open(filename.c_str(), "READ");
            tempkeyslist = (TList*)infile->GetListOfKeys();

	    nrkeys = infile->GetNkeys();

	    for(int j = 0; j < nrkeys; j++)
	    {
	       stemp[0] = string(tempkeyslist->At(j)->GetName());
	       stemp[1] = stemp[0];

	       itemp[0] = stemp[1].find("sig_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "sig_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "sig";
	       }

	       itemp[0] = stemp[1].find("back_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "back_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "back";
	       }

	       itemp[0] = stemp[1].find("data_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "data_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "data";
	       }

	       itemp[0] = stemp[1].find("datanorm_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "datanorm_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "datanorm";
	       }

	       treeType.push_back(names[0]);
	       treeObs.push_back(names[1]);
	       treeNames.push_back(stemp[0]);
	    }

	    infile->Close();
	 }

	 cout << nrbins << endl;
	 cout << nrkeys << endl;

	 for(int i = 0; i < treeNames.size(); i++)
            cout << i << ": " << treeNames[i] << endl;
	 cout << endl;

	 for(int i = 0; i < treeType.size(); i++)
            cout << i << ": " << treeType[i] << endl;
	 cout << endl;

	 for(int i = 0; i < treeObs.size(); i++)
            cout << i << ": " << treeObs[i] << endl;
	 cout << endl;

         // y-range settings
         cerr << endl << "Please set the y-axis range for all plots (comma separated): ";
         cin >> yrange[0] >> *ctemp >> yrange[1];

	 // Loop over all types (sig, back, data, normalized data) - plotting
	 for(int j = 0; j < nrkeys; j++)
	 {
	    cout << "Using tree: " << treeNames[j] << endl;
	    cout << "Type = " << treeType[j] << ", Observable = " << treeObs[j] << endl;

	    for(int i = 0; i < nrbins; i++)
               outgr[i] = new TGraphAsymmErrors();

	    // Loop over all trees
	    for(int i = 0; i < nrbins; i++)
	    {
               filename = string(argv[i+1]);
	       cout << "Opening file: " << filename << endl;
	       TFile *infile = TFile::Open(filename.c_str(), "READ");
               tempkeyslist = (TList*)infile->GetListOfKeys();

	       cout << "Tree name: " << treeNames[j+i*nrkeys] << endl;
               ingr = (TGraphAsymmErrors*)infile->Get(treeNames[j+i*nrkeys].c_str());
	       outgr[i]->Set(ingr->GetN());

	       for(int k = 0; k < ingr->GetN(); k++)
	       {
                  ingr->GetPoint(k, x[0], y[0]);
		  x[1] = ingr->GetErrorXlow(k);
		  x[2] = ingr->GetErrorXhigh(k);
		  y[1] = ingr->GetErrorYlow(k);
		  y[2] = ingr->GetErrorYhigh(k);

                  cout << k << ": " << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;

                  outgr[i]->SetPoint(k, x[0], y[0]);
                  outgr[i]->SetPointError(k, x[1], x[2], y[1], y[2]);
	       }
cout << "3" << endl;

	       SetColor(outgr[i], i, nrbins);

	       if(i == 0)
	       {
		  if(treeType[j] == "sig")
	             mystyle->SetAxisTitles(outgr[i], "FD energy [log(E/eV])", "Purity of signal simulation events");
		  else if(treeType[j] == "back")
	             mystyle->SetAxisTitles(outgr[i], "FD energy [log(E/eV])", "Purity of background simulation events");
		  else if(treeType[j] == "data")
	             mystyle->SetAxisTitles(outgr[i], "FD energy [log(E/eV])", "Signal fraction of data events");
		  else if(treeType[j] == "datanorm")
	             mystyle->SetAxisTitles(outgr[i], "FD energy [log(E/eV])", "Signal fraction of data events (normalized)");
                  outgr[i]->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);
                  outgr[i]->Draw("AP");

                  legend = new TLegend(gPad->GetLeftMargin()+0.05, gPad->GetBottomMargin()+0.05, gPad->GetLeftMargin()+.4, gPad->GetBottomMargin()+0.15);
	       }
	       else
                  outgr[i]->Draw("P; SAME");

               legend->AddEntry(outgr[i], (treeObs[j+i*nrkeys]).c_str(), "pl");
	    }

            legend->Draw("SAME");

            stemp[0] = RemoveFilename(&filename) + "/combined_fraction_plot_" + treeType[j] + ".pdf";
            c1->SaveAs(stemp[0].c_str());

	    for(int i = 0; i < nrbins; i++)
               delete outgr[i];
	 }

/*	 // Loop over all files
	 for(int i = 0; i < nrbins; i++)
	 {
            filename = string(argv[i+1]);
	    cout << "Opening file: " << filename << endl;
	    TFile *infile = TFile::Open(filename.c_str(), "READ");
            tempkeyslist = (TList*)infile->GetListOfKeys();

	    // Loop over all types (sig, back, data, normalized data) - number of trees
	    for(int j = 0; j < infile->GetNkeys(); j++)
	    {
	       stemp[0] = string(tempkeyslist->At(j)->GetName());
	       stemp[1] = stemp[0];
	       cout << "Using tree: " << stemp[0] << endl;

	       itemp[0] = stemp[1].find("sig_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "sig_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "sig";
	       }

	       itemp[0] = stemp[1].find("back_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "back_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "back";
	       }

	       itemp[0] = stemp[1].find("data_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "data_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "data";
	       }

	       itemp[0] = stemp[1].find("datanorm_");
	       if(itemp[0] != string::npos)
	       {
	          names[0] = "datanorm_";
		  names[1] = stemp[1].erase(0, names[0].length());
		  names[0] = "datanorm";
	       }

	       cout << "Type = " << names[0] << ", Observable = " << names[1] << endl;

               ingr = (TGraphAsymmErrors*)infile->Get(stemp[0].c_str());
	       outgr[j]->Set(ingr->GetN());

	       for(int k = 0; k < ingr->GetN(); k++)
	       {
                  ingr->GetPoint(i, x[0], y[0]);
		  x[1] = ingr->GetErrorXlow(i);
		  x[2] = ingr->GetErrorXhigh(i);
		  y[1] = ingr->GetErrorYlow(i);
		  y[2] = ingr->GetErrorYhigh(i);

                  cout << i << ": " << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;

                  outgr[j]->SetPoint(i, x[0], y[0]);
                  outgr[j]->SetPointError(i, x[1], x[2], y[1], y[2]);
	       }
	    }

	    if(infile->GetNkeys() < 4)
	    {
	       // Combine signal plots
	       c1->cd();
	       outgr[0]->SetMarkerColorAlpha(2, 0.75);
	       outgr[0]->SetMarkerStyle(20);
	       outgr[0]->SetMarkerSize(0.9);
	       outgr[0]->SetLineColorAlpha(2, 0.75);
	       outgr[0]->SetLineWidth(2);
	       mystyle->SetAxisTitles(outgr[0], "FD energy [log(E/eV])", "Purity of signal simulation events");
	       if(i == 0)
                  outgr[0]->Draw("AP");
	       else
                  outgr[0]->Draw("SAME");

	       // Combine background plots
	       c2->cd();
	       outgr[1]->SetMarkerColorAlpha(2, 0.75);
	       outgr[1]->SetMarkerStyle(20);
	       outgr[1]->SetMarkerSize(0.9);
	       outgr[1]->SetLineColorAlpha(2, 0.75);
	       outgr[1]->SetLineWidth(2);
	       mystyle->SetAxisTitles(outgr[1], "FD energy [log(E/eV])", "Purity of background simulation events");
	       if(i == 0)
                  outgr[1]->Draw("AP");
	       else
                  outgr[1]->Draw("SAME");

	       // Combine data plots
	       c3->cd();
	       outgr[2]->SetMarkerColorAlpha(2, 0.75);
	       outgr[2]->SetMarkerStyle(20);
	       outgr[2]->SetMarkerSize(0.9);
	       outgr[2]->SetLineColorAlpha(2, 0.75);
	       outgr[2]->SetLineWidth(2);
	       mystyle->SetAxisTitles(outgr[1], "FD energy [log(E/eV])", "Signal fraction of data events");
	       if(i == 0)
                  outgr[2]->Draw("AP");
	       else
                  outgr[2]->Draw("SAME");
	    }

	    if(infile->GetNkeys() == 4)
	    {
	       // Combine normalized data plots
	       c4->cd();
	       outgr[3]->SetMarkerColorAlpha(2, 0.75);
	       outgr[3]->SetMarkerStyle(20);
	       outgr[3]->SetMarkerSize(0.9);
	       outgr[3]->SetLineColorAlpha(2, 0.75);
	       outgr[3]->SetLineWidth(2);
	       mystyle->SetAxisTitles(outgr[1], "FD energy [log(E/eV])", "Signal fraction of data events (normalized)");
	       if(i == 0)
                  outgr[3]->Draw("AP");
	       else
                  outgr[3]->Draw("SAME");
	    }
	 }*/
      }

      // Removing allocated variables
      delete[] ftemp;
      delete ctemp;
      delete[] stemp;
      delete[] dtemp;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (application_results.txt, individual_results.txt or root files created by previous scripts)." << endl;
      return 1;
   }

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
