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

void NormalizeFrac(int nr, float **sig, float **back, float **data, double *dtemp)
{
   cerr << nr << ", Mean: x_sig = " << sig[0][nr] << ", x_back = "  << back[0][nr] << ", x_d = " << data[0][nr] << endl;
   cerr << nr << ", Neg:  x_sig = " << sig[1][nr] << ", x_back = "  << back[1][nr] << ", x_d = " << data[1][nr] << endl;
   cerr << nr << ", Pos:  x_sig = " << sig[2][nr] << ", x_back = "  << back[2][nr] << ", x_d = " << data[2][nr] << endl;

//   dtemp[0] = 1. - back[0][nr];
   // calculate the normalized fraction
   dtemp[0] = (data[0][nr] - back[0][nr])/(sig[0][nr] - back[0][nr]);
   cout << "Step 1 = " << dtemp[0] << endl;

   if( (sig != data) && (back != data) )
   {
      // calculate neg propagation error with respect to x_d (d(x_d,norm)/d(x_d))*delta(x_d)
      dtemp[1] = TMath::Power(data[1][nr]/(sig[0][nr] - back[0][nr]), 2);
      cout << "Step 2 = " << TMath::Power(data[1][nr]/(sig[0][nr] - back[0][nr]), 2) << endl;
      // calculate neg propagation error with respect to x_sig (d(x_d,norm)/d(x_sig))*delta(x_sig)
      dtemp[1] += TMath::Power(sig[1][nr]*(data[0][nr] - back[0][nr])/TMath::Power(sig[0][nr] - back[0][nr], 2), 2);
      cout << "Step 3 = " << TMath::Power(sig[1][nr]*(data[0][nr] - back[0][nr])/TMath::Power(sig[0][nr] - back[0][nr], 2), 2) << endl;
      // calculate neg propagation error with respect to x_back (d(x_d,norm)/d(x_back))*delta(x_back)
      dtemp[1] += TMath::Power(back[1][nr]*(data[0][nr] - sig[0][nr])/TMath::Power(sig[0][nr] - back[0][nr], 2), 2);
      cout << "Step 4 = " << TMath::Power(back[1][nr]*(data[0][nr] - sig[0][nr])/TMath::Power(sig[0][nr] - back[0][nr], 2), 2) << endl;
      dtemp[1] = TMath::Sqrt(dtemp[1]);
      cout << "Step 5 = " << TMath::Sqrt(dtemp[1]) << endl;
      // calculate pos propagation error with respect to x_d (d(x_d,norm)/d(x_d))*delta(x_d)
      dtemp[2] = TMath::Power(data[2][nr]/(sig[0][nr] - back[0][nr]), 2);
      // calculate pos propagation error with respect to x_sig (d(x_d,norm)/d(x_sig))*delta(x_sig)
      dtemp[2] += TMath::Power(sig[2][nr]*(data[0][nr] - back[0][nr])/TMath::Power(sig[0][nr] - back[0][nr], 2), 2);
      // calculate pos propagation error with respect to x_back (d(x_d,norm)/d(x_back))*delta(x_back)
      dtemp[2] += TMath::Power(back[2][nr]*(data[0][nr] - sig[0][nr])/TMath::Power(sig[0][nr] - back[0][nr], 2), 2);
      dtemp[2] = TMath::Sqrt(dtemp[2]);
   }
   else
   {
      cout << "These two things are the same" << endl;
      dtemp[1] = 0.;
      dtemp[2] = 0.;
   }
   
   cout << nr << ": Normfrac = " << dtemp[0] << ", Normfrac neg error = " << dtemp[1] << ", Normfrac pos error = " << dtemp[2] << endl;
}

// Read published lnA results to add them to the plot (type: 0 = EPOS, 1 = QGSJETII, 2 = SIBYLL)
int ReadLnaResults(vector<float> *val, int type)
{
   ifstream infile;
   char ctemp[1024];
   float *ftemp;
   int nrp = 0;

   string stemp;
   if(type == 0)
      stemp = string(rootdir) + "/input/lnA_moments_epos.txt";
   else if(type == 1)
      stemp = string(rootdir) + "/input/lnA_moments_qgs.txt";
   else if(type == 2)
      stemp = string(rootdir) + "/input/lnA_moments_sib.txt";
   else
      return -1;

   ftemp = new float[10];

   infile.open(stemp.c_str(), ifstream::in);

   if(infile.is_open())
   {
      infile.getline(ctemp, 1024, '\n');
      
      while(1)
      {
	 for(int i = 0; i < 10; i++)
            infile >> ftemp[i];

	 if(ftemp[0] > 18.40)
	 {
            val->push_back(ftemp[0]);
            val->push_back(ftemp[1]);
            val->push_back(ftemp[2]);

	    nrp++;
	 }

	 infile.ignore(1,' ');
	 if(infile.eof()) break;
      }
   }

   infile.close();

   delete[] ftemp;

   return nrp;
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
   bool simplNorm;
   bool xErrors;
   string plotInstr[2];

   int nrbins = -1;

   ResultRead *analRes = new ResultRead();

   string filename;
   if(argc > 1)
   {
      ftemp = new float[6];
      ctemp = new char;
      stemp = new string[3];
      itemp = new int[2];
      dtemp = new double[3];

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

	 if(fraction == -1)
	 {
            cout << "Use a simplified normalization approach, without knowing signal fraction (1 = yes, 0 = no): ";
            cin >> itemp[0];
	    while( (itemp[0] != 0) && (itemp[0] != 1) )
	    {
               cout << "Use a simplified normalization approach, without knowing signal fraction (1 = yes, 0 = no): ";
               cin >> itemp[0];
	    }
	    if(itemp[0] == 1)
               simplNorm = true;
	    else if(itemp[0] == 0)
               simplNorm = false;
	    else
               return 1;
	 }
	 else
            simplNorm = false;

         nrbins = argc-1;
         cout << "Number of bins on the plot = " << nrbins << endl;

         // Arrays for plotting
         float *xbin[3];
         float *ybinSig[3];
         float *ybinSigNorm[3];
         float *ybinBack[3];
	 float *ybinBackInv[3];
         float *ybinBackNorm[3];
         float *ybinData[3];
         float *ybinDataNorm[3];
         float *ybinDataLna[3];
         float *ybinDataLnaNorm[3];
	 float *ylowfractest[3], *yhighfractest[3];
	 float *ylowlimitfrac, *yhighlimitfrac;
	 float *ymvaCut[3];

	 vector<string> treeName;
	 vector<int> treeId;

	 double sigType;
	 double backType;

	 float *xlimit;
	 xlimit = new float[2];

         for(int i = 0; i < 3; i++)
         {
            xbin[i] = new float[nrbins];
            ybinSig[i] = new float[nrbins];
            ybinSigNorm[i] = new float[nrbins];
            ybinBack[i] = new float[nrbins];
            ybinBackInv[i] = new float[nrbins];
            ybinBackNorm[i] = new float[nrbins];
            ybinData[i] = new float[nrbins];
            ybinDataNorm[i] = new float[nrbins];
            ybinDataLna[i] = new float[nrbins];
            ybinDataLnaNorm[i] = new float[nrbins];
	    ymvaCut[i] = new float[nrbins];

	    if(autoFrac && !simplNorm)
	    {
	       ylowfractest[i] = new float[nrbins];
	       yhighfractest[i] = new float[nrbins];
	    }
         }

	 if(autoFrac && !simplNorm)
         {
            ylowlimitfrac = new float[nrbins];
            yhighlimitfrac = new float[nrbins];
	 }

         // Going through all input files
         for(int i = 0; i < nrbins; i++)
         {
            filename = string(argv[i+1]);
            analRes->ReadFile(filename);

	    if(i == 0)
	    {
               for(int j = 0; j < analRes->GetNrTrees(0); j++)
	       {
	          // Get tree name and type
                  if(analRes->GetTreeType(j) > 0)
		  {
                     treeName.push_back(analRes->GetTreeName(j));
	             treeId.push_back(analRes->GetTreeType(j));
		  }
	       }
	    }

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

	    if(i == 0)
               xlimit[0] = analRes->GetLowEnergy();

	    if(i == nrbins-1)
               xlimit[1] = analRes->GetHighEnergy();

	    // Getting MVA cut values
	    ymvaCut[0][i] = analRes->GetMvaCut(0);
	    ymvaCut[1][i] = TMath::Abs(analRes->GetMvaCut(-1) - ymvaCut[0][i]);
	    ymvaCut[2][i] = TMath::Abs(analRes->GetMvaCut(1) - ymvaCut[0][i]);

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

            // Inverted background values (get fraction instead of purity)
            ybinBackInv[0][i] = 1. - ybinBack[0][i];
            ybinBackInv[1][i] = ybinBack[1][i];
            ybinBackInv[2][i] = ybinBack[2][i];

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
               if(simplNorm)
	       {
/*                  // x_d,norm = (x_d - x_back)/(x_sig - x_back)
                  ybinSigNorm[0][i] = analRes->GetFraction(1, -1);
                  analRes->GetFractionError(ftemp);
                  ybinSigNorm[1][i] = ftemp[0];
                  ybinSigNorm[2][i] = ftemp[1];

                  ybinBackNorm[0][i] = analRes->GetFraction(0, -1);
                  analRes->GetFractionError(ftemp);
                  ybinBackNorm[1][i] = ftemp[0];
                  ybinBackNorm[2][i] = ftemp[1];

                  // Getting data ybin values (normalized)
                  ybinDataNorm[0][i] = analRes->GetFraction(2, -1);
                  analRes->GetFractionError(ftemp);
                  ybinDataNorm[1][i] = ftemp[0];
                  ybinDataNorm[2][i] = ftemp[1];*/

		  NormalizeFrac(i, ybinSig, ybinBackInv, ybinSig, dtemp);

		  ybinSigNorm[0][i] = dtemp[0];
                  ybinSigNorm[1][i] = dtemp[1];
                  ybinSigNorm[2][i] = dtemp[2];

		  NormalizeFrac(i, ybinSig, ybinBackInv, ybinBackInv, dtemp);

		  ybinBackNorm[0][i] = dtemp[0];
                  ybinBackNorm[1][i] = dtemp[1];
                  ybinBackNorm[2][i] = dtemp[2];

		  NormalizeFrac(i, ybinSig, ybinBackInv, ybinData, dtemp);

		  ybinDataNorm[0][i] = dtemp[0];
                  ybinDataNorm[1][i] = dtemp[1];
                  ybinDataNorm[2][i] = dtemp[2];
/*
	          ftemp[0] = 1. - ybinBackNorm[0][i];
                  // calculate the normalized fraction
	          ftemp[1] = (ybinDataNorm[0][i] - ftemp[0])/(ybinSigNorm[0][i] - ftemp[0]);
	          // calculate neg propagation error with respect to x_d (d(x_d,norm)/d(x_d))*delta(x_d)
	          ftemp[2] = TMath::Power(ybinDataNorm[1][i]/(ybinSigNorm[0][i] - ftemp[0]), 2);
	          // calculate neg propagation error with respect to x_sig (d(x_d,norm)/d(x_sig))*delta(x_sig)
	          ftemp[2] += TMath::Power(ybinSigNorm[1][i]*(ybinDataNorm[0][i] - ftemp[0])/TMath::Power(ybinSigNorm[0][i] - ftemp[0], 2), 2);
	          // calculate neg propagation error with respect to x_back (d(x_d,norm)/d(x_back))*delta(x_back)
	          ftemp[2] += TMath::Power(ybinBackNorm[1][i]*(ybinDataNorm[0][i] - ybinSigNorm[0][i])/TMath::Power(ybinSigNorm[0][i] - ftemp[0], 2), 2);
	          ftemp[2] = TMath::Sqrt(ftemp[2]);
	          // calculate pos propagation error with respect to x_d (d(x_d,norm)/d(x_d))*delta(x_d)
	          ftemp[3] = TMath::Power(ybinDataNorm[2][i]/(ybinSigNorm[0][i] - ftemp[0]), 2);
	          // calculate pos propagation error with respect to x_sig (d(x_d,norm)/d(x_sig))*delta(x_sig)
	          ftemp[3] += TMath::Power(ybinSigNorm[2][i]*(ybinDataNorm[0][i] - ftemp[0])/TMath::Power(ybinSigNorm[0][i] - ftemp[0], 2), 2);
	          // calculate pos propagation error with respect to x_back (d(x_d,norm)/d(x_back))*delta(x_back)
	          ftemp[3] += TMath::Power(ybinBackNorm[2][i]*(ybinDataNorm[0][i] - ybinSig[0][i])/TMath::Power(ybinSigNorm[0][i] - ftemp[0], 2), 2);
	          ftemp[3] = TMath::Sqrt(ftemp[3]);

	          cout << i << ": Back = " << ftemp[0] << ", Normfrac = " << ftemp[1] << ", Normfrac neg error = " << ftemp[2] << ", Normfrac pos error = " << ftemp[3] << endl;
*/
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
         }

	 // Set normalization for data
	 if(autoFrac && !simplNorm)
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
	    if(!simplNorm)
	    {
               // Raw data (mean + propagation of errors)
               ybinDataLna[0][i] = TMath::Log(56)*(ybinSig[0][i] - ybinData[0][i])/(ybinSig[0][i] - (1.-ybinBack[0][i]));
	       ftemp[0] = ybinSig[1][i] + ybinData[1][i];
	       ftemp[1] = ybinSig[1][i] + ybinBack[1][i];
	       ybinDataLna[1][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);
	       ftemp[0] = ybinSig[2][i] + ybinData[2][i];
	       ftemp[1] = ybinSig[2][i] + ybinBack[2][i];
	       ybinDataLna[2][i] = (TMath::Log(56)*ftemp[0]*ftemp[1])/(ftemp[0] + ftemp[1]);

//	       cerr << i << ": Raw = " << ybinDataLna[0][i] << "\t"  << ybinDataLna[1][i] << "\t" << ybinDataLna[2][i] << endl;

	       // Normalized data (mean + propagation of errors)
               ybinDataLnaNorm[0][i] = TMath::Log(56)*(1. - ybinDataNorm[0][i]);
               ybinDataLnaNorm[1][i] = TMath::Log(56)*(1. - (ybinDataNorm[0][i]-ybinDataNorm[1][i]));
               ybinDataLnaNorm[1][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[1][i]);
               ybinDataLnaNorm[2][i] = TMath::Log(56)*(1. - (ybinDataNorm[0][i]+ybinDataNorm[2][i]));
               ybinDataLnaNorm[2][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[2][i]);

//	       cerr << i << ": Norm = " << ybinDataLnaNorm[0][i] << "\t"  << ybinDataLnaNorm[1][i] << "\t" << ybinDataLnaNorm[2][i] << endl;
            }
            // New calculation (Only for a single background tree)
	    else
	    {
	       if(i == 0)
	       {
	          for(int j = 0; j < treeName.size(); j++)
	          {
	             cout << "Tree name = " << treeName[j] << endl;
	             cout << "Tree ID = " << treeId[j] << endl;

	             itemp[0] = treeName[j].find("Proton");
	             if(itemp[0] != string::npos)
                        dtemp[0] = TMath::Log(1);

	             itemp[0] = treeName[j].find("Helium");
	             if(itemp[0] != string::npos)
                        dtemp[0] = TMath::Log(4);

	             itemp[0] = treeName[j].find("Oxygen");
	             if(itemp[0] != string::npos)
                        dtemp[0] = TMath::Log(16);

	             itemp[0] = treeName[j].find("Iron");
	             if(itemp[0] != string::npos)
                        dtemp[0] = TMath::Log(56);
	             
	             if(treeId[j] == 1)
                        sigType = dtemp[0];
	             else if(treeId[j] == 2)
                        backType = dtemp[0];
	          }
	       }

               // Raw data (mean + errors): lnA = (x_d,sig)*lnA_sig + (1-x_d,sig)*lnA_back
	       ybinDataLna[0][i] = ybinData[0][i]*sigType + (1-ybinData[0][i])*backType;
	       ybinDataLna[1][i] = (ybinData[0][i]-ybinData[1][i])*sigType + (1-(ybinData[0][i]-ybinData[1][i]))*backType;
               ybinDataLna[1][i] = TMath::Abs(ybinDataLna[0][i] - ybinDataLna[1][i]);
	       ybinDataLna[2][i] = (ybinData[0][i]-ybinData[2][i])*sigType + (1-(ybinData[0][i]-ybinData[2][i]))*backType;
               ybinDataLna[2][i] = TMath::Abs(ybinDataLna[0][i] - ybinDataLna[2][i]);

	       cerr << i << ": Raw = " << ybinDataLna[0][i] << "\t"  << ybinDataLna[1][i] << "\t" << ybinDataLna[2][i] << endl;

	       // Normalized data (mean + propagation of errors): x_d,norm = (x_d - x_back)/(x_sig - x_back)
/*	       // change background purity to fraction
	       ftemp[0] = 1. - ybinBack[0][i];
               // calculate the normalized fraction
	       ftemp[1] = (ybinData[0][i] - ftemp[0])/(ybinSig[0][i] - ftemp[0]);
	       // calculate neg propagation error with respect to x_d (d(x_d,norm)/d(x_d))*delta(x_d)
	       ftemp[2] = TMath::Power(ybinData[1][i]/(ybinSig[0][i] - ftemp[0]), 2);
	       // calculate neg propagation error with respect to x_sig (d(x_d,norm)/d(x_sig))*delta(x_sig)
	       ftemp[2] += TMath::Power(ybinSig[1][i]*(ybinData[0][i] - ftemp[0])/TMath::Power(ybinSig[0][i] - ftemp[0], 2), 2);
	       // calculate neg propagation error with respect to x_back (d(x_d,norm)/d(x_back))*delta(x_back)
	       ftemp[2] += TMath::Power(ybinBack[1][i]*(ybinData[0][i] - ybinSig[0][i])/TMath::Power(ybinSig[0][i] - ftemp[0], 2), 2);
	       ftemp[2] = TMath::Sqrt(ftemp[2]);
	       // calculate pos propagation error with respect to x_d (d(x_d,norm)/d(x_d))*delta(x_d)
	       ftemp[3] = TMath::Power(ybinData[2][i]/(ybinSig[0][i] - ftemp[0]), 2);
	       // calculate pos propagation error with respect to x_sig (d(x_d,norm)/d(x_sig))*delta(x_sig)
	       ftemp[3] += TMath::Power(ybinSig[2][i]*(ybinData[0][i] - ftemp[0])/TMath::Power(ybinSig[0][i] - ftemp[0], 2), 2);
	       // calculate pos propagation error with respect to x_back (d(x_d,norm)/d(x_back))*delta(x_back)
	       ftemp[3] += TMath::Power(ybinBack[2][i]*(ybinData[0][i] - ybinSig[0][i])/TMath::Power(ybinSig[0][i] - ftemp[0], 2), 2);
	       ftemp[3] = TMath::Sqrt(ftemp[3]);

	       cout << i << ": Back = " << ftemp[0] << ", Normfrac = " << ftemp[1] << ", Normfrac neg error = " << ftemp[2] << ", Normfrac pos error = " << ftemp[3] << endl;*/

	       ybinDataLnaNorm[0][i] = ybinDataNorm[0][i]*sigType + (1-ybinDataNorm[0][i])*backType;
	       ybinDataLnaNorm[1][i] = (ybinDataNorm[0][i]-ybinDataNorm[1][i])*sigType + (1-(ybinDataNorm[0][i]-ybinDataNorm[1][i]))*backType;
               ybinDataLnaNorm[1][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[1][i]);
	       ybinDataLnaNorm[2][i] = (ybinDataNorm[0][i]+ybinDataNorm[2][i])*sigType + (1-(ybinDataNorm[0][i]+ybinDataNorm[2][i]))*backType;
               ybinDataLnaNorm[2][i] = TMath::Abs(ybinDataLnaNorm[0][i] - ybinDataLnaNorm[2][i]);

	       cerr << i << ": Norm = " << ybinDataLnaNorm[0][i] << "\t"  << ybinDataLnaNorm[1][i] << "\t" << ybinDataLnaNorm[2][i] << endl;
            }
	 }

	 cout << endl << "Point printout:" << endl;
	 ftemp[0] = 0;
	 ftemp[1] = 0;
	 ftemp[2] = 0;
	 ftemp[3] = 0;
	 ftemp[4] = 0;
	 ftemp[5] = 0;
/*	 // Skip the first point, just because
	 for(int i = 1; i < nrbins; i++)
	 {
            cerr << i << ": MVA cut    = " << ymvaCut[0][i] << "\t" << ymvaCut[1][i] << "\t" << ymvaCut[2][i] << endl;
            cerr << i << ": Signal     = " << ybinSig[0][i] << "\t" << ybinSig[1][i] << "\t" << ybinSig[2][i] << endl;
            cerr << i << ": Background = " << ybinBack[0][i] << "\t" << ybinBack[1][i] << "\t" << ybinBack[2][i] << endl;
            cerr << i << ": Raw data   = " << ybinData[0][i] << "\t" << ybinData[1][i] << "\t" << ybinData[2][i] << endl;
            cerr << i << ": Norm. data = " << ybinDataNorm[0][i] << "\t" << ybinDataNorm[1][i] << "\t" << ybinDataNorm[2][i] << endl;
	    ftemp[0] += ybinDataNorm[0][i];
	    ftemp[1] += ybinDataNorm[1][i];
	    ftemp[2] += ybinDataNorm[2][i];
	    ftemp[3] += ybinData[0][i];
	    ftemp[4] += ybinData[1][i];
	    ftemp[5] += ybinData[2][i];
	 }
	 cout << endl << "Mean data value = " << ftemp[3]/(nrbins-1.) << "\t" << ftemp[4]/(nrbins-1.) << "\t" << ftemp[5]/(nrbins-1.) << endl;
	 cout << endl << "Mean norm. data value = " << ftemp[0]/(nrbins-1.) << "\t" << ftemp[1]/(nrbins-1.) << "\t" << ftemp[2]/(nrbins-1.) << endl;
	 cout << endl;*/

         // Create directory structure for plots and delete old plots
         stemp[2] = RemoveFilename(&filename);
         itemp[0] = analRes->GetFileType();
         if(itemp[0] == 0)
            stemp[0] = RemoveFilename(&stemp[2]);
         else if(itemp[0] == 1)
            stemp[0] = stemp[2];

         stemp[1] = "mkdir -p " + RemoveFilename(&stemp[0]) + "/plots";
         system(stemp[1].c_str());

	 cout << "Is this a p/Fe (0), He/Fe (1) or O/Fe (2) treatment: ";
	 cin >> itemp[1];

         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();

	    if(itemp[1] == 0)
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fractions_" + stemp[1] + "_p-Fe_summary.root";
	    else if(itemp[1] == 1)
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fractions_" + stemp[1] + "_He-Fe_summary.root";
	    else if(itemp[1] == 2)
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fractions_" + stemp[1] + "_O-Fe_summary.root";
         }
         else if(itemp[0] == 1)
         {
            stemp[1] = analRes->GetObservableType();

	    if(itemp[1] == 0)
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fractions_" + stemp[1] + "_p-Fe_summary.root";
	    else if(itemp[1] == 1)
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fractions_" + stemp[1] + "_He-Fe_summary.root";
	    else if(itemp[1] == 2)
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fractions_" + stemp[1] + "_O-Fe_summary.root";
         }

         TFile *fractionRes = TFile::Open(stemp[2].c_str(), "RECREATE");;
         TTree *fracTree;
	 if(itemp[1] == 0)
            fracTree = new TTree("fracTree_p-Fe_treatment", "Fractions from MVA analysis");
	 else if(itemp[1] == 1)
            fracTree = new TTree("fracTree_He-Fe_treatment", "Fractions from MVA analysis");
	 else if(itemp[1] == 2)
            fracTree = new TTree("fracTree_O-Fe_treatment", "Fractions from MVA analysis");

	 vector<double> *rawFrac = new vector<double>;
	 vector<double> *normFrac = new vector<double>;
	 vector<double> *meanFrac = new vector<double>;

	 fracTree->Branch("rawFraction", rawFrac);
	 fracTree->Branch("normFraction", normFrac);
	 fracTree->Branch("meanFraction", meanFrac);

	 for(int i = 0; i < nrbins; i++)
	 {
            cerr << i << ": MVA cut    = " << ymvaCut[0][i] << "\t" << ymvaCut[1][i] << "\t" << ymvaCut[2][i] << endl;
            cerr << i << ": Signal     = " << ybinSig[0][i] << "\t" << ybinSig[1][i] << "\t" << ybinSig[2][i] << endl;
            cerr << i << ": Background = " << ybinBack[0][i] << "\t" << ybinBack[1][i] << "\t" << ybinBack[2][i] << endl;
            cerr << i << ": Raw data   = " << ybinData[0][i] << "\t" << ybinData[1][i] << "\t" << ybinData[2][i] << endl;
            cerr << i << ": Norm. data = " << ybinDataNorm[0][i] << "\t" << ybinDataNorm[1][i] << "\t" << ybinDataNorm[2][i] << endl;
	    ftemp[0] += ybinDataNorm[0][i];
	    ftemp[1] += ybinDataNorm[1][i];
	    ftemp[2] += ybinDataNorm[2][i];
	    ftemp[3] += ybinData[0][i];
	    ftemp[4] += ybinData[1][i];
	    ftemp[5] += ybinData[2][i];

	    rawFrac->push_back(ybinData[0][i]);
	    rawFrac->push_back(ybinData[1][i]);
	    rawFrac->push_back(ybinData[2][i]);
	    normFrac->push_back(ybinDataNorm[0][i]);
	    normFrac->push_back(ybinDataNorm[1][i]);
	    normFrac->push_back(ybinDataNorm[2][i]);
	 }
	 cout << endl << "Mean data value = " << ftemp[3]/(nrbins) << "\t" << ftemp[4]/(nrbins) << "\t" << ftemp[5]/(nrbins) << endl;
	 cout << endl << "Mean norm. data value = " << ftemp[0]/(nrbins) << "\t" << ftemp[1]/(nrbins) << "\t" << ftemp[2]/(nrbins) << endl;
	 cout << endl;

	 meanFrac->push_back(ftemp[3]/((double)nrbins));
	 meanFrac->push_back(ftemp[4]/((double)nrbins));
	 meanFrac->push_back(ftemp[5]/((double)nrbins));
	 meanFrac->push_back(ftemp[0]/((double)nrbins));
	 meanFrac->push_back(ftemp[1]/((double)nrbins));
	 meanFrac->push_back(ftemp[2]/((double)nrbins));

	 fracTree->Fill();

	 fractionRes->Write();
	 fractionRes->Close();
	 delete fractionRes;
	 delete rawFrac;
	 delete normFrac;
	 delete meanFrac;

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
         cerr << endl << "Please set the y-axis range for all fraction plots (comma separated): ";
         cin >> yrange[0] >> *ctemp >> yrange[1];

/*         // Create directory structure for plots and delete old plots
         stemp[2] = RemoveFilename(&filename);
         itemp[0] = analRes->GetFileType();
         if(itemp[0] == 0)
            stemp[0] = RemoveFilename(&stemp[2]);
         else if(itemp[0] == 1)
            stemp[0] = stemp[2];

         stemp[1] = "mkdir -p " + RemoveFilename(&stemp[0]) + "/plots";
         system(stemp[1].c_str());*/

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
         // MVA cut
         TGraphAsymmErrors *grMva = new TGraphAsymmErrors(nrbins, xbin[0], ymvaCut[0], xbin[1], xbin[2], ymvaCut[1], ymvaCut[2]);
         mystyle->SetGraphColor(grMva, 1);
         grMva->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
         grMva->GetYaxis()->SetRangeUser(0., 1.);

         // Signal
         TGraphAsymmErrors *grSig = new TGraphAsymmErrors(nrbins, xbin[0], ybinSig[0], xbin[1], xbin[2], ybinSig[1], ybinSig[2]);
         mystyle->SetGraphColor(grSig, 1);
         grSig->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
         grSig->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

         // Background
         TGraphAsymmErrors *grBack = new TGraphAsymmErrors(nrbins, xbin[0], ybinBack[0], xbin[1], xbin[2], ybinBack[1], ybinBack[2]);
         mystyle->SetGraphColor(grBack, 0);
         grBack->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
         grBack->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

         // Data
         TGraphAsymmErrors *grData = new TGraphAsymmErrors(nrbins, xbin[0], ybinData[0], xbin[1], xbin[2], ybinData[1], ybinData[2]);
         mystyle->SetGraphColor(grData, 2);
         grData->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
         grData->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

         TGraphAsymmErrors *grSigNorm, *grBackNorm, *grDataNorm;
         // Normalized values
         if(runnorm)
         {
            grSigNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinSigNorm[0], xbin[1], xbin[2], ybinSigNorm[1], ybinSigNorm[2]);
            mystyle->SetGraphColor(grSigNorm, 1);
            grSigNorm->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
            grSigNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

            grBackNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinBackNorm[0], xbin[1], xbin[2], ybinBackNorm[1], ybinBackNorm[2]);
            mystyle->SetGraphColor(grBackNorm, 0);
            grBackNorm->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
            grBackNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);

            grDataNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataNorm[0], xbin[1], xbin[2], ybinDataNorm[1], ybinDataNorm[2]);
            mystyle->SetGraphColor(grDataNorm, 3);
            grDataNorm->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
            grDataNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);
         }

	 // LnA graphs for data, normalized data and published data
	 TGraphAsymmErrors *grDataLna, *grDataLnaNorm, *grPubLna;
         grDataLna = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataLna[0], xbin[1], xbin[2], ybinDataLna[1], ybinDataLna[2]);
         mystyle->SetGraphColor(grDataLna, 2);
         grDataLna->GetXaxis()->SetRangeUser(xlimit[0], xlimit[1]);
         grDataLna->GetYaxis()->SetRangeUser(-0.7, 5.);
         if(runnorm)
         {
            grDataLnaNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataLnaNorm[0], xbin[1], xbin[2], ybinDataLnaNorm[1], ybinDataLnaNorm[2]);
            mystyle->SetGraphColor(grDataLnaNorm, 3);
            grDataLnaNorm->GetYaxis()->SetRangeUser(-0.7, 5.);
         }

	 cout << "Which dataset did you use (EPOS = 0, QGSJET = 1, SIBYLL = 2)? ";
	 cin >> itemp[1];

	 vector<float> returnVal;
	 itemp[1] = ReadLnaResults(&returnVal, itemp[1]);
	 cout << "Number of points = " << itemp[1] << endl;
	 float *xbinPub[3];
	 float *ybinPubLna[3];

	 for(int i = 0; i < 3; i++)
	 {
            xbinPub[i] = new float[itemp[1]];
            ybinPubLna[i] = new float[itemp[1]];
	 }

	 for(int i = 0; i < itemp[1]; i++)
	 {
            xbinPub[0][i] = returnVal[3*i];
	    xbinPub[1][i] = 0;
	    xbinPub[2][i] = 0;
	    ybinPubLna[0][i] = returnVal[3*i+1];
	    ybinPubLna[1][i] = returnVal[3*i+2];
	    ybinPubLna[2][i] = returnVal[3*i+2];

	    cout << i+1 << ", data: " << xbinPub[0][i] << "\t" << ybinPubLna[0][i] << "\t" << ybinPubLna[1][i] << endl;
	 }

         grPubLna = new TGraphAsymmErrors(itemp[1], xbinPub[0], ybinPubLna[0], xbinPub[1], xbinPub[2], ybinPubLna[1], ybinPubLna[2]);
         mystyle->SetGraphColor(grPubLna, 0);
         grPubLna->GetYaxis()->SetRangeUser(-0.7, 5.);

         // Plotting each graph separately
	 // MVA cut only
         mystyle->SetAxisTitles(grMva, "FD energy [log(E/eV)]", "MVA applied cut");
         grMva->Draw(plotInstr[0].c_str());
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_mva_cut_" + stemp[1] + ".pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            if(autoFrac)
	    {
               if(simplNorm)
	          stemp[1] = "frac-simple";
	       else
	          stemp[1] = "frac-auto";
	    }
	    else
               stemp[1] = "frac-" + ToString(fraction, 3);
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/mva_cut.pdf";
            c1->SaveAs(stemp[2].c_str());
         }

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
	    {
               if(simplNorm)
	          stemp[1] = "frac-simple";
	       else
	          stemp[1] = "frac-auto";
	    }
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
	    {
               if(simplNorm)
	          stemp[1] = "frac-simple";
	       else
	          stemp[1] = "frac-auto";
	    }
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
	    {
               if(simplNorm)
	          stemp[1] = "frac-simple";
	       else
	          stemp[1] = "frac-auto";
	    }
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
	       {
                  if(simplNorm)
	             stemp[1] = "frac-simple";
	          else
	             stemp[1] = "frac-auto";
	       }
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
	       {
                  if(simplNorm)
	             stemp[1] = "frac-simple";
	          else
	             stemp[1] = "frac-auto";
	       }
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
	       {
                  if(simplNorm)
	             stemp[1] = "frac-simple";
	          else
	             stemp[1] = "frac-auto";
	       }
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
	 legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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

/*	    // LnA data/data-norm
	    nrfigs = 3;
            mystyle->SetAxisTitles(grDataLna, "FD energy [log(E/eV)]", "<lnA> of data events");
//	    grDataLna->SetLineStyle(9);
            grDataLna->Draw(plotInstr[0].c_str());
            grDataLnaNorm->Draw(plotInstr[1].c_str());
            grPubLna->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grDataLna, "Data (raw)", "lp");
	    legend->AddEntry(grDataLnaNorm, "Data (normalized)", "lp");
	    legend->AddEntry(grPubLna, "Data (Auger published)", "lp");
            legend->SetBorderSize(1);
            legend->SetMargin(0.3);
            legend->Draw("same");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_lnA_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_lnA_data-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
//	    grDataLna->SetLineStyle(1);*/
	 }

         // Collections of three (signal/background/data, signal/background/data-norm, signal-norm/background-norm/data, signal-norm/background-norm/data-norm, signal/background/data/data-norm)
	 nrfigs = 3;
	 // Signal/background/data
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
         grSig->Draw(plotInstr[0].c_str());
         grBack->Draw(plotInstr[1].c_str());
         grData->Draw(plotInstr[1].c_str());
	 legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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

            // Signal-norm/background-norm/data-norm
            mystyle->SetAxisTitles(grSigNorm, "FD energy [log(E/eV)]", "Purity of simulation events/Signal fraction of data events");
            grSigNorm->Draw(plotInstr[0].c_str());
            grBackNorm->Draw(plotInstr[1].c_str());
            grDataNorm->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
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

	    nrfigs = 4;
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
         grPubLna->Draw(plotInstr[1].c_str());
	 legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	 legend->SetFillStyle(1001);
	 legend->SetFillColor(c_Legend);
	 legend->AddEntry(grDataLna, "Data (raw)", "lp");
	 legend->AddEntry(grPubLna, "Data (Auger published)", "lp");
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
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_lnA_data-only.pdf";
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
            grPubLna->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grDataLnaNorm, "Data (normalized)", "lp");
	    legend->AddEntry(grPubLna, "Data (Auger published)", "lp");
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
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_lnA_data-only_normalized.pdf";
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
	    nrfigs = 3;
            mystyle->SetAxisTitles(grDataLna, "FD energy [log(E/eV)]", "<lnA> of data events");
//	    grDataLna->SetLineStyle(9);
            grDataLna->Draw(plotInstr[0].c_str());
            grDataLnaNorm->Draw(plotInstr[1].c_str());
            grPubLna->Draw(plotInstr[1].c_str());
	    legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrfigs)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
	    legend->SetFillStyle(1001);
	    legend->SetFillColor(c_Legend);
	    legend->AddEntry(grDataLna, "Data (raw)", "lp");
	    legend->AddEntry(grDataLnaNorm, "Data (normalized)", "lp");
	    legend->AddEntry(grPubLna, "Data (Auger published)", "lp");
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
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_lnA_data-data_normalized.pdf";
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
            delete[] ybinBackInv[i];
            delete[] ybinBackNorm[i];
            delete[] ybinData[i];
            delete[] ybinDataNorm[i];
            delete[] ybinDataLna[i];
            delete[] ybinDataLnaNorm[i];

            delete[] xbinPub[i];
            delete[] ybinPubLna[i];

	    delete[] ymvaCut[i];

            if(autoFrac && !simplNorm)
            {
               delete[] ylowfractest[i];
               delete[] yhighfractest[i];
            }
         }      

         if(autoFrac && !simplNorm)
         {
            delete[] ylowlimitfrac;
            delete[] yhighlimitfrac;
         }

	 delete[] xlimit;
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
