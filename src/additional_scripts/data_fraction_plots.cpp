#define _STANDALONE_ 1
#include "workstation.h"
#include "separate_functions.h"
#include "root_style.h"
#include "mva_result_read.h"

using namespace std;

void SetColor(TGraphAsymmErrors *gr, int cur, int nrbins)
{
   int ci = 1738+cur;
   double colormix = (double)cur/(double)nrbins;
   TColor *color = new TColor(ci, 1.-colormix, 0, colormix, "", 1);

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

   int nrbins = -1;

   ResultRead *analRes = new ResultRead();

   string filename;
   if(argc > 1)
   {
      ftemp = new float[2];
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
         float *ybinBack[3];
         float *ybinData[3];
         float *ybinDataNorm[3];
	 float *ylowfractest[3], *yhighfractest[3];
	 float *ylowlimitfrac, *yhighlimitfrac;

         for(int i = 0; i < 3; i++)
         {
            xbin[i] = new float[nrbins];
            ybinSig[i] = new float[nrbins];
            ybinBack[i] = new float[nrbins];
            ybinData[i] = new float[nrbins];
            ybinDataNorm[i] = new float[nrbins];

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

            // Getting xbin values
            xbin[0][i] = analRes->GetEnergy();
            analRes->GetEnergyError(ftemp);
            xbin[1][i] = ftemp[0];
            xbin[2][i] = ftemp[1];

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
                  ylowlimitfrac[i] = (ylowfractest[0][i]-ylowfractest[1][i])/100.;
                  yhighlimitfrac[i] = (yhighfractest[0][i]+yhighfractest[2][i])/100.;
	       }
	       else
	       {
                  ylowlimitfrac[i] = (yhighfractest[0][i]-yhighfractest[1][i])/100.;
                  yhighlimitfrac[i] = (ylowfractest[0][i]+ylowfractest[2][i])/100.;
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

	       fraction = ybinDataNorm[0][i]/100.;
	    }

	    // Decide if we want to create plots or not
	    runnorm = true;
	 }

         // Minimum and maximum values for all bins and types
         cerr << "Minimum and maximum values for all bins (without error bars):" << endl;
         cerr << "- Signal:     \t" << TMath::MinElement(nrbins, ybinSig[0]) << "\t" << TMath::MaxElement(nrbins, ybinSig[0]) << endl;
         cerr << "- Background: \t" << TMath::MinElement(nrbins, ybinBack[0]) << "\t" << TMath::MaxElement(nrbins, ybinBack[0]) << endl;
         cerr << "- Data:       \t" << TMath::MinElement(nrbins, ybinData[0]) << "\t" << TMath::MaxElement(nrbins, ybinData[0]) << endl;
         if(runnorm)
            cerr << "- Norm. Data:\t" << TMath::MinElement(nrbins, ybinDataNorm[0]) << "\t" << TMath::MaxElement(nrbins, ybinDataNorm[0]) << endl;

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
         gStyle->SetEndErrorSize(6);

         // Preparing all graphs
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

         TGraphAsymmErrors *grDataNorm;
         // Normalized Data
         if(runnorm)
         {
            grDataNorm = new TGraphAsymmErrors(nrbins, xbin[0], ybinDataNorm[0], xbin[1], xbin[2], ybinDataNorm[1], ybinDataNorm[2]);
            mystyle->SetGraphColor(grDataNorm, 2);
            grDataNorm->GetYaxis()->SetRangeUser(yrange[0], yrange[1]);
         }

         // Plotting each graph separately
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of signal simulation events");
         grSig->Draw("AP");
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

         mystyle->SetAxisTitles(grBack, "FD energy [log(E/eV)]", "Purity of background simulation events");
         grBack->Draw("AP");
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

         mystyle->SetAxisTitles(grData, "FD energy [log(E/eV)]", "Signal fraction of data events");
         grData->Draw("AP");
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
            mystyle->SetAxisTitles(grDataNorm, "FD energy [log(E/eV)]", "Signal fraction of data events (normalized)");
            grDataNorm->Draw("AP");
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

         // Plotting combinations
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events");
         grSig->Draw("AP");
         grBack->Draw("P;SAME");
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

         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events and signal fraction of data events");
         grSig->Draw("AP");
         grData->Draw("P;SAME");
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_sig-data.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_sig-data.pdf";
            c1->SaveAs(stemp[2].c_str());
         }

         if(runnorm)
         {
            mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events and signal fraction of data events (normalized)");
            grSig->Draw("AP");
            grDataNorm->Draw("P;SAME");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_sig-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_sig-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
         }

         mystyle->SetAxisTitles(grBack, "FD energy [log(E/eV)]", "Purity of simulation events and signal fraction of data events");
         grBack->Draw("AP");
         grData->Draw("P;SAME");
         if(itemp[0] == 0)
         {
            stemp[1] = analRes->GetObservableType();
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_back-data.pdf";
            c1->SaveAs(stemp[2].c_str());
         }
         else if(itemp[0] == 1)
         {
            stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_back-data.pdf";
            c1->SaveAs(stemp[2].c_str());
         }

         if(runnorm)
         {
            mystyle->SetAxisTitles(grBack, "FD energy [log(E/eV)]", "Purity of simulation events and signal fraction of data events (normalized)");
            grBack->Draw("AP");
            grDataNorm->Draw("P;SAME");
            if(itemp[0] == 0)
            {
               stemp[1] = analRes->GetObservableType();
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/individual_fraction_plot_" + stemp[1] + "_back-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
            else if(itemp[0] == 1)
            {
               stemp[2] = RemoveFilename(&stemp[0]) + "/plots/fraction_plot_back-data_normalized.pdf";
               c1->SaveAs(stemp[2].c_str());
            }
         }

         if(runnorm)
         {
            mystyle->SetAxisTitles(grData, "FD energy [log(E/eV)]", "Signal fraction of data events (normalized)");
            grData->Draw("AP");
            grDataNorm->Draw("P;SAME");
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
         }

         // All together
         mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events and signal fraction of data events");
         grSig->Draw("AP");
         grBack->Draw("P;SAME");
         grData->Draw("P;SAME");
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
            mystyle->SetAxisTitles(grSig, "FD energy [log(E/eV)]", "Purity of simulation events and signal fraction of data events (normalized)");
            grSig->Draw("AP");
            grBack->Draw("P;SAME");
            grDataNorm->Draw("P;SAME");
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
         }

         printRes->Close();

         for(int i = 0; i < 3; i++)
         {
            delete[] xbin[i];
            delete[] ybinSig[i];
            delete[] ybinBack[i];
            delete[] ybinData[i];
            delete[] ybinDataNorm[i];

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
