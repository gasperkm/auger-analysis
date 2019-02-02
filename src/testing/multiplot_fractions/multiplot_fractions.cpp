#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "root_style.h"
#include "primary_type.h"

using namespace std;

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   string inname;
   TCanvas *c1;
   TCanvas *c2;
   TCanvas *c3;
   TPad *pads[10];
   TFile *infile;
   TGraphAsymmErrors *grindiff[15];
   TGraphAsymmErrors *grindiffOFe[5];
   TGraphAsymmErrors *tempgr[3];
   TGraphAsymmErrors *grinfinal[3];
   TGraphAsymmErrors *grinfinalSyst[3];
   TGraphAsymmErrors *grinfinalSystModel[3];
   TGraph *simpgr;
   TF1 *simpfit[2];
   TText *t;
   string *stemp = new string[6];
   int *itemp = new int[3];
   double *dtemp = new double[4];
   double maxval = 0.;

   vector<double> statmean;
   vector<double> statdev;
   vector<double> statmeanOFe;
   vector<double> statdevOFe;
   vector<double> sysdev;

   int pcount[3];

   double xpval[3];
   double ypval[3];
   double xpvalOFe[3];
   double ypvalOFe[3];

   bool noepos, noqgs, nosib;

   vector<double> energy;
   vector<double> tempenergy;
   vector<double> fraction;
   vector<double> fractionErr;

   RootStyle *mystyle;

   string plotdir;

   if(argc > 3)
   {
      // Three input files to determine systematics
      cout << "Input files:" << endl;
      cout << "- " << argv[1] << " (EPOS-LHC model)" << endl;
      cout << "- " << argv[2] << " (QGSJET-II.04 model)" << endl;
      cout << "- " << argv[3] << " (SIBYLL-2.3 model)" << endl;

      mystyle = new RootStyle();
      mystyle->SetBaseStyle();
      c2 = new TCanvas("c2","",1200,900);
      mystyle->SetSinglePlot(0, -1, c2);

      c1 = new TCanvas("c1","",1200,1600);
      mystyle->SetMultiPlot(0, 0, 4, c1);
      c1->cd();
      mystyle->SetPaddedPlot(4, c1, pads);

      gStyle->SetEndErrorSize(6);

      simpgr = new TGraph(2);
      simpgr->SetPoint(0, 18.55, 0);
      simpgr->SetPoint(1, 19.75, 0);

      inname = argv[1];
      plotdir = RemoveFilename(&inname);
      stemp[0] = "mkdir -p " + plotdir + "/multi";
      system(stemp[0].c_str());

//      cerr << "Set Y-axis title for systematic difference plots: ";
//      getline(cin, stemp[3]);

      cerr << "Select naming convention for output files: ";
      getline(cin, stemp[4]);

      cerr << "Select tree: ";
      getline(cin, stemp[5]);

      for(int j = 0; j < 3; j++)
         mystyle->CreateColorScale(j, 3);

      grindiffOFe[0] = new TGraphAsymmErrors();
      grindiffOFe[1] = new TGraphAsymmErrors();
      grindiffOFe[2] = new TGraphAsymmErrors();

      // Go through the elemental composition
      for(int j = 0; j < 4; j++)
      {
         grindiff[3*j] = new TGraphAsymmErrors();
         grindiff[3*j+1] = new TGraphAsymmErrors();
         grindiff[3*j+2] = new TGraphAsymmErrors();

         energy.clear();
         tempenergy.clear();
	 fraction.clear();
	 fractionErr.clear();

	 statmean.clear();
	 statdev.clear();
	 statmeanOFe.clear();
	 statdevOFe.clear();

         if(j == 0)
            stemp[1] = "Proton";
	 else if(j == 1)
            stemp[1] = "Helium";
	 else if(j == 2)
            stemp[1] = "Oxygen";
	 else if(j == 3)
            stemp[1] = "Iron";

	 stemp[0] = "composition_" + stemp[1] + "_" + stemp[5];

	 // Go through the open files and get the composition graphs
         for(int i = 0; i < 3; i++)
         {
            inname = argv[i+1];
            infile = TFile::Open(inname.c_str(), "READ");

            if(infile->GetListOfKeys()->Contains(stemp[0].c_str()))
               tempgr[i] = (TGraphAsymmErrors*)infile->Get(stemp[0].c_str());
            
            infile->Close();
         }

	 for(int i = 0; i < tempgr[0]->GetN(); i++)
	 {
            tempgr[0]->GetPoint(i, xpval[0], ypval[0]);
	    grindiff[3*j]->SetPoint(i, xpval[0]-0.02, ypval[0]);
	    xpval[0] = tempgr[0]->GetErrorYlow(i);
	    ypval[0] = tempgr[0]->GetErrorYhigh(i);
	    grindiff[3*j]->SetPointError(i, 0., 0., xpval[0], ypval[0]);

	    // Sum O+Fe
	    if(j == 3)
            {
               tempgr[0]->GetPoint(i, xpval[0], ypval[0]);
               grindiff[3*(j-1)]->GetPoint(i, dtemp[0], dtemp[1]);
	       if(ypval[0]+dtemp[1] < 1.)
                  grindiffOFe[0]->SetPoint(i, dtemp[0], ypval[0]+dtemp[1]);
	       else
                  grindiffOFe[0]->SetPoint(i, dtemp[0], 1.);
	       dtemp[0] = tempgr[0]->GetErrorYlow(i);
	       dtemp[1] = tempgr[0]->GetErrorYhigh(i);
	       dtemp[2] = grindiff[3*(j-1)]->GetErrorYlow(i);
	       dtemp[3] = grindiff[3*(j-1)]->GetErrorYhigh(i);
//	       grindiffOFe[0]->SetPointError(i, 0., 0., TMath::Sqrt(dtemp[0]*dtemp[0] + dtemp[2]*dtemp[2]), TMath::Sqrt(dtemp[1]*dtemp[1] + dtemp[3]*dtemp[3]));
	       grindiffOFe[0]->SetPointError(i, 0., 0., (dtemp[0] + dtemp[2])/2., (dtemp[1] + dtemp[3])/2.);
	    }
	 }

	 for(int i = 0; i < tempgr[1]->GetN(); i++)
	 {
            tempgr[1]->GetPoint(i, xpval[1], ypval[1]);
	    grindiff[3*j+1]->SetPoint(i, xpval[1], ypval[1]);
	    xpval[1] = tempgr[1]->GetErrorYlow(i);
	    ypval[1] = tempgr[1]->GetErrorYhigh(i);
	    grindiff[3*j+1]->SetPointError(i, 0., 0., xpval[1], ypval[1]);

	    // Sum O+Fe
	    if(j == 3)
            {
               tempgr[1]->GetPoint(i, xpval[1], ypval[1]);
               grindiff[3*(j-1)+1]->GetPoint(i, dtemp[0], dtemp[1]);
	       if(ypval[1]+dtemp[1] < 1.)
                  grindiffOFe[1]->SetPoint(i, dtemp[0], ypval[1]+dtemp[1]);
	       else
                  grindiffOFe[1]->SetPoint(i, dtemp[0], 1.);
	       dtemp[0] = tempgr[1]->GetErrorYlow(i);
	       dtemp[1] = tempgr[1]->GetErrorYhigh(i);
	       dtemp[2] = grindiff[3*(j-1)+1]->GetErrorYlow(i);
	       dtemp[3] = grindiff[3*(j-1)+1]->GetErrorYhigh(i);
//	       grindiffOFe[1]->SetPointError(i, 0., 0., TMath::Sqrt(dtemp[0]*dtemp[0] + dtemp[2]*dtemp[2]), TMath::Sqrt(dtemp[1]*dtemp[1] + dtemp[3]*dtemp[3]));
	       grindiffOFe[1]->SetPointError(i, 0., 0., (dtemp[0] + dtemp[2])/2., (dtemp[1] + dtemp[3])/2.);
	    }
	 }

	 for(int i = 0; i < tempgr[2]->GetN(); i++)
	 {
            tempgr[2]->GetPoint(i, xpval[2], ypval[2]);
	    grindiff[3*j+2]->SetPoint(i, xpval[2]+0.02, ypval[2]);
	    xpval[2] = tempgr[2]->GetErrorYlow(i);
	    ypval[2] = tempgr[2]->GetErrorYhigh(i);
	    grindiff[3*j+2]->SetPointError(i, 0., 0., xpval[2], ypval[2]);

	    // Sum O+Fe
	    if(j == 3)
            {
               tempgr[2]->GetPoint(i, xpval[2], ypval[2]);
               grindiff[3*(j-1)+2]->GetPoint(i, dtemp[0], dtemp[1]);
	       if(ypval[2]+dtemp[1] < 1.)
                  grindiffOFe[2]->SetPoint(i, dtemp[0], ypval[2]+dtemp[1]);
	       else
                  grindiffOFe[2]->SetPoint(i, dtemp[0], 1.);
	       dtemp[0] = tempgr[2]->GetErrorYlow(i);
	       dtemp[1] = tempgr[2]->GetErrorYhigh(i);
	       dtemp[2] = grindiff[3*(j-1)+2]->GetErrorYlow(i);
	       dtemp[3] = grindiff[3*(j-1)+2]->GetErrorYhigh(i);
//	       grindiffOFe[2]->SetPointError(i, 0., 0., TMath::Sqrt(dtemp[0]*dtemp[0] + dtemp[2]*dtemp[2]), TMath::Sqrt(dtemp[1]*dtemp[1] + dtemp[3]*dtemp[3]));
	       grindiffOFe[2]->SetPointError(i, 0., 0., (dtemp[0] + dtemp[2])/2., (dtemp[1] + dtemp[3])/2.);
	    }
	 }

	 // Get the largest number of points
	 cout << endl;
	 itemp[0] = 0;
         for(int i = 0; i < 3; i++)
	 {
            cout << i << ": Number of points = " << tempgr[i]->GetN() << endl;
	    if(tempgr[i]->GetN() > itemp[0])
	    {
               itemp[0] = tempgr[i]->GetN();
	       itemp[1] = i;
	    }
	 }

	 // Get X points
         for(int i = 0; i < itemp[0]; i++)
	 {
            tempgr[itemp[1]]->GetPoint(i, xpval[0], ypval[0]);
            tempenergy.push_back(xpval[0]);
         }

         for(int i = 0; i < itemp[0]; i++)
	 {
            // Test if some points are missing in any of the graphs
            noepos = false;
	    noqgs = false;
	    nosib = false;

	    if(!noepos)
	    {
               tempgr[0]->GetPoint(i, xpval[0], ypval[0]);
	       if(xpval[0] != tempenergy[i])
                  noepos = true;
	    }
            
	    if(!noqgs)
	    {
               tempgr[1]->GetPoint(i, xpval[1], ypval[1]);
	       if(xpval[1] != tempenergy[i])
                  noqgs = true;
	    }
               
	    if(!nosib)
	    {
               tempgr[2]->GetPoint(i, xpval[2], ypval[2]);
	       if(xpval[2] != tempenergy[i])
                  nosib = true;
	    }

            // Calculate mean of model points
	    dtemp[0] = 0;
	    dtemp[2] = 0;
	    itemp[2] = 0;
	    if(!noepos)
	    {
               dtemp[0] += ypval[0];
	       if(j == 3)
	       {
	          grindiffOFe[0]->GetPoint(i, xpvalOFe[0], ypvalOFe[0]);
		  dtemp[2] += ypvalOFe[0];
	       }
	       itemp[2]++;
	    }
	    if(!noqgs)
	    {
               dtemp[0] += ypval[1];
	       if(j == 3)
	       {
	          grindiffOFe[1]->GetPoint(i, xpvalOFe[1], ypvalOFe[1]);
		  dtemp[2] += ypvalOFe[1];
	       }
	       itemp[2]++;
	    }
	    if(!nosib)
	    {
               dtemp[0] += ypval[2];
	       if(j == 3)
	       {
	          grindiffOFe[2]->GetPoint(i, xpvalOFe[2], ypvalOFe[2]);
		  dtemp[2] += ypvalOFe[2];
	       }
	       itemp[2]++;
	    }

	    dtemp[0] = dtemp[0]/itemp[2];
            cout << "Model mean value (" << ypval[0] << ", " << ypval[1] << ", " << ypval[2] << ") = " << dtemp[0] << endl;
	    if(j == 3)
	    {
	       dtemp[2] = dtemp[2]/itemp[2];
               cout << "Model mean value, O+Fe (" << ypvalOFe[0] << ", " << ypvalOFe[1] << ", " << ypvalOFe[2] << ") = " << dtemp[2] << endl;
	    }

            // Calculate sigma of model points
	    dtemp[1] = 0;
	    dtemp[3] = 0;
	    if(!noepos)
	    {
               dtemp[1] += TMath::Power(ypval[0]-dtemp[0],2);
	       if(j == 3)
                  dtemp[3] += TMath::Power(ypvalOFe[0]-dtemp[2],2);
	    }
	    if(!noqgs)
	    {
               dtemp[1] += TMath::Power(ypval[1]-dtemp[0],2);
	       if(j == 3)
                  dtemp[3] += TMath::Power(ypvalOFe[1]-dtemp[2],2);
	    }
	    if(!nosib)
	    {
               dtemp[1] += TMath::Power(ypval[2]-dtemp[0],2);
	       if(j == 3)
                  dtemp[3] += TMath::Power(ypvalOFe[2]-dtemp[2],2);
	    }

	    dtemp[1] = TMath::Sqrt(dtemp[1]/itemp[2]);
            cout << "Model sigma value (" << ypval[0]-dtemp[0] << ", " << ypval[1]-dtemp[0] << ", " << ypval[2]-dtemp[0] << ") = " << dtemp[1] << endl;
	    if(j == 3)
	    {
	       dtemp[3] = TMath::Sqrt(dtemp[3]/itemp[2]);
               cout << "Model sigma value, O+Fe (" << ypvalOFe[0]-dtemp[2] << ", " << ypvalOFe[1]-dtemp[2] << ", " << ypvalOFe[2]-dtemp[2] << ") = " << dtemp[3] << endl;
	    }

	    statmean.push_back(dtemp[0]);
	    statdev.push_back(dtemp[1]);

	    if(j == 3)
	    {
	       statmeanOFe.push_back(dtemp[2]);
	       statdevOFe.push_back(dtemp[3]);
	    }
	 }

	 cout << endl << "Final mean and sigma results (" << stemp[1] << "):" << endl;
	 cout << "# " << stemp[1] << " (point, mean, sigma)" << endl;
	 cout << statmean.size() << endl;
	 for(int i = 0; i < statmean.size(); i++)
            cout << i << "\t" << statmean[i] << "\t" << statdev[i] << endl;
	 if(j == 3)
	 {
	    cout << endl << "Final mean and sigma results (" << "Oxygen + Iron" << "):" << endl;
	    cout << "# " << "Oxygen + Iron" << " (point, mean, sigma)" << endl;
	    cout << statmeanOFe.size() << endl;
	    for(int i = 0; i < statmeanOFe.size(); i++)
               cout << i << "\t" << statmeanOFe[i] << "\t" << statdevOFe[i] << endl;
	 }

         c2->cd();

	 mystyle->SetGraphColor(grindiff[3*j], 1);
	 mystyle->SetGraphColor(grindiff[3*j+1], 0);
	 mystyle->SetGraphColor(grindiff[3*j+2], 2);

	 mystyle->SetColorScale(grindiff[3*j], 0, 3);
	 mystyle->SetColorScale(grindiff[3*j+1], 1, 3);
	 mystyle->SetColorScale(grindiff[3*j+2], 2, 3);

	 simpgr->Draw("AP");
         simpgr->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c2));
         simpgr->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c2));
	 simpgr->GetYaxis()->SetRange(0., 1.);
	 simpgr->GetYaxis()->SetRangeUser(0., 1.);
	 mystyle->SetAxisTitles(simpgr, "FD energy [log(E/eV)]", "Elemental fractions");

	 grindiff[3*j]->Draw("P;SAME");
	 grindiff[3*j+1]->Draw("P;SAME");
	 grindiff[3*j+2]->Draw("P;SAME");

	 c2->Update();

         t = new TText();
         t->SetTextAlign(31);
         t->SetTextColor(1);
         t->SetTextSize(17);
         t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), stemp[1].c_str());

         stemp[2] = plotdir + "/multi/multiplot-fractions_" + stemp[1] + "_" + stemp[4] + "_" + stemp[5] + ".pdf";
         c2->SaveAs(stemp[2].c_str());

	 // Plot O+Fe combined
	 if(j == 3)
	 {
            c2->cd();

	    mystyle->SetGraphColor(grindiffOFe[0], 1);
	    mystyle->SetGraphColor(grindiffOFe[1], 0);
	    mystyle->SetGraphColor(grindiffOFe[2], 2);

	    mystyle->SetColorScale(grindiffOFe[0], 0, 3);
	    mystyle->SetColorScale(grindiffOFe[1], 1, 3);
	    mystyle->SetColorScale(grindiffOFe[2], 2, 3);

	    simpgr->Draw("AP");
            simpgr->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c2));
            simpgr->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c2));
	    simpgr->GetYaxis()->SetRange(0., 1.);
	    simpgr->GetYaxis()->SetRangeUser(0., 1.);
	    mystyle->SetAxisTitles(simpgr, "FD energy [log(E/eV)]", "Elemental fractions");

	    grindiffOFe[0]->Draw("P;SAME");
	    grindiffOFe[1]->Draw("P;SAME");
	    grindiffOFe[2]->Draw("P;SAME");

	    c2->Update();

            t = new TText();
            t->SetTextAlign(31);
            t->SetTextColor(1);
            t->SetTextSize(17);
            t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Oxygen + Iron");

            stemp[2] = plotdir + "/multi/multiplot-fractions_OxygenIron_" + stemp[4] + "_" + stemp[5] + ".pdf";
            c2->SaveAs(stemp[2].c_str());
	 }

	 c1->cd();
	 pads[j]->cd();

	 mystyle->SetGraphColor(grindiff[3*j], 1);
	 mystyle->SetGraphColor(grindiff[3*j+1], 0);
	 mystyle->SetGraphColor(grindiff[3*j+2], 2);

	 mystyle->SetColorScale(grindiff[3*j], 0, 3);
	 mystyle->SetColorScale(grindiff[3*j+1], 1, 3);
	 mystyle->SetColorScale(grindiff[3*j+2], 2, 3);

/*	 grindiff[3*j]->SetMarkerStyle(20);
	 grindiff[3*j+1]->SetMarkerStyle(21);
	 grindiff[3*j+2]->SetMarkerStyle(22);*/
/*         grindiff[3*j]->GetYaxis()->SetTitleOffset(mystyle->GetPaddedYoffset(4, c1));
         grindiff[3*j]->GetXaxis()->SetTitleOffset(mystyle->GetPaddedXoffset(4, c1));
	 if(j == 3)
	    mystyle->SetAxisTitles(grindiff[3*j], "FD energy [log(E/eV)]", "Elemental fractions");
	 else
	    mystyle->SetAxisTitles(grindiff[3*j], "", "Elemental fractions");*/

	 simpgr->Draw("AP");
         simpgr->GetYaxis()->SetTitleOffset(mystyle->GetPaddedYoffset(4, c1));
         simpgr->GetXaxis()->SetTitleOffset(mystyle->GetPaddedXoffset(4, c1));
	 if(j == 3)
	    mystyle->SetAxisTitles(simpgr, "FD energy [log(E/eV)]", "Elemental fractions");
	 else
	    mystyle->SetAxisTitles(simpgr, "", "Elemental fractions");

	 grindiff[3*j]->Draw("P;SAME");
	 grindiff[3*j+1]->Draw("P;SAME");
	 grindiff[3*j+2]->Draw("P;SAME");
      }

      for(int j = 0; j < 4; j++)
      {
	 pads[j]->cd();
	 simpgr->GetYaxis()->SetRange(0., 1.);
	 simpgr->GetYaxis()->SetRangeUser(0., 1.);
      }

      c1->Update();

      t = new TText();
      t->SetTextAlign(31);
      t->SetTextColor(1);
      t->SetTextSize(17);

      pads[0]->cd();
      t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Proton");
      pads[1]->cd();
      t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Helium");
      pads[2]->cd();
      t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Oxygen");
      pads[3]->cd();
      t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Iron");

      stemp[2] = plotdir + "/multi/multiplot-fractions_" + stemp[4] + "_" + stemp[5] + ".pdf";
      c1->SaveAs(stemp[2].c_str());

      // Plot combined elemental fractions: f = <f> +- sigma_stat +- sigma_sys +- sigma_model (take statistics and systematics from one model - EPOS-LHC)
      mystyle->SetBaseStyle();
      c3 = new TCanvas("c3","",1200,900);
      mystyle->SetSinglePlot(0, -1, c3);

      cout << endl << "All elemental fraction points for models and elements:" << endl;
      for(int m = 0; m < 3; m++)
      {
         if(m == 0)
            cout << "1) EPOS-LHC:" << endl;
         else if(m == 1)
            cout << "2) QGSJET-II.04:" << endl;
         else if(m == 2)
            cout << "3) Sibyll-2.3:" << endl;

         for(int j = 0; j < 5; j++)
         {
            if(j < 4)
	    {
               if(j == 0)
                  cout << "- Proton:" << endl;
               else if(j == 1)
                  cout << "- Helium:" << endl;
               else if(j == 2)
                  cout << "- Oxygen:" << endl;
               else if(j == 3)
                  cout << "- Iron:" << endl;

               for(int k = 0; k < grindiff[3*j+m]->GetN(); k++)
               {
                  grindiff[3*j+m]->GetPoint(k, dtemp[0], dtemp[1]);
	          dtemp[2] = grindiff[3*j+m]->GetErrorYlow(k);
	          dtemp[3] = grindiff[3*j+m]->GetErrorYhigh(k);
                  cout << dtemp[0] << "\t" << dtemp[1] << " (+" << dtemp[3] << ", -" << dtemp[2] << ")" << endl;
               }
	    }
	    else
	    {
               cout << "- Oxygen + Iron:" << endl;

               for(int k = 0; k < grindiffOFe[m]->GetN(); k++)
               {
                  grindiffOFe[m]->GetPoint(k, dtemp[0], dtemp[1]);
	          dtemp[2] = grindiffOFe[m]->GetErrorYlow(k);
	          dtemp[3] = grindiffOFe[m]->GetErrorYhigh(k);
                  cout << dtemp[0] << "\t" << dtemp[1] << " (+" << dtemp[3] << ", -" << dtemp[2] << ")" << endl;
               }
	    }
         }
      }

      energy.clear();
      fraction.clear();
      fractionErr.clear();
      statdev.clear();
      sysdev.clear();

      // Get energy (only the mean, not shifted energy)
      cout << endl << "Energy:" << endl;
      for(int j = 0; j < 3; j++)
      {
         for(int k = 0; k < grindiff[3*j]->GetN(); k++)
         {
            grindiff[3*j+1]->GetPoint(k, xpval[0], ypval[0]);
	    energy.push_back(xpval[0]);
	    cout << "  " << xpval[0] << " (" << j << "," << k << ")" << endl;
	 }
      }

      // Get mean elemental fractions
      // Get systematic uncertainty from hadronic interaction model (depending on energy)
      cout << endl << "Mean and HI model uncertainty:" << endl;
      for(int j = 0; j < 3; j++)
      {
	 if(j == 2)
	 {
            for(int k = 0; k < grindiffOFe[0]->GetN(); k++)
            {
               dtemp[0] = 0;
	       dtemp[1] = 0;

               grindiffOFe[0]->GetPoint(k, xpval[0], ypval[0]); // EPOS-LHC
	       dtemp[0] += ypval[0];

               grindiffOFe[1]->GetPoint(k, xpval[1], ypval[1]); // QGSJET-II.04
	       dtemp[0] += ypval[1];

               grindiffOFe[2]->GetPoint(k, xpval[2], ypval[2]); // Sibyll-2.3
	       dtemp[0] += ypval[2];

	       // Mean value
	       dtemp[0] = dtemp[0]/3.;
	       fraction.push_back(dtemp[0]);
	       cout << "  " << dtemp[0] << " (" << j << "," << k << ")" << endl;

	       // Model uncertainty
	       dtemp[1] += TMath::Power(ypval[0]-dtemp[0],2);
	       dtemp[1] += TMath::Power(ypval[1]-dtemp[0],2);
	       dtemp[1] += TMath::Power(ypval[2]-dtemp[0],2);
	       dtemp[1] = TMath::Sqrt(dtemp[1]/3.);
	       statdev.push_back(dtemp[1]);
	       cout << "  " << dtemp[1] << " (" << j << "," << k << ")" << endl;
	    }
	 }
	 else
	 {
            for(int k = 0; k < grindiff[3*j]->GetN(); k++)
            {
               dtemp[0] = 0;
	       dtemp[1] = 0;

               grindiff[3*j]->GetPoint(k, xpval[0], ypval[0]); // EPOS-LHC
	       dtemp[0] += ypval[0];

               grindiff[3*j+1]->GetPoint(k, xpval[1], ypval[1]); // QGSJET-II.04
	       dtemp[0] += ypval[1];

               grindiff[3*j+2]->GetPoint(k, xpval[2], ypval[2]); // Sibyll-2.3
	       dtemp[0] += ypval[2];

	       // Mean value
	       dtemp[0] = dtemp[0]/3.;
	       fraction.push_back(dtemp[0]);
	       cout << "  " << dtemp[0] << " (" << j << "," << k << ")" << endl;

	       // Model uncertainty
	       dtemp[1] += TMath::Power(ypval[0]-dtemp[0],2);
	       dtemp[1] += TMath::Power(ypval[1]-dtemp[0],2);
	       dtemp[1] += TMath::Power(ypval[2]-dtemp[0],2);
	       dtemp[1] = TMath::Sqrt(dtemp[1]/3.);
	       statdev.push_back(dtemp[1]);
	       cout << "  " << dtemp[1] << " (" << j << "," << k << ")" << endl;
	    }
	 }
      } 

      // Get statistic uncertainty of elemental fraction: Taking EPOS-LHC uncertainties
      cout << endl << "Statistic uncertainty:" << endl;
      for(int j = 0; j < 3; j++)
      {
	 if(j == 2)
	 {
            for(int k = 0; k < grindiffOFe[0]->GetN(); k++)
            {
	       dtemp[0] = grindiffOFe[0]->GetErrorYlow(k);
	       dtemp[1] = grindiffOFe[0]->GetErrorYhigh(k);

	       fractionErr.push_back(dtemp[0]);
	       fractionErr.push_back(dtemp[1]);
	       cout << "  " << dtemp[0] << "\t" << dtemp[1] << " (" << j << "," << k << ")" << endl;
	    }
	 }
	 else
	 {
            for(int k = 0; k < grindiff[3*j]->GetN(); k++)
            {
	       dtemp[0] = grindiff[3*j]->GetErrorYlow(k);
	       dtemp[1] = grindiff[3*j]->GetErrorYhigh(k);

	       fractionErr.push_back(dtemp[0]);
	       fractionErr.push_back(dtemp[1]);
	       cout << "  " << dtemp[0] << "\t" << dtemp[1] << " (" << j << "," << k << ")" << endl;
	    }
	 }
      } 

      // Get systematic uncertainty from other contributions: Taking EPOS-LHC uncertainties = 0.066 for proton, 0.097 for helium, (0.090+0.085)/2 = 0.088 for oxygen+iron
      cout << endl << "Systematic uncertainty:" << endl;
      for(int j = 0; j < 3; j++)
      {
         for(int k = 0; k < grindiff[3*j+1]->GetN(); k++)
         {
            if(j == 0)
               dtemp[0] = 0.066;
            else if(j == 1)
               dtemp[0] = 0.097;
            else if(j == 2)
               dtemp[0] = 0.088;

	    sysdev.push_back(dtemp[0]);
	    cout << "  " << dtemp[0] << " (" << j << "," << k << ")" << endl;
	 }
      }

      // Prepare graphs
      grinfinal[0] = new TGraphAsymmErrors();
      grinfinal[1] = new TGraphAsymmErrors();
      grinfinal[2] = new TGraphAsymmErrors();
      grinfinalSyst[0] = new TGraphAsymmErrors();
      grinfinalSyst[1] = new TGraphAsymmErrors();
      grinfinalSyst[2] = new TGraphAsymmErrors();
      grinfinalSystModel[0] = new TGraphAsymmErrors();
      grinfinalSystModel[1] = new TGraphAsymmErrors();
      grinfinalSystModel[2] = new TGraphAsymmErrors();

      itemp[0] = 0;
      for(int m = 0; m < 3; m++)
      {
         for(int k = 0; k < energy.size()/3.; k++)
	 {
            grinfinal[m]->SetPoint(k, energy[itemp[0]+k], fraction[itemp[0]+k]);
            grinfinalSyst[m]->SetPoint(k, energy[itemp[0]+k], fraction[itemp[0]+k]);
            grinfinalSystModel[m]->SetPoint(k, energy[itemp[0]+k], fraction[itemp[0]+k]);

            grinfinal[m]->SetPointError(k, 0., 0., fractionErr[2*(itemp[0]+k)], fractionErr[2*(itemp[0]+k)+1]);
            grinfinalSyst[m]->SetPointError(k, 0., 0., sysdev[itemp[0]+k], sysdev[itemp[0]+k]);
            grinfinalSystModel[m]->SetPointError(k, 0., 0., statdev[itemp[0]+k]+sysdev[itemp[0]+k], statdev[itemp[0]+k]+sysdev[itemp[0]+k]);
	 }
	 itemp[0] += energy.size()/3.;

	 // Draw each graph separately and save it
         c3->cd();

	 mystyle->SetGraphColor(grinfinalSyst[m], 0);
	 mystyle->SetGraphColor(grinfinalSystModel[m], 0);
	 grinfinalSyst[m]->SetFillColorAlpha(grinfinalSyst[m]->GetMarkerColor(),0.6);
	 grinfinalSystModel[m]->SetFillColorAlpha(grinfinalSystModel[m]->GetMarkerColor(),0.4);
	 mystyle->SetGraphColor(grinfinal[m], 2);
	 grinfinal[m]->SetMarkerSize(1.7);
	 grinfinal[m]->SetMarkerStyle(20);
	 grinfinal[m]->SetMarkerColor(kBlue+3);
	 grinfinal[m]->SetLineWidth(3);
	 grinfinal[m]->SetLineColor(kBlue);

         grinfinalSyst[m]->GetYaxis()->SetTitleOffset(mystyle->GetSingleYoffset(c3));
         grinfinalSyst[m]->GetXaxis()->SetTitleOffset(mystyle->GetSingleXoffset(c3));
	 mystyle->SetAxisTitles(grinfinalSyst[m], "FD energy [log(E/eV)]", "Elemental fractions");
	 grinfinalSyst[m]->GetYaxis()->SetRange(0., 1.);
	 grinfinalSyst[m]->GetYaxis()->SetRangeUser(0., 1.);

	 grinfinalSyst[m]->Draw("A3");
	 grinfinalSystModel[m]->Draw("3;SAME");
	 grinfinal[m]->Draw("P;SAME");

	 c3->Update();

         t = new TText();
         t->SetTextAlign(31);
         t->SetTextColor(1);
         t->SetTextSize(17);
	 if(m == 0)
	 {
            t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Proton");
            stemp[2] = plotdir + "/multi/multiplot-final_Proton_" + stemp[4] + "_" + stemp[5] + ".pdf";
	 }
	 else if(m == 1)
	 {
            t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Helium");
            stemp[2] = plotdir + "/multi/multiplot-final_Helium_" + stemp[4] + "_" + stemp[5] + ".pdf";
	 }
	 else if(m == 2)
	 {
            t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), "Oxygen + Iron");
            stemp[2] = plotdir + "/multi/multiplot-final_OxygenIron_" + stemp[4] + "_" + stemp[5] + ".pdf";
	 }

         c3->SaveAs(stemp[2].c_str());

	 if(m == 0)
            stemp[2] = plotdir + "/multi/multiplot-final_Proton_" + stemp[4] + "_" + stemp[5] + ".C";
	 else if(m == 1)
            stemp[2] = plotdir + "/multi/multiplot-final_Helium_" + stemp[4] + "_" + stemp[5] + ".C";
	 else if(m == 2)
            stemp[2] = plotdir + "/multi/multiplot-final_OxygenIron_" + stemp[4] + "_" + stemp[5] + ".C";

         c3->SaveAs(stemp[2].c_str());
      }

      delete grinfinal[0];
      delete grinfinal[1];
      delete grinfinal[2];
      delete grinfinalSyst[0];
      delete grinfinalSyst[1];
      delete grinfinalSyst[2];
      delete grinfinalSystModel[0];
      delete grinfinalSystModel[1];
      delete grinfinalSystModel[2];

      for(int j = 0; j < 4; j++)
      {
         delete grindiff[3*j];
         delete grindiff[3*j+1];
         delete grindiff[3*j+2];
      }

      delete grindiffOFe[0];
      delete grindiffOFe[1];
      delete grindiffOFe[2];

      delete c1;
      delete mystyle;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add three input files as arguments (lna_composition_results*.root)." << endl;
      return 1;
   }

   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

/*   string inname;
   TFile *ifile;
   int nrkeys;
   string *stemp = new string[3];
   int *itemp = new int[3];
   double *dtemp = new double[3];
   TTree *tempTree;
   float *inval;
   int *nentries;
   TGraph *result[300];
   RootStyle *mystyle;
   TCanvas *c1;
   TLine *line = new TLine();
   line->SetLineWidth(2);
   line->SetLineStyle(7);
   line->SetLineColor(28);
   bool combinedfile = false;
   TF1 *fitfunc[6];
   double fitparam[4];
   int angleset = 0;
   vector<string> treeName;

   TGraph *allS38[20];
   int s38count[20];
   for(int i = 0; i < 20; i++)
   {
      allS38[i] = new TGraph();
      s38count[i] = 0;
   }

   TGraph *allFitparam[60];
   int fitparamcount[20];
   vector<double> midenergy;
   for(int i = 0; i < 20; i++)
   {
      for(int j = 0; j < 3; j++)
         allFitparam[j+3*i] = new TGraph();
      fitparamcount[i] = 0;
   }

   vector<double> energy;
   vector<double> zenith;
   vector<double> s1000;
   vector<double> s38;

   ResultRead *analRes = new ResultRead();

   if(argc > 1)
   {
      cerr << "Select if zenith angle is expressed in secant of the angle (0) or in degrees (1): ";
      cin >> angleset;

      stemp[1] = "mkdir -p ./plots";
      system(stemp[1].c_str());
      stemp[1] = "rm ./plots/*.pdf ./plots/*.C";
      system(stemp[1].c_str());

      for(int k = 0; k < argc-1; k++)
      {
         mystyle = new RootStyle();
         mystyle->SetBaseStyle();
         c1 = new TCanvas("c1","",1200,900);

         inname = string(argv[k+1]);
         cout << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;
         ifile = new TFile(inname.c_str(), "READ");

         stemp[0] = RemoveFilename(&inname) + "/application_results.txt";
         itemp[0] = analRes->ReadFile(stemp[0]);

	 midenergy.push_back(analRes->GetEnergy());

         stemp[0] = "TreeA";
         combinedfile = ifile->GetListOfKeys()->Contains(stemp[0].c_str());
	 cout << "combinedfile = " << combinedfile << endl;

	 if(combinedfile)
            nrkeys = (int)ifile->GetNkeys();
	 else
            nrkeys = (int)ifile->GetNkeys()/2;

         cout << "Number of trees in file = " << nrkeys << endl;
	 inval = new float[3];
	 nentries = new int[nrkeys];

         for(int i = 0; i < nrkeys; i++)
         {
            energy.clear();
	    zenith.clear();
            s1000.clear();
	    s38.clear();

            result[i] = new TGraph();
            result[i+nrkeys] = new TGraph();
            result[i+2*nrkeys] = new TGraph();
            result[i+3*nrkeys] = new TGraph();
//            result[i+4*nrkeys] = new TGraph();

            stemp[0] = "TreeS" + ToString(i+1);
	    if(k == 0)
	       treeName.push_back(stemp[0]);
            cout << "Getting tree = " << stemp[0] << ", ";
            tempTree = (TTree*)ifile->Get(stemp[0].c_str());
            cout << "with title = " << tempTree->GetTitle() << " ---------------" << endl;

            tempTree->SetBranchAddress("shwsize", &inval[0]);
            tempTree->SetBranchAddress("zenithFD", &inval[1]);
            tempTree->SetBranchAddress("energyFD", &inval[2]);

            nentries[i] = tempTree->GetEntries();
            cout << "  Number of entries = " << nentries[i] << endl;

	    dtemp[1] = 0.;

//	    cout << "  Entries printout:" << endl;
            for(int j = 0; j < nentries[i]; j++)
            {
               tempTree->GetEntry(j);

	       if(angleset)
	       {
                  result[i]->SetPoint(j, SecTheta(inval[1],false), inval[0]);
	  	  dtemp[0] = TMath::Power(1./SecTheta(inval[1],false),2) - TMath::Power(TMath::Cos(DegToRad(38.)),2);
	       }
	       else
	       {
	          result[i]->SetPoint(j, inval[1], inval[0]);
	          dtemp[0] = TMath::Power(1./inval[1],2) - TMath::Power(TMath::Cos(DegToRad(38.)),2);
	       }
	       result[i+nrkeys]->SetPoint(j, dtemp[0], inval[0]);

               dtemp[2] += inval[0];

	       energy.push_back(inval[2]);
	       zenith.push_back(inval[1]);
	       s1000.push_back(inval[0]);
	       s38.push_back(dtemp[0]);

//	       cout << "    " << j+1 << ": s1000 = " << inval[0] << ", zenithFD = " << RadToDeg(inval[1]) << ", sec zenithFD = " << SecTheta(inval[1],false) << ", reference = " << dtemp[0] << ", s38? = " << inval[0]/(1.+0.980*dtemp[0]-1.68*dtemp[0]*dtemp[0]-1.3*dtemp[0]*dtemp[0]*dtemp[0]) << endl;
            }

	    dtemp[2] = dtemp[2]/nentries[i];
	    cout << "Mean value of S1000 = " << dtemp[2] << endl;

            mystyle->SetColorScale(result[i], i, nrkeys);
	    result[i]->SetMarkerStyle(20);
	    result[i]->SetMarkerSize(0.8);
            mystyle->SetColorScale(result[i+nrkeys], i, nrkeys);
	    result[i+nrkeys]->SetMarkerStyle(20);
	    result[i+nrkeys]->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(result[i], "FD zenith angle (sec#theta)", "S1000 (VEM)");
            mystyle->SetAxisTitles(result[i+nrkeys], "FD zenith reference [cos^{2}(#theta) - cos^{2}(38 deg)]", "S1000 (VEM)");

	    stemp[2] = "_" + ToString(analRes->GetLowEnergy(),2) + "-" + ToString(analRes->GetHighEnergy(),2);
            
	    c1->SetLogx(kFALSE);
	    c1->SetLogy(kFALSE);
	    result[i]->Draw("AP");
            fitfunc[0] = new TF1("fitfunc0", "[0]+[1]*(TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2))+[2]*TMath::Power((TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2)),2)+[3]*TMath::Power((TMath::Power(1./x,2)-TMath::Power(cos(0.663225115),2)),3)", -0.4, 0.4);
	    fitfunc[0]->SetParameters(dtemp[2],1.,-1.5,-1.3);
	    result[i]->Fit("fitfunc0");
	    c1->Update();
	    line->DrawLine(SecTheta(38.,true), gPad->GetUymin(), SecTheta(38.,true), gPad->GetUymax());
	    stemp[1] = "./plots/s1000_vs_zenithFD_" + stemp[0] + stemp[2] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());
	    stemp[1] = "./plots/s1000_vs_zenithFD_" + stemp[0] + stemp[2] + ".C";
	    c1->SaveAs(stemp[1].c_str());

	    result[i+nrkeys]->Draw("AP");
            fitfunc[1] = new TF1("fitfunc1", "[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)", -0.4, 0.4);
	    fitfunc[1]->SetParameters(dtemp[2],1.,-1.5,-1.3);
	    result[i+nrkeys]->Fit("fitfunc1","Q");
	    stemp[1] = "./plots/s1000_vs_zenithRef_" + stemp[0] + stemp[2] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());
	    stemp[1] = "./plots/s1000_vs_zenithRef_" + stemp[0] + stemp[2] + ".C";
	    c1->SaveAs(stemp[1].c_str());

	    for(int j = 0; j < 4; j++)
               fitparam[j] = fitfunc[0]->GetParameter(j);

	    cout << endl << "Fitting parameters of fCIC for " << stemp[0] <<  " are:" << endl;
	    cout << "- S1000 at 38 deg = " << fitfunc[0]->Eval(SecTheta(38.,true)) << endl;
	    cout << "- a = " << fitparam[1]/fitparam[0] << endl;
	    cout << "- b = " << fitparam[2]/fitparam[0] << endl;
	    cout << "- c = " << fitparam[3]/fitparam[0] << endl;
	    cout << endl;

            for(int j = 0; j < nentries[i]; j++)
	    {
	       s38[j] = s1000[j]/(1.+(fitparam[1]/fitparam[0])*s38[j]+(fitparam[2]/fitparam[0])*TMath::Power(s38[j],2)+(fitparam[3]/fitparam[0])*TMath::Power(s38[j],3));
	       if(angleset)
	          result[i+2*nrkeys]->SetPoint(j, SecTheta(zenith[j],false), s38[j]);
	       else
	          result[i+2*nrkeys]->SetPoint(j, zenith[j], s38[j]);
//	       result[i+3*nrkeys]->SetPoint(j, TMath::Log10(energy[j]), s38[j]);
	       result[i+3*nrkeys]->SetPoint(j, energy[j]/1e+18, s38[j]);

//	       allS38[i]->SetPoint(s38count[i], TMath::Log10(energy[j]), s38[j]);
	       allS38[i]->SetPoint(s38count[i], energy[j]/1e+18, s38[j]);
	       s38count[i]++;
	    }

	    allFitparam[3*i]->SetPoint(fitparamcount[i], midenergy[k], fitparam[1]/fitparam[0]);
	    allFitparam[3*i+1]->SetPoint(fitparamcount[i], midenergy[k], fitparam[2]/fitparam[0]);
	    allFitparam[3*i+2]->SetPoint(fitparamcount[i], midenergy[k], fitparam[3]/fitparam[0]);
	    fitparamcount[i]++;

            mystyle->SetColorScale(result[i+2*nrkeys], i, nrkeys);
	    result[i+2*nrkeys]->SetMarkerStyle(20);
	    result[i+2*nrkeys]->SetMarkerSize(0.8);
            mystyle->SetColorScale(result[i+3*nrkeys], i, nrkeys);
	    result[i+3*nrkeys]->SetMarkerStyle(20);
	    result[i+3*nrkeys]->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(result[i+2*nrkeys], "FD zenith angle (sec#theta)", "S38 (VEM)");
            mystyle->SetAxisTitles(result[i+3*nrkeys], "FD energy [EeV]", "S38 (VEM)");

	    result[i+2*nrkeys]->Draw("AP");
            fitfunc[2] = new TF1("fitfunc2", "[0]", 1., 2.);
	    fitfunc[2]->SetParameter(0,1.);
	    result[i+2*nrkeys]->Fit("fitfunc2");
	    stemp[1] = "./plots/s38_vs_zenithFD_" + stemp[0] + stemp[2] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());
	    stemp[1] = "./plots/s38_vs_zenithFD_" + stemp[0] + stemp[2] + ".C";
	    c1->SaveAs(stemp[1].c_str());

	    c1->SetLogx(kTRUE);
	    c1->SetLogy(kTRUE);
	    result[i+3*nrkeys]->GetXaxis()->SetMoreLogLabels(kTRUE);
	    result[i+3*nrkeys]->Draw("AP");
            fitfunc[3] = new TF1("fitfunc3", "[0]+[1]*x", 3., 100.);
	    fitfunc[3]->SetParameters(1.,1.);
	    result[i+3*nrkeys]->Fit("fitfunc3");
	    stemp[1] = "./plots/energy_vs_s38_" + stemp[0] + stemp[2] + ".pdf";
	    c1->SaveAs(stemp[1].c_str());
	    stemp[1] = "./plots/energy_vs_s38_" + stemp[0] + stemp[2] + ".C";
	    c1->SaveAs(stemp[1].c_str());

	    cout << endl << "Fitting parameters for linear dependence between FD energy and S38 for " << stemp[0] <<  " are:" << endl;
	    cout << "- const = " << fitfunc[3]->GetParameter(0) << " (" << fitfunc[3]->GetParError(0) << ")" << endl;
	    cout << "- slope = " << fitfunc[3]->GetParameter(1) << " (" << fitfunc[3]->GetParError(1) << ")" << endl;
	    cout << endl;
	 }

	 for(int i = 0; i < nrkeys; i++)
	 {
            delete result[i];
            delete result[i+nrkeys];
            delete result[i+2*nrkeys];
            delete result[i+3*nrkeys];
//            delete result[i+4*nrkeys];
	 }

         delete mystyle;
         delete c1;
      }

      mystyle = new RootStyle();
      mystyle->SetBaseStyle();
      c1 = new TCanvas("c1","",1200,900);

      for(int i = 0; i < nrkeys; i++)
      {
         mystyle->SetColorScale(allS38[i], i, nrkeys);
	 allS38[i]->SetMarkerStyle(20);
	 allS38[i]->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(allS38[i], "FD energy [EeV]", "S38 (VEM)");

	 c1->SetLogx(kTRUE);
	 c1->SetLogy(kTRUE);
	 allS38[i]->GetYaxis()->SetRange(4.,500.);
	 allS38[i]->GetYaxis()->SetRangeUser(4.,500.);
	 allS38[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
	 allS38[i]->Draw("AP");
//         fitfunc[4] = new TF1("fitfunc4", "[0]+[1]*x", 3., 100.);
         fitfunc[4] = new TF1("fitfunc4", "TMath::Power(x/[0],1./[1])", 3., 100.);
	 fitfunc[4]->SetParameters(0.1,1.);
	 allS38[i]->Fit("fitfunc4");
	 stemp[1] = "./plots/energy_vs_s38_" + treeName[i] + ".pdf";
	 c1->SaveAs(stemp[1].c_str());
	 stemp[1] = "./plots/energy_vs_s38_" + treeName[i] + ".C";
	 c1->SaveAs(stemp[1].c_str());

	 cout << endl << "Fitting parameters for linear dependence between FD energy and S38:" << endl;
	 cout << "- const = " << fitfunc[4]->GetParameter(0) << " (" << fitfunc[4]->GetParError(0) << ")" << endl;
	 cout << "- slope = " << fitfunc[4]->GetParameter(1) << " (" << fitfunc[4]->GetParError(1) << ")" << endl;
	 cout << endl;

	 c1->SetLogx(kFALSE);
	 c1->SetLogy(kFALSE);

         mystyle->SetColorScale(allFitparam[3*i], i, nrkeys);
	 allFitparam[3*i]->SetMarkerStyle(20);
	 allFitparam[3*i]->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(allFitparam[3*i], "FD energy [log(E/eV)]", "f_{CIC} fitting parameter a");
	 allFitparam[3*i]->Draw("AP");
	 stemp[1] = "./plots/energy_vs_fitparamA_" + treeName[i] + ".pdf";
	 c1->SaveAs(stemp[1].c_str());
	 stemp[1] = "./plots/energy_vs_fitparamA_" + treeName[i] + ".C";
	 c1->SaveAs(stemp[1].c_str());

         mystyle->SetColorScale(allFitparam[3*i+1], i, nrkeys);
	 allFitparam[3*i+1]->SetMarkerStyle(20);
	 allFitparam[3*i+1]->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(allFitparam[3*i+1], "FD energy [log(E/eV)]", "f_{CIC} fitting parameter b");
	 allFitparam[3*i+1]->Draw("AP");
	 stemp[1] = "./plots/energy_vs_fitparamB_" + treeName[i] + ".pdf";
	 c1->SaveAs(stemp[1].c_str());
	 stemp[1] = "./plots/energy_vs_fitparamB_" + treeName[i] + ".C";
	 c1->SaveAs(stemp[1].c_str());

         mystyle->SetColorScale(allFitparam[3*i+2], i, nrkeys);
	 allFitparam[3*i+2]->SetMarkerStyle(20);
	 allFitparam[3*i+2]->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(allFitparam[3*i+2], "FD energy [log(E/eV)]", "f_{CIC} fitting parameter c");
	 allFitparam[3*i+2]->Draw("AP");
	 stemp[1] = "./plots/energy_vs_fitparamC_" + treeName[i] + ".pdf";
	 c1->SaveAs(stemp[1].c_str());
	 stemp[1] = "./plots/energy_vs_fitparamC_" + treeName[i] + ".C";
	 c1->SaveAs(stemp[1].c_str());
      }

      delete c1;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (mvatree_file.root)." << endl;
      return 1;
   }

   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;
   delete analRes;
   for(int i = 0; i < 20; i++)
      delete allS38[i];

   for(int i = 0; i < 20; i++)
      for(int j = 0; j < 3; j++)
         delete allFitparam[j+3*i];*/

   return 0;
}
