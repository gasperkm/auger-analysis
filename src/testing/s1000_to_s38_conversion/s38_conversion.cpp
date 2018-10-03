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
//#include "mva_fit_histogram.h"

using namespace std;

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   string inname;
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
         delete allFitparam[j+3*i];

   return 0;
}
