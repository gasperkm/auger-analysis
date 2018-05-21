#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>

#include "TSystem.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
//#include "mvaefficiency.h"
//#include "separate_functions.h"

using namespace std;

int main(int argc, char **argv)
{
   gSystem->Load("libTree");

   string *stemp = new string[2];
   string *inname = new string;
   *inname = "temporary_mvatree_file.root";
   string *trname = new string;
   *trname = "transformation_stats.dat";
   double *dtemp = new double[4];

   if(argc > 2)
   {
      *inname = string(argv[1]);
      *trname = string(argv[2]);
   }

   string *dirname = new string;
   *dirname = "/home/gkukec/Gasper/github/auger-analysis/dbg/error_sigma";
   cout << "Directory name = " << *dirname << endl;

   TFile *ifile = TFile::Open(inname->c_str(), "READ");

   vector<string> *obser = new vector<string>;
   obser->push_back("xmax");
   obser->push_back("shwsize");
   obser->push_back("risetimerecalc");

   float *obsVals = new float[3*obser->size()];

   TTree *datatree = (TTree*)ifile->Get("TreeS5");
   cout << "Observables:" << endl;
   for(int i = 0; i < obser->size(); i++)
   {
      cout << "- " << obser->at(i) << endl;
      stemp[0] = obser->at(i);
      datatree->SetBranchAddress(stemp[0].c_str(), &obsVals[3*i]);
      stemp[0] = obser->at(i) + "_neg";
      datatree->SetBranchAddress(stemp[0].c_str(), &obsVals[3*i+1]);
      stemp[0] = obser->at(i) + "_pos";
      datatree->SetBranchAddress(stemp[0].c_str(), &obsVals[3*i+2]);
   }

   double *obslimit;
   obslimit = new double[2*obser->size()];
   char *ctemp;
   ctemp = new char[1024];
   ifstream *fstats = new ifstream;
   fstats->open(trname->c_str(), ifstream::in);
   if(fstats->is_open())
   {
      if(fstats->peek() == '#')
         fstats->getline(ctemp, 1024, '\n');

      for(int i = 0; i < obser->size(); i++)
      {
         *fstats >> dtemp[0] >> dtemp[1] >> dtemp[2] >> dtemp[3];
	 obslimit[2*i] = dtemp[2];
	 obslimit[2*i+1] = dtemp[3];
	 cout << "Observable limits: " << obslimit[2*i] << "\t" << obslimit[2*i+1] << endl;
      }
   }

   fstats->close();
   delete fstats;
   delete[] ctemp;

   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TH1F *hist, *histNeg, *histPos;

   dtemp[3] = 0;

   for(int i = 0; i < obser->size(); i++)
   {
      hist = new TH1F("h1","",100,-1.5,1.5);
      histNeg = new TH1F("h2","",100,0.,1.);
      histPos = new TH1F("h3","",100,0.,1.);

      for(int j = 0; j < datatree->GetEntries(); j++)
      {
         datatree->GetEntry(j);

	 dtemp[0] = 2*(obsVals[3*i] - obslimit[2*i])/(obslimit[2*i+1] - obslimit[2*i]) - 1;
	 dtemp[1] = dtemp[0] - (2*((obsVals[3*i]-obsVals[3*i+1]) - obslimit[2*i])/(obslimit[2*i+1] - obslimit[2*i]) - 1);
	 dtemp[2] = 2*((obsVals[3*i]+obsVals[3*i+2]) - obslimit[2*i])/(obslimit[2*i+1] - obslimit[2*i]) - 1 - dtemp[0];

         hist->Fill(dtemp[0]);
         histNeg->Fill(dtemp[1]);
         histPos->Fill(dtemp[2]);
      }

      hist->Draw();
      stemp[0] = *dirname + "/output_" + obser->at(i) + "_mean.pdf";
      c1->SaveAs(stemp[0].c_str());

      histNeg->Draw();
      stemp[0] = *dirname + "/output_" + obser->at(i) + "_neg.pdf";
      c1->SaveAs(stemp[0].c_str());

      histPos->Draw();
      stemp[0] = *dirname + "/output_" + obser->at(i) + "_pos.pdf";
      c1->SaveAs(stemp[0].c_str());

      cout << "0: Mean = " << hist->GetMean() << ", Sigma = " << hist->GetRMS() << endl;
      cout << "-: Mean = " << histNeg->GetMean() << ", Sigma = " << histNeg->GetRMS() << endl;
      cout << "+: Mean = " << histPos->GetMean() << ", Sigma = " << histPos->GetRMS() << endl;

      dtemp[3] += TMath::Power(histNeg->GetRMS(), 2);

      delete hist;
      delete histNeg;
      delete histPos;
   }

   dtemp[3] = TMath::Sqrt(dtemp[3]);
   cout << endl << "Final combined error = " << dtemp[3] << endl;

   ifile->Close();

   delete c1;
   delete[] obslimit;
   delete[] obsVals;
   delete[] stemp;
   delete[] dtemp;
   delete inname;
   delete trname;
   delete dirname;
   delete obser;

   return 0;
}
