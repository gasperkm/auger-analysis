/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAMulticlassApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained multiclass      *
 * classifiers within an analysis module                                          *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TColor.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TH1F.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace TMVA;

void SetColor(TH1F *tr, int cur, int nrbins)
{
   int ci = 1738+cur;
   double colormix = (double)cur/(double)(nrbins-1.);
   TColor *color = new TColor(ci, colormix, 0, 1.-colormix, "", 1);
   tr->SetLineColor(ci);
   tr->SetLineWidth(2);
   tr->SetFillColor(ci);
   tr->SetFillStyle(3004+cur);
}

int main(int argc, char** argv )
{
   TMVA::Tools::Instance();
   
   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 0;
   Use["BDTG"]            = 1;
   Use["FDA_GA"]          = 0;
   Use["PDEFoam"]         = 0;
   //---------------------------------------------------------------
  
   std::cout << std::endl;
   std::cout << "==> Start TMVAMulticlassApplication" << std::endl; 

   if (argc>1) {
      for (std::map<std::string,int>::iterator it = Use.begin();
           it != Use.end(); it++) {
         it->second = 0;
      }
   }
   for (int i=1; i<argc; i++) {
      std::string regMethod(argv[i]);
      if (Use.find(regMethod) == Use.end()) {
         std::cout << "Method " << regMethod << " not known in TMVA under this name. Please try one of:" << std::endl;
         for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
         std::cout << std::endl;
         return 1;
      }
      Use[regMethod] = kTRUE;
   }

   
   // create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to 
   // those given in the weight file(s) that you use
   Float_t var1, var2, var3, var4;
   reader->AddVariable( "var1", &var1 );
   reader->AddVariable( "var2", &var2 );
   reader->AddVariable( "var3", &var3 );
   reader->AddVariable( "var4", &var4 );

   // book the MVA methods
   TString dir    = "weights/";
   TString prefix = "TMVAMulticlass";
   
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
         TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }

   // book output histograms
   UInt_t nbin = 100;
   TH1F *histMLP_signal(0), *histBDTG_signal(0), *histFDAGA_signal(0), *histPDEFoam_signal(0);
   TH1F *histBDTG_bg0(0), *histBDTG_bg1(0), *histBDTG_bg2(0);
   if (Use["MLP"])    
      histMLP_signal    = new TH1F( "MVA_MLP_signal",    "MVA_MLP_signal",    nbin, 0., 1.1 );
   if (Use["BDTG"])
   {
      histBDTG_signal  = new TH1F( "MVA_BDTG_signal",   "MVA_BDTG_signal",   nbin, -0.1, 1.2 );
      histBDTG_bg0  = new TH1F( "MVA_BDTG_bg0",   "MVA_BDTG_bg0",   nbin, -0.1, 1.2 );
      histBDTG_bg1  = new TH1F( "MVA_BDTG_bg1",   "MVA_BDTG_bg1",   nbin, -0.1, 1.2 );
      histBDTG_bg2  = new TH1F( "MVA_BDTG_bg2",   "MVA_BDTG_bg2",   nbin, -0.1, 1.2 );
   }
   if (Use["FDA_GA"])
      histFDAGA_signal = new TH1F( "MVA_FDA_GA_signal", "MVA_FDA_GA_signal", nbin, 0., 1.1 );
   if (Use["PDEFoam"])
      histPDEFoam_signal = new TH1F( "MVA_PDEFoam_signal", "MVA_PDEFoam_signal", nbin, 0., 1.1 );


   TFile *input(0);
   TString fname = "./tmva_example_multiple_background.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   if (!input) {
      std::cout << "ERROR: could not open data file, please generate example data first!" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAMulticlassApp : Using input file: " << input->GetName() << std::endl;
   
   // prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   
   TCanvas *c1 = new TCanvas("c1","",1200,900);
   c1->SetGrid();
   c1->SetLogy();
   double ymax = -1;
  
   TTree* theTree;
   for(int t = 0; t < 4; t++)
   {
      if(t == 0)
         theTree = (TTree*)input->Get("TreeS");
      else if(t == 1)
         theTree = (TTree*)input->Get("TreeB0");
      else if(t == 2)
         theTree = (TTree*)input->Get("TreeB1");
      else if(t == 3)
         theTree = (TTree*)input->Get("TreeB2");

      std::cout << "--- Select signal sample" << std::endl;
      theTree->SetBranchAddress( "var1", &var1 );
      theTree->SetBranchAddress( "var2", &var2 );
      theTree->SetBranchAddress( "var3", &var3 );
      theTree->SetBranchAddress( "var4", &var4 );

      std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
      TStopwatch sw;
      sw.Start();

      for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
         if (ievt%1000 == 0){
            std::cout << "--- ... Processing event: " << ievt << std::endl;
         }
         
         theTree->GetEntry(ievt);
         if (Use["MLP"])
            histMLP_signal->Fill((reader->EvaluateMulticlass( "MLP method" ))[0]);
         if (Use["BDTG"])
	 {
            if(t == 0)
               histBDTG_signal->Fill((reader->EvaluateMulticlass( "BDTG method" ))[0]);
	    else if(t == 1)
               histBDTG_bg0->Fill((reader->EvaluateMulticlass( "BDTG method" ))[0]);
	    else if(t == 2)
               histBDTG_bg1->Fill((reader->EvaluateMulticlass( "BDTG method" ))[0]);
	    else if(t == 3)
               histBDTG_bg2->Fill((reader->EvaluateMulticlass( "BDTG method" ))[0]);
	 }
         if (Use["FDA_GA"])
            histFDAGA_signal->Fill((reader->EvaluateMulticlass( "FDA_GA method" ))[0]);
         if (Use["PDEFoam"])
            histPDEFoam_signal->Fill((reader->EvaluateMulticlass( "PDEFoam method" ))[0]);
         
      }
      
      // get elapsed time
      sw.Stop();
      std::cout << "--- End of event loop: "; sw.Print();
      
      TFile *target;
      if(t == 0)
         target  = new TFile( "TMVAMulticlassApp_histS.root","RECREATE" );
      else if(t == 1)
         target  = new TFile( "TMVAMulticlassApp_histBG0.root","RECREATE" );
      else if(t == 2)
         target  = new TFile( "TMVAMulticlassApp_histBG1.root","RECREATE" );
      else if(t == 3)
         target  = new TFile( "TMVAMulticlassApp_histBG2.root","RECREATE" );

      if (Use["MLP"])
         histMLP_signal->Write();
      if (Use["BDTG"])
      {
         if(t == 0)
	 {
            histBDTG_signal->Write(); 
	    SetColor(histBDTG_signal, t, 4);
            histBDTG_signal->Draw(); 

	    if(histBDTG_signal->GetMaximum() > ymax)
               ymax = histBDTG_signal->GetMaximum();
	 }
	 else if(t == 1)
	 {
            histBDTG_bg0->Write(); 
	    SetColor(histBDTG_bg0, t, 4);
            histBDTG_bg0->Draw("SAME"); 

	    if(histBDTG_bg0->GetMaximum() > ymax)
               ymax = histBDTG_bg0->GetMaximum();
	 }
	 else if(t == 2)
	 {
            histBDTG_bg1->Write(); 
	    SetColor(histBDTG_bg1, t, 4);
            histBDTG_bg1->Draw("SAME"); 

	    if(histBDTG_bg1->GetMaximum() > ymax)
               ymax = histBDTG_bg1->GetMaximum();
	 }
	 else if(t == 3)
	 {
            histBDTG_bg2->Write(); 
	    SetColor(histBDTG_bg2, t, 4);
            histBDTG_bg2->Draw("SAME"); 

	    if(histBDTG_bg2->GetMaximum() > ymax)
               ymax = histBDTG_bg2->GetMaximum();
	 }
      }
      if (Use["FDA_GA"])
         histFDAGA_signal->Write();
      if (Use["PDEFoam"])
         histPDEFoam_signal->Write();

      target->Close();
      std::cout << "--- Created root file: \"TMVMulticlassApp.root\" containing the MVA output histograms" << std::endl;
   }

   histBDTG_signal->GetYaxis()->SetRangeUser(0.1,(1.05*ymax));
   c1->SaveAs("mva_bdtg.pdf");

   delete reader;
   
   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}
