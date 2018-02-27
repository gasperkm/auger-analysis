#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
#include "separate_functions.h"
#include "mva_result_read.h"
#include "mva_fit_histogram.h"

using namespace std;

MvaFitHist::MvaFitHist()
{
   stemp = new string[2];
   itemp = new int[2];
   dtemp = new double[3];

   nrbins = 100;
   xlim[0] = -0.25;
   xlim[1] = 1.25;
   treeSelect = -1;

   norm = new double[40];
   nentries = new int[40];

   mystyle = new RootStyle();

   residVal = new double[nrbins];
}

MvaFitHist::MvaFitHist(int bincount, double xlow, double xhigh)
{
   stemp = new string[2];
   itemp = new int[2];
   dtemp = new double[3];

   nrbins = bincount;
   xlim[0] = xlow;
   xlim[1] = xhigh;
   treeSelect = -1;

   norm = new double[40];
   nentries = new int[40];

   mystyle = new RootStyle();

   residVal = new double[nrbins];
}

MvaFitHist::~MvaFitHist()
{
   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   delete[] norm;
   delete[] nentries;

   delete mystyle;
   delete f;

   delete[] residVal;
}

/*void MvaFitHist::fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   signorm = (double)nrdata/(double)nrsig;
   bgdnorm = (double)nrdata/(double)nrback;
   
   sim1 = (TH1F*)simsig->Clone("sim1");
   sim1->Scale(par[0]*signorm);
   sim2 = (TH1F*)simback->Clone("sim2");
   sim2->Scale((1.-par[0])*bgdnorm);
   
   sim = (TH1F*)sim1->Clone("sim");
   sim->Add(sim2);
   
   f = sim->Chi2Test(data, "WW P CHI2");
}*/

double MvaFitHist::fcn(const double *par)
{
   signorm = (double)nrdata/(double)nrsig;
   bgdnorm = (double)nrdata/(double)nrback;
   
   sim1 = (TH1F*)simsig->Clone("sim1");
   sim1->Scale(par[0]*signorm);
   sim2 = (TH1F*)simback->Clone("sim2");
   sim2->Scale((1.-par[0])*bgdnorm);
   
   sim = (TH1F*)sim1->Clone("sim");
   sim->Add(sim2);
   
   return sim->Chi2Test(data, "WU CHI2");
}

void MvaFitHist::PrepareHistograms(string *fname, int *proc, double *step)
{
   fitproc = *proc;
   ratioStep = *step;
   filename = *fname;
   mystyle->SetBaseStyle();

   cout << "Opening file: " << filename << endl;
   f = new TFile(filename.c_str(), "READ");

   float mva;

   // Check how many TTrees we have in the file
   TList *tempkeyslist = (TList*)f->GetListOfKeys();
   nrkeys = (int)f->GetNkeys()/2;
   cout << "Number of trees in file = " << nrkeys << endl;

   // Prepare variables
   for(int i = 0; i < nrkeys; i++)
   {
      stemp[0] = "result" + ToString(i+1);
      result[i] = new TH1F(stemp[0].c_str(), "", nrbins, xlim[0], xlim[1]);
   }
   TTree *tempTree;

   // Go through all TTrees in the root file
   for(int i = 0; i < nrkeys; i++)
   {
      stemp[0] = "TreeS" + ToString(i+1);
      cout << "Getting tree = " << stemp[0] << endl;
      tempTree = (TTree*)f->Get(stemp[0].c_str());

      // Extract the MVA variable
      tempTree->SetBranchAddress("MVA", &mva);
      nentries[i] = tempTree->GetEntries();
      cout << "  Number of entries = " << nentries[i] << endl;

      // Go through all events in the root file for the specific TTree and save MVA variable values to the histogram
      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);
         result[i]->Fill(mva);
      }

      // Scale the histograms by number of all events and get the maximal value
      norm[i] = result[i]->GetMaximum();
      cout << "  Highest value = " << norm[i] << endl;
   }

   // Check the application_result.txt file
   ResultRead *analRes = new ResultRead();
   stemp[0] = RemoveFilename(&filename) + "/application_results.txt";
   itemp[0] = analRes->ReadFile(stemp[0]);
   if(itemp[0] != 1)
   {
      cout << endl << "Available trees:" << endl;
      analRes->PrintVectors(0);
      cout << endl << "Use selected trees for histogram fit (0) or manually select the trees (1)? ";
      cin >> treeSelect;
      while( (treeSelect != 0) && (treeSelect != 1) )
      {
         cout << endl << "Use selected trees for histogram fit (0) or manually select the trees (1)? ";
         cin >> treeSelect;
      }

      if(treeSelect == 0)
      {
         // Saving signal tree position
         analRes->FindPos(1, 0, &trees);
         // Saving background tree position
         analRes->FindPos(2, 0, &trees);
         // Saving data tree position
         analRes->FindPos(3, 0, &trees);
      }
      else if(treeSelect == 1)
      {
         // Saving simulation trees
	 cout << "Enter tree numbers (shown in square brackets) to select simulation trees (use -1 to finish):" << endl;
	 while(itemp[1] != -1)
	 {
            cout << "  Select simulation tree: ";
   	    cin >> itemp[1];
	    if(itemp[1] != -1) 
               trees.push_back(itemp[1]);
	 }
         // Saving the data tree
	 cout << "Enter tree number (shown in square brackets) to select data tree: ";
   	 cin >> itemp[1];
         trees.push_back(itemp[1]);
      }

      cout << "The selected simulation trees are: ";
      for(int i = 0; i < trees.size()-1; i++)
      {
	 if(i > 0)
            cout << ", ";
         cout << trees[i];
      }
      cout << endl << "The selected data tree is: " << trees[trees.size()-1] << endl << endl;
   }

   // Save signal and background distributions
   simsig = (TH1F*)result[trees[0]]->Clone("simsig");
   nrsig = nentries[trees[0]];
   simback = (TH1F*)result[trees[1]]->Clone("simback");
   nrback = nentries[trees[1]];

   // Save data distributions
   data = (TH1F*)result[trees[trees.size()-1]]->Clone("data");
   nrdata = nentries[trees[2]];
   yrange[0] = data->GetMaximum();

   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/fithist";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/fithist/signal_background_data.pdf " + RemoveFilename(&filename) + "/fithist/composition_residuals*";
   system(stemp[0].c_str());
}

// Plot the distributions of signal, background and data
void MvaFitHist::PlotDistributions()
{
   TCanvas *c1 = new TCanvas("c1","",1200,900);

   mystyle->SetHistColor((TH1*)simsig, 1);
   mystyle->SetHistColor((TH1*)simback, 0);
   mystyle->SetHistColor((TH1*)data, 2);
   mystyle->SetAxisTitles((TH1*)simsig, "MVA variable", "Number of events");
   mystyle->SetAxisTitles((TH1*)simback, "MVA variable", "Number of events");
   mystyle->SetAxisTitles((TH1*)data, "MVA variable", "Number of events");

   simsig->Draw();
   simback->Draw("SAME");
   simsig->SetMaximum(1.2*TMath::MaxElement(nrkeys, norm));
   data->Draw("SAME");

   TLegend *legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-0.10, gPad->GetLeftMargin()+0.26, 1-gPad->GetTopMargin());
   int c_LegendFill = TColor::GetColor("#ffff66");
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_LegendFill);
   legend->AddEntry(simsig, "Signal events", "f");
   legend->AddEntry(simback, "Background events", "f");
   legend->AddEntry(data, "Data events", "f");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("SAME");

   stemp[0] = RemoveFilename(&filename) + "/fithist/signal_background_data.pdf";
   c1->SaveAs(stemp[0].c_str());

   delete c1;
   delete legend;
}

// Plot data, normalized signal+background and normalized residuals
void MvaFitHist::StartFitting()
{
   int nrsteps = -1;

   // Plot a range of ratios
   if(fitproc == 0)
   { 
      nrsteps = (1/ratioStep)+1;
      for(int i = 0; i < nrsteps; i++)
      {
         dtemp[0] = (double)i*ratioStep;
         PlotSumResiduals(dtemp);
      }
   }
   // Find the best fit with TMinuit
   else if(fitproc == 1)
   {
      ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
      minim->SetMaxFunctionCalls(1000000);
      minim->SetMaxIterations(10000);
      minim->SetTolerance(0.01);
      minim->SetPrintLevel(1);

      ROOT::Math::Functor fmin(this,&MvaFitHist::fcn,1);
      srand48(time(NULL));
      double step = 0.01, var = drand48();
      minim->SetFunction(fmin);
      minim->SetVariable(0, "frac", var, step);
      minim->Minimize();

      double *bestfrac = (double*)minim->X();
      double *bestfracerr = (double*)minim->Errors();
      dtemp[0] = bestfrac[0];
      dtemp[1] = bestfracerr[0];
      PlotSumResiduals(dtemp);
   }
}

// Plot data, normalized signal+background and normalized residuals
void MvaFitHist::PlotSumResiduals(double *sigFrac)
{
   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TPad *upperpad = new TPad("upperpad", "upperpad", 0.004, 0.490, 0.996, 0.996);
   TPad *lowerpad = new TPad("lowerpad", "lowerpad", 0.004, 0.000, 0.996, 0.490);
   upperpad->Draw();
   lowerpad->Draw();

   // Plot the data distribution
   upperpad->cd();
   mystyle->SetHistColor((TH1*)data, 2);
   mystyle->SetAxisTitles((TH1*)data, "MVA variable", "Number of events");
   data->GetXaxis()->SetTitleOffset(2.0);
   data->Draw();

   signorm = (double)nrdata/(double)nrsig;
   bgdnorm = (double)nrdata/(double)nrback;
   
   sim1 = (TH1F*)simsig->Clone("sim1");
   sim1->Scale(sigFrac[0]*signorm);
   sim2 = (TH1F*)simback->Clone("sim2");
   sim2->Scale((1.-sigFrac[0])*bgdnorm);
   
   sim = (TH1F*)sim1->Clone("sim");
   sim->Add(sim2);

   mystyle->SetHistColor((TH1*)sim, 0);
   mystyle->SetAxisTitles((TH1*)sim, "MVA variable", "Number of events");
   yrange[1] = sim->GetMaximum();
   data->SetMaximum(1.1*TMath::MaxElement(2, yrange));
   sim->Draw("SAME");

   chiProb = sim->Chi2TestX(data, chiVal, chiNdf, chiGood, "WU", residVal);
   cout << "Chi2 = " << chiVal << ", prob = " << chiProb << ", NDF = " << chiNdf << ", igood = " << chiGood << endl;

   TH1F *resid = new TH1F("resid", "", nrbins, xlim[0], xlim[1]);
   yresidrange[0] = 1.e+24;
   yresidrange[1] = -1.e+24;
   for(int i = 1; i <= nrbins; i++)
   {
      if( ((double)sim->GetBinContent(i) + (double)data->GetBinContent(i) == 0) )
         resid->SetBinContent(i, 0.);
      else
      {
         resid->SetBinContent(i, residVal[i-1]);

         if(residVal[i-1] < yresidrange[0])
            yresidrange[0] = residVal[i-1];
         if(residVal[i-1] > yresidrange[1])
            yresidrange[1] = residVal[i-1];
      }
   }

   TLatex chiText;
   chiText.SetTextAlign(21);
   stemp[0] = "#chi^{2}/NDF = " + ToString(chiVal, 3) + "/" + ToString(chiNdf) + " = " + ToString(chiVal/(double)chiNdf, 3);
   chiText.DrawLatex(0.5, TMath::MaxElement(2, yrange), stemp[0].c_str());
   if(fitproc == 0)
   {
      stemp[0] = "r_{p} = " + ToString(sigFrac[0], 3);
      chiText.DrawLatex(0.5, 0.94*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = "r_{p} = " + ToString(sigFrac[0], 3) + " #pm " + ToString(sigFrac[1],3);
      chiText.DrawLatex(0.5, 0.94*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }

   lowerpad->cd();
   mystyle->SetHistColor((TH1*)resid, 3);
   mystyle->SetAxisTitles((TH1*)resid, "", "Normalized residuals");
   yresidrange[0] = TMath::Abs(yresidrange[0]);
   yresidrange[1] = TMath::Abs(yresidrange[1]);
   resid->GetYaxis()->SetRangeUser(-1.1*TMath::MaxElement(2, yresidrange), 1.1*TMath::MaxElement(2, yresidrange));
   resid->Draw();

   c1->Update();
   if(fitproc == 0)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/composition_residuals_" + ToString(100*sigFrac[0], 0) + ".pdf";
      c1->SaveAs(stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/composition_residuals_minuit2.pdf";
      c1->SaveAs(stemp[0].c_str());
   }
   
   delete c1;
}

int main(int argc, char **argv)
{
   // Important variables
   int fitproc = -1;
   double ratioStep = -1;

   gSystem->Load("libTree.so");

   if(argc > 1)
   {
      int *itemp = new int;
      double *dtemp = new double;
      string *stemp = new string;

      cerr << "Perform histogram fit on a range of ratios (0) or perform a minimization (1)? ";
      cin >> *itemp;
      
      if(*itemp == 0)
      {
         cerr << "What should be the step size for ratios (between 0 and 1)? ";
	 cin >> *dtemp;

	 while( (*dtemp > 1) || (*dtemp <= 0) )
	 {
            cerr << "What should be the step size for ratios (between 0 and 1)? ";
	    cin >> *dtemp;
	 }
      }
      else if( (*itemp != 0) && (*itemp != 1) )
      {
         cerr << "Error! Wrong fitting procedure." << endl;
	 return 1;
      }

      for(int i = 0; i < argc-1; i++)
      {
         *stemp = string(argv[i+1]);

         MvaFitHist *fithist = new MvaFitHist();
         fithist->PrepareHistograms(stemp, itemp, dtemp);
         fithist->PlotDistributions();
	 fithist->StartFitting();
         delete fithist;
      }

      delete itemp;
      delete dtemp;
      delete stemp;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments." << endl;
      return 1;
   }

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
