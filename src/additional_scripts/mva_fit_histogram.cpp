#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "mva_fit_histogram.h"

using namespace std;

MvaFitHist::MvaFitHist()
{
   stemp = new string[2];
   itemp = new int[3];
   dtemp = new double[3];

   nrbins = -1;
   cout << endl;
   while(nrbins < 4)
   {
      cout << "Select number of bins (min of 4): ";
      cin >> nrbins;
   }

   cout << endl;
   PrintMethods();
   cout << "Select the used MVA type: ";
   cin >> stemp[0];

   cout << "Select constraint type (0 = function, 1 = one of the fractions): ";
   cin >> constraint;

   xlim[0] = GetMethodMin(stemp[0]);
   xlim[1] = GetMethodMax(stemp[0]);
   cout << "MVA limits (" << stemp[0] << "): " << xlim[0] << ", " << xlim[1] << endl;

   treeSelect = -1;

   norm = new double[40];
   nentries = new int[40];

   mystyle = new RootStyle();

   residVal = new double[nrbins];
   nrsim = new int[40];
   simnorm = new double[40];

   minstep = new double[40];
   minvar = new double[40];

   midLna = new double[2];
   midComposition = new vector<double>;
   midCompositionErr = new vector<double>;
   chi2value = new double;
   pvalue = new double;
   ndfvalue = new double;

   ptype = new PrimPart();

   appvar = -1;
}

/*MvaFitHist::MvaFitHist(int bincount, double xlow, double xhigh)
{
   stemp = new string[2];
   itemp = new int[3];
   dtemp = new double[3];

   nrbins = bincount;
   xlim[0] = xlow;
   xlim[1] = xhigh;
   treeSelect = -1;

   norm = new double[40];
   nentries = new int[40];

   mystyle = new RootStyle();

   residVal = new double[nrbins];

   minstep = new double[40];
   minvar = new double[40];

   appvar = -1;
}*/

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
   delete[] nrsim;
   delete[] simnorm;

   delete[] minstep;
   delete[] minvar;

   delete[] midLna;
   delete midComposition;
   delete midCompositionErr;
   delete chi2value;
   delete pvalue;
   delete ndfvalue;

   delete ptype;
}

double MvaFitHist::fcnConstrainedFunc(const double *par)
{
   dtemp[0] = 1;
   // Scale all simulation histograms with the nrdata/nrsim[i] value
   for(int i = 0; i < nrparam; i++)
   {
//      simnorm[i] = (double)nrdata/(double)nrsim[i];
      stemp[0] = "sims" + ToString(i);
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      dtemp[0] -= par[i];
      sims[i]->Scale(par[i]);
   }

   // Sum all simulation histograms together
   sim1 = (TH1F*)sims[0]->Clone("sim1");
   for(int i = 1; i < nrparam; i++)
      sim1->Add(sims[i]);

   // Define constraint histogram: sum(f_i) = 1
   sim2 = (TH1F*)sim1->Clone("sim2");

   for(int i = 0; i <= nrbins+1; i++)
   {
      if(i == 0)
         sim2->SetBinContent(i, 0);
      else if(i == nrbins+1)
         sim2->SetBinContent(i, 0);
      else
      {
         if(sim1->GetBinContent(i) == 0)
            sim2->SetBinContent(i, 0);
         else
            sim2->SetBinContent(i, dtemp[0]);
      }
   }

   // Apply constraint
   sim = (TH1F*)sim1->Clone("sim");
//   sim->Add(sim2);

   return dtemp[0];
}

double MvaFitHist::fcnFunc(const double *par)
{
   dtemp[0] = 1;
   // Scale all simulation histograms with the nrdata/nrsim[i] value
   for(int i = 0; i < nrparam; i++)
   {
//      simnorm[i] = (double)nrdata/(double)nrsim[i];
      stemp[0] = "sims" + ToString(i);
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      if(i == nrparam-1)
         sims[i]->Scale(dtemp[0]);
      else
      {
         dtemp[0] -= par[i];
         sims[i]->Scale(par[i]);
      }

/*      if(i == nrparam-1)
         sims[i]->Scale(dtemp[0]*simnorm[i]);
      else
      {
         dtemp[0] -= par[i];
         sims[i]->Scale(par[i]*simnorm[i]);
      }*/
   }

   // Sum all simulation histograms together
   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 1; i < nrparam; i++)
      sim->Add(sims[i]);
   
/*   // Go over the simulation distribution and check all the bins
   cout << "sim" << endl;
   dtemp[0] = 0;
   for(int j = 1; j < nrbins+1; j++)
   {
      dtemp[0] += sim->GetBinContent(j);
      cout << j << "\t" << sim->GetBinContent(j) << "\t" << sim->GetBinError(j) << endl;
   }
   cout << "Histogram integral = " << dtemp[0] << endl;*/

   return 0;
}

double MvaFitHist::fcn(const double *par)
{
//   cout << "New iteration -----------------------" << endl;
   if(constraint == 0)
   {
      dtemp[0] = fcnConstrainedFunc(par);
//      cout << "dtemp[0] = " << dtemp[0] << endl;
/*      cout << "dtemp[0]^2 = " << TMath::Power(dtemp[0],2) << endl;
      for(int i = 0; i < nrparam; i++)
         cout << "par[" << i << "] = " << par[i] << endl;
      cout << "chi2test = " << sim->Chi2Test(data, "WW CHI2") << endl;*/

/*      dtemp[0] = 1;
      // Scale all simulation histograms with the nrdata/nrsim[i] value
      for(int i = 0; i < nrparam; i++)
      {
         simnorm[i] = (double)nrdata/(double)nrsim[i];
         stemp[0] = "sims" + ToString(i);
         sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

         dtemp[0] -= par[i];
         sims[i]->Scale(par[i]*simnorm[i]);
      }

      // Sum all simulation histograms together
      sim1 = (TH1F*)sims[0]->Clone("sim1");
      for(int i = 1; i < nrparam; i++)
         sim1->Add(sims[i]);

      // Define constraint histogram: sum(f_i) = 1
      sim2 = (TH1F*)sim1->Clone("sim2");

      for(int i = 0; i <= nrbins+1; i++)
      {
         if(i == 0)
            sim2->SetBinContent(i, 0);
         else if(i == nrbins+1)
            sim2->SetBinContent(i, 0);
         else
         {
            if(sim1->GetBinContent(i) == 0)
               sim2->SetBinContent(i, 0);
            else
               sim2->SetBinContent(i, dtemp[0]);
         }
      }

      // Apply constraint
      sim = (TH1F*)sim1->Clone("sim");
      sim->Add(sim2);*/
   }
   else
   {
      dtemp[0] = fcnFunc(par);
//      cout << "dtemp[0] = " << dtemp[0] << endl;
/*      cout << "dtemp[0]^2 = " << TMath::Power(dtemp[0],2) << endl;
      for(int i = 0; i < nrparam; i++)
      {
         if(i < nrparam-1)
            cout << "par[" << i << "] = " << par[i] << endl;
	 else
            cout << "par[" << i << "] = " << dtemp[0] << endl;
      }
      cout << "chi2test = " << sim->Chi2Test(data, "WW CHI2") << endl;*/

/*      dtemp[0] = 1;
      // Scale all simulation histograms with the nrdata/nrsim[i] value
      for(int i = 0; i < nrparam; i++)
      {
         simnorm[i] = (double)nrdata/(double)nrsim[i];
         stemp[0] = "sims" + ToString(i);
         sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());


	 if(i == nrparam-1)
            sims[i]->Scale(dtemp[0]*simnorm[i]);
	 else
	 {
            dtemp[0] -= par[i];
            sims[i]->Scale(par[i]*simnorm[i]);
	 }
      }

      // Sum all simulation histograms together
      sim = (TH1F*)sims[0]->Clone("sim");
      for(int i = 1; i < nrparam; i++)
         sim->Add(sims[i]);*/
   }

//   return sim->Chi2Test(data, "WW CHI2") + dtemp[0];
   return sim->Chi2Test(data, "WW CHI2") + TMath::Power(dtemp[0],2);
}

void MvaFitHist::PrepareHistograms(int run, string *fname, int *proc, double *step)
{
   fitproc = proc[0];
   setopt = proc[1];
   ratioStep = *step;
   filename = *fname;
   mystyle->SetBaseStyle();

   cout << endl << "Opening file: " << filename << endl;
   f = new TFile(filename.c_str(), "READ");

   float mva;

   // Check how many TTrees we have in the file
   TList *tempkeyslist = (TList*)f->GetListOfKeys();
   nrkeys = (int)f->GetNkeys()/2;
   cout << "Number of trees in file = " << nrkeys << endl;

   // Prepare variables (histograms to hold values from each tree)
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

      vector<double> *binnentries = new vector<double>;

      // Go through all events in the root file for the specific TTree and save MVA variable values to the histogram
      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);
         result[i]->Fill(mva);
      }

      // Calculate errors
//      cout << "Err1:" << endl;
      for(int k = 1; k <= nrbins; k++)
      {
	 dtemp[0] = result[i]->GetBinContent(k);
	 if(dtemp[0] == 0)
            binnentries->push_back(1);
	 else
            binnentries->push_back(TMath::Sqrt(dtemp[0])/dtemp[0]);
//         cout << k << "\t" << result[i]->GetBinContent(k) << "\trel. err. = " << binnentries->at(k-1) << endl;
      }

      // Scale the histograms by number of all events (to get probability density functions) and get the maximal value
      result[i]->Scale(1./nentries[i]);
      norm[i] = result[i]->GetMaximum();
      cout << "  Highest value = " << norm[i] << endl;

//      cout << "Err2:" << endl;
      dtemp[1] = 0;
      for(int k = 1; k <= nrbins; k++)
      {
	 dtemp[0] = result[i]->GetBinContent(k);
	 if(dtemp[0] == 0)
            result[i]->SetBinError(k, binnentries->at(k-1));
	 else
	 {
            if((result[i]->GetBinContent(k))*binnentries->at(k-1) > dtemp[1])
               dtemp[1] = (result[i]->GetBinContent(k))*binnentries->at(k-1);
            result[i]->SetBinError(k, (result[i]->GetBinContent(k))*binnentries->at(k-1));
	 }
//         cout << k << "\t" << result[i]->GetBinContent(k) << "\t" << result[i]->GetBinError(k) << endl;
      }

      cout << "Maximum value " << norm[i] << " increased by " << dtemp[1] << endl;
      norm[i] += dtemp[1];

      delete binnentries;
   }

   ResultRead *analRes = new ResultRead();
   stemp[0] = RemoveFilename(&filename) + "/application_results.txt";
   itemp[0] = analRes->ReadFile(stemp[0]);

   if(itemp[0] != 1)
   {
      // Check the application_result.txt file
      if( ((setopt == 0) && (run == 0)) || (setopt == 1) )
      {
         if(!trees.empty())
            trees.erase(trees.begin(), trees.end());
         if(!treeNames.empty())
            treeNames.erase(treeNames.begin(), treeNames.end());

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

            for(int i = 0; i < trees.size(); i++)
               treeNames.push_back(analRes->GetTreeName(trees[i]));
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
               {
                  trees.push_back(itemp[1]);
                  treeNames.push_back(analRes->GetTreeName(itemp[1]));
               }
            }
            // Saving the data tree
            cout << "Enter tree number (shown in square brackets) to select data tree: ";
      	    cin >> itemp[1];
            trees.push_back(itemp[1]);
            treeNames.push_back(analRes->GetTreeName(itemp[1]));
         }
      }

      cout << "The selected simulation trees are: ";
      for(int i = 0; i < trees.size()-1; i++)
      {
         if(i > 0)
            cout << ", ";
         cout << trees[i] << " (" << treeNames[i] << ")";
      }
      cout << endl << "The selected data tree is: " << trees[trees.size()-1] << " (" << treeNames[trees.size()-1] << ")" << endl << endl;

      // Getting the signal fraction approximation for the data
      appvar = analRes->GetFraction(2, -1);
      mvacut = analRes->GetMvaCut(0);

      cout << "Signal fraction approximation for data: " << appvar << endl;
      cout << "Mva cut value: " << mvacut << endl;
   }

   delete analRes;

   // Save signal and background distributions
   simsig = (TH1F*)result[trees[0]]->Clone("simsig");
   nrsig = nentries[trees[0]];
   simback = (TH1F*)result[trees[trees.size()-2]]->Clone("simback");
   nrback = nentries[trees[1]];
   for(int i = 0; i < trees.size()-1; i++)
   {
      stemp[0] = "simhist" + ToString(i);
      simhist[i] = (TH1F*)result[trees[i]]->Clone(stemp[0].c_str());
//      simhist[i]->Sumw2(kTRUE);
      nrsim[i] = nentries[trees[i]];

      // Go over the simulation distribution and check all the bins
      cout << "Simhist " << i << ":" << endl;
      dtemp[0] = 0;
      for(int j = 1; j < nrbins+1; j++)
      {
	 dtemp[0] += simhist[i]->GetBinContent(j);
         cout << j << "\t" << simhist[i]->GetBinContent(j) << "\t" << simhist[i]->GetBinError(j) << endl;
      }
      cout << "Histogram integral = " << dtemp[0] << endl;
   }

   // Check if all simulation histograms have zeroes at the same place (if not, remove the large error)
   bool *iszero = new bool;
   *iszero = false;
   dtemp[0] = 0;
   for(int j = 1; j < nrbins+1; j++)
   {
      for(int i = 0; i < trees.size()-1; i++)
      {
         dtemp[0] += simhist[i]->GetBinContent(j);
	 if(simhist[i]->GetBinContent(j) == 0)
            (*iszero) = true;
      }

      if( (dtemp[0] != 0) && (*iszero) )
      {
         for(int i = 0; i < trees.size()-1; i++)
	 {
	    if(simhist[i]->GetBinContent(j) == 0)
               simhist[i]->SetBinError(j, 0);
	 }
      }
   }

   for(int i = 0; i < trees.size()-1; i++)
   {
      // Go over the simulation distribution and check all the bins
      cout << "Simhist " << i << " (edit):" << endl;
      dtemp[0] = 0;
      for(int j = 1; j < nrbins+1; j++)
      {
	 dtemp[0] += simhist[i]->GetBinContent(j);
         cout << j << "\t" << simhist[i]->GetBinContent(j) << "\t" << simhist[i]->GetBinError(j) << endl;
      }
      cout << "Histogram integral = " << dtemp[0] << endl;
   }

   // Set the number of parameters that will be needed to fit (number of simulation trees - 1)
   nrparam = trees.size() - 1;
   cout << "Number of parameters for fitting: " << nrparam << endl;

   // Save data distributions
   data = (TH1F*)result[trees[trees.size()-1]]->Clone("data");
   nrdata = nentries[trees[trees.size()-1]];
//   yrange[0] = data->GetMaximum();
   yrange[0] = norm[trees.size()-1];

/*   // Go over the data distribution and check all the bins
   cout << "Data:" << endl;
   for(int j = 1; j < nrbins+1; j++)
      cout << j << "\t" << data->GetBinContent(j) << endl;*/

   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/fithist";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/fithist/signal_background_data.pdf " + RemoveFilename(&filename) + "/fithist/composition_residuals*";
   system(stemp[0].c_str());

   // Zero empty bin errors just for plotting purposes
   newdata = (TH1F*)data->Clone("newdata");
   newsim = (TH1F*)simsig->Clone("newsim");
   newback = (TH1F*)simback->Clone("newback");
   for(int j = 1; j < nrbins+1; j++)
   {
      if(data->GetBinContent(j) == 0)
         newdata->SetBinError(j, 0);

      if(simback->GetBinContent(j) == 0)
         newback->SetBinError(j, 0);

      if(simsig->GetBinContent(j) == 0)
         newsim->SetBinError(j, 0);
   }
}

// Plot the distributions of signal, background and data
void MvaFitHist::PlotDistributions()
{
   TLegend *legend;
   TCanvas *c1 = new TCanvas("c1","",1200,900);

   mystyle->SetHistColor((TH1*)newsim, 1);
   mystyle->SetHistColor((TH1*)newback, 0);
   mystyle->SetHistColor((TH1*)newdata, 2);
   mystyle->SetAxisTitles((TH1*)newsim, "MVA variable", "Probability density");
   mystyle->SetAxisTitles((TH1*)newback, "MVA variable", "Probability density");
   mystyle->SetAxisTitles((TH1*)newdata, "MVA variable", "Probability density");

   newsim->Draw();
   newback->Draw("SAME");
   cout << "Histogram maximum value = " << TMath::MaxElement(nrkeys, norm) << endl;
   newsim->SetMaximum(1.2*TMath::MaxElement(nrkeys, norm));
   newdata->Draw("SAME");

   legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(3)), gPad->GetLeftMargin()+0.32, 1-gPad->GetTopMargin());
   int c_LegendFill = TColor::GetColor("#ffff66");
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_LegendFill);
   legend->AddEntry(newsim, "Signal events", "f");
   legend->AddEntry(newback, "Background events", "f");
   legend->AddEntry(newdata, "Data events", "f");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("SAME");

   stemp[0] = RemoveFilename(&filename) + "/fithist/signal_background_data.pdf";
   c1->SaveAs(stemp[0].c_str());

   c1->SetLogy();
   newsim->SetMaximum((TMath::Power(10,1.1))*TMath::MaxElement(nrkeys, norm));
   stemp[0] = RemoveFilename(&filename) + "/fithist/signal_background_data_logy.pdf";
   c1->SaveAs(stemp[0].c_str());

   delete c1;
   delete legend;

/*   TCanvas *c2 = new TCanvas("c2","",1200,900);
   TH1F *tempsig;
   TH1F *tempback;

   tempsig = (TH1F*)simsig->Clone("tempsig");
   tempback = (TH1F*)simback->Clone("tempback");

   tempsig->Scale((double)nrdata/(double)nrsig);
   cout << "tempsig = " << (double)nrdata/(double)nrsig << endl;
   tempback->Scale((double)nrdata/(double)nrback);
   cout << "tempback = " << (double)nrdata/(double)nrback << endl;

   mystyle->SetHistColor((TH1*)tempsig, 1);
   mystyle->SetHistColor((TH1*)tempback, 0);
   mystyle->SetHistColor((TH1*)data, 2);
   mystyle->SetAxisTitles((TH1*)tempsig, "MVA variable", "Number of events");
   mystyle->SetAxisTitles((TH1*)tempback, "MVA variable", "Number of events");
   mystyle->SetAxisTitles((TH1*)data, "MVA variable", "Number of events");

   tempsig->Draw();
   tempback->Draw("SAME");
   if(nrsig > nrback)
   {
      cout << "Histogram maximum value = " << TMath::MaxElement(nrkeys, norm)*((double)nrdata/(double)nrsig) << endl;
      tempsig->SetMaximum(1.2*TMath::MaxElement(nrkeys, norm)*((double)nrdata/(double)nrsig));
   }
   else
   {
      cout << "Histogram maximum value = " << TMath::MaxElement(nrkeys, norm)*((double)nrdata/(double)nrback) << endl;
      tempsig->SetMaximum(1.2*TMath::MaxElement(nrkeys, norm)*((double)nrdata/(double)nrback));
   }
   data->Draw("SAME");

   legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(3)), gPad->GetLeftMargin()+0.32, 1-gPad->GetTopMargin());
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_LegendFill);
   legend->AddEntry(tempsig, "Signal events", "f");
   legend->AddEntry(tempback, "Background events", "f");
   legend->AddEntry(data, "Data events", "f");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("SAME");

   stemp[0] = RemoveFilename(&filename) + "/fithist/signal_background_data_normalized.pdf";
   c2->SaveAs(stemp[0].c_str());

   delete tempsig;
   delete tempback;

   delete c2;
   delete legend;*/
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
         PlotSumResiduals(dtemp, dtemp, 0, 0);
      }
   }
   // Find the best fit with TMinuit
   else if(fitproc == 1)
   {
      ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","MIGRAD");
//      ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
      minim->SetMaxFunctionCalls(1000000);
      minim->SetMaxIterations(5000);
      minim->SetPrintLevel(0);

      if(constraint == 0)
         itemp[2] = nrparam;
      else
         itemp[2] = nrparam-1;
      ROOT::Math::Functor fmin(this,&MvaFitHist::fcn,itemp[2]);
      srand48(time(NULL));

      itemp[0] = 0;
      itemp[1] = -1;
      while(itemp[0] < 6)
      {
         for(int i = 0; i < itemp[2]; i++)
         {
            // As a rule of thumb, it should be between 5-10% of the range (0.05-0.1)
            if(itemp[0] == 0) minstep[i] = 0.05;
            if(itemp[0] == 1) minstep[i] = 0.03;
            if(itemp[0] == 2) minstep[i] = 0.01;
            if(itemp[0] == 3) minstep[i] = 0.005;
            if(itemp[0] == 4) minstep[i] = 0.003;
            if(itemp[0] == 5) minstep[i] = 0.001;

            if((appvar != -1) && (i == 0))
            {
               minvar[i] = appvar;
               cout << "Approximate value of r" << ToString(i+1) << " = " << appvar << endl;
            }
            else
	    {
               minvar[i] = drand48();
	       cout << "Random value of r" << ToString(i+1) << " = " << minvar[i] << endl;
	    }
         }
//         double step = 0.01, var = drand48();
         minim->SetFunction(fmin);
         for(int i = 0; i < itemp[2]; i++)
         {
            stemp[0] = "frac" + ToString(i);
/*	    if(constraint == 0)
               minim->SetVariable(i, stemp[0].c_str(), minvar[i], minstep[i]);
	    else*/
               minim->SetLimitedVariable(i, stemp[0].c_str(), minvar[i], minstep[i], 0.0, 1.0);
         }
         minim->Minimize();
	 itemp[0]++;
	 itemp[1] = minim->Status();

         if(itemp[1] == 0)
	 {
            cout << "Minimizer succesfully found a minimum (step = " << minstep[0] << ")." << endl;
	    break;
	 }
         else
	 {
            cout << "Minimizer did not converge and might have not found a correct minimum (step = " << minstep[0] << ")." << endl;
            cout << endl << "Minimization status = " << minim->Status() << endl;
            cout << "MinValue = " << minim->MinValue() << endl;
            cout << "Edm = " << minim->Edm() << endl;
            cout << "Number of function calls = " << minim->NCalls() << endl;
            cout << "Number of iterations = " << minim->NIterations() << endl;
	 }
         cout << endl;
      }

      cout << endl << "Minimization status = " << minim->Status() << endl;
      cout << "MinValue = " << minim->MinValue() << endl;
      cout << "Edm = " << minim->Edm() << endl;
      cout << "Number of function calls = " << minim->NCalls() << endl;
      cout << "Number of iterations = " << minim->NIterations() << endl;

      for(int i = 0; i < itemp[2]; i++)
      {
         cout << "Minos Error " << i << ": " << minim->GetMinosError(i, dtemp[0], dtemp[1]) << "\t";
	 cout << dtemp[0] << ", " << dtemp[1] << endl;
      }

      double *bestfrac = (double*)minim->X();
/*      cout << "Running Hesse = " << (int)minim->Hesse() << endl;
      cout << "Valid error: " << (int)minim->IsValidError() << endl;*/
      double *bestfracerr = (double*)minim->Errors();

      dtemp[0] = 0;
      for(int i = 0; i < itemp[2]; i++)
      {
	 dtemp[0] += bestfrac[i];
         cout << "frac" << i << "\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;
      }

      // Is this okay???
      if(constraint == 0)
      {
         cout << "Sum of all fractions: " << dtemp[0] << endl;
         for(int i = 0; i < itemp[2]; i++)
	 {
            bestfrac[i] = bestfrac[i]/dtemp[0];
            bestfracerr[i] = bestfracerr[i]/dtemp[0];
            cout << "frac" << i << "\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;
	 }
      }

      // Get the compositions
      midComposition->clear();
      midCompositionErr->clear();
      dtemp[0] = 0;
      dtemp[1] = 0;
      for(int i = 0; i < itemp[2]; i++)
      {
         midComposition->push_back(bestfrac[i]);
         midCompositionErr->push_back(bestfracerr[i]);

	 dtemp[0] += bestfrac[i];
	 dtemp[1] += bestfracerr[i];
      }

      if(constraint == 1)
      {
         midComposition->push_back(1. - dtemp[0]);
         midCompositionErr->push_back(dtemp[1]);
      }

      // Get the lnA value
      midLna[0] = 0;
      midLna[1] = 0;
     
      for(int i = 0; i < midComposition->size(); i++)
      {
         if(i == 0)
            *stemp = "Proton";
         else if(i == 1)
            *stemp = "Helium";
         else if(i == 2)
            *stemp = "Oxygen";
         else if(i == 3)
            *stemp = "Iron";

         itemp[0] = ptype->GetZ(stemp);
         midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0]));
         midLna[1] += midCompositionErr->at(i)*TMath::Log(ptype->GetA(itemp[0]));
      }

      cout << "lnA value = " << midLna[0] << " +- " << midLna[1] << endl;
/*      cout << "Composition:" << endl;
      for(int i = 0; i < midComposition->size(); i++)
         cout << "  " << midComposition->at(i) << " +- " << midCompositionErr->at(i) << endl;*/

//      for(int i = 0; i < itemp[2]; i++)
//         cout << "frac" << i << "\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;

      PlotSumResiduals(bestfrac, bestfracerr, minim->Status(), minstep[0]);
   }
   else if(fitproc == 2)
   {
      double *tempcomposition = new double[4];

      cout << "Please enter compositions for the following elements:" << endl;
      cout << "- Proton = ";
      cin >> tempcomposition[0];
      cout << "- Helium = ";
      cin >> tempcomposition[1];
      cout << "- Oxygen = ";
      cin >> tempcomposition[2];
      if(constraint == 0)
      {
         cout << "- Iron   = ";
         cin >> tempcomposition[3];
      }
      else
         tempcomposition[3] = 1. - tempcomposition[0] - tempcomposition[1] - tempcomposition[2];

      cout << endl << "Selected compositions: " << tempcomposition[0] << " p, " << tempcomposition[1] << " He, " << tempcomposition[2] << " O, " << tempcomposition[3] << " Fe" << endl;

      PlotSumResiduals(tempcomposition, tempcomposition, 0, 0);

      delete[] tempcomposition;
   }
}

double MvaFitHist::GetFinalLna(int me)
{
   return midLna[me];
}

double MvaFitHist::GetFinalComposition(int type)
{
   return midComposition->at(type);
}

double MvaFitHist::GetFinalCompositionErr(int type)
{
   return midCompositionErr->at(type);
}

// Plot data, normalized signal+background and normalized residuals
void MvaFitHist::PlotSumResiduals(double *sigFrac, double *sigFracErr, int status, double selStep)
{
   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TCanvas *c2 = new TCanvas("c2","",1200,900);
   TLatex chiText;
   c1->cd();
   TPad *upperpad = new TPad("upperpad", "upperpad", 0.004, 0.490, 0.996, 0.996);
   TPad *lowerpad = new TPad("lowerpad", "lowerpad", 0.004, 0.000, 0.996, 0.490);
   upperpad->Draw();
   lowerpad->Draw();

   // Plot the data distribution
   upperpad->cd();
   mystyle->SetHistColor((TH1*)newdata, 2);
   mystyle->SetAxisTitles((TH1*)newdata, "MVA variable", "Probability density");
   newdata->GetXaxis()->SetTitleOffset(2.0);
   newdata->Draw();

   if(constraint == 0)
      fcnConstrainedFunc(sigFrac);
   else
      fcnFunc(sigFrac);

/*   if(constraint == 0)
   {
//      dtemp[0] = -1;
      dtemp[0] = 1;
      // Scale all simulation histograms with the nrdata/nrsim[i] value
      for(int i = 0; i < nrparam; i++)
      {
         simnorm[i] = (double)nrdata/(double)nrsim[i];
         stemp[0] = "sims" + ToString(i);
         sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

//         dtemp[0] += sigFrac[i];
         dtemp[0] -= sigFrac[i];
         sims[i]->Scale(sigFrac[i]*simnorm[i]);
      }

      // Sum all simulation histograms together
      sim1 = (TH1F*)sims[0]->Clone("sim1");
      for(int i = 1; i < nrparam; i++)
         sim1->Add(sims[i]);

      // Define constraint histogram: sum(f_i) = 1
      sim2 = (TH1F*)sim1->Clone("sim2");

//      for(int i = 0; i <= nrbins+1; i++)
//      {
//         if(i == 0)
//            sim2->SetBinContent(i, 0);
//         else if(i == nrbins+1)
//            sim2->SetBinContent(i, 0);
//         else
//         {
//            if(sim1->GetBinContent(i) == 0)
//               sim2->SetBinContent(i, 0);
//            else
//               sim2->SetBinContent(i, dtemp[0]);
//         }
//      }

      // Apply constraint
      sim = (TH1F*)sim1->Clone("sim");
//      sim->Add(sim2);

      cout << "1 - sum(f_i) = " << dtemp[0] << endl;
   }
   else
   {
      dtemp[0] = 1;
      // Scale all simulation histograms with the nrdata/nrsim[i] value
      for(int i = 0; i < nrparam; i++)
      {
         simnorm[i] = (double)nrdata/(double)nrsim[i];
         stemp[0] = "sims" + ToString(i);
         sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());


	 if(i == nrparam-1)
            sims[i]->Scale(dtemp[0]*simnorm[i]);
	 else
	 {
            dtemp[0] -= sigFrac[i];
            sims[i]->Scale(sigFrac[i]*simnorm[i]);
	 }
      }

      // Sum all simulation histograms together
      sim = (TH1F*)sims[0]->Clone("sim");
      for(int i = 1; i < nrparam; i++)
         sim->Add(sims[i]);

      cout << "1 - sum(f_i) = " << dtemp[0] << endl;
   }*/

   // Zero empty bin errors just for plotting purposes
   newsim = (TH1F*)sim->Clone("newsim");
   for(int j = 1; j < nrbins+1; j++)
   {
      if(sim->GetBinContent(j) == 0)
         newsim->SetBinError(j, 0);
   }

   mystyle->SetHistColor((TH1*)newsim, 0);
   mystyle->SetAxisTitles((TH1*)newsim, "MVA variable", "Probability density");
/*   yrange[1] = 0;
   for(int i = 0; i < trees.size()-1; i++)
   {
      if(norm[i] > yrange[1])
         yrange[1] = norm[i];
   }*/
   newdata->SetMaximum(1.1*TMath::MaxElement(2, yrange));
   newsim->Draw("SAME");

   chiVal = 0;
   chiNdf = 0;
   chiGood = 0;
   chiProb = sim->Chi2TestX(data, chiVal, chiNdf, chiGood, "WW", residVal);
   cout << "Chi2 = " << chiVal << ", prob = " << chiProb << ", NDF = " << chiNdf << ", igood = " << chiGood << endl;

   *chi2value = chiVal;
   *pvalue = chiProb;
   *ndfvalue = chiNdf;

   cout << "Getting residuals" << endl;
   TH1F *resid = new TH1F("resid", "", nrbins, xlim[0], xlim[1]);
   yresidrange[0] = 1.e+24;
   yresidrange[1] = -1.e+24;
   vector<double> *residAll, *residSig, *residBack;
   residAll = new vector<double>;
   residSig = new vector<double>;
   residBack = new vector<double>;
   for(int i = 1; i <= nrbins; i++)
   {
//      cout << sim->GetBinContent(i) << endl;
/*      if( ((double)sim->GetBinContent(i) + (double)data->GetBinContent(i) == 0) )
      {
         resid->SetBinContent(i, 0.);
	 cout << "One of the two bins is empty." << endl;
      }*/
/*      if( ((double)sim->GetBinContent(i) == 0) || ((double)data->GetBinContent(i) == 0) )
      {
	 cout << "One of the two bins is empty." << endl;
      }
      else*/
      if( ((double)sim->GetBinContent(i) != 0) && ((double)data->GetBinContent(i) != 0) )
      {
         if(resid->GetBinCenter(i) > mvacut)
            residSig->push_back(residVal[i-1]);
         else
            residBack->push_back(residVal[i-1]);

	 residAll->push_back(residVal[i-1]);

         resid->SetBinContent(i, residVal[i-1]);

         if(residVal[i-1] < yresidrange[0])
            yresidrange[0] = residVal[i-1];
         if(residVal[i-1] > yresidrange[1])
            yresidrange[1] = residVal[i-1];
      }
   }

   yresidrange[0] = TMath::Abs(yresidrange[0]);
   yresidrange[1] = TMath::Abs(yresidrange[1]);

   c2->cd();
   TH1F *residDist = new TH1F("residDist", "", nrbins, -1.1*TMath::MaxElement(2, yresidrange), 1.1*TMath::MaxElement(2, yresidrange));
   double *tempd = new double[residAll->size()];
   for(int i = 0; i < residAll->size(); i++)
   {
      residDist->Fill(residAll->at(i));
   }

   mystyle->SetHistColor((TH1*)residDist, 3);
   mystyle->SetAxisTitles((TH1*)residDist, "Normalized residuals", "Nr. of residuals");
   residDist->Draw();

   c2->Update();

   chiText.SetTextAlign(31);
   stemp[0] = "Residuals (" + ToString(residAll->size()) + "):";
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.03)*(residDist->GetMaximum()), stemp[0].c_str());
   for(int i = 0; i < residAll->size(); i++)
      tempd[i] = residAll->at(i);
   cout << "Residual all: nr = " << residAll->size() << ", mean = " << TMath::Mean(residAll->size(), tempd) << ", RMS = " << TMath::RMS(residAll->size(), tempd) << endl;
   stemp[0] = "mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4);
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.06)*(residDist->GetMaximum()), stemp[0].c_str());
   stemp[0] = "RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.09)*(residDist->GetMaximum()), stemp[0].c_str());

   stemp[0] = "Residuals sig/back (" + ToString(residSig->size()) + "/" + ToString(residBack->size()) + "):";
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.15)*(residDist->GetMaximum()), stemp[0].c_str());
   for(int i = 0; i < residSig->size(); i++)
      tempd[i] = residSig->at(i);
   stemp[0] = "mean sig = " + ToString(TMath::Mean(residSig->size(), tempd), 4);
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.18)*(residDist->GetMaximum()), stemp[0].c_str());
   cout << "Residual right (sig): nr = " << residSig->size() << ", mean = " << TMath::Mean(residSig->size(), tempd) << ", RMS = " << TMath::RMS(residSig->size(), tempd) << endl;
   dtemp[0] = TMath::Mean(residSig->size(), tempd);

   for(int i = 0; i < residBack->size(); i++)
      tempd[i] = residBack->at(i);
   stemp[0] = "mean bgd = " + ToString(TMath::Mean(residBack->size(), tempd), 4);
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.21)*(residDist->GetMaximum()), stemp[0].c_str());
   cout << "Residual left (back): nr = " << residBack->size() << ", mean = " << TMath::Mean(residBack->size(), tempd) << ", RMS = " << TMath::RMS(residBack->size(), tempd) << endl;
   dtemp[1] = TMath::Mean(residBack->size(), tempd);

   stemp[0] = "mean abs(bgd-sig) = " + ToString(TMath::Abs(dtemp[0] - dtemp[1]), 4);
   chiText.DrawLatex(TMath::MaxElement(2, yresidrange), (1.-0.24)*(residDist->GetMaximum()), stemp[0].c_str());
   cout << "Residual back-sig: mean = " << TMath::Abs(dtemp[0] - dtemp[1]) << endl;

   c1->cd();
   upperpad->cd();

   for(int i = 0; i < residAll->size(); i++)
      tempd[i] = residAll->at(i);

   chiText.SetTextAlign(21);
   if(fitproc == 0)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString(chiVal, 4) + "/" + ToString(chiNdf) + " = " + ToString(chiVal/(double)chiNdf, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      stemp[0] = "r_{p} = " + ToString(sigFrac[0], 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., 0.94*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString(chiVal, 4) + "/" + ToString(chiNdf) + " = " + ToString(chiVal/(double)chiNdf, 4) + ", status = " + ToString(status) + ", step = " + ToString(selStep, 3) + ", p-value = " + ToString(chiProb, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      dtemp[0] = 1.;
      dtemp[1] = 0.;
      for(int i = 0; i < nrparam; i++)
      {
         if(constraint == 0)
	 {
            stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " #pm " + ToString(sigFracErr[i], 4) + " (" + treeNames[i] + ")";
            chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
	 }
	 else
	 {
            if(i < nrparam-1)
            {
               dtemp[0] -= sigFrac[i];
               dtemp[1] += sigFracErr[i];
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " #pm " + ToString(sigFracErr[i], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
            else
            {
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(dtemp[0], 4) + " #pm " + ToString(dtemp[1], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
	 }
      }

      stemp[0] = "Residuals: mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4) + ", RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      stemp[0] = "MVA cut = " + ToString(mvacut, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else if(fitproc == 2)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString(chiVal, 4) + "/" + ToString(chiNdf) + " = " + ToString(chiVal/(double)chiNdf, 4) + ToString(status) + ", p-value = " + ToString(chiProb, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      dtemp[0] = 1.;
      dtemp[1] = 0.;
      for(int i = 0; i < nrparam; i++)
      {
         if(constraint == 0)
	 {
            stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " #pm " + ToString(sigFracErr[i], 4) + " (" + treeNames[i] + ")";
            chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
	 }
	 else
	 {
            if(i < nrparam-1)
            {
               dtemp[0] -= sigFrac[i];
               dtemp[1] += sigFracErr[i];
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " #pm " + ToString(sigFracErr[i], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
            else
            {
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(dtemp[0], 4) + " #pm " + ToString(dtemp[1], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
	 }
      }

      stemp[0] = "Residuals: mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4) + ", RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      stemp[0] = "MVA cut = " + ToString(mvacut, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }

   delete residAll;
   delete residSig;
   delete residBack;
   delete tempd;

   lowerpad->cd();
   mystyle->SetHistColor((TH1*)resid, 3);
   mystyle->SetAxisTitles((TH1*)resid, "", "Normalized residuals");
   resid->GetYaxis()->SetRangeUser(-1.1*TMath::MaxElement(2, yresidrange), 1.1*TMath::MaxElement(2, yresidrange));
   resid->Draw();

   TLine *line = new TLine(mvacut, -1.1*TMath::MaxElement(2, yresidrange), mvacut, 1.1*TMath::MaxElement(2, yresidrange));
   line->SetLineWidth(2);
   line->SetLineStyle(1);
   line->SetLineColor(kOrange+2);
   line->Draw("SAME");

   c1->Update();

   stemp[1] = "";
   itemp[1] = 0;
   itemp[0] = Find(treeNames, "Proton");
   if((itemp[0] != -1) && (itemp[0] <= nrparam))
   {
      if(itemp[1] == 0)
      {
         stemp[1] += "p";
	 itemp[1]++;
      }
      else
         stemp[1] += "-p";
   }

   itemp[0] = Find(treeNames, "Helium");
   if((itemp[0] != -1) && (itemp[0] <= nrparam))
   {
      if(itemp[1] == 0)
      {
         stemp[1] += "He";
	 itemp[1]++;
      }
      else
         stemp[1] += "-He";
   }

   itemp[0] = Find(treeNames, "Oxygen");
   if((itemp[0] != -1) && (itemp[0] <= nrparam))
   {
      if(itemp[1] == 0)
      {
         stemp[1] += "O";
	 itemp[1]++;
      }
      else
         stemp[1] += "-O";
   }

   itemp[0] = Find(treeNames, "Iron");
   if((itemp[0] != -1) && (itemp[0] <= nrparam))
   {
      if(itemp[1] == 0)
      {
         stemp[1] += "Fe";
	 itemp[1]++;
      }
      else
         stemp[1] += "-Fe";
   }

   if(itemp[1] > 0)
      stemp[1] += "_";

   if(fitproc == 0)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_" + ToString(100*sigFrac[0], 0) + ".pdf";
      c1->SaveAs(stemp[0].c_str());
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_minuit2.pdf";
      c1->SaveAs(stemp[0].c_str());
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());
   }
   else if(fitproc == 2)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_manual.pdf";
      c1->SaveAs(stemp[0].c_str());
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());
   }
   
   delete c1;
   delete c2;
}

int MvaFitHist::GetNrElem()
{
   return midComposition->size();
}

int MvaFitHist::GetNrBins()
{
   return nrbins;
}

double MvaFitHist::GetChi2()
{
   return *chi2value;
}

double MvaFitHist::GetPvalue()
{
   return *pvalue;
}

double MvaFitHist::GetNdf()
{
   return *ndfvalue;
}

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   if(argc > 1)
   {
      int *itemp = new int[2];
      double *dtemp = new double;
      string *stemp = new string;
      vector<double> *finalLna = new vector<double>;
      vector<double> *finalChi2 = new vector<double>;
      vector<double> *finalPvalue = new vector<double>;
      vector<double> *finalNdf = new vector<double>;
      vector<double> *finalComp = new vector<double>;

      cerr << "Perform histogram fit on a range of ratios (0), perform a minimization (1) or select a composition to plot (2)? ";
      cin >> itemp[0];
      
      if(itemp[0] == 0)
      {
         cerr << "What should be the step size for ratios (between 0 and 1)? ";
	 cin >> *dtemp;

	 while( (*dtemp > 1) || (*dtemp <= 0) )
	 {
            cerr << "What should be the step size for ratios (between 0 and 1)? ";
	    cin >> *dtemp;
	 }
      }
      else if( (itemp[0] != 0) && (itemp[0] != 1) && (itemp[0] != 2) )
      {
         cerr << "Error! Wrong fitting procedure." << endl;
	 delete finalLna;
         delete finalChi2;
         delete finalPvalue;
	 delete finalNdf;
	 delete finalComp;
	 return 1;
      }

      cerr << "Keep settings the same throughout all input files (0) or set settings for each file separately (1)? ";
      cin >> itemp[1];

      MvaFitHist *fithist = new MvaFitHist();
      for(int i = 0; i < argc-1; i++)
      {
         *stemp = string(argv[i+1]);

         fithist->PrepareHistograms(i, stemp, itemp, dtemp);
         fithist->PlotDistributions();
	 fithist->StartFitting();

	 for(int j = 0; j < fithist->GetNrElem(); j++)
	 {
            *dtemp = fithist->GetFinalComposition(j);
	    finalComp->push_back(*dtemp);
            *dtemp = fithist->GetFinalCompositionErr(j);
	    finalComp->push_back(*dtemp);
	 }

         *dtemp = fithist->GetFinalLna(0);
	 finalLna->push_back(*dtemp);
         *dtemp = fithist->GetFinalLna(1);
	 finalLna->push_back(*dtemp);

         *dtemp = fithist->GetChi2();
	 finalChi2->push_back(*dtemp);
         *dtemp = fithist->GetPvalue();
	 finalPvalue->push_back(*dtemp);
         *dtemp = fithist->GetNdf();
	 finalNdf->push_back(*dtemp);
      }
      
      cout << endl << "Final results:" << endl;
      cout << "# nr. elements = " << fithist->GetNrElem() << ", nr. bins = " << fithist->GetNrBins() << endl;
      cout << "# nr\tChi2/NDF \tp-val\tlnA \t+-  \t\tp   \t+-  \tHe  \t+-  \tO   \t+-  \tFe  \t+-" << endl;
      cout.precision(4);
      for(int i = 0; i < argc-1; i++)
      {
	 *itemp = 2*(fithist->GetNrElem());
	 cout << i << fixed << "\t" << finalChi2->at(i) << "/" << (int)finalNdf->at(i) << "\t" << finalPvalue->at(i) << "\t" << finalLna->at(2*i) << "\t" << finalLna->at(2*i+1) << "\t";
         for(int j = 0; j < *itemp; j++)
            cout << "\t" << finalComp->at((*itemp)*i+j);

	 cout << endl;
      }

      delete fithist;

      delete itemp;
      delete dtemp;
      delete stemp;

      delete finalLna;
      delete finalChi2;
      delete finalPvalue;
      delete finalNdf;
      delete finalComp;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (mvatree_file.root)." << endl;
      return 1;
   }

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
