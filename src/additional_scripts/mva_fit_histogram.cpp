#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
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

   appvar = -1;
}

MvaFitHist::MvaFitHist(int bincount, double xlow, double xhigh)
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
   delete[] nrsim;
   delete[] simnorm;

   delete[] minstep;
   delete[] minvar;
}

double MvaFitHist::fcn(const double *par)
{
//   cout << "Calculating Chi2 value:" << endl;
   dtemp[0] = 1.;
   for(int i = 0; i <= nrparam; i++)
   {
      simnorm[i] = (double)nrdata/(double)nrsim[i];
//      cout << i << " (" << trees[i] << "): " << simnorm[i] << "\t";
      stemp[0] = "sims" + ToString(i);
//      cout << stemp[0] << "\t";
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      if(i < nrparam)
      {
         sims[i]->Scale(par[i]*simnorm[i]);
         dtemp[0] -= par[i];
//         cout << par[i] << "\t" << dtemp[0] << endl;
      }
      else
      {
         sims[i]->Scale(dtemp[0]*simnorm[i]);
//         cout << dtemp[0] << endl;
      }
   }

   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 1; i <= nrparam; i++)
      sim->Add(sims[i]);
   
   return sim->Chi2Test(data, "WU CHI2");
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
      appvar = analRes->GetFraction(2, -1)/100.;
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
      nrsim[i] = nentries[trees[i]];
   }

   // Set the number of parameters that will be needed to fit (number of simulation trees - 1)
   nrparam = trees.size() - 2;
   cout << "Number of parameters for fitting: " << nrparam << endl;

   // Save data distributions
   data = (TH1F*)result[trees[trees.size()-1]]->Clone("data");
   nrdata = nentries[trees[trees.size()-1]];
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
         PlotSumResiduals(dtemp, dtemp, 0, 0);
      }
   }
   // Find the best fit with TMinuit
   else if(fitproc == 1)
   {
      ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","");
      minim->SetMaxFunctionCalls(1000000);
      minim->SetMaxIterations(5000);
//      minim->SetTolerance(0.01);
      minim->SetPrintLevel(0);

      ROOT::Math::Functor fmin(this,&MvaFitHist::fcn,nrparam);
      srand48(time(NULL));

      itemp[0] = 0;
      itemp[1] = -1;
      while(itemp[0] < 6)
      {
         for(int i = 0; i < nrparam; i++)
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
               minvar[i] = drand48();
         }
//         double step = 0.01, var = drand48();
         minim->SetFunction(fmin);
         for(int i = 0; i < nrparam; i++)
         {
            stemp[0] = "frac" + ToString(i);
//            minim->SetVariable(i, stemp[0].c_str(), minvar[i], minstep[i]);
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
            cout << "Minimizer did not converge and might have not found a correct minimum (step = " << minstep[0] << ")." << endl;
         cout << endl;
      }

      double *bestfrac = (double*)minim->X();
      double *bestfracerr = (double*)minim->Errors();

      cout << endl << "Minimization status = " << minim->Status() << endl;
      cout << "MinValue = " << minim->MinValue() << endl;
      cout << "Edm = " << minim->Edm() << endl;
      cout << "Number of function calls = " << minim->NCalls() << endl;
      cout << "Number of iterations = " << minim->NIterations() << endl;
      for(int i = 0; i < nrparam; i++)
         cout << "frac" << i << "\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;

      PlotSumResiduals(bestfrac, bestfracerr, minim->Status(), minstep[0]);
   }
}

// Plot data, normalized signal+background and normalized residuals
void MvaFitHist::PlotSumResiduals(double *sigFrac, double *sigFracErr, int status, double selStep)
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

//   cout << "Calculating Chi2 value:" << endl;
   dtemp[0] = 1.;
   for(int i = 0; i <= nrparam; i++)
   {
      simnorm[i] = (double)nrdata/(double)nrsim[i];
//      cout << i << " (" << trees[i] << "): " << simnorm[i] << "\t";
      stemp[0] = "sims" + ToString(i);
//      cout << stemp[0] << "\t";
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      if(i < nrparam)
      {
         sims[i]->Scale(sigFrac[i]*simnorm[i]);
         dtemp[0] -= sigFrac[i];
//         cout << sigFrac[i] << "\t" << dtemp[0] << endl;
      }
      else
      {
         sims[i]->Scale(dtemp[0]*simnorm[i]);
//         cout << dtemp[0] << endl;
      }
   }

   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 1; i <= nrparam; i++)
      sim->Add(sims[i]);

/*   signorm = (double)nrdata/(double)nrsig;
   bgdnorm = (double)nrdata/(double)nrback;
   
   sim1 = (TH1F*)simsig->Clone("sim1");
   sim1->Scale(sigFrac[0]*signorm);
   sim2 = (TH1F*)simback->Clone("sim2");
   sim2->Scale((1.-sigFrac[0])*bgdnorm);
   
   sim = (TH1F*)sim1->Clone("sim");
   sim->Add(sim2);*/

   mystyle->SetHistColor((TH1*)sim, 0);
   mystyle->SetAxisTitles((TH1*)sim, "MVA variable", "Number of events");
   yrange[1] = sim->GetMaximum();
   data->SetMaximum(1.1*TMath::MaxElement(2, yrange));
   sim->Draw("SAME");

   chiVal = 0;
   chiNdf = 0;
   chiGood = 0;
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
   if(fitproc == 0)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString(chiVal, 4) + "/" + ToString(chiNdf) + " = " + ToString(chiVal/(double)chiNdf, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      stemp[0] = "r_{p} = " + ToString(sigFrac[0], 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., 0.94*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString(chiVal, 4) + "/" + ToString(chiNdf) + " = " + ToString(chiVal/(double)chiNdf, 4) + ", status = " + ToString(status) + ", step = " + ToString(selStep, 3);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      dtemp[0] = 1.;
      dtemp[1] = 0.;
      for(int i = 0; i <= nrparam; i++)
      {
         if(i < nrparam)
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

   lowerpad->cd();
   mystyle->SetHistColor((TH1*)resid, 3);
   mystyle->SetAxisTitles((TH1*)resid, "", "Normalized residuals");
   yresidrange[0] = TMath::Abs(yresidrange[0]);
   yresidrange[1] = TMath::Abs(yresidrange[1]);
   resid->GetYaxis()->SetRangeUser(-1.1*TMath::MaxElement(2, yresidrange), 1.1*TMath::MaxElement(2, yresidrange));
   resid->Draw();

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
   }
   else if(fitproc == 1)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "composition_residuals_minuit2.pdf";
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
      int *itemp = new int[2];
      double *dtemp = new double;
      string *stemp = new string;

      cerr << "Perform histogram fit on a range of ratios (0) or perform a minimization (1)? ";
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
      else if( (itemp[0] != 0) && (itemp[0] != 1) )
      {
         cerr << "Error! Wrong fitting procedure." << endl;
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
      }
      delete fithist;

      delete itemp;
      delete dtemp;
      delete stemp;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (mvatree_file.root)." << endl;
      return 1;
   }

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
