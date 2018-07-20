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

MvaFitHist::MvaFitHist(vector<string> *primVals)
{
   stemp = new string[3];
   itemp = new int[3];
   dtemp = new double[3];

   nrbins = -1;
   cout << endl;
   while(nrbins < 4)
   {
      cerr << "Select number of bins (min of 4): ";
      cin >> nrbins;
   }

   cout << endl;
   PrintMethods();
   cerr << "Select the used MVA type: ";
   cin >> stemp[0];

   cerr << "Select constraint type (0 = function, 1 = one of the fractions): ";
   cin >> constraint;

   xlim[0] = GetMethodMin(stemp[0]);
   xlim[1] = GetMethodMax(stemp[0]);
   cout << "MVA limits (" << stemp[0] << "): " << xlim[0] << ", " << xlim[1] << endl;

   cerr << "Transform [0,1] range of parameters to [-inf,inf] (0 = don't transform, 1 = transform) or leave parameters unlimited (2)? ";
   cin >> rangeTransform;

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
   ndfvalue = new int;
   midEnergy = new double;
   midStep = new double;

   ptype = new PrimPart();

   includePart = new vector<int>;
   nrelem = primVals->size();
   for(int i = 0; i < nrelem; i++)
      IncludePrimaryType(primVals->at(i));

   if(constraint == 0)
      nrparam = nrelem;
   else
      nrparam = nrelem-1;

   emptyhisterr = 0.;
}

MvaFitHist::~MvaFitHist()
{
   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   delete[] norm;
   delete[] nentries;

   delete mystyle;

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
   delete midEnergy;
   delete midStep;

   delete ptype;

   delete includePart;
}

void MvaFitHist::IncludePrimaryType(string type)
{
   itemp[0] = ptype->GetZ(&type);
   includePart->push_back(itemp[0]);
   cout << "Adding " << ptype->GetName(itemp[0]) << " (" << itemp[0] << ") to the composition." << endl;
}

void MvaFitHist::ResetHistograms()
{
   delete ifile;
}

double MvaFitHist::fcnATan(const double *par)
{
   return TMath::ATan(par[0])/TMath::Pi() + 1./2.;
}

double MvaFitHist::fcnInvATan(const double *par)
{
   return TMath::Tan(TMath::Pi()*(par[0] - 1./2.));
}

// step#013
double MvaFitHist::fcnConstrainedFunc(const double *par)
{
   dtemp[0] = 1;
   // Scale all simulation histograms with elemental fractions
   for(int i = 0; i < nrelem; i++)
   {
      stemp[0] = "sims" + ToString(i);
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      if(rangeTransform == 1)
      {
         dtemp[0] -= fcnATan(&par[i]);
         sims[i]->Scale(fcnATan(&par[i]));
      }
      else
      {
         dtemp[0] -= par[i];
         sims[i]->Scale(par[i]);
      }
   }

   // Sum all simulation histograms together
   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 1; i < nrelem; i++)
      sim->Add(sims[i]);

   // Go over the simulation distribution and check all the bins -> keep only bins from sim and data that are not 0
//   cout << "j   \tsim \t+-  \tdata\t+-  \tsim \t+-  \tdata\t+-" << endl;
//   cout.precision(4);
   dtemp[1] = 0;
   for(int j = 1; j < nrbins+1; j++)
   {
      dtemp[1] += sim->GetBinContent(j);
//      cout << j << fixed << "\t" << sim->GetBinContent(j) << "\t" << sim->GetBinError(j) << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j);

      if(sim->GetBinContent(j) == 0)
      {
         data->SetBinContent(j, 0);
         data->SetBinError(j, sim->GetBinError(j));
      }

      if(data->GetBinContent(j) == 0)
      {
         sim->SetBinContent(j, 0);
         sim->SetBinError(j, data->GetBinError(j));
      }

//      cout << j << fixed << "\t" << sim->GetBinContent(j) << "\t" << sim->GetBinError(j) << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j) << endl;
   }
//   cout << "Histogram integral = " << dtemp[1] << endl;

   // Return constraint: [1 - sum(f_i)]
   return dtemp[0];
}

// step#013
double MvaFitHist::fcnFunc(const double *par)
{
   dtemp[0] = 1;
   // Scale all simulation histograms with elemental fractions
   for(int i = 0; i < nrelem; i++)
   {
      stemp[0] = "sims" + ToString(i);
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      if(i == nrelem-1)
      {
         if(rangeTransform == 1)
            sims[i]->Scale(fcnATan(&dtemp[0]));
	 else
            sims[i]->Scale(dtemp[0]);
      }
      else
      {
         if(rangeTransform == 1)
         {
            dtemp[0] -= fcnATan(&par[i]);
            sims[i]->Scale(fcnATan(&par[i]));
         }
         else
         {
            dtemp[0] -= par[i];
            sims[i]->Scale(par[i]);
         }
      }
   }

   // Sum all simulation histograms together
   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 1; i < nrelem; i++)
      sim->Add(sims[i]);
   
   // Go over the simulation distribution and check all the bins
//   cout << "j   \tsim \t+-  \tdata\t+-  \tsim \t+-  \tdata\t+-" << endl;
//   cout.precision(4);
   dtemp[1] = 0;
   for(int j = 1; j < nrbins+1; j++)
   {
      dtemp[1] += sim->GetBinContent(j);
//      cout << j << fixed << "\t" << sim->GetBinContent(j) << "\t" << sim->GetBinError(j) << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j);

      if(sim->GetBinContent(j) == 0)
      {
         data->SetBinContent(j, 0);
         data->SetBinError(j, sim->GetBinError(j));
      }

      if(data->GetBinContent(j) == 0)
      {
         sim->SetBinContent(j, 0);
         sim->SetBinError(j, data->GetBinError(j));
      }

//      cout << j << fixed << "\t" << sim->GetBinContent(j) << "\t" << sim->GetBinError(j) << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j) << endl;
   }
//   cout << "Histogram integral = " << dtemp[1] << endl;

   // No returned constraint, because it is already included with the last parameter
   return 0;
}

// step#013
double MvaFitHist::fcn(const double *par)
{
//   cout << "New iteration -----------------------" << endl;
   if(constraint == 0)
      dtemp[0] = fcnConstrainedFunc(par);
   else
      dtemp[0] = fcnFunc(par);

//   return sim->Chi2Test(data, "WW CHI2") + dtemp[0];
   return sim->Chi2Test(data, "WW CHI2") + TMath::Power(dtemp[0],2);
}

void MvaFitHist::PrepareHistograms(int run, string *fname, int *proc, double *step)
{
   fitproc = proc[0];
//   setopt = proc[1];
   ratioStep = *step;
   filename = *fname;
   mystyle->SetBaseStyle();

// step#002
   cout << endl << "Opening file: " << filename << endl;
   ifile = new TFile(filename.c_str(), "READ");

   float mva;

   // Check how many TTrees we have in the file
   nrkeys = (int)ifile->GetNkeys()/2;
   cout << "Number of trees in file = " << nrkeys << endl;

   // Prepare variables (histograms to hold values from each tree)
// step#003
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
      cout << "Getting tree = " << stemp[0] << ", ";
      tempTree = (TTree*)ifile->Get(stemp[0].c_str());
      cout << "with title = " << tempTree->GetTitle() << endl;

      // Extract the MVA variable
      tempTree->SetBranchAddress("MVA", &mva);
      nentries[i] = tempTree->GetEntries();
      cout << "  Number of entries = " << nentries[i] << endl;

      vector<double> *binnentries = new vector<double>;

// step#004
      // Go through all events in the root file for the specific TTree and save MVA variable values to the histogram
      for(int j = 0; j < tempTree->GetEntries(); j++)
      {
         tempTree->GetEntry(j);
         result[i]->Fill(mva);
      }

// step#005
      // Calculate errors
      cout << "result" << i << ":" << endl;
      for(int k = 1; k <= nrbins; k++)
      {
	 dtemp[0] = result[i]->GetBinContent(k);
	 cout << k << "\t" << dtemp[0] << "\t" << result[i]->GetBinError(k) << endl;
	 if(dtemp[0] == 0)
            binnentries->push_back(emptyhisterr);
	 else
            binnentries->push_back(TMath::Sqrt(dtemp[0])/dtemp[0]);
      }

// step#006
      // Scale the histograms by number of all events (to get probability density functions) and get the maximal value
      result[i]->Scale(1./nentries[i]);
      norm[i] = result[i]->GetMaximum();
      cout << "  Highest value = " << norm[i] << endl;

      dtemp[1] = 0;
// step#007
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
      if( /*((setopt == 0) &&*/ (run == 0)/*) || (setopt == 1)*/ )
      {
         if(!trees.empty())
            trees.erase(trees.begin(), trees.end());
         if(!treeNames.empty())
            treeNames.erase(treeNames.begin(), treeNames.end());

	 itemp[1] = 0;
         for(int i = 0; i < nrkeys; i++)
         {
            // Saving simulation trees
//            cout << "Analysis = " << analRes->GetTreeName(i) << endl;
            for(int j = 0; j < nrelem; j++)
            {
//               cout << "Element = " << ptype->GetName(includePart->at(j)) << endl;
               stemp[1] = analRes->GetTreeName(i);
               stemp[2] = ptype->GetName(includePart->at(j));
               if( stemp[1].find(stemp[2]) != string::npos )
               {
		  // Make sure to only select pure samples for simulations, not mixed mock datasets
                  if( (stemp[1].find("[") == string::npos) && (stemp[1].find("]") == string::npos) )
		  {
                     trees.push_back(i);
                     treeNames.push_back(stemp[2]);
//                     cout << "They are equal (" << i << "," << j << ")" << endl;
		     itemp[1]++;
		  }
               }
            }

	    if(itemp[1] >= nrelem)
               break;
         }

         // Saving the data tree
         cerr << endl << "Available trees:" << endl;
         analRes->PrintVectors(0, 0);
         cerr << "Enter tree number (shown in square brackets) to select data tree: ";
      	 cin >> itemp[1];
	 cout << "Selected data tree: [" << itemp[1] << "], " << analRes->GetTreeName(itemp[1]) << endl;
         trees.push_back(itemp[1]);
         treeNames.push_back(analRes->GetTreeName(itemp[1]));
      }

      cout << "The selected simulation trees are: ";
      for(int i = 0; i < nrelem; i++)
      {
         if(i > 0)
            cout << ", ";
         cout << trees[i] << " (" << treeNames[i] << ")";
      }
      cout << endl << "The selected data tree is: " << trees[nrelem] << " (" << treeNames[nrelem] << ")" << endl << endl;

      // Getting the signal fraction approximation for the data
      mvacut = analRes->GetMvaCut(0);
      *midEnergy = analRes->GetEnergy();

      cout << "Mva cut value: " << mvacut << endl;
      cout << "Energy value: " << *midEnergy << endl;
   }

// step#008
   // Save signal and background distributions (lightest and heaviest elements)
   for(int i = 0; i < nrelem; i++)
   {
      if(i == 0)
      {
         simsig = (TH1F*)result[trees[0]]->Clone("simsig");
         nrsig = nentries[trees[0]];
      }
      else if(i == nrelem-1)
      {
         simback = (TH1F*)result[trees[nrelem-1]]->Clone("simback");
         nrback = nentries[trees[nrelem-1]];
      }
   }

   // Save all simulation histograms
   for(int i = 0; i < nrelem; i++)
   {
      stemp[0] = "simhist" + ToString(i);
      simhist[i] = (TH1F*)result[trees[i]]->Clone(stemp[0].c_str());
      nrsim[i] = nentries[trees[i]];

      // Go over the simulation distribution and check all the bins
      cout << "Simhist " << i << " (" << treeNames[i] << "):" << endl;
      dtemp[0] = 0;
      for(int j = 1; j < nrbins+1; j++)
      {
	 dtemp[0] += simhist[i]->GetBinContent(j);
         cout << j << "\t" << simhist[i]->GetBinContent(j) << "\t" << simhist[i]->GetBinError(j) << endl;
      }
      cout << "Histogram integral = " << dtemp[0] << endl;
   }

// step#009
   // Check if all simulation histograms have zeroes at the same place (if not, remove the large error)
   iszero = new bool;
   *iszero = false;
   dtemp[0] = 0;
   for(int j = 1; j < nrbins+1; j++)
   {
      *iszero = false;
      dtemp[0] = 0;

      for(int i = 0; i < nrelem; i++)
      {
         dtemp[0] += simhist[i]->GetBinContent(j);
	 if(simhist[i]->GetBinContent(j) == 0)
            (*iszero) = true;
      }

      if( (dtemp[0] != 0) && (*iszero) )
      {
         for(int i = 0; i < nrelem; i++)
	 {
	    if(simhist[i]->GetBinContent(j) == 0)
               simhist[i]->SetBinError(j, 0);
	 }
      }
   }

   delete iszero;

   for(int i = 0; i < nrelem; i++)
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
//   nrparam = trees.size() - 1;
   cout << "Number of parameters for fitting: " << nrparam << endl;

   // Save data distributions
   data = (TH1F*)result[trees[nrelem]]->Clone("data");
   nrdata = nentries[trees[nrelem]];
   yrange[0] = norm[nrelem];

   // Go over the data distribution and check all the bins
   cout << "Data:" << endl;
   for(int j = 1; j < nrbins+1; j++)
      cout << j << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j) << endl;

   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(&filename) + "/fithist";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(&filename) + "/fithist/signal_background_data.pdf " + RemoveFilename(&filename) + "/fithist/composition_residuals*";
   system(stemp[0].c_str());

// step#010
   // Zero empty bin errors just for plotting purposes
   newdata = (TH1F*)data->Clone("newdata");
   newsig = (TH1F*)simsig->Clone("newsig");
   newback = (TH1F*)simback->Clone("newback");

   for(int j = 1; j < nrbins+1; j++)
   {
      if(data->GetBinContent(j) == 0)
         newdata->SetBinError(j, 0);

      if(simback->GetBinContent(j) == 0)
         newback->SetBinError(j, 0);

      if(simsig->GetBinContent(j) == 0)
         newsig->SetBinError(j, 0);
   }

   delete analRes;
}

// Plot the distributions of signal, background and data
void MvaFitHist::PlotDistributions()
{
   TLegend *legend;
   TCanvas *c1 = new TCanvas("c1","",1200,900);

   mystyle->SetHistColor((TH1*)newsig, 1);
   mystyle->SetHistColor((TH1*)newback, 0);
   mystyle->SetHistColor((TH1*)newdata, 2);
   mystyle->SetAxisTitles((TH1*)newsig, "MVA variable", "Probability density");
   mystyle->SetAxisTitles((TH1*)newback, "MVA variable", "Probability density");
   mystyle->SetAxisTitles((TH1*)newdata, "MVA variable", "Probability density");

   newsig->Draw();
   newback->Draw("SAME");
   cout << "Histogram maximum value = " << TMath::MaxElement(nrkeys, norm) << endl;
   newsig->SetMaximum(1.2*TMath::MaxElement(nrkeys, norm));
   newdata->Draw("SAME");

   legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(3)), gPad->GetLeftMargin()+0.32, 1-gPad->GetTopMargin());
   int c_LegendFill = TColor::GetColor("#ffff66");
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_LegendFill);
   for(int i = 0; i < nrelem; i++)
   {
      if(i == 0)
      {
         stemp[0] = "Lightest element (" + ptype->GetShortName(includePart->at(i)) + ")";
         legend->AddEntry(newsig, stemp[0].c_str(), "l");
      }
      else if(i == nrelem-1)
      {
         stemp[0] = "Heaviest element (" + ptype->GetShortName(includePart->at(i)) + ")";
         legend->AddEntry(newback, stemp[0].c_str(), "l");
      }
   }
/*   legend->AddEntry(newsig, "Signal events", "f");
   legend->AddEntry(newback, "Background events", "f");*/
   legend->AddEntry(newdata, "Data", "l");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("SAME");

   stemp[0] = RemoveFilename(&filename) + "/fithist/signal_background_data.pdf";
   c1->SaveAs(stemp[0].c_str());

   c1->SetLogy();
   newsig->SetMaximum((TMath::Power(10,1.1))*TMath::MaxElement(nrkeys, norm));
   stemp[0] = RemoveFilename(&filename) + "/fithist/signal_background_data_logy.pdf";
   c1->SaveAs(stemp[0].c_str());

   delete c1;
   delete legend;
}

// Plot data, normalized signal+background and normalized residuals
int MvaFitHist::StartFitting()
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
// step#012
      ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2","MIGRAD");
      minim->SetMaxFunctionCalls(1000000);
      minim->SetMaxIterations(5000);
      minim->SetPrintLevel(0);

// step#013
      ROOT::Math::Functor fmin(this,&MvaFitHist::fcn, nrparam);
      srand48(time(NULL));

      itemp[0] = 0;
      itemp[1] = -1;
      while(itemp[0] < 6)
      {
         for(int i = 0; i < nrparam; i++)
         {
// step#014
            // As a rule of thumb, it should be between 5-10% of the range (for [0,1] MVA -> 0.05-0.1, for [-1,1] MVA -> 0.1-0.2)
            if(itemp[0] == 0) minstep[i] = 0.05;
            if(itemp[0] == 1) minstep[i] = 0.03;
            if(itemp[0] == 2) minstep[i] = 0.01;
            if(itemp[0] == 3) minstep[i] = 0.005;
            if(itemp[0] == 4) minstep[i] = 0.003;
            if(itemp[0] == 5) minstep[i] = 0.001;
	 }

         cout << "New minimization run ------------------------" << endl;
         dtemp[0] = 1.;
         for(int i = 0; i < nrparam; i++)
         {
            minvar[i] = drand48()*dtemp[0];
	    if(rangeTransform == 1)
               minvar[i] = fcnInvATan(&minvar[i]);
//               minvar[i] = TMath::Tan(TMath::Pi()*(minvar[i] - 1./2.));
	    dtemp[0] -= minvar[i];
	    cout << "Random value of r" << ToString(i+1) << " = " << minvar[i] << endl;
         }

	 cout << "Remaining value = " << dtemp[0] << endl;

         minim->SetFunction(fmin);
         for(int i = 0; i < nrparam; i++)
         {
            stemp[0] = "frac" + ToString(i);
	    if(rangeTransform == 0)
               minim->SetLimitedVariable(i, stemp[0].c_str(), minvar[i], minstep[i], 0.0, 1.0);
	    else
               minim->SetVariable(i, stemp[0].c_str(), minvar[i], minstep[i]);
         }
         minim->Minimize();
	 minim->PrintResults();

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

      if(itemp[1] != 0)
      {
         cout << "Minimizer failed to converge! Dismissing this file from further analysis." << endl;
         cerr << "Minimizer failed to converge! Dismissing this file from further analysis." << endl;
	 return -1;
      }

      cout << endl << "Minimization status = " << minim->Status() << endl;
      cout << "MinValue = " << minim->MinValue() << endl;
      cout << "Edm = " << minim->Edm() << endl;
      cout << "Number of function calls = " << minim->NCalls() << endl;
      cout << "Number of iterations = " << minim->NIterations() << endl;

      for(int i = 0; i < nrparam; i++)
      {
         cout << "Minos Error " << i << ": " << minim->GetMinosError(i, dtemp[0], dtemp[1]) << "\t";
	 cout << dtemp[0] << ", " << dtemp[1] << endl;
      }

      double *bestfrac = (double*)minim->X();
      double *bestfracerr = (double*)minim->Errors();

      dtemp[0] = 0;
      for(int i = 0; i < nrparam; i++)
      {
	 if(rangeTransform == 1)
	 {
            bestfrac[i] = fcnATan(&bestfrac[i]);/*TMath::ATan(bestfrac[i])/TMath::Pi() + 1./2.;*/
            bestfracerr[i] = fcnATan(&bestfracerr[i]);/*TMath::ATan(bestfracerr[i])/TMath::Pi() + 1./2.;*/
         }

	 dtemp[0] += bestfrac[i];
         cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;
      }

      // Is this okay???
      if(constraint == 0)
      {
         cout << "Sum of all fractions: " << dtemp[0] << endl;
         for(int i = 0; i < nrparam; i++)
	 {
            bestfrac[i] = bestfrac[i]/dtemp[0];
            bestfracerr[i] = bestfracerr[i]/dtemp[0];
            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;
	 }
      }

      // Get the compositions
      midComposition->clear();
      midCompositionErr->clear();
      dtemp[0] = 0;
      dtemp[1] = 0;
      for(int i = 0; i < nrparam; i++)
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
     
      for(int i = 0; i < nrelem; i++)
      {
         itemp[0] = includePart->at(i);
         midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0]));
         midLna[1] += midCompositionErr->at(i)*TMath::Log(ptype->GetA(itemp[0]));
      }

      cout << "lnA value = " << midLna[0] << " +- " << midLna[1] << endl;

      // Get the parameter step value
      *midStep = minstep[0];

      PlotSumResiduals(bestfrac, bestfracerr, minim->Status(), minstep[0]);
   }
   else if(fitproc == 2)
   {
      double *tempcomposition = new double[nrelem];

      cerr << "Please enter compositions for the following elements:" << endl;
      dtemp[0] = 0.;
      for(int i = 0; i < nrelem; i++)
      {
         if(constraint == 0)
	 {
            cerr << "- " << ptype->GetName(includePart->at(i)) << " = ";
            cin >> tempcomposition[i];
	 }
	 else
	 {
            if(i < nrelem-1)
	    {
               cerr << "- " << ptype->GetName(includePart->at(i)) << " = ";
               cin >> tempcomposition[i];
	       dtemp[0] += tempcomposition[i];
	    }
	    else
	    {
               tempcomposition[i] = 1. - dtemp[0];
	    }
	 }
      }

      cout << endl << "Selected compositions:" << endl;
      for(int i = 0; i < nrelem; i++)
         cout << "- " << ptype->GetName(includePart->at(i)) << " = " << tempcomposition[i] << endl;

      PlotSumResiduals(tempcomposition, tempcomposition, 0, 0);

      delete[] tempcomposition;
   }
   // Find the best fit with TMinuit (new interface)
   else if(fitproc == 3)
   {
      fitproc = 1;

      double *params = new double[nrparam];

      ROOT::Fit::Fitter fitter;
      fitter.Config().SetParamsSettings(nrparam, params);
      srand48(time(NULL));

      ROOT::Math::MinimizerOptions opt;
      opt.SetMinimizerType("Minuit2");
      opt.SetMinimizerAlgorithm("MIGRAD");
      opt.SetTolerance(1);
      opt.SetMaxFunctionCalls(10000);
      opt.SetPrintLevel(0);
      opt.Print();

      fitter.Config().SetMinimizerOptions(opt);
      ROOT::Math::Functor fmin(this,&MvaFitHist::fcn, nrparam);
      fitter.SetFCN(fmin);

      ROOT::Fit::FitResult result;

      itemp[0] = 0;
      itemp[1] = -1;
      while(itemp[0] < 6)
      {
         for(int i = 0; i < nrparam; i++)
         {
            // As a rule of thumb, it should be between 5-10% of the range (for [0,1] MVA -> 0.05-0.1, for [-1,1] MVA -> 0.1-0.2)
            if(itemp[0] == 0) minstep[i] = 0.05;
            if(itemp[0] == 1) minstep[i] = 0.03;
            if(itemp[0] == 2) minstep[i] = 0.01;
            if(itemp[0] == 3) minstep[i] = 0.005;
            if(itemp[0] == 4) minstep[i] = 0.003;
            if(itemp[0] == 5) minstep[i] = 0.001;
	 }

         cout << "New minimization run ------------------------" << endl;
         dtemp[0] = 1.;
         for(int i = 0; i < nrparam; i++)
         {
            if( (constraint == 0) && (i == nrparam-1) )
               minvar[i] = dtemp[0];
	    else
               minvar[i] = drand48()*dtemp[0];

	    if(rangeTransform == 1)
               minvar[i] = fcnInvATan(&minvar[i]);
//               minvar[i] = TMath::Tan(TMath::Pi()*(minvar[i] - 1./2.));
            dtemp[0] -= minvar[i];
            cout << "Random value of r" << ToString(i+1) << " = " << minvar[i] << endl;

            fitter.Config().ParSettings(i).SetValue(minvar[i]);
	    if(rangeTransform == 0)
	    {
               fitter.Config().ParSettings(i).SetLimits(0.,1.);
               fitter.Config().ParSettings(i).SetStepSize(minstep[i]);
	    }
            else
               fitter.Config().ParSettings(i).SetStepSize(minstep[i]);
         }

/*	 cout << "# Fit 1 (no helium and oxygen): ---------------------------------" << endl;
	 fitter.Config().ParSettings(1).Fix();
	 fitter.Config().ParSettings(2).Fix();*/
         fitter.FitFCN();
         cout << "Minos errors: " << fitter.CalculateMinosErrors() << endl;
         result = fitter.Result();
         result.Print(cout);
/*	 cout << "# Fit 2 (no helium and iron): -----------------------------------" << endl;
	 fitter.Config().ParSettings(2).Release();
	 fitter.Config().ParSettings(3).Fix();
         fitter.FitFCN();
         cout << "Minos errors: " << fitter.CalculateMinosErrors() << endl;
         result = fitter.Result();
         result.Print(cout);
	 cout << "# Fit 3 (no proton and iron): -----------------------------------" << endl;
	 fitter.Config().ParSettings(1).Release();
	 fitter.Config().ParSettings(0).Fix();
         fitter.FitFCN();
         cout << "Minos errors: " << fitter.CalculateMinosErrors() << endl;
         result = fitter.Result();
         result.Print(cout);
	 cout << "# Fit 4 (all parameters): ---------------------------------------" << endl;
	 fitter.Config().ParSettings(0).Release();
	 fitter.Config().ParSettings(3).Release();*/
/*         fitter.FitFCN();
//         cout << "Minos errors: " << fitter.CalculateMinosErrors() << endl;
         result = fitter.Result();
         result.Print(cout);*/

	 itemp[0]++;
         itemp[1] = result.Status();

         if(itemp[1] == 0)
	 {
            cout << "Minimizer succesfully found a minimum (step = " << minstep[0] << ")." << endl;
	    break;
	 }
         else
            cout << "Minimizer did not converge and might have not found a correct minimum (step = " << minstep[0] << ")." << endl;
         cout << endl;

         cout << "---------------------------------------------" << endl;
      }

      if(itemp[1] != 0)
      {
         cout << "Minimizer failed to converge! Dismissing this file from further analysis." << endl;
         cerr << "Minimizer failed to converge! Dismissing this file from further analysis." << endl;
	 return -1;
      }

      cout << endl << "Minimization status = " << result.Status() << endl;
      cout << "MinValue = " << result.MinFcnValue() << endl;
      cout << "Edm = " << result.Edm() << endl;
      cout << "Number of function calls = " << result.NCalls() << endl;

      for(int i = 0; i < nrparam; i++)
      {
         if(result.HasMinosError(i))
            cout << "Parameter " << i << " has Minos errors." << endl;
	 else
            cout << "Parameter " << i << " does not have Minos errors." << endl;
      }

      double *bestfrac = (double*)result.GetParams();
      double *bestfracerr = (double*)result.GetErrors();

      dtemp[0] = 0;
      for(int i = 0; i < nrparam; i++)
      {
	 if(rangeTransform == 1)
	 {
            bestfrac[i] = fcnATan(&bestfrac[i]);/*TMath::ATan(bestfrac[i])/TMath::Pi() + 1./2.;*/
            bestfracerr[i] = fcnATan(&bestfracerr[i]);/*TMath::ATan(bestfracerr[i])/TMath::Pi() + 1./2.;*/
         }

	 dtemp[0] += bestfrac[i];
         cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;
      }

      // Is this okay???
      if(constraint == 0)
      {
         cout << "Sum of all fractions: " << dtemp[0] << endl;
         for(int i = 0; i < nrparam; i++)
	 {
            bestfrac[i] = bestfrac[i]/dtemp[0];
            bestfracerr[i] = bestfracerr[i]/dtemp[0];
            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[i], 5) << endl;
	 }
      }

      // Get the compositions
      midComposition->clear();
      midCompositionErr->clear();
      dtemp[0] = 0;
      dtemp[1] = 0;
      for(int i = 0; i < nrparam; i++)
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
         itemp[0] = includePart->at(i);
         midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0]));
         midLna[1] += midCompositionErr->at(i)*TMath::Log(ptype->GetA(itemp[0]));
      }

      cout << "lnA value = " << midLna[0] << " +- " << midLna[1] << endl;

      // Get the parameter step value
      *midStep = minstep[0];

      PlotSumResiduals(bestfrac, bestfracerr, result.Status(), 0.05);

      delete[] params;
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

double MvaFitHist::GetFinalEnergy()
{
   return *midEnergy;
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
      cout << "Constraint = " << fcnConstrainedFunc(sigFrac) << endl;
   else
      cout << "Constraint = " << fcnFunc(sigFrac) << endl;

   // Zero empty bin errors just for plotting purposes
   newsig = (TH1F*)sim->Clone("newsig");
   for(int j = 1; j < nrbins+1; j++)
   {
      if(sim->GetBinContent(j) == 0)
         newsig->SetBinError(j, 0);
   }

   mystyle->SetHistColor((TH1*)newsig, 0);
   mystyle->SetAxisTitles((TH1*)newsig, "MVA variable", "Probability density");
   newdata->SetMaximum(1.1*TMath::MaxElement(2, yrange));
   newsig->Draw("SAME");

   *pvalue = sim->Chi2TestX(data, *chi2value, *ndfvalue, chiGood, "WW", residVal);
   cout << "Chi2 = " << *chi2value << ", p-value = " << *pvalue << ", NDF = " << *ndfvalue << ", igood = " << chiGood << endl;

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
      // Only set residuals, if any of the two histogram bins are not zero
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
      stemp[0] = "#chi^{2}/NDF = " + ToString((*chi2value), 4) + "/" + ToString((*ndfvalue)) + " = " + ToString((*chi2value)/(double)(*ndfvalue), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      stemp[0] = "r_{p} = " + ToString(sigFrac[0], 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., 0.94*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString((*chi2value), 4) + "/" + ToString((*ndfvalue)) + " = " + ToString((*chi2value)/(double)(*ndfvalue), 4) + ", status = " + ToString(status) + ", step = " + ToString(selStep, 3) + ", p-value = " + ToString((*pvalue), 4);
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
      stemp[0] = "#chi^{2}/NDF = " + ToString((*chi2value), 4) + "/" + ToString((*ndfvalue)) + " = " + ToString((*chi2value)/(double)(*ndfvalue), 4) + ToString(status) + ", p-value = " + ToString((*pvalue), 4);
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
   for(int i = 0; i < nrelem; i++)
   {
      if(i == 0)
         stemp[1] += ptype->GetShortName(includePart->at(i));
      else
         stemp[1] += "-" + ptype->GetShortName(includePart->at(i));
   }

/*   stemp[1] = "";
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
      stemp[1] += "_";*/

   if(fitproc == 0)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_" + ToString(100*sigFrac[0], 0) + ".pdf";
      c1->SaveAs(stemp[0].c_str());
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());
   }
   else if(fitproc == 1)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_minuit2.pdf";
      c1->SaveAs(stemp[0].c_str());
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());
   }
   else if(fitproc == 2)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_manual.pdf";
      c1->SaveAs(stemp[0].c_str());
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());
   }
   
   delete c1;
   delete c2;
}

int MvaFitHist::GetNrElem()
{
   return nrelem;/*midComposition->size();*/
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

int MvaFitHist::GetNdf()
{
   return *ndfvalue;
}

double MvaFitHist::GetStep()
{
   return *midStep;
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

int MvaFitHist::GetElemType(int type)
{
   return includePart->at(type);
}

void ExtractElements(string instring, vector<string> *outstring)
{
   int smth = 1;
   size_t pos[2] = {0,0};
   string tempstring = instring;

   outstring->clear();

   while(smth == 1)
   {
      pos[1] = tempstring.find_first_of(',', pos[0]);
//      cerr << "Position of ',' is: " << (int)pos[1] << endl;
//      cerr << "Information is located between [" << pos[0] << "," << pos[1] << "]" << endl;

      if(pos[1] == string::npos)
      {
         outstring->push_back(tempstring.substr(pos[0]));
         cerr << "Found element: " << tempstring.substr(pos[0]) << endl;
         break;
      }
      else
      {
         outstring->push_back(tempstring.substr(pos[0], pos[1]-pos[0]));
         cerr << "Found element: " << tempstring.substr(pos[0], pos[1]-pos[0]) << endl;
      }

      pos[0] = pos[1]+1;
   }
}

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   if(argc > 1)
   {
      int *itemp = new int[3];
      double *dtemp = new double[5];
      string *stemp = new string[3];
      vector<double> *finalLna = new vector<double>;
      vector<double> *finalChi2 = new vector<double>;
      vector<double> *finalPvalue = new vector<double>;
      vector<double> *finalNdf = new vector<double>;
      vector<double> *finalComp = new vector<double>;
      vector<double> *finalEnergy = new vector<double>;
      vector<double> *finalStep = new vector<double>;

      string inname;

      cerr << "Perform histogram fit on a range of ratios (0), perform a minimization (1), select a composition to plot (2) or perform a minimization with a new approach (3)? ";
      cin >> itemp[0];
      
      if(itemp[0] == 0)
      {
         cerr << "What should be the step size for ratios (between 0 and 1)? ";
	 cin >> dtemp[0];

	 while( (dtemp[0] > 1) || (dtemp[0] <= 0) )
	 {
            cerr << "What should be the step size for ratios (between 0 and 1)? ";
	    cin >> dtemp[0];
	 }
      }
      else if( (itemp[0] != 0) && (itemp[0] != 1) && (itemp[0] != 2) && (itemp[0] != 3) )
      {
         cerr << "Error! Wrong fitting procedure." << endl;
	 delete finalLna;
         delete finalChi2;
         delete finalPvalue;
	 delete finalNdf;
	 delete finalComp;
	 delete finalEnergy;
	 delete finalStep;
	 delete[] itemp;
	 delete[] stemp;
	 delete[] dtemp;
	 return 1;
      }

/*      cerr << "Keep settings the same throughout all input files (0) or set settings for each file separately (1)? ";
      cin >> itemp[1];*/

      vector<string> *types = new vector<string>;
      cerr << "Select the elemental composition (enter element names separated by commas): ";
      cin >> stemp[0];
//      stemp[0] = "Proton,Helium,Oxygen,Iron";
      ExtractElements(stemp[0], types);
/*      types->push_back("Proton");
      types->push_back("Helium");
      types->push_back("Oxygen");
      types->push_back("Iron");*/

      MvaFitHist *fithist = new MvaFitHist(types);
      int nrp = 0;
      for(int i = 0; i < argc-1; i++)
      {
         inname = string(argv[i+1]);

// step#001
         fithist->PrepareHistograms(i, &inname, &itemp[0], &dtemp[0]);
// step#011
         fithist->PlotDistributions();
	 itemp[2] = fithist->StartFitting();
	 fithist->ResetHistograms();

	 if(itemp[2] != -1)
	 {
	    for(int j = 0; j < fithist->GetNrElem(); j++)
	    {
               dtemp[0] = fithist->GetFinalComposition(j);
	       finalComp->push_back(dtemp[0]);
               dtemp[0] = fithist->GetFinalCompositionErr(j);
	       finalComp->push_back(dtemp[0]);
	    }

            dtemp[0] = fithist->GetFinalLna(0);
	    finalLna->push_back(dtemp[0]);
            dtemp[0] = fithist->GetFinalLna(1);
	    finalLna->push_back(dtemp[0]);

	    dtemp[0] = fithist->GetFinalEnergy();
	    finalEnergy->push_back(dtemp[0]);

            dtemp[0] = fithist->GetChi2();
	    finalChi2->push_back(dtemp[0]);
            dtemp[0] = fithist->GetPvalue();
	    finalPvalue->push_back(dtemp[0]);
            dtemp[0] = (double)fithist->GetNdf();
	    finalNdf->push_back(dtemp[0]);
            dtemp[0] = fithist->GetStep();
	    finalStep->push_back(dtemp[0]);

	    nrp++;
	 }
      }

      // Printout finar results for each file
      double *xval = new double[nrp];
      double *xerr = new double[nrp];
      double *yval = new double[nrp];
      double *yerr = new double[nrp];

      PrimPart *ptype = new PrimPart();

      cout << endl << "Final results:" << endl;
      cout << "# nr. elements = " << fithist->GetNrElem() << ", nr. bins = " << fithist->GetNrBins() << endl;
      cout << "# E   \tChi2/NDF \tp-val\tstep\tlnA \t+-  \t";
      for(int i = 0; i < fithist->GetNrElem(); i++)
         cout << "\t" << ptype->GetShortName(fithist->GetElemType(i)) << "  \t+-  ";
      cout << endl;
      cout.precision(4);
      for(int i = 0; i < nrp; i++)
      {
	 itemp[0] = 2*(fithist->GetNrElem());
	 cout << fixed << finalEnergy->at(i) << "\t" << finalChi2->at(i) << "/" << (int)finalNdf->at(i) << "\t" << finalPvalue->at(i) << "\t" << finalStep->at(i) << "\t" << finalLna->at(2*i) << "\t" << finalLna->at(2*i+1) << "\t";

         xval[i] = finalEnergy->at(i);
	 xerr[i] = 0.;
	 yval[i] = finalLna->at(2*i);
	 yerr[i] = finalLna->at(2*i+1);

         for(int j = 0; j < itemp[0]; j++)
            cout << "\t" << finalComp->at(itemp[0]*i+j);

	 cout << endl;
      }

      cout << endl;

      RootStyle *mystyle = new RootStyle();
      mystyle->SetBaseStyle();
      TCanvas *c1 = new TCanvas("c1","",1200,900);

      // Prepare lnA plot, with appended published values
      TGraphErrors *grLna;
      TGraphAsymmErrors *grPubLna;
      TLegend *legend;
      int c_Legend = TColor::GetColor("#ffff66");

      grLna = new TGraphErrors(nrp, xval, yval, xerr, yerr);
      mystyle->SetGraphColor(grLna, 2);
      grLna->GetXaxis()->SetRangeUser(18.5, 20.0);
//      grLna->GetXaxis()->SetLimits(18.5, 20.0);
      grLna->GetYaxis()->SetRangeUser(-0.7, 5.);

      cerr << "Which dataset did you use (EPOS = 0, QGSJET = 1, SIBYLL = 2)? ";
      cin >> itemp[2];

      vector<float> returnVal;
      itemp[2] = ReadLnaResults(&returnVal, itemp[2]);
      cout << "Number of points = " << itemp[2] << endl;
      float *xbinPub[3];
      float *ybinPubLna[3];
     
      for(int i = 0; i < 3; i++)
      {
         xbinPub[i] = new float[itemp[2]];
         ybinPubLna[i] = new float[itemp[2]];
      }
     
      for(int i = 0; i < itemp[2]; i++)
      {
         xbinPub[0][i] = returnVal[3*i];
         xbinPub[1][i] = 0;
         xbinPub[2][i] = 0;
         ybinPubLna[0][i] = returnVal[3*i+1];
         ybinPubLna[1][i] = returnVal[3*i+2];
         ybinPubLna[2][i] = returnVal[3*i+2];
     
//         cout << i+1 << ", data: " << xbinPub[0][i] << "\t" << ybinPubLna[0][i] << "\t" << ybinPubLna[1][i] << endl;
      }
     
      grPubLna = new TGraphAsymmErrors(itemp[2], xbinPub[0], ybinPubLna[0], xbinPub[1], xbinPub[2], ybinPubLna[1], ybinPubLna[2]);
      mystyle->SetGraphColor(grPubLna, 0);
      grPubLna->GetYaxis()->SetRangeUser(-0.7, 5.);

      mystyle->SetAxisTitles(grLna, "FD energy [log(E/eV)]", "<lnA> of data events");
      grLna->Draw("APL");
      grPubLna->Draw("PL;SAME");

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(2)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
      legend->SetFillStyle(1001);
      legend->SetFillColor(c_Legend);
      legend->AddEntry(grLna, "Data (Histogram fit)", "lp");
      legend->AddEntry(grPubLna, "Data (Auger published)", "lp");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      TText *t = new TText();
      t->SetTextAlign(12);
      t->SetTextColor(28);
      t->SetTextSize(24);
      TLine *l = new TLine();
      l->SetLineWidth(2);
      l->SetLineStyle(7);
      l->SetLineColor(28);

      // Right side markings and particle type lines
      c1->Update();
      stemp[1] = "Proton";
      dtemp[0] = ptype->GetA(&stemp[1]);
      t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(dtemp[0]),  "p");
      l->DrawLine(gPad->GetUxmin(), TMath::Log(dtemp[0]), gPad->GetUxmax(), TMath::Log(dtemp[0]));
      stemp[1] = "Helium";
      dtemp[0] = ptype->GetA(&stemp[1]);
      t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(dtemp[0]),  "He");
      l->DrawLine(gPad->GetUxmin(), TMath::Log(dtemp[0]), gPad->GetUxmax(), TMath::Log(dtemp[0]));
      stemp[1] = "Nitrogen";
      dtemp[0] = ptype->GetA(&stemp[1]);
      t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(dtemp[0]), "N");
      l->DrawLine(gPad->GetUxmin(), TMath::Log(dtemp[0]), gPad->GetUxmax(), TMath::Log(dtemp[0]));
      stemp[1] = "Oxygen";
      dtemp[0] = ptype->GetA(&stemp[1]);
      t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(dtemp[0]), "O");
      l->DrawLine(gPad->GetUxmin(), TMath::Log(dtemp[0]), gPad->GetUxmax(), TMath::Log(dtemp[0]));
      stemp[1] = "Iron";
      dtemp[0] = ptype->GetA(&stemp[1]);
      t->DrawText(gPad->GetUxmax() + (gPad->GetUxmax()-gPad->GetUxmin())/100., TMath::Log(dtemp[0]), "Fe");
      l->DrawLine(gPad->GetUxmin(), TMath::Log(dtemp[0]), gPad->GetUxmax(), TMath::Log(dtemp[0]));

      stemp[1] = RemoveFilename(&inname);
      stemp[2] = "mkdir -p " + RemoveFilename(&stemp[1]) + "/plots";
      system(stemp[2].c_str());

      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/histogram_fit_lnA.pdf";
      c1->SaveAs(stemp[2].c_str());
     
      for(int i = 0; i < 3; i++)
      {
         delete[] xbinPub[i];
         delete[] ybinPubLna[i];
      }

      delete c1;

      cout << endl;

      double *ybinComp[40];
      double *ybinCompErr[40];

      // Plot composition plots (all together and separate)
      c1 = new TCanvas("c1","",1200,900);

      TGraphErrors *grComp[40];

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         ybinComp[j] = new double[nrp];
         ybinCompErr[j] = new double[nrp];

         for(int i = 0; i < nrp; i++)
         {
	    itemp[0] = 2*(fithist->GetNrElem());
	    ybinComp[j][i] = finalComp->at(2*j+itemp[0]*i);
	    ybinCompErr[j][i] = finalComp->at(2*j+1+itemp[0]*i);
//	    cout << "i = " << i << ", j = " << j << ": " << ybinComp[j][i] << "\t" << ybinCompErr[j][i] << endl;
         }
      }

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         grComp[j] = new TGraphErrors(nrp, xval, ybinComp[j], xerr, ybinCompErr[j]);
//         mystyle->SetGraphColor(grComp[j], 2);
         grComp[j]->SetMarkerStyle(20);
         grComp[j]->SetMarkerSize(1.4);
         grComp[j]->SetLineWidth(2);
         grComp[j]->GetYaxis()->SetRangeUser(0., 1.);
	 mystyle->SetColorScale(grComp[j], j, fithist->GetNrElem());

         if(j == 0)
	 {
            mystyle->SetAxisTitles(grComp[j], "FD energy [log(E/eV)]", "Elemental fractions");
	    grComp[j]->Draw("ALP");
	 }
	 else
            grComp[j]->Draw("LP;SAME");
      }

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(fithist->GetNrElem())), gPad->GetLeftMargin()+.20, 1-gPad->GetTopMargin());
      legend->SetFillStyle(1001);
      legend->SetFillColor(c_Legend);
      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
	 stemp[1] = ptype->GetName(fithist->GetElemType(j));
         legend->AddEntry(grComp[j], stemp[1].c_str(), "lp");
      }
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      stemp[1] = RemoveFilename(&inname);
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot.pdf";
      c1->SaveAs(stemp[2].c_str());

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         grComp[j]->Draw("ALP");
         stemp[1] = RemoveFilename(&inname);
         stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_" + ptype->GetName(fithist->GetElemType(j)) + ".pdf";
         c1->SaveAs(stemp[2].c_str());
      }

      // Plot composition plot with light and heavy compositions together
      double compMassBreak = 7.4;

      for(int i = 0; i < nrp; i++)
      {
	 dtemp[1] = 0.;
	 dtemp[2] = 0.;
	 dtemp[3] = 0.;
	 dtemp[4] = 0.;
	 if(i == 0)
	 {
	    stemp[1] = "";
	    stemp[2] = "";
	 }
	 itemp[0] = 0;
	 itemp[1] = 0;
         for(int j = 0; j < fithist->GetNrElem(); j++)
	 {
            dtemp[0] = ptype->GetA(fithist->GetElemType(j));

            if(dtemp[0] < compMassBreak)
	    {
               dtemp[1] += ybinComp[j][i];
               dtemp[2] += ybinCompErr[j][i];

	       if(i == 0)
	       {
		  if(itemp[0] == 0)
		  {
                     stemp[1] += ptype->GetShortName(fithist->GetElemType(j));
		     itemp[0]++;
		  }
		  else
		  {
                     stemp[1] += "+";
                     stemp[1] += ptype->GetShortName(fithist->GetElemType(j));
		  }
	       }
	    }
            else
	    {
               dtemp[3] += ybinComp[j][i];
               dtemp[4] += ybinCompErr[j][i];

	       if(i == 0)
	       {
		  if(itemp[1] == 0)
		  {
                     stemp[2] += ptype->GetShortName(fithist->GetElemType(j));
		     itemp[1]++;
		  }
		  else
		  {
                     stemp[2] += "+";
                     stemp[2] += ptype->GetShortName(fithist->GetElemType(j));
		  }
	       }
	    }
	 }

         ybinComp[0][i] = dtemp[1];
         ybinCompErr[0][i] = dtemp[2];
         ybinComp[1][i] = dtemp[3];
         ybinCompErr[1][i] = dtemp[4];
      }

      grComp[0] = new TGraphErrors(nrp, xval, ybinComp[0], xerr, ybinCompErr[0]);
      grComp[0]->SetMarkerStyle(20);
      grComp[0]->SetMarkerSize(1.4);
      grComp[0]->SetLineWidth(2);
      grComp[0]->GetYaxis()->SetRangeUser(0., 1.);
      mystyle->SetColorScale(grComp[0], 0, fithist->GetNrElem());
      mystyle->SetAxisTitles(grComp[0], "FD energy [log(E/eV)]", "Elemental fractions");
      grComp[0]->Draw("ALP");
      grComp[1] = new TGraphErrors(nrp, xval, ybinComp[1], xerr, ybinCompErr[1]);
      grComp[1]->SetMarkerStyle(20);
      grComp[1]->SetMarkerSize(1.4);
      grComp[1]->SetLineWidth(2);
      grComp[1]->GetYaxis()->SetRangeUser(0., 1.);
      mystyle->SetColorScale(grComp[1], fithist->GetNrElem()-1, fithist->GetNrElem());
      grComp[1]->Draw("LP;SAME");

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(2)), gPad->GetLeftMargin()+.26, 1-gPad->GetTopMargin());
      legend->SetFillStyle(1001);
      legend->SetFillColor(c_Legend);
      stemp[0] = "Light (" + stemp[1] + ")";
      legend->AddEntry(grComp[0], stemp[0].c_str(), "lp");
      stemp[0] = "Heavy (" + stemp[2] + ")";
      legend->AddEntry(grComp[1], stemp[0].c_str(), "lp");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      stemp[1] = RemoveFilename(&inname);
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_light-heavy.pdf";
      c1->SaveAs(stemp[2].c_str());


/*      for(int i = 0; i < nrp; i++)
      {
         ybinComp[0][i] += ybinComp[1][i];
         ybinCompErr[0][i] += ybinCompErr[1][i];
         ybinComp[2][i] += ybinComp[3][i];
         ybinCompErr[2][i] += ybinCompErr[3][i];
      }

      for(int j = 0; j < fithist->GetNrElem(); j++)
         delete grComp[j];

      grComp[0] = new TGraphErrors(nrp, xval, ybinComp[0], xerr, ybinCompErr[0]);
      grComp[0]->SetMarkerStyle(20);
      grComp[0]->SetMarkerSize(1.4);
      grComp[0]->SetLineWidth(2);
      grComp[0]->GetYaxis()->SetRangeUser(0., 1.);
      mystyle->SetColorScale(grComp[0], 0, fithist->GetNrElem());
      mystyle->SetAxisTitles(grComp[0], "FD energy [log(E/eV)]", "Elemental fractions");
      grComp[0]->Draw("ALP");

      grComp[1] = new TGraphErrors(nrp, xval, ybinComp[1], xerr, ybinCompErr[1]);
      grComp[1]->SetMarkerStyle(20);
      grComp[1]->SetMarkerSize(1.4);
      grComp[1]->SetLineWidth(2);
      grComp[1]->GetYaxis()->SetRangeUser(0., 1.);
      mystyle->SetColorScale(grComp[1], fithist->GetNrElem()-1, fithist->GetNrElem());
      grComp[1]->Draw("LP;SAME");

      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(2)), gPad->GetLeftMargin()+.26, 1-gPad->GetTopMargin());
      legend->SetFillStyle(1001);
      legend->SetFillColor(c_Legend);
      legend->AddEntry(grComp[0], "proton+helium", "lp");
      legend->AddEntry(grComp[1], "oxygen+iron", "lp");
      legend->SetBorderSize(1);
      legend->SetMargin(0.3);
      legend->Draw("same");

      stemp[1] = RemoveFilename(&stemp[0]);
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_light-heavy.pdf";
      c1->SaveAs(stemp[2].c_str());*/

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         delete[] ybinComp[j];
         delete[] ybinCompErr[j];
      }

      delete grComp[0];
      delete grComp[1];

      delete c1;

      delete mystyle;

      delete[] xval;
      delete[] xerr;
      delete[] yval;
      delete[] yerr;

      delete fithist;

      delete[] itemp;
      delete[] dtemp;
      delete[] stemp;

      delete finalLna;
      delete finalChi2;
      delete finalPvalue;
      delete finalNdf;
      delete finalComp;
      delete finalEnergy;
      delete finalStep;

      delete types;
      delete ptype;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (mvatree_file.root)." << endl;
      return 1;
   }

   cerr << "Plotting program finished correctly." << endl;
   return 0;
}
