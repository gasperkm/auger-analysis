#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "mva_fit_histogram.h"

using namespace std;

MvaFitHist::MvaFitHist(vector<string> *primVals)
{
   stemp = new string[3];
   itemp = new int[4];
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

   midLna = new double[3];
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

   sim = new TH1F("sim", "", nrbins, xlim[0], xlim[1]);
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

//   delete sim;
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

/*   // Scale all simulation histograms with elemental fractions
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
   }*/

   // Sum all simulation histograms together (and scale them depending on the fraction)
   for(int j = 1; j < nrbins+1; j++)
   {
      sim->SetBinContent(j, 0);
      sim->SetBinError(j, 0);
   }

//   sim = new TH1F("sim", "", nrbins, xlim[0], xlim[1]);
//   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 0; i < nrelem; i++)
   {
      if(rangeTransform == 1)
      {
         sim->Add(simhist[i], fcnATan(&par[i]));
         dtemp[0] -= fcnATan(&par[i]);
      }
      else
      {
         sim->Add(simhist[i], par[i]);
         dtemp[0] -= par[i];
      }
   }

   cout.precision(7);
/*   cout << "TEST |";
   for(int i = 0; i < nrelem; i++)
      cout << par[i] << "|";
   cout << endl;
   cout << "j   \tsim \t+-  \tdata\t+-" << endl;
   for(int j = 1; j < nrbins+1; j++)
   {
      cout << j << fixed << "\t" << sim->GetBinContent(j) << "\t" << sim->GetBinError(j) << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j) << endl;
   }*/

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
         newdata->SetBinContent(j, 0);
         newdata->SetBinError(j, sim->GetBinError(j));
         data->SetBinContent(j, 0);
         data->SetBinError(j, sim->GetBinError(j));
      }

      if( (newdata->GetBinContent(j) == 0) || (data->GetBinContent(j) == 0) )
      {
         sim->SetBinContent(j, 0);
         sim->SetBinError(j, newdata->GetBinError(j));
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

/*   cout << "Parameters: ";
   for(int i = 0; i < nrelem; i++)
      cout << par[i] << ", ";
   cout << "Chi2 = " << sim->Chi2Test(data, "WW CHI2") << ", Constraint = " << dtemp[0] << endl;*/

//   return sim->Chi2Test(data, "WW CHI2") + dtemp[0];
   return sim->Chi2Test(data, "WW CHI2") + TMath::Power(dtemp[0],2);
}

double MvaFitHist::PearsonChi2(TH1 *sim, TH1 *data, int chi2dist)
{
   double chi2ret = 0;

   for(int j = 1; j < nrbins+1; j++)
   {
      if( (sim->GetBinContent(j) != 0) && (data->GetBinContent(j) != 0) )
      {
         if(chi2dist == 0)
            chi2ret += (TMath::Power((sim->GetBinContent(j)) - (data->GetBinContent(j)), 2))/((sim->GetBinContent(j)));
	 else if(chi2dist == 1)
            chi2ret += (TMath::Power((sim->GetBinContent(j)) - (data->GetBinContent(j)), 2))/((sim->GetBinContent(j)) + (data->GetBinContent(j)));
	 else
            return -1;
      }
   }

   if(chi2dist == 1)
      chi2ret = chi2ret/2.;

   return chi2ret;
}

void MvaFitHist::PrepareHistograms(int run, string *fname, int *proc, double *step)
{
   fitproc = proc[0];
//   setopt = proc[1];
   ratioStep = *step;
   filename = *fname;
   mystyle->SetBaseStyle();

   if(fitproc == 5)
      scaleHist = false;
   else
      scaleHist = true;

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

   // In case if not scaling histograms, the simulation histograms need to have the same number of events
/*   equateVal = 1000000;	// EQUATING
   if(fitproc == 5)
   {
      for(int i = 0; i < nrkeys; i++)
      {
         stemp[0] = "TreeS" + ToString(i+1);
         tempTree = (TTree*)ifile->Get(stemp[0].c_str());
         stemp[1] = (string)tempTree->GetTitle();
         // Make sure to only select pure samples for simulations, not mixed mock datasets
         if( (stemp[1].find("[") == string::npos) && (stemp[1].find("]") == string::npos) && (stemp[1].find("Data") == string::npos) )
	 {
	    if(tempTree->GetEntries() < equateVal)
               equateVal = tempTree->GetEntries();
	 }
      }

      cout << "Minimum number of events for simulation trees = " << equateVal << endl;

      equateMcEvents = true;
   }
   else
      equateMcEvents = false;*/

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
/*      if(equateMcEvents && (equateVal < nentries[i]))		// EQUATING
         nentries[i] = equateVal;
      cout << "  Number of used entries = " << nentries[i] << endl;*/

// step#004
      // Go through all events in the root file for the specific TTree and save MVA variable values to the histogram
      for(int j = 0; j < nentries[i]; j++)
      {
         tempTree->GetEntry(j);
         result[i]->Fill(mva);
      }

// step#005
      // Calculate errors
      if(scaleHist)
      {
         vector<double> *binnentries = new vector<double>;

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
         // Reevaluate the bin errors on the scaled histogram
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
      else
      {
         norm[i] = result[i]->GetMaximum();
         cout << "  Highest value = " << norm[i] << endl;

         dtemp[1] = 0;
         for(int k = 1; k <= nrbins; k++)
	 {
            dtemp[0] = result[i]->GetBinContent(k);

            if(result[i]->GetBinError(k) > dtemp[1])
               dtemp[1] = result[i]->GetBinError(k);
	 }

         cout << "Maximum value " << norm[i] << " increased by " << dtemp[1] << endl;
         norm[i] += dtemp[1];
      }
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
		     simnorm[j] = norm[i];
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
	 simnorm[nrelem] = norm[itemp[1]];

/*	 cout << "simnorm:" << endl;
	 for(int j = 0; j < nrelem+1; j++)
	 {
            cout << j << ":\t" << simnorm[j] << endl;
	 }
	 cout << endl;*/
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
   if(scaleHist)
   {
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
   cout << "Histogram maximum value = " << TMath::MaxElement(nrelem+1, simnorm) << endl;

   newsig->SetMaximum(1.2*TMath::MaxElement(nrelem, simnorm));
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
   newsig->SetMaximum((TMath::Power(10,1.1))*TMath::MaxElement(nrelem, simnorm));
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

      return -1;
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

      stepsel = 0;
      itemp[1] = -1;
      while(stepsel < 6)
      {
// step#014
         // As a rule of thumb, it should be between 5-10% of the range (for [0,1] MVA -> 0.05-0.1, for [-1,1] MVA -> 0.1-0.2)
         if(stepsel == 0) minstep[stepsel] = 0.05;
         if(stepsel == 1) minstep[stepsel] = 0.025;
         if(stepsel == 2) minstep[stepsel] = 0.01;
         if(stepsel == 3) minstep[stepsel] = 0.005;
         if(stepsel == 4) minstep[stepsel] = 0.0025;
         if(stepsel == 5) minstep[stepsel] = 0.001;

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
               minim->SetLimitedVariable(i, stemp[0].c_str(), minvar[i], minstep[stepsel], 0.0, 1.0);
	    else
               minim->SetVariable(i, stemp[0].c_str(), minvar[i], minstep[stepsel]);
         }
         minim->Minimize();
	 minim->PrintResults();

	 itemp[1] = minim->Status();

         if(itemp[1] == 0)
	 {
            cout << "Minimizer succesfully found a minimum (step = " << minstep[stepsel] << ")." << endl;
	    break;
	 }
         else
            cout << "Minimizer did not converge and might have not found a correct minimum (step = " << minstep[stepsel] << ")." << endl;
         cout << endl;

	 stepsel++;
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
      double *dtempres = (double*)minim->Errors();
      double *bestfracerr = new double[2*nrparam];

      dtemp[0] = 0;
      for(int i = 0; i < nrparam; i++)
      {
         bestfracerr[2*i] = dtempres[i];	// for now, saving same errors for + and -
         bestfracerr[2*i+1] = dtempres[i];	// for now, saving same errors for + and -

	 if(rangeTransform == 1)
	 {
            bestfrac[i] = fcnATan(&bestfrac[i]);/*TMath::ATan(bestfrac[i])/TMath::Pi() + 1./2.;*/
            bestfracerr[2*i] = fcnATan(&bestfracerr[2*i]);/*TMath::ATan(bestfracerr[i])/TMath::Pi() + 1./2.;*/
            bestfracerr[2*i+1] = fcnATan(&bestfracerr[2*i+1]);/*TMath::ATan(bestfracerr[i])/TMath::Pi() + 1./2.;*/
         }

	 dtemp[0] += bestfrac[i];
         cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
      }

      // Is this okay???
      if(constraint == 0)
      {
         cout << "Sum of all fractions: " << dtemp[0] << endl;
         for(int i = 0; i < nrparam; i++)
	 {
            bestfrac[i] = bestfrac[i]/dtemp[0];
            bestfracerr[2*i] = bestfracerr[2*i]/dtemp[0];
            bestfracerr[2*i+1] = bestfracerr[2*i+1]/dtemp[0];
            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
	 }
      }

      // Get the compositions
      midComposition->clear();
      midCompositionErr->clear();
      dtemp[0] = 0;
      dtemp[1] = 0;
      dtemp[2] = 0;
      for(int i = 0; i < nrparam; i++)
      {
         midComposition->push_back(bestfrac[i]);
         midCompositionErr->push_back(bestfracerr[2*i]);
         midCompositionErr->push_back(bestfracerr[2*i+1]);

	 dtemp[0] += bestfrac[i];
	 dtemp[1] += bestfracerr[2*i]*bestfracerr[2*i];
	 dtemp[2] += bestfracerr[2*i+1]*bestfracerr[2*i+1];
      }

      if(constraint == 1)
      {
         midComposition->push_back(1. - dtemp[0]);
         midCompositionErr->push_back(TMath::Sqrt(dtemp[1]));
         midCompositionErr->push_back(TMath::Sqrt(dtemp[2]));
      }

      // Get the lnA value
      midLna[0] = 0;
      midLna[1] = 0;
      midLna[2] = 0;
     
      for(int i = 0; i < nrelem; i++)
      {
         itemp[0] = includePart->at(i);
	 cout << midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0])) << endl;
	 cout << midCompositionErr->at(i)*TMath::Log(ptype->GetA(itemp[0])) << endl;
         midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0]));
         midLna[1] += (midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[0])))*(midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[0])));
         midLna[2] += (midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[0])))*(midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[0])));
      }

      midLna[1] = TMath::Sqrt(midLna[1]);
      midLna[2] = TMath::Sqrt(midLna[2]);

      cout << "lnA value = " << midLna[0] << " +/- " << midLna[1] << "/" << midLna[2] << endl;

      // Get the parameter step value
      *midStep = minstep[stepsel];
      cout << "Parameter step = " << minstep[stepsel] << " (" << stepsel << ")" << endl;

      PlotSumResiduals(bestfrac, bestfracerr, minim->Status(), minstep[stepsel]);

      delete[] bestfracerr;
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

      return -1;
   }
   // Find the best fit with TMinuit (new interface)
   else if(fitproc == 3)
   {
      fitproc = 1;

      double *params = new double[nrparam];

      ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
      fitter->Config().SetParamsSettings(nrparam, params);
      srand48(time(NULL));

      ROOT::Math::MinimizerOptions opt;
//      opt.SetMinimizerType("Fumili");
      opt.SetMinimizerType("Minuit2");
      opt.SetMinimizerAlgorithm("MIGRAD");
      opt.SetTolerance(1);
      opt.SetMaxFunctionCalls(10000);
      opt.SetPrintLevel(0);
      opt.Print();

      fitter->Config().SetMinimizerOptions(opt);
      ROOT::Math::Functor fmin(this, &MvaFitHist::fcn, nrparam);
      fitter->SetFCN(fmin);

      ROOT::Fit::FitResult *result = new ROOT::Fit::FitResult();

      stepsel = 0;
      itemp[1] = -1;
      while(stepsel < 6)
      {
         // As a rule of thumb, it should be between 5-10% of the range (for [0,1] MVA -> 0.05-0.1, for [-1,1] MVA -> 0.1-0.2)
         if(stepsel == 0) minstep[stepsel] = 0.05;
         if(stepsel == 1) minstep[stepsel] = 0.025;
         if(stepsel == 2) minstep[stepsel] = 0.01;
         if(stepsel == 3) minstep[stepsel] = 0.005;
         if(stepsel == 4) minstep[stepsel] = 0.0025;
         if(stepsel == 5) minstep[stepsel] = 0.001;

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

            fitter->Config().ParSettings(i).SetValue(minvar[i]);
	    if(rangeTransform == 0)
	    {
               fitter->Config().ParSettings(i).SetLimits(0.,1.);
               fitter->Config().ParSettings(i).SetStepSize(minstep[stepsel]);
	    }
            else
               fitter->Config().ParSettings(i).SetStepSize(minstep[stepsel]);
         }

/*	 cout << "# Fit 1 (no helium and oxygen): ---------------------------------" << endl;
	 fitter->Config().ParSettings(1).Fix();
	 fitter->Config().ParSettings(2).Fix();*/
         fitter->FitFCN();
         cout << "Minos errors: " << fitter->CalculateMinosErrors() << endl;
         *result = fitter->Result();
         result->Print(cout);

         dFitError = new double[2];
         iFitError = new int[2];

	 // Check if sum is close to 1
	 dFitError[0] = 0.;
	 for(int i = 0; i < nrparam; i++)
            dFitError[0] += result->Parameter(i);

	 if(rangeTransform == 1)
            dFitError[0] = fcnInvATan(&dFitError[0]);

	 if( (dFitError[0] > 1.25) || (dFitError[0] < 0.75) )
            cout << "Redoing the minimization, because parameter sum is only " << dFitError[0] << " (above 1.25 or below 0.75, ideally should be 1)." << endl;
         else
	 {
/* New stuff -------------------------------------- */
	    // Set minimal limit close to the limits (0 or 1)
	    dFitError[1] = 0.001;

	    if(nrelem > 2)
	    {
	       // Set selection numbers
	       iFitError[0] = -1;
	       iFitError[1] = -1;
	       // Fixing first parameter with smallest error
	       dFitError[0] = 20.;
	       cout << endl << "Parameters:" << endl;
	       for(int i = 0; i < nrparam; i++)
	       {
                  cout << "- par" << i << ": " << result->Parameter(i) << "\t" << result->ParError(i) << endl;
	          if(result->ParError(i) != 0)
	          {
	             if(result->ParError(i) < dFitError[0])
                     {
                        dFitError[0] = result->ParError(i);
	                iFitError[0] = i;
	             }

	             if(result->Parameter(i) < dFitError[1])
	             {
                        iFitError[1] = i;
	             }
	          }
	       }
	       cout << endl;
/*	       if(iFitError[0] == iFitError[1])
	       {
	          cout << "Parameter " << iFitError[0] << " has smallest error and is close to limiting value, so fixing it." << endl;
	          fitter->Config().ParSettings(iFitError[0]).Fix();
	       }
	       else
	       {
                  if(iFitError[1] != -1)
	          {
	             cout << "Parameter " << iFitError[1] << " is close to limiting value, so fixing it." << endl;
	             fitter->Config().ParSettings(iFitError[1]).Fix();
	          }
	          else
	          {*/
	             cout << "Parameter " << iFitError[0] << " has smallest error, so fixing it." << endl;
                     fitter->Config().ParSettings(iFitError[0]).Fix();
/*	          }
	       }*/
               fitter->FitFCN();
               cout << "Minos errors: " << fitter->CalculateMinosErrors() << endl;
               *result = fitter->Result();
               result->Print(cout);
	    }

	    if(nrelem > 3)
	    {
	       // Set selection numbers
	       iFitError[0] = -1;
	       iFitError[1] = -1;
	       // Fixing second parameter with smallest error
	       dFitError[0] = 20.;
	       cout << endl << "Parameters:" << endl;
	       for(int i = 0; i < nrparam; i++)
	       {
                  cout << "- par" << i << ": " << result->Parameter(i) << "\t" << result->ParError(i) << endl;
	          if(result->ParError(i) != 0)
	          {
	             if(result->ParError(i) < dFitError[0])
                     {
                        dFitError[0] = result->ParError(i);
	                iFitError[0] = i;
	             }

	             if(result->Parameter(i) < dFitError[1])
	             {
                        iFitError[1] = i;
	             }
	          }
	       }
	       cout << endl;
/*	       if(iFitError[0] == iFitError[1])
	       {
	          cout << "Parameter " << iFitError[0] << " has smallest error and is close to limiting value, so fixing it." << endl;
	          fitter->Config().ParSettings(iFitError[0]).Fix();
	       }
	       else
	       {
                  if(iFitError[1] != -1)
	          {
	             cout << "Parameter " << iFitError[1] << " is close to limiting value, so fixing it." << endl;
	             fitter->Config().ParSettings(iFitError[1]).Fix();
	          }
	          else
	          {*/
	             cout << "Parameter " << iFitError[0] << " has smallest error, so fixing it." << endl;
                     fitter->Config().ParSettings(iFitError[0]).Fix();
/*	          }
	       }*/
               fitter->FitFCN();
               cout << "Minos errors: " << fitter->CalculateMinosErrors() << endl;
               *result = fitter->Result();
               result->Print(cout);

	       cout << endl << "Parameters:" << endl;
	       for(int i = 0; i < nrparam; i++)
                  cout << "- par" << i << ": " << result->Parameter(i) << "\t" << result->ParError(i) << endl;
	       cout << endl;
	    }

	    if(nrelem > 2)
	    {
	       // Unfix all parameters
	       cout << "Unfixing all parameters" << endl;
	       for(int i = 0; i < nrparam; i++)
	          fitter->Config().ParSettings(i).Release();
               fitter->FitFCN();
               cout << "Minos errors: " << fitter->CalculateMinosErrors() << endl;
               *result = fitter->Result();
               result->Print(cout);

	       cout << endl << "Parameters:" << endl;
	       for(int i = 0; i < nrparam; i++)
                  cout << "- par" << i << ": " << result->Parameter(i) << "\t" << result->ParError(i) << endl;
	       cout << endl;
	    }

/* New stuff -------------------------------------- */

            itemp[1] = result->Status();

	    iFitError[1] = 0;

            if(itemp[1] == 0)
	    {
               cout << "Minimizer succesfully found a minimum (step = " << minstep[stepsel] << ")." << endl;
	       cerr << "Redo the minimization (1 = yes, 0 = no)? ";
	       cin >> iFitError[1];
	       if(iFitError[1] == 0)
	       {
                  delete[] iFitError;
	          break;
	       }
	    }
            else
               cout << "Minimizer did not converge and might have not found a correct minimum (step = " << minstep[stepsel] << ")." << endl;
            cout << endl;

	    if(iFitError[1] == 0)
	       stepsel++;
	 }

	 delete[] dFitError;
	 delete[] iFitError;

         cout << "---------------------------------------------" << endl;
      }

      if(itemp[1] != 0)
      {
         cout << "Minimizer failed to converge! Dismissing this file from further analysis." << endl;
         cerr << "Minimizer failed to converge! Dismissing this file from further analysis." << endl;
	 return -1;
      }

      cout << endl << "Minimization status = " << result->Status() << endl;
      cout << "MinValue = " << result->MinFcnValue() << endl;
      cout << "Edm = " << result->Edm() << endl;
      cout << "Number of function calls = " << result->NCalls() << endl;

      for(int i = 0; i < nrparam; i++)
      {
         if(result->HasMinosError(i))
            cout << "Parameter " << i << " has Minos errors." << endl;
	 else
            cout << "Parameter " << i << " does not have Minos errors." << endl;
      }

      double *bestfrac = (double*)result->GetParams();
      double *dtempres = (double*)result->GetErrors();
      double *bestfracerr = new double[2*nrparam];

      dtemp[0] = 0;
      for(int i = 0; i < nrparam; i++)
      {
         bestfracerr[2*i] = dtempres[i];	// for now, saving same errors for + and -
         bestfracerr[2*i+1] = dtempres[i];	// for now, saving same errors for + and -

	 if(rangeTransform == 1)
	 {
            bestfrac[i] = fcnATan(&bestfrac[i]);/*TMath::ATan(bestfrac[i])/TMath::Pi() + 1./2.;*/
            bestfracerr[2*i] = fcnATan(&bestfracerr[2*i]);/*TMath::ATan(bestfracerr[i])/TMath::Pi() + 1./2.;*/
            bestfracerr[2*i+1] = fcnATan(&bestfracerr[2*i+1]);/*TMath::ATan(bestfracerr[i])/TMath::Pi() + 1./2.;*/
         }

	 dtemp[0] += bestfrac[i];
         cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
      }

      // Is this okay???
      if(constraint == 0)
      {
         cout << "Sum of all fractions: " << dtemp[0] << endl;
         for(int i = 0; i < nrparam; i++)
	 {
            bestfrac[i] = bestfrac[i]/dtemp[0];
            bestfracerr[2*i] = bestfracerr[2*i]/dtemp[0];
            bestfracerr[2*i+1] = bestfracerr[2*i+1]/dtemp[0];
            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
	 }
      }

      // Get the compositions
      midComposition->clear();
      midCompositionErr->clear();
      dtemp[0] = 0;
      dtemp[1] = 0;
      dtemp[2] = 0;
      for(int i = 0; i < nrparam; i++)
      {
         midComposition->push_back(bestfrac[i]);
         midCompositionErr->push_back(bestfracerr[2*i]);
         midCompositionErr->push_back(bestfracerr[2*i+1]);

	 dtemp[0] += bestfrac[i];
	 dtemp[1] += bestfracerr[2*i]*bestfracerr[2*i];
	 dtemp[2] += bestfracerr[2*i+1]*bestfracerr[2*i+1];
      }

      if(constraint == 1)
      {
         midComposition->push_back(1. - dtemp[0]);
         midCompositionErr->push_back(TMath::Sqrt(dtemp[1]));
         midCompositionErr->push_back(TMath::Sqrt(dtemp[2]));
      }

      // Get the lnA value
      midLna[0] = 0;
      midLna[1] = 0;
      midLna[2] = 0;
     
      for(int i = 0; i < midComposition->size(); i++)
      {
         itemp[0] = includePart->at(i);
	 cout << midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0])) << endl;
	 cout << midCompositionErr->at(i)*TMath::Log(ptype->GetA(itemp[0])) << endl;
         midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[0]));
         midLna[1] += (midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[0])))*(midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[0])));
         midLna[2] += (midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[0])))*(midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[0])));
      }

      midLna[1] = TMath::Sqrt(midLna[1]);
      midLna[2] = TMath::Sqrt(midLna[2]);

      cout << "lnA value = " << midLna[0] << " +/- " << midLna[1] << "/" << midLna[2] << endl;

      // Get the parameter step value
      *midStep = minstep[stepsel];
      cout << "Parameter step = " << minstep[stepsel] << " (" << stepsel << ")" << endl;

      PlotSumResiduals(bestfrac, bestfracerr, result->Status(), minstep[stepsel]);

      delete[] bestfracerr;
      delete[] params;
      delete fitter;
      delete result;
   }
   else if(fitproc == 4)
   {
      // Only reasonable, if we have 2 elements that are free and the others are fixed to some value
      if(nrelem > 1)
      {
         int *elemsteps = new int;
	 int *nrp = new int;
         vector<double> *elemval[2];
//	 double *midelemval = new double[2];
         double *elemSub = new double[5];
	 double *valFrac;
         vector<double> *chi2vect = new vector<double>;
         vector<double> *pvalvect = new vector<double>;
	 int *fixVal = new int[40];

	 vector<double> *gridval[2];

	 for(int i = 0; i < nrelem; i++)
	 {
  	    fixVal[i] = -1;
	    elemval[i] = new vector<double>;
	    gridval[i] = new vector<double>;
	 }
         elemSub[4] = 0;

	 if(nrelem > 2)
	 {
            // If we have 3 elements, fix one of them
            if(nrelem == 3)
            {
               cerr << endl << "Please select 1 element to keep fixed (use number in brackets):" << endl;
               for(int i = 0; i < nrelem; i++)
	       {
		  cerr << "[" << i << "]\t" << ptype->GetName(includePart->at(i)) << endl;
	       }
	       cerr << "Select element: ";
	       cin >> itemp[0];
	       fixVal[itemp[0]] = 1;

	       cerr << "Select its fixed value: ";
               cin >> dtemp[0];
	       elemSub[itemp[0]] = dtemp[0];
	       elemSub[4] += dtemp[0];

	       cout << endl << "Fixing element " << ptype->GetName(includePart->at(itemp[0])) << " to " << elemSub[itemp[0]] << endl;
	    }

	    // If we have 4 elements, fix two of them
            if(nrelem == 4)
            {
               cerr << endl << "Please select 2 elements to keep fixed (use number in brackets):" << endl;
               for(int i = 0; i < nrelem; i++)
	       {
		  cerr << "[" << i << "]\t" << ptype->GetName(includePart->at(i)) << endl;
	       }
	       cerr << "Select first element: ";
	       cin >> itemp[0];
	       fixVal[itemp[0]] = 1;

	       cerr << "Select its fixed value: ";
               cin >> dtemp[0];
	       elemSub[itemp[0]] = dtemp[0];
	       elemSub[4] += dtemp[0];

	       cout << endl << "Fixing element " << ptype->GetName(includePart->at(itemp[0])) << " to " << elemSub[itemp[0]] << endl;

	       cerr << "Select second element: ";
	       cin >> itemp[0];
	       fixVal[itemp[0]] = 1;

	       cerr << "Select its fixed value: ";
               cin >> dtemp[0];
	       elemSub[itemp[0]] = dtemp[0];
	       elemSub[4] += dtemp[0];

	       cout << endl << "Fixing element " << ptype->GetName(includePart->at(itemp[0])) << " to " << elemSub[itemp[0]] << endl;
	    }
	 }

	 // Printout fixed parameter values
	 cout << "Fixed values:" << endl;
	 for(int i = 0; i < nrelem; i++)
            cout << i << ": " << fixVal[i] << endl;
	 cout << endl;

	 // Limiting values, where chi2 is calculated
	 double limitchi[2] = {-0.25,1.25};

	 // Number of steps to be used when evaluating chi2
	 *elemsteps = ((limitchi[1] - limitchi[0])/ratioStep)+1;

         cout << "1.-elemSub[4] = " << 1.-elemSub[4] << endl;
         cout << "ratioStep = " << ratioStep << endl;
	 cout << "elemsteps = " << *elemsteps << endl;
	 cout << "limitchi = " << limitchi[0] << ", " << limitchi[1] << endl;

	 // Go through all steps for both free parameters
	 for(int i = 0; i < (*elemsteps); i++)
	 {
            for(int j = 0; j < (*elemsteps); j++)
	    {
               // If 2 elements, only select steps that give sum of parameters equal to 1
	       if(nrelem == 2)
	       {
                  if( ((double)i*ratioStep+limitchi[0]) + ((double)j*ratioStep+limitchi[0]) == 1. )
	          {
                     elemval[0]->push_back(((double)i*ratioStep+limitchi[0]));
                     elemval[1]->push_back(((double)j*ratioStep+limitchi[0]));
	          }
	       }
	       // If more than 2 elements, correctly select the parameters and subtract the fixed parameters from 1
               else
	       {
                  itemp[0] = 0;
		  dtemp[0] = (1. - elemSub[4] - ((double)i*ratioStep+limitchi[0]) - ((double)j*ratioStep+limitchi[0]));
		  // If remaining free parameter space is sufficiently small, treat it as 0
		  if( (dtemp[0] < (double)ratioStep/100.) && (dtemp[0] > (double)ratioStep/(-100.)) )
                     dtemp[0] = 0.;
                  if( dtemp[0] == 0. )
	          {
		     for(int k = 0; k < nrelem; k++)
		     {
	                if(fixVal[k] == -1)
		        {
			   if(itemp[0] == 0)
			   {
                              elemval[k]->push_back(((double)i*ratioStep+limitchi[0]));
			      itemp[0]++;
			   }
                           else if(itemp[0] == 1)
			   {
                              elemval[k]->push_back(((double)j*ratioStep+limitchi[0]));
			      itemp[0]++;
			   }
	                }
		        else
		        {
		           elemval[k]->push_back(elemSub[k]);
		        }
		     }
		  }
	       }
	    }
	 }

	 // Get the number of points that were saved (1D)
	 *nrp = elemval[0]->size();

	 // Prepare free parameters and calculate chi2 for them (1D)
	 for(int i = 0; i < *nrp; i++)
	 {
            valFrac = new double[nrelem];

	    itemp[0] = 0;
	    for(int j = 0; j < nrelem; j++)
               valFrac[j] = elemval[j]->at(i);

            if(constraint == 0)
               cout << "Constraint = " << fcnConstrainedFunc(valFrac) << endl;
            else
               cout << "Constraint = " << fcnFunc(valFrac) << endl;

            *pvalue = sim->Chi2TestX(data, *chi2value, *ndfvalue, chiGood, "WW");
            cout << "Chi2 = " << *chi2value << ", p-value = " << *pvalue << ", NDF = " << *ndfvalue << ", igood = " << chiGood << endl;
//	    *chi2value = PearsonChi2(sim, data, 1);
//	    cout << "Pearson Chi2 = " << *chi2value << endl;
	    chi2vect->push_back(*chi2value);
	    pvalvect->push_back(*pvalue);

	    delete[] valFrac;
	 }

	 // Find the minimal value for each parameter (1D)
	 itemp[0] = FindMinElement(chi2vect);
	 cout << endl << "Minimal value: " << itemp[0] << ", chi2 = " << chi2vect->at(itemp[0]) << ", p-value = " << pvalvect->at(itemp[0]) << endl;
	 for(int j = 0; j < nrelem; j++)
            cout << ptype->GetName(includePart->at(j)) << " = " << elemval[j]->at(itemp[0]) << endl;
	 cout << endl;

	 // Plot chi2 and p-value with respect to elemental fraction (1D)
         TCanvas *c1 = new TCanvas("c1","",1200,800);
         TCanvas *c2 = new TCanvas("c2","",1200,800);
	 TGraph *plotchi2dep[40];
	 TGraph *plotpdep[40];

	 itemp[0] = 0;
	 itemp[1] = 0;
	 for(int j = 0; j < nrelem; j++)
	 {
	    plotchi2dep[j] = new TGraph();
	    plotpdep[j] = new TGraph();

	    if(fixVal[j] == -1)
	    {
               for(int i = 0; i < *nrp; i++)
               {
                  plotchi2dep[j]->SetPoint(i, elemval[j]->at(i), chi2vect->at(i));
                  plotpdep[j]->SetPoint(i, elemval[j]->at(i), pvalvect->at(i));
	       }

               c1->cd();
               plotchi2dep[j]->SetMarkerStyle(20);
               plotchi2dep[j]->SetMarkerSize(1.2);
               plotchi2dep[j]->SetLineWidth(2);
	       plotchi2dep[j]->GetXaxis()->SetRangeUser(limitchi[0],limitchi[1]-elemSub[4]);
	       mystyle->SetColorScale(plotchi2dep[j], j, nrelem);

               if(itemp[0] == 0)
	       {
                  mystyle->SetAxisTitles(plotchi2dep[j], "Elemental fractions", "Chi-Square value");
	          plotchi2dep[j]->Draw("AL");
	          itemp[0]++;
	       }
	       else
	       {
                  mystyle->SetAxisTitles(plotchi2dep[j], "Elemental fractions", "Chi-Square value");
                  plotchi2dep[j]->Draw("L;SAME");
	       }

	       c2->cd();
               plotpdep[j]->SetMarkerStyle(20);
               plotpdep[j]->SetMarkerSize(1.2);
               plotpdep[j]->SetLineWidth(2);
	       plotpdep[j]->GetXaxis()->SetRangeUser(limitchi[0],limitchi[1]-elemSub[4]);
	       mystyle->SetColorScale(plotpdep[j], j, nrelem);

               if(itemp[1] == 0)
	       {
                  mystyle->SetAxisTitles(plotpdep[j], "Elemental fractions", "p-value");
	          plotpdep[j]->Draw("AL");
	          itemp[1]++;
	       }
	       else
	       {
                  mystyle->SetAxisTitles(plotpdep[j], "Elemental fractions", "p-value");
                  plotpdep[j]->Draw("L;SAME");
	       }
	    }
	    else
	    {
               plotchi2dep[j]->SetMarkerStyle(20);
               plotchi2dep[j]->SetMarkerSize(1.2);
               plotchi2dep[j]->SetLineWidth(2);
	       mystyle->SetColorScale(plotchi2dep[j], j, nrelem);

               plotpdep[j]->SetMarkerStyle(20);
               plotpdep[j]->SetMarkerSize(1.2);
               plotpdep[j]->SetLineWidth(2);
	       mystyle->SetColorScale(plotpdep[j], j, nrelem);
	    }
	 }

	 c1->cd();
	 c1->Update();
	 for(int j = 0; j < nrelem; j++)
	 {
            dtemp[0] = (gPad->GetUymax()) + 0.025*nrelem*((gPad->GetUymax()) - (gPad->GetUymin()));
            plotchi2dep[j]->GetYaxis()->SetRangeUser((gPad->GetUymin()), dtemp[0]);
	 }
	 c2->cd();
	 c2->Update();
	 for(int j = 0; j < nrelem; j++)
	 {
            dtemp[0] = (gPad->GetUymax()) + 0.025*nrelem*((gPad->GetUymax()) - (gPad->GetUymin()));
            plotpdep[j]->GetYaxis()->SetRangeUser((gPad->GetUymin()), dtemp[0]);
	 }

         TLegend *legend;
         int c_Legend = TColor::GetColor("#ffff66");

	 c1->cd();
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrelem)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
         legend->SetFillStyle(1001);
         legend->SetFillColor(c_Legend);
	 for(int j = 0; j < nrelem; j++)
	 {
            if(fixVal[j] == -1)
	    {
               stemp[0] = ptype->GetName(includePart->at(j)) + " (free, min = " + ToString(elemval[j]->at(FindMinElement(chi2vect)), 4) + ")";
               legend->AddEntry(plotchi2dep[j], stemp[0].c_str(), "lp");
	    }
	    else
	    {
               stemp[0] = ptype->GetName(includePart->at(j)) + " (fixed, " + ToString(elemval[j]->at(0), 4) + ")";
               legend->AddEntry(plotchi2dep[j], stemp[0].c_str(), "lp");
	    }
	 }
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");

	 stemp[0] = RemoveFilename(&filename) + "/fithist/chi2_dependence_" + ToString(nrelem) + "var";
	 for(int j = 0; j < nrelem; j++)
	 {
            if(fixVal[j] == -1)
               stemp[0] += "_" + ptype->GetShortName(includePart->at(j));
	 }
	 stemp[0] += ".pdf";
         c1->SaveAs(stemp[0].c_str());

	 c2->cd();
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(nrelem)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
         legend->SetFillStyle(1001);
         legend->SetFillColor(c_Legend);
	 for(int j = 0; j < nrelem; j++)
	 {
            if(fixVal[j] == -1)
	    {
               stemp[0] = ptype->GetName(includePart->at(j)) + " (free, max = " + ToString(elemval[j]->at(FindMinElement(chi2vect)), 4) + ")";
               legend->AddEntry(plotpdep[j], stemp[0].c_str(), "lp");
	    }
	    else
	    {
               stemp[0] = ptype->GetName(includePart->at(j)) + " (fixed, " + ToString(elemval[j]->at(0), 4) + ")";
               legend->AddEntry(plotpdep[j], stemp[0].c_str(), "lp");
	    }
	 }
         legend->SetBorderSize(1);
         legend->SetMargin(0.3);
         legend->Draw("same");

	 stemp[0] = RemoveFilename(&filename) + "/fithist/p-value_dependence_" + ToString(nrelem) + "var";
	 for(int j = 0; j < nrelem; j++)
	 {
            if(fixVal[j] == -1)
               stemp[0] += "_" + ptype->GetShortName(includePart->at(j));
	 }
	 stemp[0] += ".pdf";
         c2->SaveAs(stemp[0].c_str());

	 delete c1;
	 delete c2;

	 delete nrp;
	 delete[] elemSub;
	 delete pvalvect;
	 for(int i = 0; i < nrelem; i++)
	    delete elemval[i];

         chi2vect->clear();

	 // Select grid points that will be used for determining chi2 on a 2D surface
	 for(int i = 0; i < (*elemsteps); i++)
	 {
            gridval[0]->push_back(((double)i*ratioStep+limitchi[0]));
            gridval[1]->push_back(((double)i*ratioStep+limitchi[0]));
	 }

	 cout << "Number of points (X) = " << gridval[0]->size() << endl;
	 cout << "Number of points (Y) = " << gridval[1]->size() << endl;

	 // Prepare free parameters and calculate chi2 for them (2D)
	 for(int i = 0; i < gridval[0]->size(); i++)
	 {
            for(int j = 0; j < gridval[1]->size(); j++)
	    {
               valFrac = new double[nrelem];

	       itemp[0] = 0;
	       valFrac[0] = gridval[0]->at(i);
	       valFrac[1] = gridval[1]->at(j);

	       if((gridval[0]->at(i) != 0) || (gridval[1]->at(j) != 0))
	       {
/*                  if(constraint == 0)
                     cout << "Constraint = " << fcnConstrainedFunc(valFrac) << endl;
                  else
                     cout << "Constraint = " << fcnFunc(valFrac) << endl;*/
	          if(constraint == 0)
                     fcnConstrainedFunc(valFrac);
	          else
                     fcnFunc(valFrac);

                  *pvalue = sim->Chi2TestX(data, *chi2value, *ndfvalue, chiGood, "WW");
//                  cout << "Chi2 = " << *chi2value << ", p-value = " << *pvalue << ", NDF = " << *ndfvalue << ", igood = " << chiGood << endl;
//	          *chi2value = PearsonChi2(sim, data, 1);
//	          cout << "Pearson Chi2 = " << *chi2value << endl;
	          chi2vect->push_back(*chi2value);
	       }

	       delete[] valFrac;
	    }
	 }

	 // Plot chi2 on a surface plot (2D)
         c1 = new TCanvas("c1","",1000,800);
	 gPad->SetRightMargin(0.18);
 	 TGraph2D *plotchi2surf = new TGraph2D();

	 cout << "Number of points (chi2) = " << chi2vect->size() << endl;

	 itemp[0] = 0;

	 for(int i = 0; i < gridval[0]->size(); i++)
	 {
            for(int j = 0; j < gridval[1]->size(); j++)
	    {
	       if((gridval[0]->at(i) != 0) || (gridval[1]->at(j) != 0))
	       {
                  cout << itemp[0] << "\t" << gridval[0]->at(i) << "\t" << gridval[1]->at(j) << "\t" << chi2vect->at(itemp[0]) << endl;
                  plotchi2surf->SetPoint(itemp[0], gridval[0]->at(i), gridval[1]->at(j), chi2vect->at(itemp[0]));
		  itemp[0]++;
	       }
	    }
         }

	 for(int j = 0; j < nrelem; j++)
	 {
            if(fixVal[j] == -1)
	    {
	       stemp[j] = "Elemental fraction (" + ptype->GetShortName(includePart->at(j)) + ")";
	    }
	 }

	 itemp[0] = FindMinElement(chi2vect);
	 cout << endl << "Minimal value: " << itemp[0] << ", chi2 = " << chi2vect->at(itemp[0]) << endl;

	 plotchi2surf->SetNpx(100);
	 plotchi2surf->SetNpy(100);
         mystyle->SetAxisTitles(plotchi2surf, stemp[0].c_str(), stemp[1].c_str(), "Chi-Square value");
	 plotchi2surf->GetHistogram()->GetZaxis()->SetRangeUser(chi2vect->at(itemp[0]),5*chi2vect->at(itemp[0]));
	 plotchi2surf->Draw("colz");

         TLine *line = new TLine(limitchi[0], limitchi[1], limitchi[1], limitchi[0]);
         line->SetLineWidth(1);
         line->SetLineStyle(1);
         line->SetLineColor(1);
         line->Draw("SAME");

	 stemp[0] = RemoveFilename(&filename) + "/fithist/chi2_surface_plot_" + ToString(nrelem) + "var";
	 for(int j = 0; j < nrelem; j++)
	 {
            if(fixVal[j] == -1)
               stemp[0] += "_" + ptype->GetShortName(includePart->at(j));
	 }
	 stemp[0] += ".pdf";
         c1->SaveAs(stemp[0].c_str());

	 delete c1;

	 delete elemsteps;
	 delete[] fixVal;
	 delete chi2vect;
	 for(int i = 0; i < nrelem; i++)
            delete gridval[i];
      }

      return -1;
   }
   // Find the best fit with TFractionFitter
   else if(fitproc == 5)
   {
      // Prepare variables to hold best parameter fits and their errors
      double *bestfrac = new double[nrelem];
      double *bestfracerr = new double[2*nrelem];

      // Prepare the simulation histograms
      mc =  new TObjArray(nrelem);
      for(int i = 0; i < nrelem; i++)
         mc->Add(simhist[i]);

      // Prepare the FractionFitter class and set constraints
      fracFitter = new TFractionFitter(data, mc);

      if(rangeTransform == 0)
      {
         for(int i = 0; i < nrelem; i++)
            fracFitter->Constrain(i, 0., 1.);
      }
      
      // Exclude any bins that are 0 in simulations or data
      for(int i = 1; i < nrbins+1; i++)
      {
         dtemp[0] = 0;
         for(int j = 0; j < nrelem; j++)
            dtemp[0] += simhist[j]->GetBinContent(i);

	 if( (data->GetBinContent(i) == 0) || (dtemp[0] == 0) )
	 {
            cout << "Excluding bin " << i-1 << " from further analysis." << endl;
	    fracFitter->ExcludeBin(i);
	 }
      }

      itemp[0] = fracFitter->Fit();
      fracFitter->ErrorAnalysis(1);
      cout << "Fitting status = " << itemp[0] << endl;

      if(itemp[0] == 0)
      {
/*         cout << "Stats " << (fracFitter->GetFitter())->GetStats(dtempres[0],dtempres[1],dtempres[2],itempres[0],itempres[1]) << ": ";
	 cout << "amin = " << dtempres[0] << ", edm = " << dtempres[1] << ", errdef = " << dtempres[2] << ", nvpar = " << itempres[0] << ", nparx = " << itempres[1] << endl;

         cout << "Results:" << endl;
         (fracFitter->GetFitter())->PrintResults(3,dtempres[0]);

	 cout << "Error simple (par0) = " << (fracFitter->GetFitter())->GetParError(0) << endl;
	 cout << "Error info " << (fracFitter->GetFitter())->GetErrors(0,dtempres[0],dtempres[1],dtempres[2],dtempres[3]) << ": ";
	 cout << "E+ = " << dtempres[0] << ", E- = " << dtempres[1] << ", E2 = " << dtempres[2] << ", globcc = " << dtempres[3] << endl;

	 dtemp[0] = (fracFitter->GetFitter())->GetParError(0)/dtempres[2];
	 if((dtempres[0] == 0) || (dtempres[1] == 0))
            cout << "Wrong parabolic error estimation (" << dtemp[0] << ", " << dtempres[0] << ", " << dtempres[1] << ")!" << endl;
	 else
            cout << "Correct parabolic error estimation (" << dtemp[0] << ")!" << endl;
	 cout << endl;

	 cout << "Error simple (par1) = " << (fracFitter->GetFitter())->GetParError(1) << endl;
	 cout << "Error info " << (fracFitter->GetFitter())->GetErrors(1,dtempres[0],dtempres[1],dtempres[2],dtempres[3]) << ": ";
	 cout << "E+ = " << dtempres[0] << ", E- = " << dtempres[1] << ", E2 = " << dtempres[2] << ", globcc = " << dtempres[3] << endl;

	 dtemp[0] = (fracFitter->GetFitter())->GetParError(1)/dtempres[2];
	 if((dtempres[0] == 0) || (dtempres[1] == 0))
            cout << "Wrong parabolic error estimation (" << dtemp[0] << ", " << dtempres[0] << ", " << dtempres[1] << ")!" << endl;
	 else
            cout << "Correct parabolic error estimation (" << dtemp[0] << ")!" << endl;
	 cout << endl;*/

/*         sim = (TH1F*) fracFitter->GetPlot();
	 *chi2value = fracFitter->GetChisquare();
	 *ndfvalue = fracFitter->GetNDF();
	 *pvalue = fracFitter->GetProb();
         cout << "Chi2/NDF = " << *chi2value << "/" << *ndfvalue << endl;
         cout << "Probability = " << *pvalue << endl;*/

	 for(int i = 0; i < nrelem; i++)
	 {
//            fracFitter->GetResult(i, bestfrac[i], bestfracerr[i]);
            bestfrac[i] = (fracFitter->GetFitter())->GetParameter(i);
	    (fracFitter->GetFitter())->GetErrors(i, bestfracerr[2*i], bestfracerr[2*i+1], dtemp[1], dtemp[0]);
            dtemp[0] = (fracFitter->GetFitter())->GetParError(i)/dtemp[1];

	    cout << i << ": Parabolic error: " << (fracFitter->GetFitter())->GetParError(i) << "/" << dtemp[1] << " = " << dtemp[0] << "\t" << "\t";
            if( (dtemp[0] > 0.9) && (dtemp[0] < 1.1) )
               cout << "(GOOD)" << endl;
	    else
               cout << "(BAD)" << endl;

            if(rangeTransform == 0)
	    {
	       // Both errors missing
	       if((bestfracerr[2*i] == 0.) && (bestfracerr[2*i+1] == 0.))
	       {
                  cout << "Both errors are zero!" << endl;
                  bestfracerr[2*i] = TMath::Abs(bestfracerr[2*i]);
                  bestfracerr[2*i+1] = TMath::Abs(bestfracerr[2*i+1]);
	       }
	       // Positive error is missing
	       else if(bestfracerr[2*i] == 0.)
	       {
                  bestfracerr[2*i] = TMath::Abs(bestfracerr[2*i]);
	          cout << "Zero positive error!" << endl;

		  // when errors are missing and parabolic errors are calculated, use those
		  if( (dtemp[0] > 0.9) && (dtemp[0] < 1.1) )
		  {
                     bestfracerr[2*i] = TMath::Abs(dtemp[1]);
		     cout << "Using parabolic error instead (" << dtemp[1] << ")." << endl;
		  }
/*                  // positive error must extend up to the max value
		  if((1. - bestfrac[i]) < bestfrac[i])
 	             bestfracerr[2*i] = TMath::Abs(1. - bestfrac[i]);*/

	          bestfracerr[2*i+1] = TMath::Abs(bestfracerr[2*i+1]);
	       }
	       // Negative error
	       else if(bestfracerr[2*i+1] == 0.)
	       {
	          bestfracerr[2*i+1] = TMath::Abs(bestfracerr[2*i+1]);
	          cout << "Zero negative error!" << endl;

		  // when errors are missing and parabolic errors are calculated, use those
		  if( (dtemp[0] > 0.9) && (dtemp[0] < 1.1) )
		  {
                     bestfracerr[2*i+1] = TMath::Abs(dtemp[1]);
		     cout << "Using parabolic error instead (" << dtemp[1] << ")." << endl;
		  }
/*                  // negative error must extend down to the min value
		  if((1. - bestfrac[i]) > bestfrac[i])
	             bestfracerr[2*i+1] = TMath::Abs(bestfrac[i]);*/

                  bestfracerr[2*i] = TMath::Abs(bestfracerr[2*i]);
	       }
	       else
	       {
                  bestfracerr[2*i] = TMath::Abs(bestfracerr[2*i]);
	          bestfracerr[2*i+1] = TMath::Abs(bestfracerr[2*i+1]);
	          cout << "Both errors defined!" << endl;
	       }
	    }
	    else
	    {
               bestfracerr[2*i] = TMath::Abs(bestfracerr[2*i]);
	       bestfracerr[2*i+1] = TMath::Abs(bestfracerr[2*i+1]);
	       cout << "Both errors defined!" << endl;
	    }

            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
         }

         // Get the compositions
         midComposition->clear();
         midCompositionErr->clear();
         for(int i = 0; i < nrparam; i++)
         {
            midComposition->push_back(bestfrac[i]);
            midCompositionErr->push_back(bestfracerr[2*i]);
            midCompositionErr->push_back(bestfracerr[2*i+1]);
         }

         // Get the lnA value
         midLna[0] = 0;
         midLna[1] = 0;
         midLna[2] = 0;
         
         for(int i = 0; i < midComposition->size(); i++)
         {
            itemp[1] = includePart->at(i);
/*	    cout << midComposition->at(i)*TMath::Log(ptype->GetA(itemp[1])) << endl;
	    cout << midCompositionErr->at(i)*TMath::Log(ptype->GetA(itemp[1])) << endl;*/
            midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[1]));
            midLna[1] += (midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[1])))*(midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[1])));
            midLna[2] += (midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[1])))*(midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[1])));
         }

         midLna[1] = TMath::Sqrt(midLna[1]);
         midLna[2] = TMath::Sqrt(midLna[2]);

         cout << "lnA value = " << midLna[0] << " +/- " << midLna[1] << "/" << midLna[2] << endl;

         // Get the parameter step value (not applicable in this case)
         *midStep = -1.;
//         cout << "Parameter step = " << minstep[stepsel] << " (" << stepsel << ")" << endl;

//	 delete c1;

         PlotSumResiduals(bestfrac, bestfracerr, itemp[0], -1.);

         delete mc;
         delete fracFitter;

         delete[] bestfrac;
         delete[] bestfracerr;
      }
      else
      {
         cout << "Error! Minimization failed." << endl;

	 for(int i = 0; i < nrelem; i++)
	 {
//            fracFitter->GetResult(i, bestfrac[i], bestfracerr[i]);
            bestfrac[i] = (fracFitter->GetFitter())->GetParameter(i);
	    (fracFitter->GetFitter())->GetErrors(i, bestfracerr[2*i], bestfracerr[2*i+1], dtemp[1], dtemp[0]);
            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
	 }

         delete mc;
         delete fracFitter;

         delete[] bestfrac;
         delete[] bestfracerr;

	 return -1;
      }
   }

   return 0;
}

double MvaFitHist::GetFinalLna(int me)
{
   return midLna[me];
}

double MvaFitHist::GetFinalComposition(int type)
{
   return midComposition->at(type);
}

double MvaFitHist::GetFinalCompositionErr(int type, int hilo)
{
   return midCompositionErr->at(2*type+hilo);
}

double MvaFitHist::GetFinalEnergy()
{
   return *midEnergy;
}

// Plot data, normalized signal+background and normalized residuals
void MvaFitHist::PlotSumResiduals(double *sigFrac, double *sigFracErr, int status, double selStep)
{
   mystyle->SetBaseStyle();
   TCanvas *c1 = new TCanvas("c1","",1200,900);
//   TCanvas *c2 = new TCanvas("c2","",1200,900);
   TLatex chiText;
   c1->cd();
   TPad *upperpad = new TPad("upperpad", "upperpad", 0.004, 0.490, 0.996, 0.996);
   TPad *lowerpad = new TPad("lowerpad", "lowerpad", 0.004, 0.000, 0.996, 0.490);
   upperpad->Draw();
   lowerpad->Draw();

   // Plot the data distribution
   upperpad->cd();

   if(fitproc != 5)
   {
      mystyle->SetHistColor((TH1*)newdata, 2);
      mystyle->SetAxisTitles((TH1*)newdata, "MVA variable", "Probability density");
      newdata->GetXaxis()->SetTitleOffset(2.0);
      newdata->Draw();

      if(constraint == 0)
//         cout << "Constraint = " << fcnConstrainedFunc(sigFrac) << endl;
         cout << "Chi2+constraint^2 = " << fcn(sigFrac) << endl;
      else
         cout << "Constraint = " << fcnFunc(sigFrac) << endl;

      // Zero empty bin errors just for plotting purposes
      newsig = (TH1F*)sim->Clone("newsig");
      for(int j = 1; j < nrbins+1; j++)
      {
         if(sim->GetBinContent(j) == 0)
            newsig->SetBinError(j, 0);
      }

      yrange[0] = newdata->GetMaximum();
      yrange[1] = newsig->GetMaximum();

      mystyle->SetHistColor((TH1*)newsig, 0);
      mystyle->SetAxisTitles((TH1*)newsig, "MVA variable", "Probability density");
      newsig->SetMaximum(1.1*TMath::MaxElement(2, yrange));
      newsig->Draw("SAME");

      *pvalue = sim->Chi2TestX(data, *chi2value, *ndfvalue, chiGood, "WW", residVal);
      cout << "Chi2 = " << *chi2value << ", p-value = " << *pvalue << ", NDF = " << *ndfvalue << ", igood = " << chiGood << endl;

      cout << "Pearson's Chi2 = " << PearsonChi2(sim, data, 0) << endl;
      cout << "Chi2 distance = " << PearsonChi2(sim, data, 1) << endl;
   }
   else
   {
      sim = (TH1F*) fracFitter->GetPlot();
      *chi2value = fracFitter->GetChisquare();
      *ndfvalue = fracFitter->GetNDF();
      *pvalue = fracFitter->GetProb();
      cout << "Chi2 = " << *chi2value << ", p-value = " << *pvalue << ", NDF = " << *ndfvalue << endl;

      yrange[0] = sim->GetMaximum();
      yrange[1] = data->GetMaximum();

      mystyle->SetHistColor((TH1*)sim, 0);
      mystyle->SetAxisTitles((TH1*)sim, "MVA variable", "Number of events");
      sim->Draw();
      
      mystyle->SetHistColor((TH1*)data, 2);
      mystyle->SetAxisTitles((TH1*)data, "MVA variable", "Number of events");
      sim->SetMaximum(1.1*TMath::MaxElement(2, yrange));
      data->Draw("Ep;same");

      for(int i = 1; i <= nrbins; i++)
      {
         if(sim->GetBinContent(i) == 0)
            residVal[i-1] = 0;
	 else
	 {
            residVal[i-1] = (data->GetBinContent(i) - sim->GetBinContent(i))/TMath::Sqrt(sim->GetBinContent(i));

//            cout << i-1 << ": Values = " << sim->GetBinContent(i) << ", " << data->GetBinContent(i) << ", " << residVal[i-1] << endl;
	 }
      }

/*      dtemp[0] = TMath::RMS(nrbins, residVal);
      cout << "RMS = " << dtemp[0] << endl;
      for(int i = 0; i < nrbins; i++)
      {
         residVal[i] = residVal[i]/dtemp[0];
	 cout << i << ": " << residVal[i] << endl;
      }*/
   }

   cout << "Getting residuals" << endl;
   TH1F *resid = new TH1F("resid", "", nrbins, xlim[0], xlim[1]);
   yresidrange[0] = 1.e+24;
   yresidrange[1] = -1.e+24;
//   vector<double> *residAll, *residSig, *residBack;
   vector<double> *residAll;
   residAll = new vector<double>;
/*   residSig = new vector<double>;
   residBack = new vector<double>;*/
   for(int i = 1; i <= nrbins; i++)
   {
      // Only set residuals, if any of the two histogram bins are not zero
      if( ((double)sim->GetBinContent(i) != 0) && ((double)data->GetBinContent(i) != 0) )
      {
/*         if(resid->GetBinCenter(i) > mvacut)
            residSig->push_back(residVal[i-1]);
         else
            residBack->push_back(residVal[i-1]);*/

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

   double *tempd = new double[residAll->size()];
/*   c2->cd();
   TH1F *residDist = new TH1F("residDist", "", nrbins, -1.1*TMath::MaxElement(2, yresidrange), 1.1*TMath::MaxElement(2, yresidrange));
   for(int i = 0; i < residAll->size(); i++)
   {
      residDist->Fill(residAll->at(i));
   }

   mystyle->SetHistColor((TH1*)residDist, 3);
   mystyle->SetAxisTitles((TH1*)residDist, "Adjusted residuals", "Nr. of residuals");
   residDist->Draw();

   c2->Update();*/

/*   chiText.SetTextAlign(31);
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
   cout << "Residual back-sig: mean = " << TMath::Abs(dtemp[0] - dtemp[1]) << endl;*/

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
      dtemp[2] = 0.;
      for(int i = 0; i < nrparam; i++)
      {
         if(constraint == 0)
	 {
            stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames[i] + ")";
            chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
	 }
	 else
	 {
            if(i < nrparam-1)
            {
               dtemp[0] -= sigFrac[i];
               dtemp[1] += sigFracErr[2*i];
               dtemp[2] += sigFracErr[2*i+1];
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
            else
            {
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(dtemp[0], 4) + " +" + ToString(dtemp[1], 4) + "/-" + ToString(dtemp[2], 4) + " (" + treeNames[i] + ")";
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
      dtemp[2] = 0.;
      for(int i = 0; i < nrparam; i++)
      {
         if(constraint == 0)
	 {
            stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames[i] + ")";
            chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
	 }
	 else
	 {
            if(i < nrparam-1)
            {
               dtemp[0] -= sigFrac[i];
               dtemp[1] += sigFracErr[2*i];
               dtemp[2] += sigFracErr[2*i+1];
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
            else
            {
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(dtemp[0], 4) + " +" + ToString(dtemp[1], 4) + "/-" + ToString(dtemp[2], 4) + " (" + treeNames[i] + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
	 }
      }

      stemp[0] = "Residuals: mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4) + ", RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      stemp[0] = "MVA cut = " + ToString(mvacut, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else if(fitproc == 5)
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString((*chi2value), 4) + "/" + ToString((*ndfvalue)) + " = " + ToString((*chi2value)/(double)(*ndfvalue), 4) + ", status = " + ToString(status) + ", p-value = " + ToString((*pvalue), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      dtemp[0] = 1.;
      dtemp[1] = 0.;
      dtemp[2] = 0.;
      for(int i = 0; i < nrparam; i++)
      {
         stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames[i] + ")";
         chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      }

      stemp[0] = "Residuals: mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4) + ", RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      stemp[0] = "MVA cut = " + ToString(mvacut, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)nrparam+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }

   delete residAll;
/*   delete residSig;
   delete residBack;*/
   delete tempd;

   lowerpad->cd();
   mystyle->SetHistColor((TH1*)resid, 3);
   if(fitproc != 5)
     mystyle->SetAxisTitles((TH1*)resid, "", "Adjusted residuals");
   else
     mystyle->SetAxisTitles((TH1*)resid, "", "Standardized residuals");
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

   if(fitproc == 0)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_" + ToString(100*sigFrac[0], 0) + ".pdf";
      c1->SaveAs(stemp[0].c_str());
/*      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());*/
   }
   else if(fitproc == 1)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_minuit2.pdf";
      c1->SaveAs(stemp[0].c_str());
/*      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());*/
   }
   else if(fitproc == 2)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_manual.pdf";
      c1->SaveAs(stemp[0].c_str());
/*      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());*/
   }
   else if(fitproc == 5)
   {
      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_fracfitter.pdf";
      c1->SaveAs(stemp[0].c_str());
/*      stemp[0] = RemoveFilename(&filename) + "/fithist/" + stemp[1] + "_composition_residuals_distribution.pdf";
      c2->SaveAs(stemp[0].c_str());*/
   }
   
   delete c1;
//   delete c2;
}

int MvaFitHist::GetNrElem()
{
   return nrelem;
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
      double *dtemp = new double[7];
      string *stemp = new string[3];
      vector<double> *finalLna = new vector<double>;
      vector<double> *finalChi2 = new vector<double>;
      vector<double> *finalPvalue = new vector<double>;
      vector<double> *finalNdf = new vector<double>;
      vector<double> *finalComp = new vector<double>;
      vector<double> *finalEnergy = new vector<double>;
      vector<double> *finalStep = new vector<double>;

      string inname;

      cerr << "Perform histogram fit on a range of ratios (0), perform a minimization (1), select a composition to plot (2), perform a minimization with the new approach (3), plot parameters versus Chi2 value (4) or perform a minimization with TFractionFitter (5)? ";
      cin >> itemp[0];
      
      if( (itemp[0] == 0) || (itemp[0] == 4) )
      {
         cerr << "What should be the step size (between 0 and 1)? ";
	 cin >> dtemp[0];

	 while( (dtemp[0] > 1) || (dtemp[0] <= 0) )
	 {
            cerr << "What should be the step size (between 0 and 1)? ";
	    cin >> dtemp[0];
	 }
      }
      else if( (itemp[0] != 0) && (itemp[0] != 1) && (itemp[0] != 2) && (itemp[0] != 3) && (itemp[0] != 4) && (itemp[0] != 5) )
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
	       cout << dtemp[0] << "\t";
	       finalComp->push_back(dtemp[0]);
               dtemp[0] = fithist->GetFinalCompositionErr(j,0);		// 1 is low error
	       cout << dtemp[0] << "\t";
	       finalComp->push_back(dtemp[0]);
               dtemp[0] = fithist->GetFinalCompositionErr(j,1);		// 0 is high error
	       cout << dtemp[0] << endl;
	       finalComp->push_back(dtemp[0]);
	    }

            dtemp[0] = fithist->GetFinalLna(0);
	    cout << dtemp[0] << "\t";
	    finalLna->push_back(dtemp[0]);
            dtemp[0] = fithist->GetFinalLna(1);
	    cout << dtemp[0] << "\t";
	    finalLna->push_back(dtemp[0]);
            dtemp[0] = fithist->GetFinalLna(2);
	    cout << dtemp[0] << endl;
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

      if( (itemp[0] == 2) || (itemp[0] == 0) || (itemp[0] == 4) )
         return 0;

      // Printout final results for each file
      double *xval = new double[nrp];
      double *xerr = new double[nrp];
      double *yval = new double[nrp];
      double *yerrlo = new double[nrp];
      double *yerrhi = new double[nrp];

      PrimPart *ptype = new PrimPart();

      cout << endl << "Final results:" << endl;
      cout << "# nr. elements = " << fithist->GetNrElem() << ", nr. bins = " << fithist->GetNrBins() << endl;
      cout << "# E   \tChi2/NDF \tp-val\tstep\tlnA \t+   \t-   \t";
      for(int i = 0; i < fithist->GetNrElem(); i++)
         cout << "\t" << ptype->GetShortName(fithist->GetElemType(i)) << "  \t+   \t-   ";
      cout << endl;
      cout.precision(4);
      for(int i = 0; i < nrp; i++)
      {
	 itemp[0] = 3*(fithist->GetNrElem());
	 cout << fixed << finalEnergy->at(i) << "\t" << finalChi2->at(i) << "/" << (int)finalNdf->at(i) << "\t" << finalPvalue->at(i) << "\t" << finalStep->at(i) << "\t" << finalLna->at(3*i) << "\t" << finalLna->at(3*i+1) << "\t" << finalLna->at(3*i+2) << "\t";

         xval[i] = finalEnergy->at(i);
	 xerr[i] = 0.;
	 yval[i] = finalLna->at(3*i);
	 yerrhi[i] = finalLna->at(3*i+1);
	 yerrlo[i] = finalLna->at(3*i+2);

         for(int j = 0; j < itemp[0]; j++)
            cout << "\t" << finalComp->at(itemp[0]*i+j);

	 cout << endl;
      }

      cout << endl;

      RootStyle *mystyle = new RootStyle();
      mystyle->SetBaseStyle();
      TCanvas *c1 = new TCanvas("c1","",1200,900);

      // Prepare lnA plot, with appended published values
      TGraphAsymmErrors *grLna;
      TGraphAsymmErrors *grPubLna;
      TLegend *legend;
      int c_Legend = TColor::GetColor("#ffff66");

      grLna = new TGraphAsymmErrors(nrp, xval, yval, xerr, xerr, yerrlo, yerrhi);
      mystyle->SetGraphColor(grLna, 2);
      grLna->GetXaxis()->SetRangeUser(18.5, 20.0);
//      grLna->GetXaxis()->SetLimits(18.5, 20.0);
      grLna->GetYaxis()->SetRangeUser(-0.7, 5.);

      cerr << "Which dataset did you use (EPOS = 0, QGSJET = 1, SIBYLL = 2, NONE = -1)? ";
      cin >> itemp[2];

      bool addPublished;

      float *xbinPub[3];
      float *ybinPubLna[3];

      if(itemp[2] != -1)
      {
         addPublished = true;
         vector<float> returnVal;
         itemp[2] = ReadLnaResults(&returnVal, itemp[2]);
//         cout << "Number of points = " << itemp[2] << endl;
         
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
         
//            cout << i+1 << ", data: " << xbinPub[0][i] << "\t" << ybinPubLna[0][i] << "\t" << ybinPubLna[1][i] << endl;
         }
         
         grPubLna = new TGraphAsymmErrors(itemp[2], xbinPub[0], ybinPubLna[0], xbinPub[1], xbinPub[2], ybinPubLna[1], ybinPubLna[2]);
         mystyle->SetGraphColor(grPubLna, 0);
         grPubLna->GetYaxis()->SetRangeUser(-0.7, 5.);
      }
      else
         addPublished = false;

      mystyle->SetAxisTitles(grLna, "FD energy [log(E/eV)]", "<lnA> of data events");
      grLna->Draw("APL");
      if(addPublished)
         grPubLna->Draw("PL;SAME");

      if(addPublished)
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(2)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
      else
         legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(1)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
      legend->SetFillStyle(1001);
      legend->SetFillColor(c_Legend);
      legend->AddEntry(grLna, "Data (Histogram fit)", "lp");
      if(addPublished)
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
     
      if(addPublished)
      {
         for(int i = 0; i < 3; i++)
         {
            delete[] xbinPub[i];
            delete[] ybinPubLna[i];
         }
      }

      delete c1;

      cout << endl;

      double *ybinComp[40];
      double *ybinCompErrlo[40];
      double *ybinCompErrhi[40];

      // Plot composition plots (all together and separate)
      c1 = new TCanvas("c1","",1200,900);

      TGraphAsymmErrors *grComp[40];

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         ybinComp[j] = new double[nrp];
         ybinCompErrlo[j] = new double[nrp];
         ybinCompErrhi[j] = new double[nrp];

         for(int i = 0; i < nrp; i++)
         {
	    itemp[0] = 3*(fithist->GetNrElem());
	    ybinComp[j][i] = finalComp->at(3*j+itemp[0]*i);
	    ybinCompErrhi[j][i] = finalComp->at(3*j+1+itemp[0]*i);
	    ybinCompErrlo[j][i] = finalComp->at(3*j+2+itemp[0]*i);
//	    cout << "i = " << i << ", j = " << j << ": " << ybinComp[j][i] << "\t" << ybinCompErrlo[j][i] << "\t" << ybinCompErrhi[j][i] << endl;
         }
      }

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         grComp[j] = new TGraphAsymmErrors(nrp, xval, ybinComp[j], xerr, xerr, ybinCompErrlo[j], ybinCompErrhi[j]);
//         mystyle->SetGraphColor(grComp[j], 2);
         grComp[j]->SetMarkerStyle(20);
         grComp[j]->SetMarkerSize(1.4);
         grComp[j]->SetLineWidth(2);
         grComp[j]->GetYaxis()->SetRangeUser(0., 1.2);
	 mystyle->SetColorScale(grComp[j], j, fithist->GetNrElem());

         if(j == 0)
	 {
            mystyle->SetAxisTitles(grComp[j], "FD energy [log(E/eV)]", "Elemental fractions");
	    grComp[j]->Draw("ALP");
	 }
	 else
	 {
            mystyle->SetAxisTitles(grComp[j], "FD energy [log(E/eV)]", "Elemental fractions");
            grComp[j]->Draw("LP;SAME");
	 }
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
         grComp[j]->GetYaxis()->SetRangeUser(0., 1.);
         grComp[j]->Draw("ALP");
         stemp[1] = RemoveFilename(&inname);
         stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_" + ptype->GetName(fithist->GetElemType(j)) + ".pdf";
         c1->SaveAs(stemp[2].c_str());
      }

      // Plot all four on separate pads
      TCanvas *c2 = new TCanvas("c2", "", 1200, 400.*(fithist->GetNrElem()));
      c2->Divide(1,fithist->GetNrElem());

      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
	 c2->cd(j+1);
         mystyle->SetAxisTitles(grComp[j], "", "Elemental fractions");
	 grComp[j]->GetYaxis()->SetTitleOffset(1.9 + (fithist->GetNrElem()-2)*0.45);
         grComp[j]->Draw("ALP");
	 grComp[j]->GetYaxis()->SetRangeUser(0.,1.);
//	 grComp[j]->GetYaxis()->SetRangeUser(-0.2,1.2);	// TFractionFitter limits
      }

      c2->Update();

      t->SetTextAlign(31);
      t->SetTextColor(1);
      t->SetTextSize(17);
      for(int j = 0; j < fithist->GetNrElem(); j++)
      {
         stemp[1] = ptype->GetName(fithist->GetElemType(j));
	 c2->cd(j+1);
         t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), stemp[1].c_str());
      }

      stemp[1] = RemoveFilename(&inname);
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_combine.pdf";
      c2->SaveAs(stemp[2].c_str());

      delete c2;

      // Plot composition plot with light and heavy compositions together
      c1->cd();
      double compMassBreak = 7.4;

      for(int i = 0; i < nrp; i++)
      {
	 dtemp[1] = 0.;
	 dtemp[2] = 0.;
	 dtemp[3] = 0.;
	 dtemp[4] = 0.;
	 dtemp[5] = 0.;
	 dtemp[6] = 0.;
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
               dtemp[2] += (ybinCompErrlo[j][i])*(ybinCompErrlo[j][i]);
               dtemp[3] += (ybinCompErrhi[j][i])*(ybinCompErrhi[j][i]);

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
               dtemp[4] += ybinComp[j][i];
               dtemp[5] += (ybinCompErrlo[j][i])*(ybinCompErrlo[j][i]);
               dtemp[6] += (ybinCompErrhi[j][i])*(ybinCompErrhi[j][i]);

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
         ybinCompErrlo[0][i] = TMath::Sqrt(dtemp[2]);
         ybinCompErrhi[0][i] = TMath::Sqrt(dtemp[3]);
         ybinComp[1][i] = dtemp[4];
         ybinCompErrlo[1][i] = TMath::Sqrt(dtemp[5]);
         ybinCompErrhi[1][i] = TMath::Sqrt(dtemp[6]);
      }

      grComp[0] = new TGraphAsymmErrors(nrp, xval, ybinComp[0], xerr, xerr, ybinCompErrlo[0], ybinCompErrhi[0]);
      grComp[0]->SetMarkerStyle(20);
      grComp[0]->SetMarkerSize(1.4);
      grComp[0]->SetLineWidth(2);
      grComp[0]->GetYaxis()->SetRangeUser(0., 1.2);
//      grComp[0]->GetYaxis()->SetRangeUser(-0.2,1.2);	// TFractionFitter limits
      mystyle->SetColorScale(grComp[0], 0, fithist->GetNrElem());
      mystyle->SetAxisTitles(grComp[0], "FD energy [log(E/eV)]", "Elemental fractions");
      grComp[0]->Draw("ALP");
      grComp[1] = new TGraphAsymmErrors(nrp, xval, ybinComp[1], xerr, xerr, ybinCompErrlo[1], ybinCompErrhi[1]);
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
         delete[] ybinCompErrlo[j];
         delete[] ybinCompErrhi[j];
      }

      delete grComp[0];
      delete grComp[1];

      delete c1;

      delete mystyle;

      delete[] xval;
      delete[] xerr;
      delete[] yval;
      delete[] yerrlo;
      delete[] yerrhi;

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
