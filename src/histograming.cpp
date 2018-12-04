/*#define _STANDALONE_ 1*/
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "histograming.h"

using namespace std;

ScatHist::ScatHist(int *inType)
{
   stemp = new string[3];
   itemp = new int[4];
   dtemp = new double[3];

   yrange = new double[2];

   uflowcount = new int[MAXINFILES];
   oflowcount = new int[MAXINFILES];
   totcount[0] = new int[MAXINFILES];
   totcount[1] = new int[MAXINFILES];
   totcount[2] = new int[MAXINFILES];
   totcount[3] = new int[MAXINFILES];

   type = *inType;

   mystyle = new RootStyle();

/*   xlim = new double[2];
   ylim = new double[2];

   yrange = new double[2];
   yresidrange = new double[2];

   nrkeys = new int;
   norm = new double[MAXINFILES];
   nentries = new int[MAXINFILES];

   filename = new string;

   trees = new vector<int>;
   treeNames = new vector<string>;

   mystyle = new RootStyle();
   ptype = new PrimPart();

   mvacut = new double;
   nrsim = new int[MAXINFILES];

   midLna = new double[3];
   midComposition = new vector<double>;
   midCompositionErr = new vector<double>;
   chi2value = new double;
   pvalue = new double;
   ndfvalue = new int;
   midEnergy = new double;
   midStep = new double;

   chiGood = new int;*/
}

ScatHist::~ScatHist()
{
   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   delete[] xlim;
   delete[] ylim;

   delete[] yrange;

   delete[] uflowcount;
   delete[] oflowcount;
   delete[] totcount[0];
   delete[] totcount[1];
   delete[] totcount[2];
   delete[] totcount[3];

   delete filename;
   delete obser;
   delete trees;

   delete nrbins;
   delete[] otherSettings;

   delete mystyle;

/*   delete[] xlim;
   delete[] ylim;

   delete[] yrange;
   delete[] yresidrange;

   delete nrkeys;
   delete[] norm;
   delete[] nentries;

   delete filename;

   delete trees;
   delete treeNames;

   delete mystyle;
   delete ptype;

   delete mvacut;
   delete[] nrsim;

   delete[] midLna;
   delete midComposition;
   delete midCompositionErr;
   delete chi2value;
   delete pvalue;
   delete ndfvalue;
   delete midEnergy;
   delete midStep;

   delete chiGood;

   delete includePart;
   delete nrelem;
   delete method;
   delete[] mvalim;
   delete himodel;
   delete simProd;
   delete dataNum;
   delete nrbins;
   delete[] otherSettings;
   delete nrparam;

   delete[] residVal;

   delete compnrp;
   delete[] xval;
   delete[] xerr;
   delete[] yval;
   delete[] yerrlo;
   delete[] yerrhi;
   delete composition;*/
}

void ScatHist::SetFilenames(vector<string> *inFiles)
{
   filename = new vector<string>;
   for(int i = 0; i < inFiles->size(); i++)
      filename->push_back(inFiles->at(i));
}

void ScatHist::SetObservables(vector<string> *inObs)
{
   obser = new vector<string>;
   for(int i = 0; i < inObs->size(); i++)
      obser->push_back(inObs->at(i));
}

void ScatHist::SetTrees(vector<string> *inTrees)
{
   trees = new vector<string>;
   for(int i = 0; i < inTrees->size(); i++)
      trees->push_back(inTrees->at(i));
}

void ScatHist::SetNrBins(int *inNrBins)
{
   nrbins = new int;
   *nrbins = *inNrBins; 
}

void ScatHist::SetXaxisLimits(bool use, double *inLimits)
{
   if(use)
   {
      xlim[0] = inLimits[0]; 
      xlim[1] = inLimits[1]; 

      cout << "Custom X-axis limits: " << xlim[0] << ", " << xlim[1] << endl;
   }
}

void ScatHist::SetYaxisLimits(bool use, double *inLimits)
{
   if(use)
   {
      ylim[0] = inLimits[0]; 
      ylim[1] = inLimits[1]; 

      cout << "Custom Y-axis limits: " << ylim[0] << ", " << ylim[1] << endl;
   }
}

void ScatHist::SetOtherSettings(bool *inConst)
{
   otherSettings = new bool[2];
   otherSettings[0] = inConst[0];
   otherSettings[1] = inConst[1];
}
// PlotSumResiduals
int ScatHist::StartPlotting()
{
   mystyle->SetBaseStyle();
   TCanvas *c1 = new TCanvas("c1","",1200,900);
//   TLatex chiText;

   // Comparison of observables on the same plot
   if(otherSettings[1])
   {
      cout << "TODO" << endl;
   }
   // Don't mix observables on the same plot
   else
   {
      // Merge trees to the same plot
      if(otherSettings[0])
      {
         for(int i = 0; i < filename->size(); i++)
	 {
            ifile = new TFile((filename->at(i)).c_str(), "READ");
	    itemp[0] = GetRootKeys(ifile);

	    delete ifile;
	 }
         cout << "TODO" << endl;
      }
      // Run plotting on single trees
      else
      {
         cout << "TODO" << endl;
      }
   }

   delete c1;
}

/*void MvaFitHist::SetPrimaries(vector<int> *primVals)
{
   includePart = new vector<int>;
   nrelem = new int;
   *nrelem = primVals->size();
   for(int i = 0; i < *nrelem; i++)
      includePart->push_back(primVals->at(i));
}

void MvaFitHist::SetMethod(string *inMethod)
{
   method = new string;
   mvalim = new double[2];
   mvalim[0] = GetMethodMin(*inMethod);
   mvalim[1] = GetMethodMax(*inMethod);
}

void MvaFitHist::SetHImodel(int *inModel)
{
   himodel = new int;
   *himodel = *inModel;
}

void MvaFitHist::SetSimProduction(int *inProd)
{
   simProd = new int;
   *simProd = *inProd;
}

void MvaFitHist::SetData(int *inTree)
{
   dataNum = new int;
   *dataNum = *inTree;
}

void MvaFitHist::SetNrBins(int *inNrBins)
{
   nrbins = new int;
   *nrbins = *inNrBins; 

   residVal = new double[*nrbins];
}

void MvaFitHist::SetXaxisLimits(bool use, double *inLimits)
{
   if(use)
   {
      xlim[0] = inLimits[0]; 
      xlim[1] = inLimits[1]; 

      cout << "Custom X-axis limits: " << xlim[0] << ", " << xlim[1] << endl;
   }
}

void MvaFitHist::SetYaxisLimits(bool use, double *inLimits)
{
   if(use)
   {
      ylim[0] = inLimits[0]; 
      ylim[1] = inLimits[1]; 

      cout << "Custom Y-axis limits: " << ylim[0] << ", " << ylim[1] << endl;
   }
}

void MvaFitHist::SetOtherSettings(bool *inConst)
{
   otherSettings = new bool[4];
   otherSettings[0] = inConst[0];
   otherSettings[1] = inConst[1];
   otherSettings[2] = inConst[2];
   otherSettings[3] = inConst[3];

   // Running an Xmax distribution fit
   if(otherSettings[1])
   {
      xlim[0] = 560.;
      xlim[1] = 1040.;
      cout << "Xmax limits: " << xlim[0] << ", " << xlim[1] << endl;
   }
   // Running a MVA distribution fit
   else
   {
      xlim[0] = GetMethodMin(stemp[0]);
      xlim[1] = GetMethodMax(stemp[0]);
      cout << "MVA limits: " << xlim[0] << ", " << xlim[1] << endl;
   }

   nrparam = new int;
   if(otherSettings[0])
      *nrparam = (*nrelem)-1;
   else
      *nrparam = (*nrelem);
}

void MvaFitHist::ResetHistograms()
{
   for(int i = 0; i < (*nrkeys); i++)
      delete result[i];
   delete ifile;
}

double MvaFitHist::fcnConstrainedFunc(const double *par)
{
   dtemp[0] = 1;

   // Sum all simulation histograms together (and scale them depending on the fraction)
   for(int j = 1; j < (*nrbins)+1; j++)
   {
      sim->SetBinContent(j, 0);
      sim->SetBinError(j, 0);
   }

   for(int i = 0; i < (*nrelem); i++)
   {
      sim->Add(simhist[i], par[i]);
      dtemp[0] -= par[i];
   }

   // Go over the simulation distribution and check all the bins -> keep only bins from sim and data that are not 0
   dtemp[1] = 0;
   for(int j = 1; j < (*nrbins)+1; j++)
   {
      dtemp[1] += sim->GetBinContent(j);

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
   }

   // Return constraint: [1 - sum(f_i)]
   return dtemp[0];
}

// step#013
double MvaFitHist::fcnFunc(const double *par)
{
   dtemp[0] = 1;
   // Scale all simulation histograms with elemental fractions
   for(int i = 0; i < (*nrelem); i++)
   {
      stemp[0] = "sims" + ToString(i);
      sims[i] = (TH1F*)simhist[i]->Clone(stemp[0].c_str());

      if(i == (*nrelem)-1)
         sims[i]->Scale(dtemp[0]);
      else
      {
         dtemp[0] -= par[i];
         sims[i]->Scale(par[i]);
      }
   }

   // Sum all simulation histograms together
   sim = (TH1F*)sims[0]->Clone("sim");
   for(int i = 1; i < (*nrelem); i++)
      sim->Add(sims[i]);
   
   // Go over the simulation distribution and check all the bins
   dtemp[1] = 0;
   for(int j = 1; j < (*nrbins)+1; j++)
   {
      dtemp[1] += sim->GetBinContent(j);

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
   }

   // No returned constraint, because it is already included with the last parameter
   return 0;
}

// step#013
double MvaFitHist::fcn(const double *par)
{
   if(otherSettings[0])
      dtemp[0] = fcnConstrainedFunc(par);
   else
      dtemp[0] = fcnFunc(par);

   return sim->Chi2Test(data, "WW CHI2") + TMath::Power(dtemp[0],2);
}

void MvaFitHist::SetInputFile(string *inFile)
{
   *filename = *inFile;
}

void MvaFitHist::PrepareHistograms(int run)
{
   mystyle->SetBaseStyle();

   cout << endl << "Opening file: " << *filename << endl;
   ifile = new TFile(filename->c_str(), "READ");

   // Check how many TTrees we have in the file
   *nrkeys = GetRootKeys(ifile);

   // Prepare variable histograms to hold values from each tree
// step#003
   for(int i = 0; i < *nrkeys; i++)
   {
      stemp[0] = "result" + ToString(i+1);
      result[i] = new TH1F(stemp[0].c_str(), "", *nrbins, xlim[0], xlim[1]);
   }
   TTree *tempTree;

   float *mva = new float;

   // Go through all TTrees in the root file
   for(int i = 0; i < (*nrkeys); i++)
   {
      stemp[0] = "TreeS" + ToString(i+1);
      cout << "Getting tree = " << stemp[0] << ", ";
      tempTree = (TTree*)ifile->Get(stemp[0].c_str());
      cout << "with title = " << tempTree->GetTitle() << endl;

      // Use xmax instead of MVA
      if(otherSettings[1])
         tempTree->SetBranchAddress("xmax", mva);
      // Use MVA variable
      else
         tempTree->SetBranchAddress("MVA", mva);
      
      nentries[i] = tempTree->GetEntries();
      cout << "  Number of entries = " << nentries[i] << endl;

// step#004
      // Go through all events in the root file for the specific TTree and save MVA variable values to the histogram
      for(int j = 0; j < nentries[i]; j++)
      {
         tempTree->GetEntry(j);
         result[i]->Fill(*mva);
      }

// step#005
      // Calculate errors
      if(!otherSettings[3]) // only when not using FractionFitter
      {
         vector<double> *binnentries = new vector<double>;

         cout << "result" << i << ":" << endl;
         for(int k = 1; k <= (*nrbins); k++)
         {
            dtemp[0] = result[i]->GetBinContent(k);
            cout << k << "\t" << dtemp[0] << "\t" << result[i]->GetBinError(k) << endl;
	    // If bin value is N = 0, save value for error as 0,
	    //   otherwise, use a Poissonian error -> sqrt(N)/N
            if(dtemp[0] == 0)
               binnentries->push_back(0.);
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
         // Reevaluate the bin errors on scaled histograms
         for(int k = 1; k <= (*nrbins); k++)
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

	 // Modify the maximal value with the error
         cout << "Maximum value " << norm[i] << " increased by " << dtemp[1] << endl;
         norm[i] += dtemp[1];

         delete binnentries;
      }
      else
      {
         // Get the maximal value
         norm[i] = result[i]->GetMaximum();
         cout << "  Highest value = " << norm[i] << endl;

	 // Modify the maximal value with the error
         dtemp[1] = 0;
         for(int k = 1; k <= (*nrbins); k++)
	 {
            dtemp[0] = result[i]->GetBinContent(k);

            if(result[i]->GetBinError(k) > dtemp[1])
               dtemp[1] = result[i]->GetBinError(k);
	 }

         cout << "Maximum value " << norm[i] << " increased by " << dtemp[1] << endl;
         norm[i] += dtemp[1];
      }
   }

   delete mva;

   // Read the application_result.txt file in order to get information on energy bins
   ResultRead *analRes = new ResultRead();
   stemp[0] = RemoveFilename(filename) + "/application_results.txt";
   itemp[0] = analRes->ReadFile(stemp[0]);

   if(itemp[0] != 1)
   {
      // Check the application_result.txt file (only once)
      if(run == 0)
      {
         trees->clear();
         treeNames->clear();

	 itemp[1] = 0;
         for(int i = 0; i < (*nrkeys); i++)
         {
            // Saving simulation trees
            for(int j = 0; j < (*nrelem); j++)
            {
               stemp[1] = analRes->GetTreeName(i);
               stemp[2] = ptype->GetName(includePart->at(j));

	       // Make sure to only select pure samples for simulations, not mixed mock datasets
	       if(stemp[1].compare(stemp[2]) == 0)
	       {
                  trees->push_back(i);
                  treeNames->push_back(stemp[2]);
	          itemp[1]++;
	       }
            }

	    if(itemp[1] >= (*nrelem))
               break;
         }

         // Saving the data tree
         trees->push_back(*dataNum);
         treeNames->push_back(analRes->GetTreeName(*dataNum));
      }

      cout << "The selected simulation trees are: ";
      for(int i = 0; i < (*nrelem); i++)
      {
         if(i > 0)
            cout << ", ";
         cout << trees->at(i) << " (" << treeNames->at(i) << ")";
      }
      cout << endl << "The selected data tree is: " << trees->at(*nrelem) << " (" << treeNames->at(*nrelem) << ")" << endl << endl;

      // Getting the signal fraction approximation for the data
      *mvacut = analRes->GetMvaCut(0);
      *midEnergy = analRes->GetEnergy();

      cout << "Mva cut value: " << *mvacut << endl;
      cout << "Energy value: " << *midEnergy << endl;
   }

   // Save all simulation histograms (save them from result histograms to simhist histograms)
   for(int i = 0; i < (*nrelem); i++)
   {
      stemp[0] = "simhist" + ToString(i);
      simhist[i] = (TH1F*)result[trees->at(i)]->Clone(stemp[0].c_str());
      nrsim[i] = nentries[trees->at(i)];

      // Go over the simulation distribution and check all the bins
      cout << "Simhist " << i << " (" << treeNames->at(i) << "):" << endl;
      dtemp[0] = 0;
      for(int j = 1; j < (*nrbins)+1; j++)
      {
	 dtemp[0] += simhist[i]->GetBinContent(j);
         cout << j << "\t" << simhist[i]->GetBinContent(j) << "\t" << simhist[i]->GetBinError(j) << endl;
      }
      cout << "Histogram integral = " << dtemp[0] << endl;
   }

// step#009
   // Check if all simulation histograms have zeroes at the same place (if not, remove the large error)
   if(!otherSettings[3])
   {
      iszero = new bool;
      *iszero = false;
      dtemp[0] = 0;
      for(int j = 1; j < (*nrbins)+1; j++)
      {
         *iszero = false;
         dtemp[0] = 0;

         for(int i = 0; i < (*nrelem); i++)
         {
            dtemp[0] += simhist[i]->GetBinContent(j);
            if(simhist[i]->GetBinContent(j) == 0)
               (*iszero) = true;
         }

         if( (dtemp[0] != 0) && (*iszero) )
         {
            for(int i = 0; i < (*nrelem); i++)
            {
               if(simhist[i]->GetBinContent(j) == 0)
                  simhist[i]->SetBinError(j, 0);
            }
         }
      }

      delete iszero;

      for(int i = 0; i < (*nrelem); i++)
      {
         // Go over the simulation distribution and check all the bins
         cout << "Simhist " << i << " (edit):" << endl;
         dtemp[0] = 0;
         for(int j = 1; j < (*nrbins)+1; j++)
         {
            dtemp[0] += simhist[i]->GetBinContent(j);
            cout << j << "\t" << simhist[i]->GetBinContent(j) << "\t" << simhist[i]->GetBinError(j) << endl;
         }
         cout << "Histogram integral = " << dtemp[0] << endl;
      }
   }

   cout << "Number of parameters for fitting: " << *nrparam << endl;

   // Save data distributions
   data = (TH1F*)result[trees->at(*nrelem)]->Clone("data");
   yrange[0] = norm[*nrelem];

   // Go over the data distribution and check all the bins
   cout << "Data:" << endl;
   for(int j = 1; j < (*nrbins)+1; j++)
      cout << j << "\t" << data->GetBinContent(j) << "\t" << data->GetBinError(j) << endl;

   // Create directory structure for plots and delete old plots
   stemp[0] = "mkdir -p " + RemoveFilename(filename) + "/fithist";
   system(stemp[0].c_str());
   stemp[0] = "rm -fr " + RemoveFilename(filename) + "/fithist/signal_background_data.pdf " + RemoveFilename(filename) + "/fithist/composition_residuals*";
   system(stemp[0].c_str());

// step#010
   // Zero empty bin errors just for plotting purposes
   newdata = (TH1F*)data->Clone("newdata");

   for(int j = 1; j < (*nrbins)+1; j++)
   {
      if(data->GetBinContent(j) == 0)
         newdata->SetBinError(j, 0);
   }

   delete analRes;
}

// Plot data, normalized signal+background and normalized residuals
int MvaFitHist::StartFitting()
{
   // Find the best fit with TMinuit (new interface)
   if(!otherSettings[3])
   {
*//*      fitproc = 1;

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
	       cout << "Parameter " << iFitError[0] << " has smallest error, so fixing it." << endl;
               fitter->Config().ParSettings(iFitError[0]).Fix();
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
	       cout << "Parameter " << iFitError[0] << " has smallest error, so fixing it." << endl;
               fitter->Config().ParSettings(iFitError[0]).Fix();
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
            bestfrac[i] = fcnATan(&bestfrac[i]);
            bestfracerr[2*i] = fcnATan(&bestfracerr[2*i]);
            bestfracerr[2*i+1] = fcnATan(&bestfracerr[2*i+1]);
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
      delete result;*//*
   }
   // Find the best fit with FracFitter
   else
   {
      // Prepare variables to hold best parameter fits and their errors
      double *bestfrac = new double[*nrelem];
      double *bestfracerr = new double[2*(*nrelem)];

      // Prepare the simulation histograms
      mc =  new TObjArray(*nrelem);
      for(int i = 0; i < *nrelem; i++)
         mc->Add(simhist[i]);

      // Prepare the FractionFitter class and set constraints
      fracFitter = new TFractionFitter(data, mc);

      // Constrain fitting parameters between 0 and 1
      for(int i = 0; i < *nrelem; i++)
         fracFitter->Constrain(i, 0., 1.);

      fracVirtFitter = (TVirtualFitter*)fracFitter->GetFitter();

      // Remove zero bins
      if(otherSettings[2])
      {
         // Exclude any bins that are 0 in simulations or data
         for(int i = 1; i < (*nrbins)+1; i++)
         {
            dtemp[0] = 0;
            for(int j = 0; j < (*nrelem); j++)
               dtemp[0] += simhist[j]->GetBinContent(i);

            if( (data->GetBinContent(i) == 0) || (dtemp[0] == 0) )
            {
               cout << "Excluding bin " << i-1 << " from further analysis." << endl;
               fracFitter->ExcludeBin(i);
            }
         }
      }

      // Run fitting
      itemp[0] = fracFitter->Fit();
      fracFitter->ErrorAnalysis(1);
      cout << "Fitting status = " << itemp[0] << endl;

      // Retry to fit with initial values from a previous fit (if the previous returns an error)
      if(itemp[0] != 0)
      {
	 cout << "Retrying a fit" << endl;
         itemp[0] = fracFitter->Fit();
         fracFitter->ErrorAnalysis(1);
         cout << "Fitting status = " << itemp[0] << endl;
      }

      // If fitting status is okay, continue with saving values and plotting
      if(itemp[0] == 0)
      {
	 for(int i = 0; i < (*nrelem); i++)
	 {
            bestfrac[i] = fracVirtFitter->GetParameter(i);
	    fracVirtFitter->GetErrors(i, bestfracerr[2*i], bestfracerr[2*i+1], dtemp[1], dtemp[0]);
            dtemp[0] = fracVirtFitter->GetParError(i)/dtemp[1];

	    // bestfrac[2*i]   -> positive error
	    // bestfrac[2*i+1] -> negative error
	    // dtemp[1]        -> parabolic error

	    // Getting fitting uncertainties (separate values from error analysis, not the parabolic error)
	    cout << i << ": Parabolic error ratio: " << fracVirtFitter->GetParError(i) << "/" << dtemp[1] << " = " << dtemp[0] << "\t" << "\t";

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
	       else
                  bestfracerr[2*i] = 0;

	       bestfracerr[2*i+1] = TMath::Abs(bestfracerr[2*i+1]);
	    }
	    // Negative error is missing
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
	       else
                  bestfracerr[2*i+1] = 0;

               bestfracerr[2*i] = TMath::Abs(bestfracerr[2*i]);
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
         for(int i = 0; i < (*nrparam); i++)
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
            midLna[0] += midComposition->at(i)*TMath::Log(ptype->GetA(itemp[1]));
            midLna[1] += (midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[1])))*(midCompositionErr->at(2*i)*TMath::Log(ptype->GetA(itemp[1])));
            midLna[2] += (midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[1])))*(midCompositionErr->at(2*i+1)*TMath::Log(ptype->GetA(itemp[1])));
         }

         midLna[1] = TMath::Sqrt(midLna[1]);
         midLna[2] = TMath::Sqrt(midLna[2]);

         cout << "lnA value = " << midLna[0] << " +/- " << midLna[1] << "/" << midLna[2] << endl;

         // Get the parameter step value (not applicable in this case)
         *midStep = -1.;

         PlotSumResiduals(bestfrac, bestfracerr, itemp[0], -1.);

         delete mc;
         delete fracFitter;

         delete[] bestfrac;
         delete[] bestfracerr;

	 return 0;
      }
      else
      {
         cout << "Error! Minimization failed." << endl;

	 // Although minimization failed, still return the parameter values
	 for(int i = 0; i < (*nrelem); i++)
	 {
            bestfrac[i] = fracVirtFitter->GetParameter(i);
	    fracVirtFitter->GetErrors(i, bestfracerr[2*i], bestfracerr[2*i+1], dtemp[1], dtemp[0]);
            cout << "frac" << i << " (" << ptype->GetName(includePart->at(i)) << ")\t = " << ToString(bestfrac[i], 5) << "\t+/-\t" << ToString(bestfracerr[2*i], 5) << "/" << ToString(bestfracerr[2*i+1], 5) << endl;
	 }

         delete mc;
         delete fracFitter;

         delete[] bestfrac;
         delete[] bestfracerr;

	 return -1;
      }
   }
}

// Plot data, normalized signal+background and normalized residuals
void MvaFitHist::PlotSumResiduals(double *sigFrac, double *sigFracErr, int status, double selStep)
{
   mystyle->SetBaseStyle();
   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TLatex chiText;
   c1->cd();
   TPad *upperpad = new TPad("upperpad", "upperpad", 0.004, 0.490, 0.996, 0.996);
   TPad *lowerpad = new TPad("lowerpad", "lowerpad", 0.004, 0.000, 0.996, 0.490);
   upperpad->Draw();
   lowerpad->Draw();

   // Plot the data distribution
   upperpad->cd();

   if(!otherSettings[3])
   {
      mystyle->SetHistColor((TH1*)newdata, 2);
      if(otherSettings[1])
         mystyle->SetAxisTitles((TH1*)newdata, "X_{max} (g/cm^{2})", "Probability density");
      else
         mystyle->SetAxisTitles((TH1*)newdata, "MVA variable", "Probability density");
      newdata->GetXaxis()->SetTitleOffset(2.0);
      newdata->Draw();

      if(!otherSettings[0])
         cout << "Chi2+constraint^2 = " << fcn(sigFrac) << endl;
      else
         cout << "Constraint = " << fcnFunc(sigFrac) << endl;

      // Zero empty bin errors just for plotting purposes
      newsig = (TH1F*)sim->Clone("newsig");
      for(int j = 1; j < (*nrbins)+1; j++)
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

      *pvalue = sim->Chi2TestX(data, *chi2value, *ndfvalue, *chiGood, "WW", residVal);
      cout << "Chi2 = " << *chi2value << ", p-value = " << *pvalue << ", NDF = " << *ndfvalue << ", igood = " << *chiGood << endl;
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
      if(otherSettings[1])
         mystyle->SetAxisTitles((TH1*)sim, "X_{max} (g/cm^{2})", "Number of events");
      else
         mystyle->SetAxisTitles((TH1*)sim, "MVA variable", "Number of events");
      sim->Draw();
      
      mystyle->SetHistColor((TH1*)data, 2);
      if(otherSettings[1])
         mystyle->SetAxisTitles((TH1*)data, "X_{max} (g/cm^{2})", "Number of events");
      else
         mystyle->SetAxisTitles((TH1*)data, "MVA variable", "Number of events");
      sim->SetMaximum(1.1*TMath::MaxElement(2, yrange));
      data->Draw("Ep;same");

      for(int i = 1; i <= (*nrbins); i++)
      {
         if(sim->GetBinContent(i) == 0)
            residVal[i-1] = 0;
	 else
	 {
            residVal[i-1] = (data->GetBinContent(i) - sim->GetBinContent(i))/TMath::Sqrt(sim->GetBinContent(i));
	 }
      }
   }

   cout << "Getting residuals" << endl;
   TH1F *resid = new TH1F("resid", "", *nrbins, xlim[0], xlim[1]);
   yresidrange[0] = 1.e+24;
   yresidrange[1] = -1.e+24;
   vector<double> *residAll = new vector<double>;
   for(int i = 1; i <= (*nrbins); i++)
   {
      // Only set residuals, if any of the two histogram bins are not zero
      if( ((double)sim->GetBinContent(i) != 0) && ((double)data->GetBinContent(i) != 0) )
      {
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

   c1->cd();
   upperpad->cd();

   for(int i = 0; i < residAll->size(); i++)
      tempd[i] = residAll->at(i);

   // Use xmax instead of MVA:
   if(otherSettings[1])
      chiText.SetTextAlign(31);
   else
      chiText.SetTextAlign(21);

   // Print information regarding the fit on the plot
   if(!otherSettings[3])
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString((*chi2value), 4) + "/" + ToString((*ndfvalue)) + " = " + ToString((*chi2value)/(double)(*ndfvalue), 4) + ", status = " + ToString(status) + ", step = " + ToString(selStep, 3) + ", p-value = " + ToString((*pvalue), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      dtemp[0] = 1.;
      dtemp[1] = 0.;
      dtemp[2] = 0.;
      for(int i = 0; i < (*nrparam); i++)
      {
         if(!otherSettings[1])
	 {
            stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames->at(i) + ")";
            chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
	 }
	 else
	 {
            if(i < (*nrparam)-1)
            {
               dtemp[0] -= sigFrac[i];
               dtemp[1] += sigFracErr[2*i];
               dtemp[2] += sigFracErr[2*i+1];
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames->at(i) + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
            else
            {
               stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(dtemp[0], 4) + " +" + ToString(dtemp[1], 4) + "/-" + ToString(dtemp[2], 4) + " (" + treeNames->at(i) + ")";
               chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
            }
	 }
      }

      stemp[0] = "Residuals: mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4) + ", RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)(*nrparam)+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      stemp[0] = "MVA cut = " + ToString(*mvacut, 4);
      chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)(*nrparam)+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }
   else
   {
      stemp[0] = "#chi^{2}/NDF = " + ToString((*chi2value), 4) + "/" + ToString((*ndfvalue)) + " = " + ToString((*chi2value)/(double)(*ndfvalue), 4) + ", status = " + ToString(status) + ", p-value = " + ToString((*pvalue), 4);
      // Use xmax instead of MVA:
      if(otherSettings[1])
         chiText.DrawLatex(xlim[1]*0.98, TMath::MaxElement(2, yrange), stemp[0].c_str());
      else
         chiText.DrawLatex((xlim[0]+xlim[1])/2., TMath::MaxElement(2, yrange), stemp[0].c_str());

      dtemp[0] = 1.;
      dtemp[1] = 0.;
      dtemp[2] = 0.;
      for(int i = 0; i < (*nrparam); i++)
      {
         stemp[0] = "r_{" + ToString(i+1) + "} = " + ToString(sigFrac[i], 4) + " +" + ToString(sigFracErr[2*i], 4) + "/-" + ToString(sigFracErr[2*i+1], 4) + " (" + treeNames->at(i) + ")";
         // Use xmax instead of MVA:
	 if(otherSettings[1])
            chiText.DrawLatex(xlim[1]*0.98, (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
	 // Use MVA
	 else
            chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)i+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      }

      stemp[0] = "Residuals: mean = " + ToString(TMath::Mean(residAll->size(), tempd), 4) + ", RMS = " + ToString(TMath::RMS(residAll->size(), tempd), 4);
      // Use xmax instead of MVA:
      if(otherSettings[1])
         chiText.DrawLatex(xlim[1]*0.98, (1.-0.06*((double)(*nrparam)+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      // Use MVA
      else
         chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)(*nrparam)+1.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      stemp[0] = "MVA cut = " + ToString(*mvacut, 4);
      // Use xmax instead of MVA:
      if(otherSettings[1])
         chiText.DrawLatex(xlim[1]*0.98, (1.-0.06*((double)(*nrparam)+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
      // Use MVA
      else
         chiText.DrawLatex((xlim[0]+xlim[1])/2., (1.-0.06*((double)(*nrparam)+2.))*TMath::MaxElement(2, yrange), stemp[0].c_str());
   }

   delete residAll;
   delete tempd;

   lowerpad->cd();
   mystyle->SetHistColor((TH1*)resid, 3);
   if(!otherSettings[3])
     mystyle->SetAxisTitles((TH1*)resid, "", "Adjusted residuals");
   else
     mystyle->SetAxisTitles((TH1*)resid, "", "Standardized residuals");
   resid->GetYaxis()->SetRangeUser(-1.1*TMath::MaxElement(2, yresidrange), 1.1*TMath::MaxElement(2, yresidrange));
   resid->Draw();

   TLine *line = new TLine(*mvacut, -1.1*TMath::MaxElement(2, yresidrange), *mvacut, 1.1*TMath::MaxElement(2, yresidrange));
   line->SetLineWidth(2);
   line->SetLineStyle(1);
   line->SetLineColor(kOrange+2);
   line->Draw("SAME");

   c1->Update();

   stemp[1] = "";
   for(int i = 0; i < (*nrelem); i++)
   {
      if(i == 0)
         stemp[1] += ptype->GetShortName(includePart->at(i));
      else
         stemp[1] += "-" + ptype->GetShortName(includePart->at(i));
   }

   if(!otherSettings[3])
   {
      stemp[0] = RemoveFilename(filename) + "/fithist/" + stemp[1] + "_composition_residuals_minuit2.pdf";
      c1->SaveAs(stemp[0].c_str());
   }
   else
   {
      if(otherSettings[1])
         stemp[0] = RemoveFilename(filename) + "/fithist/" + stemp[1] + "_composition_residuals_fracfitter-xmax.pdf";
      else
         stemp[0] = RemoveFilename(filename) + "/fithist/" + stemp[1] + "_composition_residuals_fracfitter.pdf";
      c1->SaveAs(stemp[0].c_str());
   }

   delete c1;
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

void MvaFitHist::PrintoutFinalResults(int *nrp, vector<double> *lnas, vector<double> *chi2, vector<double> *pvalues, vector<double> *ndfs, vector<double> *comps, vector<double> *energies, vector<double> *steps)
{
   compnrp = new int;
   *compnrp = *nrp;
   xval = new double[*compnrp];
   xerr = new double[*compnrp];
   yval = new double[*compnrp];
   yerrlo = new double[*compnrp];
   yerrhi = new double[*compnrp];
   composition = new vector<double>;

   cout << endl << "Final results:" << endl;
   cout << "# nr. elements = " << (*nrelem) << ", nr. bins = " << (*nrbins) << endl;
   cout << "# E   \tChi2/NDF \tp-val\tstep\tlnA \t+   \t-   \t";
   for(int i = 0; i < (*nrelem); i++)
      cout << "\t" << ptype->GetShortName(includePart->at(i)) << "  \t+   \t-   ";
   cout << endl;
   cout.precision(4);
   for(int i = 0; i < *compnrp; i++)
   {
      itemp[0] = 3*(*nrelem);
      cout << fixed << energies->at(i) << "\t" << chi2->at(i) << "/" << (int)ndfs->at(i) << "\t" << pvalues->at(i) << "\t" << steps->at(i) << "\t" << lnas->at(3*i) << "\t" << lnas->at(3*i+1) << "\t" << lnas->at(3*i+2) << "\t";

      xval[i] = energies->at(i);
      xerr[i] = 0.;
      yval[i] = lnas->at(3*i);
      yerrhi[i] = lnas->at(3*i+1);
      yerrlo[i] = lnas->at(3*i+2);

      for(int j = 0; j < itemp[0]; j++)
         cout << "\t" << comps->at(itemp[0]*i+j);

      cout << endl;
   }

   for(int i = 0; i < comps->size(); i++)
      composition->push_back(comps->at(i));

   cout << endl;
}

// Plot Lna and composition plots
void MvaFitHist::PlotLnaComposition()
{
   mystyle->SetBaseStyle();
   TCanvas *c1 = new TCanvas("c1","",1200,900);

   // Prepare lnA plot, with appended published values
   TGraphAsymmErrors *grLna;
   TGraphAsymmErrors *grPubLnaFD, *grPubLnaSD, *grPubLnaFDFrac;
   TLegend *legend;
   int c_Legend = TColor::GetColor("#ffff66");

   grLna = new TGraphAsymmErrors(*compnrp, xval, yval, xerr, xerr, yerrlo, yerrhi);
   mystyle->SetGraphColor(grLna, 2);
   grLna->GetXaxis()->SetRange(18.4, 20.0);
   grLna->GetXaxis()->SetRangeUser(18.4, 20.0);
   grLna->GetXaxis()->SetLimits(18.4, 20.0);
   grLna->GetYaxis()->SetRangeUser(-0.7, 5.);

   cout << "Used dataset type: " << *himodel << endl;

   bool *addPublishedFD = new bool;
   bool *addPublishedSD = new bool;
   bool *addPublishedFrac = new bool;
   bool *pubPRD = new bool;

   if(*simProd == 0)
      *pubPRD = true;
   else
      *pubPRD = false;

   float *xbinPubFD[3], *xbinPubSD[3], *xbinPubFDFrac[3];
   float *ybinPubLnaFD[3], *ybinPubLnaSD[3], *ybinPubLnaFDFrac[3];

   vector<float> *returnFrac = new vector<float>;

   // Read the published results, if we have not selected None in the high-energy hadron interaction model
   if(*himodel != 3)
   {
      *addPublishedFD = true;
      vector<float> *returnVal = new vector<float>;
      cout << "Reading FD published Xmax results" << endl;
      itemp[0] = ReadLnaResultsFD(returnVal);

      if(itemp[0] != -1)
      {
         cout << "Number of FD points = " << itemp[0] << endl;
         
         for(int i = 0; i < 3; i++)
         {
            xbinPubFD[i] = new float[itemp[0]];
            ybinPubLnaFD[i] = new float[itemp[0]];
         }
         
         for(int i = 0; i < itemp[0]; i++)
         {
            xbinPubFD[0][i] = returnVal->at(3*i);
            xbinPubFD[1][i] = 0;
            xbinPubFD[2][i] = 0;
            ybinPubLnaFD[0][i] = returnVal->at(3*i+1);
            ybinPubLnaFD[1][i] = returnVal->at(3*i+2);
            ybinPubLnaFD[2][i] = returnVal->at(3*i+2);
         }
         
         grPubLnaFD = new TGraphAsymmErrors(itemp[0], xbinPubFD[0], ybinPubLnaFD[0], xbinPubFD[1], xbinPubFD[2], ybinPubLnaFD[1], ybinPubLnaFD[2]);
         mystyle->SetGraphColor(grPubLnaFD, 0);
         grPubLnaFD->GetXaxis()->SetRange(18.5, 20.0);
         grPubLnaFD->GetXaxis()->SetRangeUser(18.5, 20.0);
         grPubLnaFD->GetXaxis()->SetLimits(18.5, 20.0);
         grPubLnaFD->GetYaxis()->SetRangeUser(-0.7, 5.);

         if(*himodel != 2)
         {
            *addPublishedSD = true;
            returnVal->clear();
            cout << "Reading SD published results" << endl;
            itemp[0] = ReadLnaResultsSD(returnVal, 1);
            if(itemp[0] != -1)
            {
               cout << "Number of SD points = " << itemp[0] << endl;
               
               for(int i = 0; i < 3; i++)
               {
                  xbinPubSD[i] = new float[itemp[0]];
                  ybinPubLnaSD[i] = new float[itemp[0]];
               }
               
               for(int i = 0; i < itemp[0]; i++)
               {
                  xbinPubSD[0][i] = returnVal->at(3*i);
                  xbinPubSD[1][i] = 0;
                  xbinPubSD[2][i] = 0;
                  ybinPubLnaSD[0][i] = returnVal->at(3*i+1);
                  ybinPubLnaSD[1][i] = returnVal->at(3*i+2);
                  ybinPubLnaSD[2][i] = returnVal->at(3*i+2);
               }
               
               grPubLnaSD = new TGraphAsymmErrors(itemp[0], xbinPubSD[0], ybinPubLnaSD[0], xbinPubSD[1], xbinPubSD[2], ybinPubLnaSD[1], ybinPubLnaSD[2]);
               mystyle->SetGraphColor(grPubLnaSD, 1);
               grPubLnaSD->GetXaxis()->SetRange(18.5, 20.0);
               grPubLnaSD->GetXaxis()->SetRangeUser(18.5, 20.0);
               grPubLnaSD->GetXaxis()->SetLimits(18.5, 20.0);
               grPubLnaSD->GetYaxis()->SetRangeUser(-0.7, 5.);
            }
            else
               *addPublishedSD = false;
         }
         else
            *addPublishedSD = false;
      }
      else
         *addPublishedFD = false;

      *addPublishedFrac = true;
      cout << "Reading FD lnA from published fraction results" << endl;
      itemp[0] = ReadFracResultsFD(returnFrac, true);
      if(itemp[0] != -1)
      {
         cout << "Number of FD points = " << itemp[0] << endl;
         
         for(int i = 0; i < 3; i++)
         {
            xbinPubFDFrac[i] = new float[itemp[0]];
            ybinPubLnaFDFrac[i] = new float[itemp[0]];
         }
         
         for(int i = 0; i < itemp[0]; i++)
         {
            xbinPubFDFrac[0][i] = returnFrac->at(4*i);
            xbinPubFDFrac[1][i] = 0;
            xbinPubFDFrac[2][i] = 0;
            ybinPubLnaFDFrac[0][i] = returnFrac->at(4*i+1);
            ybinPubLnaFDFrac[1][i] = returnFrac->at(4*i+2);
            ybinPubLnaFDFrac[2][i] = returnFrac->at(4*i+3);
         }
         
         grPubLnaFDFrac = new TGraphAsymmErrors(itemp[0], xbinPubFDFrac[0], ybinPubLnaFDFrac[0], xbinPubFDFrac[1], xbinPubFDFrac[2], ybinPubLnaFDFrac[1], ybinPubLnaFDFrac[2]);
         mystyle->SetGraphColor(grPubLnaFDFrac, 0);
         grPubLnaFDFrac->GetXaxis()->SetRange(18.5, 20.0);
         grPubLnaFDFrac->GetXaxis()->SetRangeUser(18.5, 20.0);
         grPubLnaFDFrac->GetXaxis()->SetLimits(18.5, 20.0);
         grPubLnaFDFrac->GetYaxis()->SetRangeUser(-0.7, 5.);
      }
      else
         *addPublishedFrac = false;

      delete returnVal;
   }

   // Plot lnA graph, where we add Xmax and Delta analysis (both from ICRC 2017)
   mystyle->SetAxisTitles(grLna, "FD energy [log(E/eV)]", "<lnA> of data events");
   grLna->Draw("APL");
   if(*addPublishedFD)
      grPubLnaFD->Draw("PL;SAME");
   if(*addPublishedSD)
      grPubLnaSD->Draw("PL;SAME");

   if(*addPublishedFD && *addPublishedSD)
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(3)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
   else if(*addPublishedFD || *addPublishedSD)
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(2)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
   else
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(1)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_Legend);
   if(*himodel == 3)
      legend->AddEntry(grLna, "Mock data (Histogram fit)", "lp");
   else
      legend->AddEntry(grLna, "Data (Histogram fit)", "lp");
   if(*addPublishedFD && !(*pubPRD))
      legend->AddEntry(grPubLnaFD, "Data (Xmax analysis ICRC2017)", "lp");
   if(*addPublishedFD && *pubPRD)
      legend->AddEntry(grPubLnaFD, "Data (Xmax analysis PRD2014)", "lp");
   if(*addPublishedSD)
      legend->AddEntry(grPubLnaSD, "Data (Delta analysis ICRC2017)", "lp");
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

   stemp[1] = RemoveFilename(filename);
   stemp[2] = "mkdir -p " + RemoveFilename(&stemp[1]) + "/plots";
   system(stemp[2].c_str());

   if(otherSettings[1])
   {
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/histogram_fit_lnA_xmax.pdf";
      c1->SaveAs(stemp[2].c_str());
   }
   else
   {
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/histogram_fit_lnA.pdf";
      c1->SaveAs(stemp[2].c_str());
   }

   // Plot lnA graph, where we add Xmax analysis from histogram fitting approach (longitudinal distribution fit)
   mystyle->SetAxisTitles(grLna, "FD energy [log(E/eV)]", "<lnA> of data events");
   grLna->Draw("APL");
   if(*addPublishedFrac)
      grPubLnaFDFrac->Draw("PL;SAME");
   if(*addPublishedSD)
      grPubLnaSD->Draw("PL;SAME");

   if(*addPublishedFrac && *addPublishedSD)
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(3)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
   else if(*addPublishedFrac || *addPublishedSD)
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(2)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
   else
      legend = new TLegend(gPad->GetLeftMargin(), 1-gPad->GetTopMargin()-(mystyle->SetLegendHeight(1)), gPad->GetLeftMargin()+.32, 1-gPad->GetTopMargin());
   legend->SetFillStyle(1001);
   legend->SetFillColor(c_Legend);
   if(*himodel == 3)
      legend->AddEntry(grLna, "Mock data (Histogram fit)", "lp");
   else
      legend->AddEntry(grLna, "Data (Histogram fit)", "lp");
   if(*addPublishedFrac)
      legend->AddEntry(grPubLnaFDFrac, "Data (Longitudinal fit analysis ICRC2017)", "lp");
   if(*addPublishedSD)
      legend->AddEntry(grPubLnaSD, "Data (Delta analysis ICRC2017)", "lp");
   legend->SetBorderSize(1);
   legend->SetMargin(0.3);
   legend->Draw("same");

   t = new TText();
   t->SetTextAlign(12);
   t->SetTextColor(28);
   t->SetTextSize(24);
   l = new TLine();
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

   stemp[1] = RemoveFilename(filename);
   stemp[2] = "mkdir -p " + RemoveFilename(&stemp[1]) + "/plots";
   system(stemp[2].c_str());

   if(otherSettings[1])
   {
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/histogram_fit_lnA-frac_xmax.pdf";
      c1->SaveAs(stemp[2].c_str());
   }
   else
   {
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/histogram_fit_lnA-frac.pdf";
      c1->SaveAs(stemp[2].c_str());
   }
     
   if(*addPublishedFD)
   {
      for(int i = 0; i < 3; i++)
      {
         delete[] xbinPubFD[i];
         delete[] ybinPubLnaFD[i];
      }
   }
   
   if(addPublishedFrac)
   {
      for(int i = 0; i < 3; i++)
      {
         delete[] xbinPubFDFrac[i];
         delete[] ybinPubLnaFDFrac[i];
      }
   }
   
   if(addPublishedSD)
   {
      for(int i = 0; i < 3; i++)
      {
         delete[] xbinPubSD[i];
         delete[] ybinPubLnaSD[i];
      }
   }

   delete addPublishedFD;
   delete addPublishedSD;
   delete pubPRD;

   delete c1;

   cout << endl;

   double *ybinComp[40];
   double *ybinCompErrlo[40];
   double *ybinCompErrhi[40];

   double *xbinCompPub;
   double *xbinCompPubErr;
   double *ybinCompPub[40];
   double *ybinCompPubErrlo[40];
   double *ybinCompPubErrhi[40];

   // Plot composition plots (all together)
   c1 = new TCanvas("c1","",1200,400.*(*nrelem));

   TGraphAsymmErrors *grComp[40], *grCompPub[40];

   // Get published composition results
   if(*himodel != 3)
   {
      *addPublishedFrac = true;
      returnFrac->clear();
      cout << "Reading FD composition from published fraction results" << endl;
      itemp[2] = ReadFracResultsFD(returnFrac, false);
      if(itemp[2] != -1)
      {
         cout << "Number of FD points = " << itemp[2] << endl;

         xbinCompPub = new double[itemp[2]];
         xbinCompPubErr = new double[itemp[2]];
         for(int j = 0; j < (*nrelem); j++)
         {
            ybinCompPub[j] = new double[itemp[2]];
            ybinCompPubErrlo[j] = new double[itemp[2]];
            ybinCompPubErrhi[j] = new double[itemp[2]];
         }

         cout << "Points:" << endl;
         for(int i = 0; i < itemp[2]; i++)
         {
            itemp[0] = (3*(*nrelem)+1)*i;
            cout << i << ", energy = " << returnFrac->at(itemp[0]);
            xbinCompPub[i] = returnFrac->at(itemp[0]);
            xbinCompPubErr[i] = 0.;

            for(int j = 0; j < (*nrelem); j++)
            {
               ybinCompPub[j][i] = returnFrac->at(3*j+itemp[0]+1);
               ybinCompPubErrlo[j][i] = TMath::Abs(returnFrac->at(3*j+itemp[0]+2));
               ybinCompPubErrhi[j][i] = TMath::Abs(returnFrac->at(3*j+itemp[0]+3));

     	  cout << ", fraction = " << ybinCompPub[j][i] << ", +" << ybinCompPubErrhi[j][i] << "/-" << ybinCompPubErrlo[j][i];
            }

            cout << endl;
         }
      }
      else
         *addPublishedFrac = false;
   }
   else
      *addPublishedFrac = false;

   for(int j = 0; j < (*nrelem); j++)
   {
      ybinComp[j] = new double[*compnrp];
      ybinCompErrlo[j] = new double[*compnrp];
      ybinCompErrhi[j] = new double[*compnrp];

      dtemp[0] = 0.;

      for(int i = 0; i < (*compnrp); i++)
      {
         itemp[0] = 3*(*nrelem);
         ybinComp[j][i] = composition->at(3*j+itemp[0]*i);
         ybinCompErrhi[j][i] = composition->at(3*j+1+itemp[0]*i);
         ybinCompErrlo[j][i] = composition->at(3*j+2+itemp[0]*i);
      }
   }

   for(int j = 0; j < (*nrelem); j++)
   {
      // Composition from histogram fit
      grComp[j] = new TGraphAsymmErrors(*compnrp, xval, ybinComp[j], xerr, xerr, ybinCompErrlo[j], ybinCompErrhi[j]);
      grComp[j]->SetMarkerStyle(20);
      grComp[j]->SetMarkerSize(1.4);
      grComp[j]->SetLineWidth(2);
      grComp[j]->GetXaxis()->SetRange(18.5, 20.0);
      grComp[j]->GetXaxis()->SetRangeUser(18.5, 20.0);
      grComp[j]->GetXaxis()->SetLimits(18.5, 20.0);
      grComp[j]->GetYaxis()->SetRangeUser(0., 1.2);
      mystyle->SetColorScale(grComp[j], j, *nrelem);

      // Composition from histogram fit
      if(*addPublishedFrac)
      {
         grCompPub[j] = new TGraphAsymmErrors(itemp[2], xbinCompPub, ybinCompPub[j], xbinCompPubErr, xbinCompPubErr, ybinCompPubErrlo[j], ybinCompPubErrhi[j]);

         cout << "New element (" << j+1 << "):" << endl;
         for(int i = 0; i < itemp[2]; i++)
            cout << "Point" << i << " = " << xbinCompPub[i] << ", " << ybinCompPub[j][i] << ", +" << ybinCompPubErrhi[j][i] << "/-" << ybinCompPubErrlo[j][i] << endl;

         grCompPub[j]->SetMarkerStyle(21);
         grCompPub[j]->SetMarkerSize(1.4);
         grCompPub[j]->SetMarkerColor(12);
         grCompPub[j]->SetLineWidth(2);
         grCompPub[j]->SetLineColor(12);
         grCompPub[j]->GetXaxis()->SetRange(18.5, 20.0);
         grCompPub[j]->GetXaxis()->SetRangeUser(18.5, 20.0);
         grCompPub[j]->GetXaxis()->SetLimits(18.5, 20.0);
         grCompPub[j]->GetYaxis()->SetRangeUser(0., 1.2);
      }

*//*      if(j == 0)
      {
         mystyle->SetAxisTitles(grComp[j], "FD energy [log(E/eV)]", "Elemental fractions");
         grComp[j]->Draw("ALP");
      }
      else
      {
         mystyle->SetAxisTitles(grComp[j], "FD energy [log(E/eV)]", "Elemental fractions");
         grComp[j]->Draw("LP;SAME");
      }*//*
   }

   // Plot all four on separate pads
   c1->Divide(1,*nrelem);

   for(int j = 0; j < (*nrelem); j++)
   {
      c1->cd(j+1);
      if(*addPublishedFrac)
      {
         mystyle->SetAxisTitles(grCompPub[j], "", "Elemental fractions");
         grCompPub[j]->GetXaxis()->SetRange(18.5, 20.0);
         grCompPub[j]->GetXaxis()->SetRangeUser(18.5, 20.0);
         grCompPub[j]->GetXaxis()->SetLimits(18.5, 20.0);
         grCompPub[j]->GetYaxis()->SetTitleOffset(1.9 + ((*nrelem)-2)*0.45);
         grCompPub[j]->Draw("ALP");
         grCompPub[j]->GetYaxis()->SetRangeUser(0.,1.);
         mystyle->SetAxisTitles(grComp[j], "", "Elemental fractions");
         grComp[j]->GetXaxis()->SetRange(18.5, 20.0);
         grComp[j]->GetXaxis()->SetRangeUser(18.5, 20.0);
         grComp[j]->GetXaxis()->SetLimits(18.5, 20.0);
         grComp[j]->GetYaxis()->SetTitleOffset(1.9 + ((*nrelem)-2)*0.45);
         grComp[j]->Draw("LP;SAME");
         grComp[j]->GetYaxis()->SetRangeUser(0.,1.);
      }
      else
      {
         mystyle->SetAxisTitles(grComp[j], "", "Elemental fractions");
         grComp[j]->GetXaxis()->SetRange(18.5, 20.0);
         grComp[j]->GetXaxis()->SetRangeUser(18.5, 20.0);
         grComp[j]->GetXaxis()->SetLimits(18.5, 20.0);
         grComp[j]->GetYaxis()->SetTitleOffset(1.9 + ((*nrelem)-2)*0.45);
         grComp[j]->Draw("ALP");
         grComp[j]->GetYaxis()->SetRangeUser(0.,1.);
      }
   }

   c1->Update();

   t->SetTextAlign(31);
   t->SetTextColor(1);
   t->SetTextSize(17);
   for(int j = 0; j < (*nrelem); j++)
   {
      stemp[1] = ptype->GetName(includePart->at(j));
      c1->cd(j+1);
      if(*addPublishedFrac)
      {
         if(stemp[1] == "Oxygen")
         {
            stemp[1] += ", Nitrogen (published)";
            t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), stemp[1].c_str());
         }
         else
            t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), stemp[1].c_str());
      }
      else
         t->DrawText(gPad->GetUxmax(), (gPad->GetUymax())+(0.01*(gPad->GetUymax()-gPad->GetUymin())), stemp[1].c_str());
   }

   stemp[1] = RemoveFilename(filename);

   if(otherSettings[1])
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_xmax_combine.pdf";
   else
      stemp[2] = RemoveFilename(&stemp[1]) + "/plots/mass_composition_plot_combine.pdf";
   c1->SaveAs(stemp[2].c_str());

   if(*addPublishedFrac)
   {
      delete[] xbinCompPub;
      delete[] xbinCompPubErr;

      for(int j = 0; j < (*nrelem); j++)
      {
         delete[] ybinCompPub[j];
         delete[] ybinCompPubErrlo[j];
         delete[] ybinCompPubErrhi[j];
      }
   }

   for(int j = 0; j < (*nrelem); j++)
   {
      delete[] ybinComp[j];
      delete[] ybinCompErrlo[j];
      delete[] ybinCompErrhi[j];
   }

   delete c1;

   delete addPublishedFrac;
   delete returnFrac;
}

// Read published lnA FD results to add them to the plot (type: 0 = EPOS, 1 = QGSJETII, 2 = SIBYLL)
int MvaFitHist::ReadLnaResultsFD(vector<float> *val)
{
   ifstream *infile;
   char *ctemp;
   float *ftemp;
   ret = 0;

   // v2r9p5 simulation production
   if(*simProd == 0)
   {
      itemp[0] = 8;

      // EPOS_LHC
      if(*himodel == 0)
         stemp[0] = string(rootdir) + "/input/lnA_moments_prd2014_epos.txt";
      // QGSJET
      else if(*himodel == 1)
         stemp[0] = string(rootdir) + "/input/lnA_moments_prd2014_qgs.txt";
      // Sibyll
      else if(*himodel == 2)
         stemp[0] = string(rootdir) + "/input/lnA_moments_prd2014_sib.txt";
      else
	 return -1;
   }
   // v3r3p4 simulation production
   else if(*simProd == 1)
   {
      itemp[0] = 10;

      // EPOS_LHC
      if(*himodel == 0)
         stemp[0] = string(rootdir) + "/input/lnA_moments_icrc2017_epos.txt";
      // QGSJET
      else if(*himodel == 1)
         stemp[0] = string(rootdir) + "/input/lnA_moments_icrc2017_qgs.txt";
      // Sibyll
      else if(*himodel == 2)
         stemp[0] = string(rootdir) + "/input/lnA_moments_icrc2017_sib.txt";
      else
	 return -1;
   }
   else
      return -1;

   infile = new ifstream;
   ctemp = new char[1024];
   ftemp = new float[itemp[0]];

   infile->open(stemp[0].c_str(), ifstream::in);

   if(infile->is_open())
   {
      infile->getline(ctemp, 1024, '\n');
      
      while(1)
      {
	 for(int i = 0; i < itemp[0]; i++)
            *infile >> ftemp[i];

	 if(ftemp[0] > 18.40)
	 {
            val->push_back(ftemp[0]);
            val->push_back(ftemp[1]);
            val->push_back(ftemp[2]);

	    cout << ftemp[0] << ", " << ftemp[1] << ", " << ftemp[2] << endl;

	    ret++;
	 }

	 infile->ignore(1,' ');
	 if(infile->eof()) break;
      }
   }

   infile->close();

   delete infile;
   delete[] ftemp;
   delete[] ctemp;

   return ret;
}

// Read published lnA SD results to add them to the plot (type: 0 = EPOS, 1 = QGSJETII; array: 0 = 750m, 1 = 1500m)
int MvaFitHist::ReadLnaResultsSD(vector<float> *val, int array)
{
   ifstream *infile;
   char *ctemp;
   float *ftemp;
   ret = 0;

   // EPOS_LHC
   if(*himodel == 0)
   {
      if(array == 0)
         stemp[0] = string(rootdir) + "/input/lnA_moments_sd750_epos.txt";
      else if(array == 1)
         stemp[0] = string(rootdir) + "/input/lnA_moments_sd1500_epos.txt";
      else
         return -1;
   }
   // QGSJET
   else if(*himodel == 1)
   {
      if(array == 0)
         stemp[0] = string(rootdir) + "/input/lnA_moments_sd750_qgs.txt";
      else if(array == 1)
         stemp[0] = string(rootdir) + "/input/lnA_moments_sd1500_qgs.txt";
      else
         return -1;
   }
   else
      return -1;

   infile = new ifstream;
   ctemp = new char[1024];
   ftemp = new float[4];

   infile->open(stemp[0].c_str(), ifstream::in);

   if(infile->is_open())
   {
      infile->getline(ctemp, 1024, '\n');
      
      while(1)
      {
	 for(int i = 0; i < 4; i++)
            *infile >> ftemp[i];

	 if(ftemp[0] > 18.40)
	 {
            val->push_back(ftemp[0]);
            val->push_back(ftemp[1]);
            val->push_back(ftemp[2]);

	    ret++;
	 }

	 infile->ignore(1,' ');
	 if(infile->eof()) break;
      }
   }

   infile->close();

   delete infile;
   delete[] ftemp;
   delete[] ctemp;

   return ret;
}

// Read published fraction FD results to add them to the plots (type: 0 = EPOS, 1 = QGSJETII, 2 = SIBYLL)
int MvaFitHist::ReadFracResultsFD(vector<float> *val, bool uselna)
{
   ifstream *infile;
   char *ctemp;
   float *ftemp, *lna;
   ret = 0;

   stemp[0] = string(rootdir) + "/input/fractions_icrc17.txt";

   infile = new ifstream;
   ftemp = new float[10];
   lna = new float[3];
   ctemp = new char[1024];

   infile->open(stemp[0].c_str(), ifstream::in);

   if(infile->is_open())
   {
      while(1)
      {
         lna[0] = 0.;
         lna[1] = 0.;
         lna[2] = 0.;

         // Get energy
         *infile >> ftemp[0];
	 if(infile->eof()) break;

	 if(ftemp[0] > 18.40)
	 {
	    val->push_back(ftemp[0]);

            // Read any values from other models (before the one we are using)
            for(int j = 0; j < 20*(*himodel); j++)
               *infile >> ftemp[0];

	    // Get elements
	    itemp[1] = 0;
	    for(int i = 0; i < (*nrelem); i++)
	    {
               itemp[0] = includePart->at(i);
               // Get mean and errors
               for(int j = 0; j < 5; j++)
                  *infile >> ftemp[j];

	       if(!uselna)
	       {
                  val->push_back(ftemp[0]);
                  val->push_back(ftemp[1]);
                  val->push_back(ftemp[2]);
	       }

	       // Using Nitrogen instead of Oxygen for published values
	       if(ptype->GetName(itemp[0]) == "Oxygen")
                  itemp[0]--;

	       if(uselna)
	       {
	          lna[0] += ftemp[0]*TMath::Log(ptype->GetA(itemp[0]));
	          lna[1] += (ftemp[1]*TMath::Log(ptype->GetA(itemp[0])))*(ftemp[1]*TMath::Log(ptype->GetA(itemp[0])));
	          lna[2] += (ftemp[0]*TMath::Log(ptype->GetA(itemp[0])))*(ftemp[2]*TMath::Log(ptype->GetA(itemp[0])));
	       }

	       itemp[1]++;
	    }

            // Read any values from other models (after the one we are using)
            for(int j = 0; j < 20*(2-(*himodel)); j++)
               *infile >> ftemp[0];

	    if(uselna)
	    {
	       // Save lnA values
	       val->push_back(lna[0]);

	       lna[1] = TMath::Sqrt(lna[1]);
	       lna[2] = TMath::Sqrt(lna[2]);

	       val->push_back(lna[1]);
	       val->push_back(lna[2]);
	    }

            ret++;
	 }
	 else
	 {
            for(int j = 0; j < 60; j++)
               *infile >> ftemp[0];
	 }
      }
   }

   infile->close();

   delete infile;
   delete[] ftemp;
   delete[] lna;
   delete[] ctemp;

   return ret;
}*/
