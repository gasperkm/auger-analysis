#include "mva_methods.h"
#include "workstation.h"

void PrintMethods()
{
   cout << "Available MVA methods:" << endl
        << "- Cuts" << endl
        << "- CutsD" << endl
        << "- CutsPCA" << endl
        << "- CutsSA" << endl
        << "- Likelihood" << endl
        << "- LikelihoodD" << endl
        << "- LikelihoodPCA" << endl
        << "- LikelihoodKDE" << endl
        << "- PDERS" << endl
        << "- PDERSD" << endl
        << "- PDERSPCA" << endl
        << "- PDEFoam" << endl
        << "- PDEFoamBoost" << endl
        << "- KNN" << endl
        << "- LD" << endl
        << "- Fisher" << endl
        << "- FisherG" << endl
        << "- BoostedFisher" << endl
        << "- HMatrix" << endl
        << "- FDA_GA" << endl
        << "- FDA_SA" << endl
        << "- FDA_MC" << endl
        << "- FDA_MT" << endl
        << "- FDA_GAMT" << endl
        << "- FDA_MCMT" << endl
        << "- MLP" << endl
        << "- MLPBFGS" << endl
        << "- MLPBNN" << endl
        << "- CFMlpANN" << endl
        << "- TMlpANN" << endl
        << "- SVM" << endl
        << "- BDT" << endl
        << "- BDTG" << endl
        << "- BDTB" << endl
        << "- BDTD" << endl
        << "- BDTF" << endl
        << "- RuleFit" << endl;
}

double GetMethodMin(TString name)
{
   if(name.Contains("CutsD"))
      return 0.;
   else if(name.Contains("Cuts"))
      return 0.;
   else if(name.Contains("LikelihoodPCA"))
      return 0.;
   else if(name.Contains("Likelihood"))
      return 0.;
   else if(name.Contains("PDERS"))
      return 0.;
   else if(name.Contains("KNN"))
      return 0.;
   else if(name.Contains("LD"))
      return -0.1;
   else if(name.Contains("BoostedFisher"))
      return -1.;
   else if(name.Contains("Fisher"))
      return -3.;
   else if(name.Contains("FDA_GA"))
      return 0.;
   else if(name.Contains("MLPBNN"))
      return -0.1;
   else if(name.Contains("MLPBFGS"))
      return -0.1;
   else if(name.Contains("SVM"))
      return 0.;
   else if(name.Contains("BDT"))
      return -1.;
//      return -1.1;
   else if(name.Contains("RuleFit"))
      return -1.6;
   else
      return -1.;
}

double GetMethodMax(TString name)
{
   if(name.Contains("CutsD"))
      return 1.;
   else if(name.Contains("Cuts"))
      return 1.;
   else if(name.Contains("LikelihoodPCA"))
      return 1.;
   else if(name.Contains("Likelihood"))
      return 1.;
   else if(name.Contains("PDERS"))
      return 1.;
   else if(name.Contains("KNN"))
      return 1.;
   else if(name.Contains("LD"))
      return 4.;
   else if(name.Contains("BoostedFisher"))
      return 1.4;
   else if(name.Contains("Fisher"))
      return 5.;
   else if(name.Contains("FDA_GA"))
      return 1.;
   else if(name.Contains("MLPBNN"))
      return 1.1;
   else if(name.Contains("MLPBFGS"))
      return 1.1;
   else if(name.Contains("SVM"))
      return 1.;
   else if(name.Contains("BDT"))
      return 1.;
//      return 1.1;
   else if(name.Contains("RuleFit"))
      return 1.6;
   else
      return 1.;
}

bool GetMethodDist(TString name)
{
   if(name.Contains("CutsD"))
      return 0;
   else if(name.Contains("Cuts"))
      return 0;
   else if(name.Contains("LikelihoodPCA"))
      return 0;
   else if(name.Contains("Likelihood"))
      return 0;
   else if(name.Contains("PDERS"))
      return 0;
   else if(name.Contains("KNN"))
      return 0;
   else if(name.Contains("LD"))
      return 1;
   else if(name.Contains("BoostedFisher"))
      return 0;
   else if(name.Contains("Fisher"))
      return 1;
   else if(name.Contains("FDA_GA"))
      return 0;
   else if(name.Contains("MLPBNN"))
      return 0;
   else if(name.Contains("MLPBFGS"))
      return 0;
   else if(name.Contains("SVM"))
      return 0;
   else if(name.Contains("BDT"))
      return 0;
   else if(name.Contains("RuleFit"))
      return 0;
   else
      return 0;
}
