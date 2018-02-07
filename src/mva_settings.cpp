#include "frame.h"
#include <fstream>
#include <algorithm>

// Compare method names from values on the drop-down list
string MyFrame::GetMethodName(string name)
{
   vector<string>::iterator it;
   it = find(methodsDesc.begin(), methodsDesc.end(), name);
//   cout << "Finding 1 = " << *it << ", " << methodsDesc[distance(methodsDesc.begin(), it)] << endl;
   if(it != methodsDesc.end())
      return methods[distance(methodsDesc.begin(), it)];
   it = find(methods.begin(), methods.end(), name);
//   cout << "Finding 2 = " << *it << ", " << methods[distance(methods.begin(), it)] << endl;
   if(it != methods.end())
      return methodsDesc[distance(methods.begin(), it)];

/*   for(int i = 0; i < nrmethods; i++)
   {
      if(name == methodsDesc[i])
         return methods[i];

      if(name == methods[i])
         return methodsDesc[i];
   }*/
/*   if(name == "Test all methods")
      return "All";
   else if(name == "Cut optimization (Cuts)")
      return "Cuts";
   else if(name == "Cut optimization (CutsD)")
      return "CutsD";
   else if(name == "1D likelihood (Likelihood)")
      return "Likelihood";
   else if(name == "1D likelihood (LikelihoodPCA)")
      return "LikelihoodPCA";
   else if(name == "Multidimensional likelihood (PDERS)")
      return "PDERS";
   else if(name == "Nearest neighbours (KNN)")
      return "KNN";
   else if(name == "Fisher discriminants (Fisher)")
      return "Fisher";
   else if(name == "Linear discriminant (LD)")
      return "LD";
   else if(name == "Functional discriminant (FDA_GA)")
      return "FDA_GA";
   else if(name == "Neural network (MLPBNN)")
      return "MLPBNN";
   else if(name == "Support vector machine (SVM)")
      return "SVM";
   else if(name == "Boosted decision trees (BDT)")
      return "BDT";
   else if(name == "Friedman's rulefit (RuleFit)")
      return "RuleFit";*/

/*   if(name == "All")
      return "Test all methods";
   else if(name == "Cuts")
      return "Cut optimization (Cuts)";
   else if(name == "CutsD")
      return "Cut optimization (CutsD)";
   else if(name == "Likelihood")
      return "1D likelihood (Likelihood)";
   else if(name == "LikelihoodPCA")
      return "1D likelihood (LikelihoodPCA)";
   else if(name == "PDERS")
      return "Multidimensional likelihood (PDERS)";
   else if(name == "KNN")
      return "Nearest neighbours (KNN)";
   else if(name == "Fisher")
      return "Fisher discriminants (Fisher)";
   else if(name == "LD")
      return "Linear discriminant (LD)";
   else if(name == "FDA_GA")
      return "Functional discriminant (FDA_GA)";
   else if(name == "MLPBNN")
      return "Neural network (MLPBNN)";
   else if(name == "SVM")
      return "Support vector machine (SVM)";
   else if(name == "BDT")
      return "Boosted decision trees (BDT)";
   else if(name == "RuleFit")
      return "Friedman's rulefit (RuleFit)";*/
}
