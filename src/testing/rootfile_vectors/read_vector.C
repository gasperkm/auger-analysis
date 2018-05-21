#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TSystem.h"

using namespace std;

void read_vector()
{
   TFile *f = TFile::Open("xval_output.root","READ");
   TTree *t;
   f->GetObject("vect",t);

   int nrentries = t->GetEntries();
   cout << "Number of entries in the tree = " << nrentries << endl;

   vector<float> *xval = 0;
   int simpval;
   float anotherval;

   t->SetBranchAddress("simpval",&simpval);
   t->SetBranchAddress("anotherval",&anotherval);

   TBranch *bxval = 0;
   t->SetBranchAddress("xval",&xval,&bxval);

   for(int i = 0; i < nrentries; i++)
   {
      t->GetEntry(i);
      cout << "i = " << i << endl;
      cout << "  simpval = " << simpval << endl;
      cout << "  anotherval = " << anotherval << endl;

      Long64_t tentry = t->LoadTree(i);
      bxval->GetEntry(tentry);

      for(int j = 0; j < xval->size(); ++j)
      {
         cout << "  xval[" << j << "] = " << xval->at(j) << endl;
      }
   }
}
