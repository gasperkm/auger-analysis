#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TSystem.h"

using namespace std;

void write_vector()
{
   TFile *f = TFile::Open("xval_output.root","RECREATE");

   vector<float> xval;
   int simpval;
   float anotherval;

   TTree *t = new TTree("vect","Tree with a vector");
   t->Branch("xval",&xval);
   t->Branch("simpval",&simpval,"simpval/I");
   t->Branch("anotherval",&anotherval,"anotherval/F");

   gRandom->SetSeed();

   float p[3];

   for(int i = 0; i < 10; i++)
   {
      int npx = (int)(gRandom->Rndm(1)*15)+1;

      xval.clear();

      simpval = npx;
      anotherval = gRandom->Gaus();

      cout << "i = " << i << endl;
      cout << "  simpval = " << simpval << endl;
      cout << "  anotherval = " << anotherval << endl;

      for(int j = 0; j < npx; ++j)
      {
         xval.push_back(gRandom->Gaus());
         cout << "  xval[" << j << "] = " << xval[j] << endl;
      }

      t->Fill();
   }
   f->Write();

   delete f;
}
