#define _STANDALONE_ 1
#include "workstation.h"
#include "adst_mva.h"
#include "separate_functions.h"

#include "src/adst_mva.cpp"
#include "src/calc_observables.cpp"

using namespace std;

int main(int argc, char **argv)
{
   int nrfiles = argc;
   if(nrfiles < 2)
   {
      cout << "Error! No input files supplied." << endl << "Usage instructions:" << endl << "  ./start.sh rewrite [ADST version] [input files]" << endl << "  ./bin/batch_rewrite-[ADST version] [input files]" << endl;
      return 1;
   }

   // Initialize observables
   ifstream ifs;
   string *stemp = new string[3];
   int *itemp = new int;
   float *ftemp = new float[2];
   char *ctemp = new char[1024];
   stemp[0] = string(rootdir) + "/input/observables.txt";
   ifs.open(stemp[0].c_str(), ifstream::in);

   Observables *generalObservables;
   int nrobs;
   double mvalimit[2];
   vector<string> observables;
   vector<bool> obssel;
   vector<bool> obsorigsel;

   vector<float> tempmin;
   vector<float> tempmax;
   vector<string> tempdesc;

   if(ifs.is_open())
   {
      if(!observables.empty())
         observables.erase(observables.begin(), observables.end());
      if(!obssel.empty())
         obssel.erase(obssel.begin(), obssel.end());
      if(!obsorigsel.empty())
         obsorigsel.erase(obsorigsel.begin(), obsorigsel.end());
      if(!tempmin.empty())
         tempmin.erase(tempmin.begin(), tempmin.end());
      if(!tempmax.empty())
         tempmax.erase(tempmax.begin(), tempmax.end());
      if(!tempdesc.empty())
         tempdesc.erase(tempdesc.begin(), tempdesc.end());

      while(1)
      {
         if(ifs.peek() == '#')
         {
            ifs.getline(ctemp, 1024);
         }
         else
         {
            ifs >> stemp[1] >> *itemp >> ftemp[0] >> ftemp[1];
            if(ifs.eof()) break;
            observables.push_back(stemp[1]);
            obssel.push_back((bool)*itemp);
            obsorigsel.push_back((bool)*itemp);
            tempmin.push_back(ftemp[0]);
            tempmax.push_back(ftemp[1]);

            ifs.getline(ctemp, 1024);
            stemp[2] = string(ctemp);
            tempdesc.push_back(stemp[2]);

cout << "Values: " << stemp[1] << "\t" << *itemp << "\t" << ftemp[0] << "\t" << ftemp[1] << "\t" << stemp[2] << endl;
         }
      }
   }

   nrobs = observables.size();
   cout << "Number of observables = " << nrobs << endl;

   generalObservables = new Observables(observables);
   for(int i = 0; i < nrobs; i++)
   {
      generalObservables->SetMin(i, tempmin[i]);
      generalObservables->SetMax(i, tempmax[i]);
      generalObservables->SetLabel(i, tempdesc[i]);
   }

   mvalimit[0] = -2.;
   mvalimit[1] = 2.;

   ifs.close();

   cout << endl << "----------------------------------------------" << endl;
   int ret;

   for(int i = 1; i < nrfiles; i++)
   {
      stemp[0] = string(argv[i]);
      cout << "Rewriting file: " << stemp[0] << endl;

      // Prepare writeout holders
      AdstMva *mvatool = new AdstMva;
      stemp[1] = RemoveFilename(&stemp[0]);
      stemp[2] = RemovePath(&stemp[0]);
      mvatool->outname = stemp[1] + "/rewrite_" + stemp[2];

      cout << "To file: " << mvatool->outname << endl;

      Observables *obssig[3];
      Observables *obsall[3];

      for(int j = 0; j < 3; j++)
      {
         obssig[j] = new Observables(observables);
         obsall[j] = new Observables(observables);
      }

      mvatool->outfile = TFile::Open((mvatool->outname).c_str(), "RECREATE");
      mvatool->all_tree = new TTree("TreeA", "Background tree with all events, including signal events.");
      mvatool->rewritecode = 0;

      mvatool->inname = stemp[0];
      ret = mvatool->RewriteObservables(1, 0, obssig, obsall);
      if(ret == -1)
      {
         for (int j = 0; j < 3; j++)
         {
            delete obssig[j];
            delete obsall[j];
         }

	 cout << "Error! Observables not rewritten." << endl;

         (mvatool->outfile)->Close();
         stemp[0] = "rm -fr " + (mvatool->outname);
         system(stemp[0].c_str());

         delete mvatool;

/*         delete[] stemp;
         delete itemp;
         delete[] ftemp;
         delete[] ctemp;
         return -1;*/
      }
      else
      {
         mvatool->PrepareOtherTrees(1, 0, observables);

         (mvatool->all_tree)->Write();
         (mvatool->outfile)->Close();

         for (int j = 0; j < 3; j++)
         {
            delete obssig[j];
            delete obsall[j];
         }

         delete mvatool;
      }

      cout << endl;
   }

   delete[] stemp;
   delete itemp;
   delete[] ftemp;
   delete[] ctemp;

   return 0;
}
