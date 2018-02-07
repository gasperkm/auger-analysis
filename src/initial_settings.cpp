#include "frame.h"
#include "popups.h"
#include <fstream>

void MyFrame::InitObservables()
{
   ifstream ifs;
   string *stemp = new string[3];
   int *itemp = new int;
   float *ftemp = new float[2];
   char *ctemp = new char[1024];
   stemp[0] = string(rootdir) + "/input/observables.txt";
   ifs.open(stemp[0].c_str(), ifstream::in);

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
	    if(DBGSIG > 1)
               cout << "Line: " << ctemp << endl;
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

	    if(DBGSIG > 0)
               cout << "Values: " << stemp[1] << "\t" << *itemp << "\t" << ftemp[0] << "\t" << ftemp[1] << "\t" << stemp[2] << endl;
	 }
      }
   }
   else
   {
      AlertPopup("Observables list error", "Can't open the observables list file " + stemp[0] + ". Please make sure it exists. If not, create it.");
   }

   nrobs = observables.size();
   if(DBGSIG > 0)
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
   delete[] stemp;
   delete itemp;
   delete[] ftemp;
   delete[] ctemp;
}

void MyFrame::InitMethods()
{
   ifstream ifs;
   string *stemp = new string[4];
   char *ctemp = new char[1024];
   stemp[0] = string(rootdir) + "/input/mva_options.txt";
   ifs.open(stemp[0].c_str(), ifstream::in);

   if(ifs.is_open())
   {
      if(!methods.empty())
         methods.erase(methods.begin(), methods.end());
      if(!methodsOpt.empty())
         methodsOpt.erase(methodsOpt.begin(), methodsOpt.end());
      if(!methodsDesc.empty())
         methodsDesc.erase(methodsDesc.begin(), methodsDesc.end());

      methods.push_back("All");
      methodsDesc.push_back("Test all methods");
      methodsOpt.push_back("AllMethods");

      while(1)
      {
         if(ifs.peek() == '#')
         {
            ifs.getline(ctemp, 1024);
            if(DBGSIG > 1)
               cout << "Line: " << ctemp << endl;
         }
         else
         {
            ifs >> stemp[1] >> stemp[2] >> stemp[3];
//	    cout << "Method:" << endl << "  " << stemp[1] << endl << "  " << stemp[2] << endl << "  " << stemp[3] << endl;
            methods.push_back(stemp[1]);
            methodsDesc.push_back(stemp[2]);
            methodsOpt.push_back(stemp[3]);
	    ifs.ignore(1,' ');
            if(ifs.eof()) break;
         }
      }

      nrmethods = methods.size();
   }

   ifs.close();
   delete[] stemp;
   delete[] ctemp;
/*}
   if(!methods.empty())
      methods.erase(methods.begin(), methods.end());

   methods.push_back("All");
   methods.push_back("Cuts");
   methods.push_back("CutsD");
   methods.push_back("Likelihood");
   methods.push_back("LikelihoodPCA");
   methods.push_back("PDERS");
   methods.push_back("KNN");
   methods.push_back("Fisher");
   methods.push_back("LD");
   methods.push_back("FDA_GA");
   methods.push_back("MLPBNN");
   methods.push_back("SVM");
   methods.push_back("BDT");
   methods.push_back("RuleFit");

   nrmethods = methods.size();*/
}

void MyFrame::InitVariables()
{
   oldselect[0] = 0;
   oldselect[1] = 2;
   oldselect[2] = 3;

   currentMvaDir = new string;
   *currentMvaDir = string(rootdir) + "/results";
   currentRewriteDir = new string;
   *currentRewriteDir = string(rootdir) + "/results";
   currentAnalysisDir = new string;
   *currentAnalysisDir = string(rootdir) + "/results";
   currentCutsDir = new string;
   *currentCutsDir = string(rootdir) + "/results";
   currentCutsInputDir = new string;
   *currentCutsInputDir = string(rootdir) + "/results";

   freshAnalysis = false;

   nrlists = 0;
}
