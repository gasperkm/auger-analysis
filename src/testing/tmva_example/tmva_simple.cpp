#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "root_include.h"
#include "mvaefficiency.h"
#include "separate_functions.h"

using namespace std;

template <class T>
string ToString(T value)
{
   ostringstream ss;
   ss << setprecision(5) << value;
   return ss.str();
}

// Add observables to the MVA analysis
int MvaNoteObservables(int count)
{
   // Write out the number of different selected observables
   ofstream fobser;
   fobser.open("/home/gkukec/Gasper/github/auger-analysis/dbg/tmva_example/observables_nr.dat", ofstream::out | ofstream::trunc);
   if(fobser.is_open())
      fobser << count << endl;
   fobser.close();

   return count;
}

int main(int argc, char **argv)
{
   string stemp[2];
   string inname = "example.root";
   double cut = -1;
   double dtemp[4];

   if(argc > 1)
   {
      inname = string(argv[1]);

      if(argc > 2)
         cut = atof(argv[2]);
   }

   string dirname = "/home/gkukec/Gasper/github/auger-analysis/dbg/tmva_example";
   cout << "Directory name = " << dirname << endl;

   TFile *ifile = TFile::Open(inname.c_str(), "READ");
   TFile *ofile = TFile::Open("tmva_output.root", "RECREATE");

   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",ofile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
   (TMVA::gConfig().GetIONames()).fWeightFileDir = "./weights";

   vector<string> obser;
   obser.push_back("xmax");
   obser.push_back("shwsize");
   obser.push_back("risetimerecalc");

   MvaNoteObservables(obser.size());

   size_t foundneg, foundpos;

   for(int i = 0; i < obser.size(); i++)
   {
      stemp[1] = obser[i];

      foundneg = stemp[1].find("_neg");
      foundpos = stemp[1].find("_pos");

      if(foundneg != string::npos)
      {
//	 cout << "Found position of _neg = " << foundneg << endl;
	 stemp[1] = "myvar" + ToString(i+1) + " := " + stemp[1].substr(0, foundneg) + "-" + obser[i];
//	 cout << "stemp[1] = " << stemp[1] << endl;
         factory->AddVariable(stemp[1], 'F');
      }
      else if(foundpos != string::npos)
      {
//	 cout << "Found position of _pos = " << foundpos << endl;
	 stemp[1] = "myvar" + ToString(i+1) + " := " + stemp[1].substr(0, foundpos) + "+" + obser[i];
//	 cout << "stemp[1] = " << stemp[1] << endl;
         factory->AddVariable(stemp[1], 'F');
      }
      else
         factory->AddVariable(obser[i].c_str(), 'F');
   }

   TTree *signalTree = (TTree*)ifile->Get("TreeS1");
   TTree *backgroundTree = (TTree*)ifile->Get("TreeS4");

   int sigcount, backcount;
   sigcount = signalTree->GetEntries();
   backcount = backgroundTree->GetEntries();

   factory->AddSignalTree(signalTree, 1.0);
   factory->AddBackgroundTree(backgroundTree, 1.0);
   factory->PrepareTrainingAndTestTree("", "", "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

   factory->BookMethod(TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator:CreateMVAPdfs");

   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   ifile->Close();
   delete factory;
   ofile->Close();

   double *obslimit;
   obslimit = new double[2*obser.size()];
   char *ctemp;
   ctemp = new char[1024];
   stemp[1] = dirname + "/transformation_stats.dat";
   ifstream fstats;
   fstats.open(stemp[1].c_str(), ifstream::in);
   if(fstats.is_open())
   {
      if(fstats.peek() == '#')
         fstats.getline(ctemp, 1024, '\n');

      for(int i = 0; i < obser.size(); i++)
      {
         fstats >> dtemp[0] >> dtemp[1] >> dtemp[2] >> dtemp[3];
	 obslimit[2*i] = dtemp[2];
	 obslimit[2*i+1] = dtemp[3];
	 cout << "Observable limits: " << obslimit[2*i] << "\t" << obslimit[2*i+1] << endl;
      }
   }

   fstats.close();
   delete[] ctemp;

//   system("./tmvagui tmva_output.root");
   TString *instrname = new TString;
   *instrname = (TString)"tmva_output.root";
   MvaEfficiency *effplot = new MvaEfficiency(sigcount, backcount, &dirname);
   effplot->RunMvaEfficiency(instrname);

   if(cut == -1)
   {
      cout << "Cut value = " << effplot->sigpurCut << endl;
      cut = effplot->sigpurCut;
   }
   else
      cout << "Set cut value = " << cut << endl;

   TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
   float *obsvars;
   obsvars = new float[obser.size()];
   for(int i = 0; i < obser.size(); i++)
   {
      stemp[1] = obser[i];

      foundneg = stemp[1].find("_neg");
      foundpos = stemp[1].find("_pos");

      if(foundneg != string::npos)
      {
//	 cout << "Found position of _neg = " << foundneg << endl;
	 stemp[1] = "myvar" + ToString(i+1) + " := " + stemp[1].substr(0, foundneg) + "-" + obser[i];
//	 cout << "stemp[1] = " << stemp[1] << endl;
         reader->AddVariable(stemp[1], &obsvars[i]);
      }
      else if(foundpos != string::npos)
      {
//	 cout << "Found position of _pos = " << foundpos << endl;
	 stemp[1] = "myvar" + ToString(i+1) + " := " + stemp[1].substr(0, foundpos) + "+" + obser[i];
//	 cout << "stemp[1] = " << stemp[1] << endl;
         reader->AddVariable(stemp[1], &obsvars[i]);
      }
      else
         reader->AddVariable(obser[i].c_str(), &obsvars[i]);
   }
/*   for(int i = 0; i < obser.size(); i++)
      reader->AddVariable(obser[i].c_str(), &obsvars[i]);*/
/*   reader->AddVariable("xmax", &obsvars[0]);
   reader->AddVariable("shwsize", &obsvars[1]);
   reader->AddVariable("risetimerecalc", &obsvars[2]);*/

   reader->BookMVA("MLPBNN method", "./weights/TMVAClassification_MLPBNN.weights.xml");
   ifile = TFile::Open(inname.c_str(), "READ");

   int evcnt[2];

   vector<float> *event;
   double **norm;
   double *mean;
   double **diff;
   double finalerr;
   TMatrixD *cov = new TMatrixD(obser.size(),obser.size());
   TMatrixD *coverr = new TMatrixD(obser.size(),obser.size());

   float *addobsvars;
   addobsvars = new float[2*obser.size()];

   for(int j = 1; j <= ifile->GetNkeys(); j++)
   {
      evcnt[0] = 0;
      evcnt[1] = 0;

      stemp[0] = "TreeS" + ToString(j);
      TTree *signalapp = (TTree*)ifile->Get(stemp[0].c_str());

      for(int i = 0; i < obser.size(); i++)
      {
         stemp[1] = obser[i];

         foundneg = stemp[1].find("_neg");
         foundpos = stemp[1].find("_pos");

         if(foundneg != string::npos)
         {
	    stemp[1] = stemp[1].substr(0, foundneg);
//	    cout << "stemp[1] = " << stemp[1] << endl;
            signalapp->SetBranchAddress(stemp[1].c_str(), &addobsvars[2*i]);
            signalapp->SetBranchAddress(obser[i].c_str(), &addobsvars[2*i+1]);
	 }
         else if(foundpos != string::npos)
	 {
	    stemp[1] = stemp[1].substr(0, foundpos);
//	    cout << "stemp[1] = " << stemp[1] << endl;
            signalapp->SetBranchAddress(stemp[1].c_str(), &addobsvars[2*i]);
            signalapp->SetBranchAddress(obser[i].c_str(), &addobsvars[2*i+1]);
	 }
	 else
            signalapp->SetBranchAddress(obser[i].c_str(), &obsvars[i]);
      }
/*      signalapp->SetBranchAddress("xmax", &obsvars[0]);
      signalapp->SetBranchAddress("shwsize", &obsvars[1]);
      signalapp->SetBranchAddress("risetimerecalc", &obsvars[2]);*/

      event = new vector<float>[obser.size()];

      cout << "Current tree (" << stemp[0] << ") entries: " << signalapp->GetEntries() << endl;
      for(int ievt = 0; ievt < signalapp->GetEntries(); ievt++)
      {
         signalapp->GetEntry(ievt);

         for(int i = 0; i < obser.size(); i++)
         {
            stemp[1] = obser[i];

            foundneg = stemp[1].find("_neg");
            foundpos = stemp[1].find("_pos");

            if(foundneg != string::npos)
	    {
               obsvars[i] = addobsvars[2*i] - addobsvars[2*i+1];
	       if(ievt < 10)
                  cout << "obsvars = " << obsvars[i] << ", addobsvars" << 2*i << " = " << addobsvars[2*i] << ", addobsvars" << 2*i+1 << " = " << addobsvars[2*i+1] << endl;
	    }
            else if(foundpos != string::npos)
	    {
               obsvars[i] = addobsvars[2*i] + addobsvars[2*i+1];
	       if(ievt < 10)
                  cout << "obsvars = " << obsvars[i] << ", addobsvars" << 2*i << " = " << addobsvars[2*i] << ", addobsvars" << 2*i+1 << " = " << addobsvars[2*i+1] << endl;
	    }
	    else
	    {
	       if(ievt < 10)
                  cout << "obsvars = " << obsvars[i] << endl;
	    }
	 }

	 if(reader->EvaluateMVA("MLPBNN method") >= cut)
	 {
//            cout << "     Signal " << ievt << ": MVA value = " << reader->EvaluateMVA("MLPBNN method") << endl;
	    evcnt[0]++;
	 }
	 else
	 {
//            cout << " Background " << ievt << ": MVA value = " << reader->EvaluateMVA("MLPBNN method") << endl;
	    evcnt[1]++;
	 }

         for(int i = 0; i < obser.size(); i++)
            event[i].push_back(obsvars[i]);

	 if(ievt < 10)
	    cout << ievt << ": " << event[0][ievt] << "\t" << event[1][ievt] << "\t" << event[2][ievt] << endl;
      }

      cout << stemp[0] << endl << "   signal events = " << evcnt[0] << endl
	                    << "   background events = " << evcnt[1] << endl;

      if(stemp[0].compare("TreeS5") == 0)
      {
	 cout << "Setting data tree" << endl;
         norm = new double*[obser.size()];
         mean = new double[obser.size()];
         diff = new double*[obser.size()];
         for(int i = 0; i < obser.size(); i++)
         {
            norm[i] = new double[event[0].size()];
            diff[i] = new double[event[0].size()];

	    mean[i] = 0;
         }
         
         // Get normalized values
         for(int j = 0; j < event[0].size(); j++)
         {
            if(j < 10)
               cout << j << ": Norm = ";

            for(int i = 0; i < obser.size(); i++)
	    {
               norm[i][j] = 2*(event[i][j] - obslimit[2*i])/(obslimit[2*i+1] - obslimit[2*i]) - 1.;
	       mean[i] += norm[i][j];

               if(j < 10)
                  cout << norm[i][j] << "\t";
	    }

            if(j < 10)
               cout << endl;
         }

	 cout << "Mean = ";
         for(int i = 0; i < obser.size(); i++)
	 {
            mean[i] = mean[i]/event[0].size();
	    cout << mean[i] << "\t";
	 }
	 cout << endl;

         // Get differences
         for(int j = 0; j < event[0].size(); j++)
         {
            if(j < 10)
               cout << j << ": Diff = ";

            for(int i = 0; i < obser.size(); i++)
	    {
               diff[i][j] = norm[i][j] - mean[i];

               if(j < 10)
                  cout << diff[i][j] << "\t";
	    }

            if(j < 10)
               cout << endl;
         }

         cout << "Matrix:" << endl;
         for(int k = 0; k < obser.size(); k++)
         {
            for(int i = 0; i < obser.size(); i++)
            {
               dtemp[0] = 0;
               for(int j = 0; j < event[0].size(); j++)
               {
                  dtemp[1] = diff[i][j];
                  dtemp[2] = diff[k][j];
                  dtemp[0] += (dtemp[1]*dtemp[2]);
               }
               dtemp[0] = dtemp[0]/(event[0].size()-1.);
               (*cov)(i,k) = dtemp[0];

	       cout << (*cov)(i,k) << " | ";
            }
            
            cout << endl;
         }

         TMatrixDEigen *eigencov = new TMatrixDEigen((const TMatrixD)*cov);
         TMatrixD *dcov = new TMatrixD(obser.size(),obser.size());
         (*dcov) = eigencov->GetEigenValues();
         
         cout << "Diagonal values: ";
         finalerr = 0;
         for(int j = 0; j < obser.size(); j++)
         {
            finalerr += TMath::Power((*dcov)(j,j), 2);
            cout << (*dcov)(j,j) << " ";
         }
         cout << endl;
         
         cout << "Final error = " << TMath::Sqrt(finalerr) << endl;
	 
      }

      delete[] event;
   }

   ifile->Close();
   delete reader;

   delete effplot;
   delete instrname;

   return 0;
}
