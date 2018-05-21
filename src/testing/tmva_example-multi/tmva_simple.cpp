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

void SetColor(TH1F *tr, int cur, int nrbins)
{
   int ci = 1738+cur;
   double colormix = (double)cur/(double)(nrbins-1.);
   TColor *color = new TColor(ci, colormix, 0, 1.-colormix, "", 1);
   tr->SetLineColor(ci);
   tr->SetLineWidth(2);
   tr->SetFillColor(ci);
   tr->SetFillStyle(3004+cur);
}


int main(int argc, char **argv)
{
   string stemp[3];
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
   TFile *ofile = TFile::Open("tmva_multiclass_output.root", "RECREATE");

   TMVA::Factory *factory = new TMVA::Factory("TMVAMulticlass",ofile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass");
   (TMVA::gConfig().GetIONames()).fWeightFileDir = "./weights";

   vector<string> obser;
   obser.push_back("xmax");
   obser.push_back("shwsize");
//   obser.push_back("risetimerecalc");

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

   vector<string> treeNames;
   treeNames.push_back("TreeS1");
   treeNames.push_back("TreeS2");
   treeNames.push_back("TreeS3");
   treeNames.push_back("TreeS4");

   TTree *fileTree[10];

   for(int t = 0; t < treeNames.size(); t++)
   {
      fileTree[t] = (TTree*)ifile->Get(treeNames[t].c_str());
      factory->AddTree(fileTree[t], treeNames[t]);
   }

/*   backgroundTree[0] = (TTree*)ifile->Get(treeNames[0].c_str());
   backgroundTree[1] = (TTree*)ifile->Get(treeNames[1].c_str());
   backgroundTree[2] = (TTree*)ifile->Get(treeNames[2].c_str());
   backgroundTree[3] = (TTree*)ifile->Get(treeNames[3].c_str());

   factory->AddTree(backgroundTree[0], treeNames[0]);
   factory->AddTree(backgroundTree[1], treeNames[1]);
   factory->AddTree(backgroundTree[2], treeNames[2]);
   factory->AddTree(backgroundTree[3], treeNames[3]);*/

   factory->PrepareTrainingAndTestTree("", "SplitMode=Random:NormMode=NumEvents:!V");

//   factory->BookMethod(TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator:CreateMVAPdfs");
   factory->BookMethod(TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   ifile->Close();
   delete factory;
   ofile->Close();

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

   reader->BookMVA("BDTG method", "./weights/TMVAMulticlass_BDTG.weights.xml");
   ifile = TFile::Open(inname.c_str(), "READ");

   vector<float> *event;
   double **norm;
   double *mean;
   double **diff;
   double finalerr;

   float *addobsvars;
   addobsvars = new float[2*obser.size()];

   TH1F *hist[10], *mvahist[10];
   TCanvas *c1 = new TCanvas("c1","",1200,900);
   TCanvas *c2 = new TCanvas("c2","",1200,900);
   TCanvas *c3 = new TCanvas("c3","",1200,900);
   TCanvas *c4 = new TCanvas("c4","",1200,900);
   c1->SetGrid();
   c2->SetGrid();
   c3->SetGrid();
   c4->SetGrid();

   double ymax[2];
   ymax[0] = -1;
   ymax[1] = -1;
   double mvaValue;
   double multiValue;

   for(int t = 0; t < treeNames.size(); t++)
   {
      stemp[2] = "MVA_BDTG_MULTI_" + treeNames[t];//ToString(t+1);
      hist[t] = new TH1F(stemp[2].c_str(), stemp[2].c_str(), 80, -0.1, 1.2);
      stemp[2] = "MVA_BDTG_" + treeNames[t];//ToString(t+1);
      mvahist[t] = new TH1F(stemp[2].c_str(), stemp[2].c_str(), 80, -1.2, 1.2);
   }
   stemp[2] = "rm -fr mva_bdtg_multi_Tree*.pdf mva_bdtg_Tree*.pdf";
   system(stemp[2].c_str());
/*   hist[0] = new TH1F("MVA_BDTG_proton", "MVA_BDTG_proton", 80, 0., 1.1);
   hist[1] = new TH1F("MVA_BDTG_helium", "MVA_BDTG_helium", 80, 0., 1.1);
   hist[2] = new TH1F("MVA_BDTG_oxygen", "MVA_BDTG_oxygen", 80, 0., 1.1);
   hist[3] = new TH1F("MVA_BDTG_iron",   "MVA_BDTG_iron",   80, 0., 1.1);*/

//   for(int j = 1; j <= ifile->GetNkeys(); j++)
   for(int t = 0; t < treeNames.size(); t++)
   {
/*      stemp[0] = "TreeS" + ToString(j);
      TTree *signalapp = (TTree*)ifile->Get(stemp[0].c_str());*/
      TTree *signalapp = (TTree*)ifile->Get(treeNames[t].c_str());

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

/*	 if(reader->EvaluateMVA("MLPBNN method") >= cut)
	 {
//            cout << "     Signal " << ievt << ": MVA value = " << reader->EvaluateMVA("MLPBNN method") << endl;
	    evcnt[0]++;
	 }
	 else
	 {
//            cout << " Background " << ievt << ": MVA value = " << reader->EvaluateMVA("MLPBNN method") << endl;
	    evcnt[1]++;
	 }*/

	 mvaValue = reader->EvaluateMVA("BDTG method");
	 multiValue = reader->EvaluateMulticlass(0,"BDTG method");

         for(int i = 0; i < obser.size(); i++)
            event[i].push_back(obsvars[i]);

	 if(ievt < 10)
	 {
	    cout << ievt << ": " << event[0][ievt] << "\t" << event[1][ievt] << endl;
	    if(t < treeNames.size())
	    {
	       cout << "  Multiclass evaluate (" << t << "): " << multiValue << endl;
	       cout << "  MVA evaluate (" << t << "): " << mvaValue << endl;
	    }
	 }

         hist[t]->Fill(multiValue);
         mvahist[t]->Fill(mvaValue);
      }

      c1->cd();
      SetColor(hist[t], t, treeNames.size());
      if(t == 0)
         hist[t]->Draw();
      else
         hist[t]->Draw("SAME");

      if(hist[t]->GetMaximum() > ymax[0])
         ymax[0] = hist[t]->GetMaximum();

      c2->cd();
      SetColor(mvahist[t], t, treeNames.size());
      if(t == 0)
         mvahist[t]->Draw();
      else
         mvahist[t]->Draw("SAME");

      if(mvahist[t]->GetMaximum() > ymax[1])
         ymax[1] = mvahist[t]->GetMaximum();

      c3->cd();
      hist[t]->Draw();
      stemp[2] = "mva_bdtg_multi_" + treeNames[t] + ".pdf";
      c3->SaveAs(stemp[2].c_str());

      c4->cd();
      mvahist[t]->Draw();
      stemp[2] = "mva_bdtg_" + treeNames[t] + ".pdf";
      c4->SaveAs(stemp[2].c_str());

/*      cout << stemp[0] << endl << "   signal events = " << evcnt[0] << endl
	                    << "   background events = " << evcnt[1] << endl;*/

      delete[] event;
   }

   hist[0]->GetYaxis()->SetRangeUser(0,(1.05*ymax[0]));
   stemp[2] = "mva_bdtg_multi.pdf";
   c1->SaveAs(stemp[2].c_str());

   mvahist[0]->GetYaxis()->SetRangeUser(0,(1.05*ymax[1]));
   stemp[2] = "mva_bdtg.pdf";
   c2->SaveAs(stemp[2].c_str());

   for(int t = 0; t < treeNames.size(); t++)
   {
      delete hist[t];
      delete mvahist[t];
   }
/*   delete hist[0];
   delete hist[1];
   delete hist[2];
   delete hist[3];*/

/*   double *obslimit;
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
   delete instrname;*/

   return 0;
}
