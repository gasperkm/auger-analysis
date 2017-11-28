#include "combine.h"

void CombineRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames )
{
   cout << "Combining root files." << endl;
   for(int i = 0; i < filenames.size(); i++) cout << "  " << filenames[i] << endl;

   int sigKeys = 0;

   for(int i = 0; i < nrkeys.size(); i++)
      sigKeys += (nrkeys[i]-1)/2;

   cout << "Number of signal keys in all files = " << sigKeys << endl;

   int nrfiles = sourcelist->GetSize();
   string strOldS, strNewS;
   TTree *treeOldS, *treeNewS, *treeAll;
   TTree *outSig[sigKeys], *outAll;
   int selectCount[2];
   int fileCount = 0;

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove(0,2);

   TFile *nextsource = (TFile*)sourcelist->First();
   nextsource->cd(path);
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   cout << "Get number of files = " << nrfiles << endl;

   fileCount = 0;
   selectCount[0] = 0;
   selectCount[1] = 0;

   TList *listAll, *listSig[sigKeys];
   target->cd();
   listAll = new TList;
   for(int j = 0; j < sigKeys; j++)
      listSig[j] = new TList;

   while(nextsource)
   {
      cout << "Number of keys in file " << nextsource->GetName() << ": " << nrkeys[fileCount] << "(" << (nrkeys[fileCount]-1)/2 << ")" << endl;
      
      for(int j = 1; j <= (nrkeys[fileCount]-1)/2; j++)
      {
         strOldS = "TreeOldS" + ToString(j);
         strNewS = "TreeNewS" + ToString(j);

         nextsource->cd(path);

         treeOldS = (TTree*)nextsource->Get(strOldS.c_str());
         treeNewS = (TTree*)nextsource->Get(strNewS.c_str());
         treeAll = (TTree*)nextsource->Get("TreeA");

         // Merge the signal trees
         printf("#- %d %d a ---------------------------------------------\n", fileCount, j);

         if(treeOldS->GetEntries() > 0)	// Mean values
         {
            printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeOldS->GetEntries(), strOldS.c_str(), nextsource->GetKey(strOldS.c_str())->GetTitle());
            listSig[selectCount[0]]->Add(treeOldS);
            selectCount[0]++;
         }

         if(treeNewS->GetEntries() > 0)	// Mean values
         {
            printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeNewS->GetEntries(), strNewS.c_str(), nextsource->GetKey(strNewS.c_str())->GetTitle());
            listSig[selectCount[0]]->Add(treeNewS);
            selectCount[0]++;
         }
      }

      cout << "selectCount = " << selectCount[0] << ", " << selectCount[1] <</* ", " << selectCount[2] <<*/ endl;

      // Merge all events
      printf("#- %d c ---------------------------------------------\n", fileCount);
      if(treeAll->GetEntries() > 0)	// Mean values
      {
         printf("Found %d entries in tree \"TreeA\" and title \"%s\"\n", (int)treeAll->GetEntries(), nextsource->GetKey("TreeA")->GetTitle());
         listAll->Add(treeAll);
      }

      // Move to the next file
      printf("#- %d d -----------------------------------------------\n", fileCount);
      fileCount++;
      nextsource = (TFile*)sourcelist->After(nextsource);
   }

   printf("#- e saving the file --------------------------------------\n");

   string treeTitle[sigKeys];

   target->cd();
   printf("Saving signal trees\n");
   for(int i = 0; i < sigKeys; i++)
   {
      outSig[i] = TTree::MergeTrees(listSig[i]);
      cout << "Events in output signal tree = " << outSig[i]->GetEntries() << endl;
      strOldS = "TreeS" + ToString(i+1);
      outSig[i]->SetName(strOldS.c_str());
      treeTitle[i] = outSig[i]->GetTitle();
      outSig[i]->Write();
   }
   printf("Saving all events tree\n");
   outAll = TTree::MergeTrees(listAll);
   cout << "Events in output total tree = " << outAll->GetEntries() << endl;
   outAll->Write();

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);

   for(int j = 0; j < sigKeys; j++)
      delete listSig[j];
   delete listAll;
}

void MergeRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames, int opt )
{
   cout << "Merging root files." << endl;
   for(int i = 0; i < filenames.size(); i++) cout << "  " << filenames[i] << endl;

   int sigKeys = 0;

   for(int i = 0; i < nrkeys.size(); i++)
      sigKeys += (nrkeys[i]-1)/2;
   cout << "Number of signal keys in all files = " << sigKeys << endl;

   int nrfiles = sourcelist->GetSize();
   string strOldS, strNewS;
   TTree *treeOldS, *treeNewS, *treeAll;
   TTree *outOld, *outNew, *outAll;
   int selectCount[2];
   int fileCount = 0;

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove(0,2);

   TFile *nextsource = (TFile*)sourcelist->First();
   nextsource->cd(path);
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   cout << "Get number of files = " << nrfiles << endl;

   fileCount = 0;
   selectCount[0] = 0;
   selectCount[1] = 0;
    
   TList *listAll, *listOld, *listNew;
   target->cd();
   listAll = new TList;
   listOld = new TList;
   listNew = new TList;

   while(nextsource)
   {
      for(int j = 1; j <= (nrkeys[fileCount]-1)/2; j++)
      {
         strOldS = "TreeOldS" + ToString(j);
         strNewS = "TreeNewS" + ToString(j);

         nextsource->cd(path);

         treeOldS = (TTree*)nextsource->Get(strOldS.c_str());
         treeNewS = (TTree*)nextsource->Get(strNewS.c_str());
         treeAll = (TTree*)nextsource->Get("TreeA");

         // Merge the signal trees
         printf("#- %d %d a ---------------------------------------------\n", fileCount, j);

         if(treeOldS->GetEntries() > 0)	// All Old trees
         {
            printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeOldS->GetEntries(), strOldS.c_str(), nextsource->GetKey(strOldS.c_str())->GetTitle());
            listOld->Add(treeOldS);
            selectCount[0]++;
         }

         if(treeNewS->GetEntries() > 0)	// All New trees
         {
            printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeNewS->GetEntries(), strNewS.c_str(), nextsource->GetKey(strNewS.c_str())->GetTitle());
            listNew->Add(treeNewS);
            selectCount[1]++;
         }
      }

      cout << "selectCount = " << selectCount[0] << ", " << selectCount[1] << ", " << selectCount[2] << endl;

      // Merge all events
      printf("#- %d c ---------------------------------------------\n", fileCount);
      if(treeAll->GetEntries() > 0)	// Mean values
      {
         printf("Found %d entries in tree \"TreeA\" and title \"%s\"\n", (int)treeAll->GetEntries(), nextsource->GetKey("TreeA")->GetTitle());
         listAll->Add(treeAll);
      }

      // Move to the next file
      printf("#- %d d -----------------------------------------------\n", fileCount);
      fileCount++;
      nextsource = (TFile*)sourcelist->After(nextsource);
   }

   cout << "Entries in old list: " << listOld->GetEntries() << endl;
   cout << "Entries in new list: " << listNew->GetEntries() << endl;
   cout << "Entries in all list: " << listAll->GetEntries() << endl;

   printf("#- e saving the file --------------------------------------\n");

   string treeTitleOld, treeTitleNew;

   target->cd();
   printf("Saving old signal tree\n");

   if(listOld->GetEntries() > 0)
   {
      outOld = TTree::MergeTrees(listOld);
      strOldS = "TreeOldS1";
      outOld->SetName(strOldS.c_str());
   }
   else
   {
      outOld = (TTree*)((TFile*)sourcelist->First())->Get("TreeOldS1");
   }

   treeTitleOld = outOld->GetTitle();
   outOld->Write();
   
   printf("Saving new signal tree\n");
   if(listNew->GetEntries() > 0)
   {
      outNew = TTree::MergeTrees(listNew);
      strOldS = "TreeNewS1";
      outNew->SetName(strOldS.c_str());
   }
   else
   {
      outNew = (TTree*)((TFile*)sourcelist->First())->Get("TreeNewS1");
   }
   treeTitleNew = outNew->GetTitle();
   outNew->Write();

   printf("Saving all events tree\n");
   outAll = TTree::MergeTrees(listAll);
   outAll->Write();

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);

   delete listAll;
   delete listOld;
   delete listNew;
}

void CheckKeys( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts )
{
   int nrfiles = sourcelist->GetSize();
   int fileCount = 0;
   string strOldS, strNewS;
   TTree *treeOldS, *treeNewS, *treeAll;

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove(0,2);

   TFile *nextsource = (TFile*)sourcelist->First();
   nextsource->cd(path);
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   while(nextsource)
   {
      cout << fileCount << ": Number of keys = " << nextsource->GetNkeys() << endl;
      nrkeys.push_back(nextsource->GetNkeys());

      treeAll = (TTree*)nextsource->Get("TreeA");
cout << "Number of events in tree TreeA = " << treeAll->GetEntries() << endl;
      nrevts.push_back(treeAll->GetEntries());

      for(int j = 1; j <= (nrkeys[fileCount]-1)/2; j++)
      {
         strOldS = "TreeOldS" + ToString(j);
         strNewS = "TreeNewS" + ToString(j);

         treeOldS = (TTree*)nextsource->Get(strOldS.c_str());
         treeNewS = (TTree*)nextsource->Get(strNewS.c_str());

cout << "Number of events in tree " << strOldS << " = " << treeOldS->GetEntries() << endl;
cout << "Number of events in tree " << strNewS << " = " << treeNewS->GetEntries() << endl;
         nrevts.push_back(treeOldS->GetEntries());
         nrevts.push_back(treeNewS->GetEntries());
      }

      nextsource = (TFile*)sourcelist->After(nextsource);
      fileCount++;
   }
}

void hmerge(int nrfiles, string *files, string *outname)
{
   cout << "Number of files: " << nrfiles << endl;
   cout << "Saving into file: " << *outname << endl;

   TList *FileList;
   TFile *Target;

   char ctemp[1024];
   vector<int> nrkeys;
   vector<int> nrevts;
   vector<string> filenames;

   // Prepare the files to be merged
   sprintf(ctemp, "%s/results/hsimple1.root", rootdir);
   if(gSystem->AccessPathName(ctemp))
   {
     for(int i = 0; i < nrfiles; i++)
     {
       sprintf(ctemp, "%s/results/hsimple%d.root", rootdir, i+1);
       gSystem->CopyFile(files[i].c_str(), ctemp);
       filenames.push_back(files[i]);
     }
   }
   else
     return;

   Target = TFile::Open(outname->c_str(), "RECREATE");

   FileList = new TList();
   for(int i = 0; i < nrfiles; i++)
   {
     sprintf(ctemp, "%s/results/hsimple%d.root", rootdir, i+1);
     FileList->Add(TFile::Open(ctemp));   
   }

   CheckKeys(Target, FileList, nrkeys, nrevts);
   for(int i = 0; i < nrkeys.size(); i++)
      printf("Keys in file %d = %d\n", i, nrkeys[i]);
   for(int i = 0; i < nrevts.size(); i++)
      printf("Events in key %d = %d\n", i, nrevts[i]);

   if(nrfiles > 1)
      MergeRootfile(Target, FileList, nrkeys, nrevts, filenames, 0);
   else
      MergeRootfile(Target, FileList, nrkeys, nrevts, filenames, 1);

   Target->Close();

   for(int i = 0; i < nrfiles; i++)
   {
     sprintf(ctemp, "rm -v %s/results/hsimple%d.root", rootdir, i+1);
     system(ctemp);
   }

   delete FileList;
}

void hadd(int nrfiles, string *files, string *outname)
{
   cout << "Number of files: " << nrfiles << endl;
   cout << "Saving into file: " << *outname << endl;

   TList *FileList;
   TFile *Target;

   char ctemp[1024];
   vector<int> nrkeys;
   vector<int> nrevts;
   vector<string> filenames;

   // Prepare the files to be merged
   sprintf(ctemp, "%s/results/hsimple1.root", rootdir);
   if(gSystem->AccessPathName(ctemp))
   {
     for(int i = 0; i < nrfiles; i++)
     {
       sprintf(ctemp, "%s/results/hsimple%d.root", rootdir, i+1);
       gSystem->CopyFile(files[i].c_str(), ctemp);
       filenames.push_back(files[i]);
     }
   }

   Target = TFile::Open(outname->c_str(), "RECREATE");

   FileList = new TList();
   for(int i = 0; i < nrfiles; i++)
   {
     sprintf(ctemp, "%s/results/hsimple%d.root", rootdir, i+1);
     FileList->Add(TFile::Open(ctemp));   
   }

   CheckKeys(Target, FileList, nrkeys, nrevts);
   for(int i = 0; i < nrkeys.size(); i++)
      printf("Keys in file %d = %d\n", i, nrkeys[i]);
   for(int i = 0; i < nrevts.size(); i++)
      printf("Events in key %d = %d\n", i, nrevts[i]);

   if(nrfiles > 1)
      CombineRootfile(Target, FileList, nrkeys, nrevts, filenames);

   Target->Close();

   for(int i = 0; i < nrfiles; i++)
   {
     sprintf(ctemp, "rm -v %s/results/hsimple%d.root", rootdir, i+1);
     system(ctemp);
   }

   delete FileList;
}
