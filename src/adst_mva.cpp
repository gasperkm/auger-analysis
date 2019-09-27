#if _STANDALONE_ == 0
   #include "adst_mva.h"
   #include "separate_functions.h"
#endif

AdstMva::AdstMva()
{
#if ROOTADST == 1
   // ADST objects
   fRecEvent = new RecEvent();
   fDetGeo = new DetectorGeometry();
   sdrecshw = new SdRecShower();
   genshw = new GenShower();

   // Initialise number of eyes
   nreyes = 0;

   // Risetime calculation settings
   fRTWeights = new TFormula("RiseTimeWeights", "(80.0+(5.071e-7+6.48e-4*y-3.051e-4*y*y)*x*x)/z-16.46*y+36.16");

   // Initializing some variables
   shfootlimit = 0.1;
#endif
}

AdstMva::~AdstMva()
{
#if ROOTADST == 1
   delete fRTWeights;
   delete genshw;
   delete sdrecshw;
   delete fDetGeo;
   delete fRecEvent;
#endif
}

#if ROOTADST == 1
void AdstMva::PrepareOtherTrees(unsigned int nrfiles, int innr, vector<string> obs)
{
   Observables *othersig = new Observables(obs);
   Observables *othersig_neg = new Observables(obs);
   Observables *othersig_pos = new Observables(obs);

   string *stemp = new string[2];

   for(int i = 0; i < 3; i++)
   {
      stationDistance[i].clear();
      stationRisetime[i].clear();
   }
   stationHSat.clear();

   TTree *other_sig_tree;
#if OFFVER == 0
   stemp[0] = "TreeNewS" + ToString(innr+1);
   other_sig_tree = new TTree(stemp[0].c_str(), "Signal tree from new file.");
#elif OFFVER == 1
   stemp[0] = "TreeOldS" + ToString(innr+1);
   other_sig_tree = new TTree(stemp[0].c_str(), "Signal tree from old file.");
#endif

   // Determine branches for all observables
   other_sig_tree->Branch("rewritecode", &rewritecode, "rewritecode/I");
   for(int i = 0; i < othersig->Count(); i++)
   { 
      stemp[0] = othersig->GetName(i);
      stemp[1] = othersig->GetName(i) + "/F";
      other_sig_tree->Branch(stemp[0].c_str(), &(othersig->obsstruct[i].value), stemp[1].c_str());
      stemp[0] = othersig->GetName(i) + "_neg";
      stemp[1] = othersig->GetName(i) + "_neg/F";
      other_sig_tree->Branch(stemp[0].c_str(), &(othersig_neg->obsstruct[i].value), stemp[1].c_str());
      stemp[0] = othersig->GetName(i) + "_pos";
      stemp[1] = othersig->GetName(i) + "_pos/F";
      other_sig_tree->Branch(stemp[0].c_str(), &(othersig_neg->obsstruct[i].value), stemp[1].c_str());
   }
   other_sig_tree->Branch("stationdistance", &stationDistance[0]);
   other_sig_tree->Branch("stationrisetime", &stationRisetime[0]);
   other_sig_tree->Branch("stationhighsat", &stationHSat);
   other_sig_tree->Branch("stationdistance_neg", &stationDistance[1]);
   other_sig_tree->Branch("stationrisetime_neg", &stationRisetime[1]);
   other_sig_tree->Branch("stationdistance_pos", &stationDistance[2]);
   other_sig_tree->Branch("stationrisetime_pos", &stationRisetime[2]);

   // Write an empty tree
   other_sig_tree->Write();

   delete other_sig_tree;
   delete othersig_neg;
   delete othersig_pos;
   delete othersig;
   delete[] stemp;
}

#if _STANDALONE_ == 0
int AdstMva::RewriteObservables(int nrfiles, int innr, Observables **sig, Observables **all, wxProgressDialog *progress)
#elif _STANDALONE_ == 1
int AdstMva::RewriteObservables(int nrfiles, int innr, Observables **sig, Observables **all)
#endif
{
   string *stemp;
   int *itemp;

   stemp = new string[4];
   itemp = new int[2];

   cout << "# RewriteObservables    #: " << "# New input file (" << inname << ") ---------------------------------" << endl;

   // Set the names for the tree (old are simulations, new are data events)
#if OFFVER == 0
   stemp[0] = "TreeOldS" + ToString(innr+1);
#elif OFFVER == 1
   stemp[0] = "TreeNewS" + ToString(innr+1);
#endif

   // Create a tree to hold signal observable values
   stemp[1] = "Signal tree from file " + inname + ".";
   sig_tree = new TTree(stemp[0].c_str(), stemp[1].c_str());

   // Run over two loops to set branches for each observable -> loop in i is for all observables, loop in j is for mean, negative and positive values
   for(int j = 0; j < 3; j++)
   {
      for(int i = 0; i < sig[j]->Count(); i++)
      {
         if(j == 0)
	 {
            stemp[2] = sig[j]->GetName(i);
            stemp[3] = sig[j]->GetName(i) + "/F";
	 }
	 else if(j == 1)
	 {
            stemp[2] = sig[j]->GetName(i) + "_neg";
            stemp[3] = sig[j]->GetName(i) + "_neg/F";
	 }
	 else if(j == 2)
	 {
            stemp[2] = sig[j]->GetName(i) + "_pos";
            stemp[3] = sig[j]->GetName(i) + "_pos/F";
	 }

         // Set branch address for the rewrite code (defines if rewrite has failed at some point)
	 if( (i == 0) && (j == 0) )
	 {
	    sig_tree->Branch("rewritecode", &rewritecode, "rewritecode/I");
	    all_tree->Branch("rewritecode", &rewritecode, "rewritecode/I");
	 }

         // Set branch addresses for signal + sum observables
         sig_tree->Branch(stemp[2].c_str(), &(sig[j]->obsstruct[i].value), stemp[3].c_str());
         all_tree->Branch(stemp[2].c_str(), &(all[j]->obsstruct[i].value), stemp[3].c_str());
      }

      // Add station distance and station risetime vectors
      if(j == 0)
      {
         sig_tree->Branch("stationdistance", &stationDistance[j]);
         all_tree->Branch("stationdistance", &stationDistance[j]);
         sig_tree->Branch("stationrisetime", &stationRisetime[j]);
         all_tree->Branch("stationrisetime", &stationRisetime[j]);
         sig_tree->Branch("stationhighsat", &stationHSat);
         all_tree->Branch("stationhighsat", &stationHSat);
      }
      else if(j == 1)
      {
         sig_tree->Branch("stationdistance_neg", &stationDistance[j]);
         all_tree->Branch("stationdistance_neg", &stationDistance[j]);
         sig_tree->Branch("stationrisetime_neg", &stationRisetime[j]);
         all_tree->Branch("stationrisetime_neg", &stationRisetime[j]);
      }
      else if(j == 2)
      {
         sig_tree->Branch("stationdistance_pos", &stationDistance[j]);
         all_tree->Branch("stationdistance_pos", &stationDistance[j]);
         sig_tree->Branch("stationrisetime_pos", &stationRisetime[j]);
         all_tree->Branch("stationrisetime_pos", &stationRisetime[j]);
      }
   }

   // Open and prepare the ADST files for reading
   fFile = new RecEventFile(inname.c_str(), RecEventFile::eRead);
   itemp[0] = fFile->GetNEvents();
   cout << "# RewriteObservables    #: " << "Number of events: " << itemp[0] << endl;
   if(itemp[0] == 0)
   {
#if _STANDALONE_ == 0
      progress->Update(100);
#endif
      delete[] stemp;
      delete[] itemp;
      return -1;
   }
   fFile->SetBuffers(&fRecEvent);
   fFile->ReadDetectorGeometry(*fDetGeo);

#if _STANDALONE_ == 0
   progress->SetRange(itemp[0]);	// Set progress maximum for the progress dialog
#endif

   // Prepare vector to determine if reconstruction has failed at some stage
   vector<int> *recfail = new vector<int>;
   for(int i = 0; i < 10; i++)
      recfail->push_back(0);

   // Go over all events in the ADST file and write them out to the output file
   for(int j = 0; j < itemp[0]; j++)
   {
      // Update the progress bar
      if(itemp[0] < 20)
      {
#if _STANDALONE_ == 0
         progress->Update(j);
#elif _STANDALONE_ == 1
	 cerr << "Currently at " << j << "/" << itemp[0] << endl;
#endif
      }
      else if(itemp[0] >= 50)
      {
         if(j%((int)(itemp[0]*0.05)) == 0)
	 {
#if _STANDALONE_ == 0
            progress->Update(j);
#elif _STANDALONE_ == 1
	    cerr << "Currently at " << j << "/" << itemp[0] << endl;
#endif
	 }
      }
      else
      {
         if(j%((int)(itemp[0]*0.1)) == 0)
	 {
#if _STANDALONE_ == 0
            progress->Update(j);
#elif _STANDALONE_ == 1
	    cerr << "Currently at " << j << "/" << itemp[0] << endl;
#endif
	 }
      }

      // goodrec variable determines if FD and SD reconstructions have failed
      goodrec = true;
      goodrecsim = true;
      goodrecsd = true;
      goodrecfd = true;
      rewritecode = 0;

      fFile->ReadEvent(j);

      cout << "# RewriteObservables    #: " << "# New event (" << j+1 << ", ID = " << fRecEvent->GetEventId() << ", Time = [" << fRecEvent->GetYYMMDD() << "," << fRecEvent->GetHHMMSS() << "]) ---------------------------------" << endl;

      if(DBGSIG > 1)
         cout << "# RewriteObservables    #: " << "rewritecode (init) = " << rewritecode << endl;

      // Prepare simulated events --------------------------------------------
      *genshw = fRecEvent->GetGenShower();
      // Check if we have a simulated event or not (in order to set the name of the output tree)
      stemp[2] = string(genshw->GetPrimaryName());
      ReplaceAllCharacters(&stemp[2], string(" "), string(""));
      stemp[3] = string(sig_tree->GetTitle());

      if( (stemp[2] != "Photon") && (stemp[2] != "Proton") && (stemp[2] != "Iron") && (stemp[2] != "Carbon") && (stemp[2] != "Helium") && (stemp[2] != "Silicon") && (stemp[2] != "Nitrogen") && (stemp[2] != "Oxygen") )
         stemp[2] = "Data";

      if(stemp[2].compare(stemp[3]) != 0)
      {
         sig_tree->SetObject(stemp[0].c_str(), stemp[2].c_str() );

         if(DBGSIG > 0)
            cout << "# RewriteObservables    #: " << "Renaming tree to " << stemp[2] << endl;
      }

      // (Muon number at ground level) - only if we have simulations!
      itemp[1] = SetSimObservables(sig);
      itemp[1] = SetSimObservables(all);
      if(DBGSIG > 1)
         cout << "# RewriteObservables    #: " << "rewritecode (Sim) = " << rewritecode << endl;

      // Prepare FD eye events ----------------------------------------------
      if(DBGSIG > 0)
         cout << "# RewriteObservables    #: " << "   Number of eyes = " << fRecEvent->GetNEyes() << endl;
      if(fRecEvent->GetNEyes() == 0)
      {
         // If event has no FD eyes, disregard the event
         if(DBGSIG > 0)
            cout << "# RewriteObservables    #: " << "   Error! No reconstructed eyes for this event." << endl;
         if(goodrecfd)
	    (recfail->at(0))++;
         goodrecfd = false;
      }
      else
      {
         // Save all active FD eyes
	 acteyes.clear();
         acteyes = fRecEvent->GetFDEvents();
	 nreyes = acteyes.size();
         if(DBGSIG > 0)
	    cout << "# RewriteObservables    #: " << "   Number of active eyes = " << nreyes << endl;

         // Check if the saved event is considered to be hybrid (SD+FD) or not
         for(int i = 0; i < nreyes; i++)
         {
            if(!acteyes[i].IsHybridEvent())
            {
               if(DBGSIG > 0)
                  cout << "# RewriteObservables    #: " << "   Error! Eye " << acteyes[i].GetEyeId() << " not a hybrid event." << endl;
               goodrecfd = false;
            }
            else
            {
               if(DBGSIG > 0)
                  cout << "# RewriteObservables    #: " << "   Event is a hybrid event." << endl;
               goodrecfd = true;
               break;
            }
         }

         if(!goodrecfd)
	    (recfail->at(1))++;

         if(goodrecfd)
	 {
            itemp[1] = SetFdObservables(sig);
            itemp[1] = SetFdObservables(all);
	 }

         if(DBGSIG > 1)
            cout << "# RewriteObservables    #: " << "rewritecode (FD) = " << rewritecode << endl;
      }

      // Prepare SD station events -------------------------------------------
      *sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();

      // Check if there are triggered SD stations
      if(!(fRecEvent->GetSDEvent().HasTriggeredStations()))
      {
         if(DBGSIG > 0)
            cout << "# RewriteObservables    #: " << "Error! No triggered stations in SD reconstruction." << endl;
         if(goodrecsd) (recfail->at(3))++;
         goodrecsd = false;
      }

      // Check if there are any SD stations in the event
      if(!(fRecEvent->GetSDEvent().HasStations()))
      {
         if(DBGSIG > 0)
            cout << "# RewriteObservables    #: " << "Error! No stations in SD reconstruction." << endl;
         if(goodrecsd) (recfail->at(4))++;
         goodrecsd = false;
      }

      // Check if SD stations have a VEM trace
      if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
      {
         if(DBGSIG > 0)
            cout << "# RewriteObservables    #: " << "Error! No VEM traces in SD tanks." << endl;
         if(goodrecsd) (recfail->at(5))++;
         goodrecsd = false;
      }

      InitRisetimeVariables();
      itemp[1] = SetStationValues();
      if(itemp[1] == -1)
         goodrecsd = false;
      itemp[1] = SetSdObservables(sig);
      itemp[1] = SetSdObservables(all);
      if(DBGSIG > 1)
         cout << "# RewriteObservables    #: " << "rewritecode (SD) = " << rewritecode << endl;

      goodrec = (goodrecsim && goodrecsd && goodrecfd);

      if(goodrec)
      {
         sig_tree->Fill();
	 all_tree->Fill();
      }
      else
      {
         cout << "# RewriteObservables    #: " << "rewrite with problems = " << rewritecode << endl;
	 if(!goodrecsim)
	 {
            ZeroSimObservables(sig);
            ZeroSimObservables(all);
	 }

	 if(!goodrecfd)
	 {
            ZeroFdObservables(sig);
            ZeroFdObservables(all);
	 }

	 if(!goodrecsd)
	 {
            ZeroSdObservables(sig);
            ZeroSdObservables(all);

	    for(int i = 0; i < 3; i++)
	    {
               stationDistance[i].clear();
               stationRisetime[i].clear();
	    }
	    stationHSat.clear();
	 }

         sig_tree->Fill();
	 all_tree->Fill();
      }
   }

#if _STANDALONE_ == 0
   progress->Update(itemp[0]);
#elif _STANDALONE_ == 1
   cerr << "Finished rewriting" << endl;
#endif

   cout << "# RewriteObservables    #: " << (recfail->at(0))+(recfail->at(1))+(recfail->at(2))+(recfail->at(3))+(recfail->at(4))+(recfail->at(5))+(recfail->at(6))+(recfail->at(7)) << " events have been removed:" << endl;
   cout << "# RewriteObservables    #: " << " - No reconstructed FD eyes:      " << recfail->at(0) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - Not a hybrid event:            " << recfail->at(1) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - No valid FD reconstructions:   " << recfail->at(2) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - No triggered SD stations:      " << recfail->at(3) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - No reconstructed SD stations:  " << recfail->at(4) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - No reconstructed VEM traces:   " << recfail->at(5) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - No reconstructed risetime:     " << recfail->at(6) << " events" << endl;
   cout << "# RewriteObservables    #: " << " - No Area-over-Peak calculation: " << recfail->at(7) << " events" << endl;

   sig_tree->Write();
   fFile->Close();

   delete recfail;
   delete fFile;
   delete sig_tree;
   delete[] stemp;
   delete[] itemp;

   return 0;
}
#endif // ROOTADST
