#if _STANDALONE_ == 0
   #include "adst_mva.h"
   #include "separate_functions.h"
#endif

// TODO

// Functions to hold rules for saving observables ------------------------------
/*
     Return values determine if the selected observable was valid or not -> these values will set the rewrite code
*/
int AdstMva::SetSimObservables(Observables **cursig)
{
   cout << "# SetSimObservables     #: " << "Entering function AdstMva::SetSimObservables()" << endl;
   SetBinary(0, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values
   for(int j = 0; j < 3; j++)
   {
      cursig[j]->SetValue("nrmu", GetNrMuons(j));
   }
   return 0;
}

void AdstMva::ZeroSimObservables(Observables **cursig)
{
   cout << "# ZeroSimObservables     #: " << "Entering function AdstMva::ZeroSimObservables()" << endl;
   SetBinary(0, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values
   for(int j = 0; j < 3; j++)
   {
      cursig[j]->SetValue("nrmu", -1.);
   }
}

int AdstMva::SetSdObservables(Observables **cursig)
{
   cout << "# SetSdObservables      #: " << "Entering function AdstMva::SetSdObservables()" << endl;
   SetBinary(1, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values
   for(int j = 0; j < 3; j++)
   {
      cursig[j]->SetValue("nrstations", GetNrStations(j));
      cursig[j]->SetValue("shwsize", GetShowerSize(j));
      cursig[j]->SetValue("deltas38", 0.);
      cursig[j]->SetValue("energySD", GetSdEnergy(j));
      cursig[j]->SetValue("ldfbeta", GetBeta(j));
      cursig[j]->SetValue("curvature", GetCurvature(j));
      cursig[j]->SetValue("risetime", GetRisetime(j, false));
      InitRisetimeVariables();
      CalculateRisetime();
      cursig[j]->SetValue("risetimerecalc", GetRisetime(j, true));
      cursig[j]->SetValue("deltarisetime", 0.);
      cursig[j]->SetValue("zenithSD", GetSdZenith(j));
      cursig[j]->SetValue("azimuthSD", GetSdAzimuth(j));
      cursig[j]->SetValue("latitudeSD", GetSdLatitude(j));
      cursig[j]->SetValue("longitudeSD", GetSdLongitude(j));
      CalculateAoP();
      cursig[j]->SetValue("aop", GetAoP(j));
   }
   return 0;
}

void AdstMva::ZeroSdObservables(Observables **cursig)
{
   cout << "# ZeroSdObservables      #: " << "Entering function AdstMva::ZeroSdObservables()" << endl;
   SetBinary(1, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values
   for(int j = 0; j < 3; j++)
   {
      cursig[j]->SetValue("nrstations", 0.);
      cursig[j]->SetValue("shwsize", -1.);
      cursig[j]->SetValue("deltas38", -1.);
      cursig[j]->SetValue("energySD", -1.);
      cursig[j]->SetValue("ldfbeta", -1.);
      cursig[j]->SetValue("curvature", -1.);
      cursig[j]->SetValue("risetime", -1.);
      cursig[j]->SetValue("risetimerecalc", -1.);
      cursig[j]->SetValue("deltarisetime", -1.);
      cursig[j]->SetValue("zenithSD", -1.);
      cursig[j]->SetValue("azimuthSD", -1.);
      cursig[j]->SetValue("latitudeSD", -1.);
      cursig[j]->SetValue("longitudeSD", -1.);
      cursig[j]->SetValue("aop", -1.);
   }
}

int AdstMva::SetFdObservables(Observables **cursig)
{
   cout << "# SetFdObservables      #: " << "# Entering function AdstMva::SetFdObservables()" << endl;
   int *count = new int;
   double *dtemp = new double[3];
   bool *fdobservable = new bool;
   *fdobservable = true;
   int *realnr = new int;
   string *stemp = new string;
   isheco = false;

   // Zero all variables used for combining stereo events
   wQuantity = new double[2*((*cursig)->Count())];
   quantitySum = new double[2*((*cursig)->Count())];
   wQuantitySum = new double[2*((*cursig)->Count())];
   quantitySumErr = new double[2*((*cursig)->Count())];
   for(int i = 0; i < 2*((*cursig)->Count()); i++)
   {
      wQuantity[i] = 0;
      quantitySum[i] = 0;
      quantitySumErr[i] = 0;
      wQuantitySum[i] = 0;
   }

   SetBinary(2, 1, &rewritecode);

   // Check if we have HeCo and use it instead of CO or HE, and evaluate if nreyes was defined correctly
   *realnr = 0;
   for(int j = 0; j < nreyes; j++)
   {
      if(acteyes[j].GetEyeId() == 6)
      {
         isheco = true;
	 cout << "# SetFdObservables      #: " << "   Selected eye is HeCo." << endl;
      }

      if(acteyes[j].IsHybridEvent())
         (*realnr)++;
   }

   if(*realnr != nreyes)
      cout << "Nr. eyes vs. real nr. = " << nreyes << "/" << *realnr << endl;

   // Go over all observables and all eyes, and calculate the combined value from active eyes
   for(int j = 0; j < (*cursig)->Count(); j++)
   {
      *count = 0;
      for(int i = 0; i < nreyes; i++)
      {
         cout << "i = " << i << ", nreyes = " << nreyes << ", eyeID-1 = " << acteyes[i].GetEyeId()-1 << endl;

         if(isheco && (acteyes[i].GetEyeId()-1 == 3) && (nreyes > 1))
	 {
            i++;
	    cout << "new i = " << i << ", nreyes = " << nreyes << ", eyeID-1 = " << acteyes[i].GetEyeId()-1 << endl;
	 }
         else if(isheco && (acteyes[i].GetEyeId()-1 == 3) && (nreyes > 1))
	 {
            i++;
	    cout << "new i = " << i << ", nreyes = " << nreyes << ", eyeID-1 = " << acteyes[i].GetEyeId()-1 << endl;
	 }

//         if(DBGSIG > 1)
            cout << "Rewritting eye " << i+1 << "/" << nreyes << endl;

         CalculateShowerFoot(i);

	 for(int k = 0; k < 3; k++)
	 {
            if(j == 0)
               dtemp[k] = GetXmax(i, k);
            else if(j == 1)
               dtemp[k] = GetX0(i, k);
            else if(j == 2)
               dtemp[k] = GetLambda(i, k);
            else if(j == 3)
               dtemp[k] = GetFdEnergy(i, k);
            else if(j == 4)
               dtemp[k] = GetFdZenith(i, k);
            else if(j == 5)
               dtemp[k] = GetFdAzimuth(i, k);
            else if(j == 6)
               dtemp[k] = GetFdLatitude(i, k);
            else if(j == 7)
               dtemp[k] = GetFdLongitude(i, k);
            else if(j == 8)
               dtemp[k] = GetShowerFoot(i, k);
	    else
	    {
               *fdobservable = false;
               break;
	    }
	 }

         if(!(*fdobservable))
            break;

	 cout << "Input values: " << dtemp[0] << ", " << dtemp[1] << ", " << dtemp[2] << endl;

         wQuantity[2*j]   = 1./TMath::Power(dtemp[1],2);
         wQuantity[2*j+1] = 1./TMath::Power(dtemp[2],2);
         quantitySum[2*j]   += dtemp[0]*wQuantity[2*j];
         quantitySum[2*j+1] += dtemp[0]*wQuantity[2*j+1];
         wQuantitySum[2*j]   += wQuantity[2*j];
         wQuantitySum[2*j+1] += wQuantity[2*j+1];

         (*count)++;
      }

      if(*fdobservable)
      {
         if(j == 0)
            *stemp = "xmax";
         else if(j == 1)
            *stemp = "x0";
         else if(j == 2)
            *stemp = "lambda";
         else if(j == 3)
            *stemp = "energyFD";
         else if(j == 4)
            *stemp = "zenithFD";
         else if(j == 5)
            *stemp = "azimuthFD";
         else if(j == 6)
            *stemp = "latitudeFD";
         else if(j == 7)
            *stemp = "longitudeFD";
         else if(j == 8)
            *stemp = "shfoot";

         if(nreyes > 0)
	 {
            quantitySumErr[2*j]   = TMath::Sqrt(1./wQuantitySum[2*j]);
            quantitySumErr[2*j+1] = TMath::Sqrt(1./wQuantitySum[2*j+1]);
            quantitySum[2*j]   = quantitySum[2*j]/wQuantitySum[2*j];
            quantitySum[2*j+1] = quantitySum[2*j+1]/wQuantitySum[2*j+1];

	    if(quantitySum[2*j] != quantitySum[2*j+1])
	    {
               quantitySum[2*j] = (quantitySum[2*j] + quantitySum[2*j+1])/2.;
               quantitySumErr[2*j] += TMath::Abs(quantitySum[2*j+1] - quantitySum[2*j]);
               quantitySumErr[2*j+1] += TMath::Abs(quantitySum[2*j+1] - quantitySum[2*j]);
	    }

            cout << "Writing out FD observable: qSum = " << quantitySum[2*j] << ", qSumErr = " << quantitySumErr[2*j] << ", " << quantitySumErr[2*j+1] << endl;

	    cursig[0]->SetValue(*stemp, (float)quantitySum[2*j]);
	    cursig[1]->SetValue(*stemp, (float)quantitySumErr[2*j]);
	    cursig[2]->SetValue(*stemp, (float)quantitySumErr[2*j+1]);
	 }
	 else
	 {

            cout << "FD event has no valid eyes!" << endl;

	    cursig[0]->SetValue(*stemp, -1.);
	    cursig[1]->SetValue(*stemp, -1.);
	    cursig[2]->SetValue(*stemp, -1.);
	 }
      }
      else
         break;
   }

   delete[] wQuantity;
   delete[] quantitySum;
   delete[] quantitySumErr;
   delete[] wQuantitySum;

   delete fdobservable;
   delete realnr;
   delete stemp;
   delete[] dtemp;
   delete count;
   return 0;
}

void AdstMva::ZeroFdObservables(Observables **cursig)
{
   cout << "# ZeroFdObservables      #: " << "Entering function AdstMva::ZeroFdObservables()" << endl;
   SetBinary(1, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values
   for(int j = 0; j < 3; j++)
   {
      cursig[j]->SetValue("xmax", -1.);
      cursig[j]->SetValue("x0", -1.);
      cursig[j]->SetValue("lambda", -1.);
      cursig[j]->SetValue("energyFD", -1.);
      cursig[j]->SetValue("zenithFD", -1.);
      cursig[j]->SetValue("azimuthFD", -1.);
      cursig[j]->SetValue("latitudeFD", -1.);
      cursig[j]->SetValue("longitudeFD", -1.);
      cursig[j]->SetValue("shfoot", -1.);
   }
}

float AdstMva::TotalSignalFromPMT(SdRecStation *station, bool silent)
{
   int *itemp = new int;
   float retVal;
   vector<float> *vemtrace = new vector<float>;

   *itemp = 0;
   retVal = 0;

   for(int iPMT = 1; iPMT <= 3; iPMT++)
   {
      *vemtrace = station->GetVEMTrace(iPMT);
      if(vemtrace->size() > 0)
      {
         retVal += station->GetPMTTraces((ETraceType)0, iPMT).GetVEMSignal();
	 (*itemp)++;
      }
   }

   retVal = retVal/(*itemp);
   if(!silent)
      cout << "Total signal vs. calculated signal from " << *itemp << " PMTs: " << station->GetTotalSignal() << " vs. " << retVal << endl;

   delete vemtrace;
   delete itemp;

   return retVal;
}

int AdstMva::SetStationValues()
{
   cout << "# SetStationValues      #: " << "# Entering function AdstMva::SetStationValues()" << endl;

   string *stemp = new string[3];
   int *itemp = new int[3];
   float *ftemp = new float[4];
   bool *btemp = new bool;

   int *start_bin = new int;
   int *stop_bin = new int;

   double *xp = new double;
   double *yp = new double;
   double *maxval = new double;

   double *byrange = new double[2];
   double *bzrange = new double[2];
   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;

   vector<float> *time = new vector<float>;
   vector<float> *vemtrace = new vector<float>;
   vector<float> *yvalue = new vector<float>;
   vector<float> *tempVect = new vector<float>;

   float *eventThetaRec = new float;
   float *secZenith = new float;
   float *alpha = new float;
   float *gamma = new float;
   float *g = new float;
   float *zeta = new float;

   for(int i = 0; i < 3; i++)
   {
      stationDistance[i].clear();
      stationRisetime[i].clear();
   }
   stationHSat.clear();

   tempVect->clear();

   actstations = fRecEvent->GetSDEvent().GetStationVector();

   ftemp[0] = 0;
   itemp[0] = 0;
   itemp[2] = 0;

   if(TMath::Log10((fRecEvent->GetSDEvent().GetSdRecShower()).GetEnergy()) > 19.6)
      limitTankDistance[1] = 2000.;
   else
      limitTankDistance[1] = 1400.;

   for(int i = 0; i < actstations.size(); i++)
   {
      // Use only the stations that are valid candidates
      if(actstations[i].IsCandidate())
      {
/*-----------------------------------*/
         cout << ">>> New station (" << actstations[i].GetId() << ")" << endl;	// DEBUG
         *start_bin = actstations[i].GetSignalStartSlot() - 4;
         *stop_bin = actstations[i].GetSignalEndSlot();
   
         if( (*start_bin >= *stop_bin) || (*start_bin < 0) || (*start_bin > 5000) )
            *start_bin = 0;

         ftemp[1] = 0;
         itemp[1] = 0;

	 // Check all PMTs
	 for(int iPMT = 1; iPMT <= 3; iPMT++)
	 {
            time->clear();
            yvalue->clear();
            vemtrace->clear();
   
            *yp = 0;
            *maxval = -1.e40;
   
            *vemtrace = actstations[i].GetVEMTrace(iPMT);
   
            //cout << "PMT " << iPMT << ": Number of points in the VEM trace: " << vemtrace->size() << " --------------------------------------------------------" << endl;	// DEBUG
   
            // Continue if there is a VEM trace
            if(vemtrace->size() > 0)
            {
               itemp[0]++;
               itemp[1]++;
   
               ftemp[2] = 0;
   
               // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
               for(int iVEM = 0; iVEM < vemtrace->size(); iVEM++)
               {
                  if( (iVEM >= *start_bin) && (iVEM <= *stop_bin) )
                  {
                     time->push_back((float)iVEM*25.);
   
                     *yp += vemtrace->at(iVEM);
                     if(*yp > *maxval)
                        *maxval = *yp;
                  
                     yvalue->push_back(*yp);
                     ftemp[2] += *yp;
                  }
               }
   
               //cout << "Number of points in the signal slot: " << yvalue.size() << ", Maxval: " << *maxval << endl;	// DEBUG
   
               if(ftemp[2] < 0)
               {
                  cout << "Rejected PMT " << iPMT << " in tank " << actstations[i].GetId() << ": Negative signal integral value = " << ftemp[2] << endl;
                  itemp[0]--;
                  itemp[1]--;
               }
               else
               {
                  for(int iy = 0; iy < yvalue->size(); iy++)
                  {
                     //cout << time->at(iy)/25. << "\t" << yvalue->at(iy)/(*maxval) << endl;	// DEBUG
   
                     if(yvalue->at(iy)/(*maxval) > 0.95)
                        break;
   
                     if(yvalue->at(iy)/(*maxval) <= 0.10)
                     {
                        byrange[0] = yvalue->at(iy)/(*maxval);
                        byrange[1] = yvalue->at(iy+1)/(*maxval);
   
                        *yp = 0.1;
                        //cout << "yp = " << *yp << ", byrange = " << byrange[0] << ", " << byrange[1] << ", time = " << time->at(iy) << ", " << time->at(iy+1) << endl;	// DEBUG
                           // Find the x value of point with y value = yp = 0.1, that lies on a line between two points
                        // y = k*x + a
                        //    k = (y2 - y1)/(x2 - x1)
                        //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                        // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                        *xp = ((time->at(iy+1) - time->at(iy))*(*yp - byrange[1]))/(byrange[1] - byrange[0]) + time->at(iy+1);
   
                        byrange[0] = *xp;
                        byrange[1] = *yp;
                     }
   
                     if(yvalue->at(iy)/(*maxval) <= 0.50)
                     {
                        bzrange[0] = yvalue->at(iy)/(*maxval);
                        bzrange[1] = yvalue->at(iy+1)/(*maxval);
   
                        *yp = 0.5;
                        // Find the x value of point with y value = yp = 0.5, that lies on a line between two points
                        // y = k*x + a
                        //    k = (y2 - y1)/(x2 - x1)
                        //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                        // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                        *xp = ((time->at(iy+1) - time->at(iy))*(*yp - bzrange[1]))/(bzrange[1] - bzrange[0]) + time->at(iy+1);
   
                        bzrange[0] = *xp;
                        bzrange[1] = *yp;
                     }
                  }
   
                  //cout << "Calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;	// DEBUG
   
                  ftemp[0] += bzrange[0] - byrange[0];
                  ftemp[1] += bzrange[0] - byrange[0];
               }
            }
	 }
   
         ftemp[1] = ftemp[1]/itemp[1];

	 // Some ADST files don't calculate the VEM for an event
         ftemp[3] = actstations[i].GetTotalSignal();
	 if(ftemp[3] == 0)
            ftemp[3] = TotalSignalFromPMT(&actstations[i], false);

         // Asymmetry correction
         *eventThetaRec = fRecEvent->GetSDEvent().GetSdRecShower().GetZenith();
         *secZenith = 1./cos(*eventThetaRec);
         *alpha = 96.73 + (*secZenith)*(-282.40 + (*secZenith)*(241.80 - 62.61*(*secZenith)));
         *gamma = -0.0009572 + (*secZenith)*(0.002068 + (*secZenith)*(-0.001362 + 0.0002861*(*secZenith)));
         *g = (*alpha) + (*gamma) * actstations[i].GetSPDistance()*actstations[i].GetSPDistance();
         *zeta = actstations[i].GetAzimuthSP();
   
         risemean = ftemp[1] - (*g)*cos(*zeta);
         riseerr = fRTWeights->Eval(actstations[i].GetSPDistance(), (*secZenith), ftemp[3]);

	 *btemp = true;

	 // Check if station signal is below the minimum VEM signal
         if( (ftemp[3] < minSignal) )
	 {
            *btemp = false;
            cout << "Rejected: Station signal " << ftemp[3] << " is below the minimum accepted (" << minSignal << " VEM)." << endl;
	 }

	 // Check if station is low gain saturated (making it impossible to use)
         if( (actstations[i].IsLowGainSaturated()) || (includeSaturated) )
	 {
            *btemp = false;
            cout << "Rejected: Station signal is low gain saturated." << endl;
	 }

	 // Check if station is outside the distance limits: (300m,1400m) for low energy and (300m,2000m) for high energy
         if( (actstations[i].GetSPDistance() < limitTankDistance[0]) || (actstations[i].GetSPDistance() > limitTankDistance[1]) )
	 {
            *btemp = false;
            cout << "Rejected: Station distance " << actstations[i].GetSPDistance() << " is outside the limits (" << limitTankDistance[0] << "m, " << limitTankDistance[1] << "m)." << endl;
	 }

	 // Check if asymmetry corrected station risetime is negative
         if(risemean < 0.)
	 {
            *btemp = false;
            cout << "Rejected: Asymmetry corrected station risetime is negative (" << risemean << ")." << endl;
	 }
   
         if(*btemp)
         {
            tempVect->push_back((double)actstations[i].IsHighGainSaturated());
            tempVect->push_back(actstations[i].GetSPDistance());
            tempVect->push_back(actstations[i].GetSPDistanceError());
            tempVect->push_back(risemean);
            tempVect->push_back(riseerr);
   
            itemp[2]++;
         }
      }
   }

   if(itemp[2] < minPoints)
   {
      cout << "Rejected: Only " << itemp[2] << " valid tanks." << endl;
      cout << "Event reconstruction has failed." << endl;

      delete eventThetaRec;
      delete secZenith;
      delete alpha;
      delete gamma;
      delete g;
      delete zeta;

      delete xp;
      delete yp;
      delete maxval;

      delete[] byrange;

      delete time;
      delete vemtrace;
      delete yvalue;
      delete tempVect;

      delete start_bin;
      delete stop_bin;

      delete[] stemp;
      delete[] itemp;
      delete[] ftemp;
      delete btemp;
      return -1;
   }
   else
   {
      for(int i = 0; i < itemp[2]; i++)
      {
         // Save if station is high gain saturated
         stationHSat.push_back((bool)tempVect->at(5*i));
         // Save distance of station to shower axis
         stationDistance[0].push_back(tempVect->at(5*i+1));
         stationDistance[1].push_back(tempVect->at(5*i+2));
         stationDistance[2].push_back(tempVect->at(5*i+2));
         // Save risetime for each station
         stationRisetime[0].push_back(tempVect->at(5*i+3));
         stationRisetime[1].push_back(tempVect->at(5*i+4));
         stationRisetime[2].push_back(tempVect->at(5*i+4));
      }
   }

   nractstations = stationHSat.size();
/*-----------------------------------*/

   delete eventThetaRec;
   delete secZenith;
   delete alpha;
   delete gamma;
   delete g;
   delete zeta;

   delete xp;
   delete yp;
   delete maxval;

   delete[] byrange;

   delete time;
   delete vemtrace;
   delete yvalue;
   delete tempVect;

   delete start_bin;
   delete stop_bin;

   delete[] stemp;
   delete[] itemp;
   delete[] ftemp;
   delete btemp;
   return 0;
}
// Functions to hold rules for saving observables ------------------------------

// User defined functions for getting observable values ------------------------
float AdstMva::GetXmax(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetXmax();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetXmaxError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetXmaxError();
   else
      return -1.0;
}

float AdstMva::GetX0(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetX0();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetX0Error();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetX0Error();
   else
      return -1.0;
}

float AdstMva::GetLambda(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetLambda();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetLambdaError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetLambdaError();
   else
      return -1.0;
}

float AdstMva::GetFdEnergy(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetEnergy();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetEnergyError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetEnergyError();
   else
      return -1.0;
}

float AdstMva::GetFdZenith(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetZenith();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetZenithError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetZenithError();
   else
      return -1.0;
}

float AdstMva::GetFdAzimuth(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetAzimuth();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetAzimuthError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetAzimuthError();
   else
      return -1.0;
}

float AdstMva::GetFdLatitude(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetGalacticLatitude();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetGalacticLatitudeError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetGalacticLatitudeError();
   else
      return -1.0;
}

float AdstMva::GetFdLongitude(int eye, int type)
{
   // Mean value
   if(type == 0)
      return (acteyes[eye].GetFdRecShower()).GetGalacticLongitude();
   // Mean - neg value
   else if(type == 1)
      return (acteyes[eye].GetFdRecShower()).GetGalacticLongitudeError();
   // Mean + pos value
   else if(type == 2)
      return (acteyes[eye].GetFdRecShower()).GetGalacticLongitudeError();
   else
      return -1.0;
}

void AdstMva::CalculateShowerFoot(int eye)
{
//   cout << "# Entering function AdstMva::CalculateShowerFoot()" << endl;
   double *x, *xerr;
   int *itemp;
   double *dtemp;

   x = new double;
   xerr = new double;
   itemp = new int;
   dtemp = new double[3];

   *x = 0;
   *xerr = 0;

   vector<double> *slantDepth = new vector<double>;
   vector<double> *profiledEdX = new vector<double>;
   vector<double> *profiledEdXerr = new vector<double>;

   *slantDepth = (acteyes[eye].GetFdRecShower()).GetDepth();
   *profiledEdX = (acteyes[eye].GetFdRecShower()).GetEnergyDeposit();
   *profiledEdXerr = (acteyes[eye].GetFdRecShower()).GetEnergyDepositError();

   vector<double> xfoot;
   vector<double> yfoot;
   vector<double> yerrfoot;

   xfoot.clear();
   yfoot.clear();
   yerrfoot.clear();

   for(int i = 0; i < slantDepth->size(); i++)
   {
      *x += profiledEdX->at(i);
      *xerr += profiledEdXerr->at(i);
   
      xfoot.push_back(slantDepth->at(i));
      yfoot.push_back(*x);
      yerrfoot.push_back(*xerr);
   }

   *itemp = 0;

   for(int i = 0; i < yfoot.size(); i++)
   {
      if( (yfoot[i] >= shfootlimit*(yfoot[yfoot.size()-1])) && (*itemp == 0) )
      {
	 if(xfoot[i] != xfoot[i-1])	// In case a value is duplicated in the xfoot vector
	 {
            *itemp = 1;

            // Find the x value of point with y value that is a fraction of the maximum, that lies on a line between two points
            // y = k*x + a
            //    k = (y2 - y1)/(x2 - x1)
            //    a = y2 - (y2 - y1)/(x2 - x1)*x2
            // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
            dtemp[0] = (xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]); // 1/k = (x2 - x1)/(y2 - y1)
            dtemp[1] = shfootlimit*(yfoot[yfoot.size()-1]) - yfoot[i]; // y - y2
            *x = (dtemp[0])*(dtemp[1]) + xfoot[i]; // x = (1/k)*(y - y2) + x2

            dtemp[0] = (xfoot[i] - xfoot[i-1])/(yfoot[i]+yerrfoot[i] - (yfoot[i-1]+yerrfoot[i-1])); // 1/kerr = (x2 - x1)/(y2err - y1err)
            dtemp[1] = (yfoot[i]+yerrfoot[i]) - (1/(dtemp[0]))*(xfoot[i]); // aerr = y2err - (y2err - y1err)/(x2 - x1)*x2
            dtemp[2] = (1/(dtemp[0]))*(*x) + (dtemp[1]); // yerr = kerr*x + aerr
            *xerr = (dtemp[2]) - (((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i]))); // Dy = yerr - y

            dtemp[0] = ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i])); // y = k*x + a

            shfootmean = *x;

            for(int j = i; ; j++)	// Calculate positive error
            {
               if( yfoot[j] >= (dtemp[0])+(*xerr) )
               {
                  dtemp[1] = TMath::Abs(xfoot[j] - shfootmean);
                  break;
               }
            }

            for(int j = i; ; j--)	// Calculate negative error
            {
               if(j == 0)
               {
                  dtemp[2] = TMath::Abs(xfoot[j] - shfootmean);
                  break;
               }

               if( yfoot[j] <= (dtemp[0])-(*xerr) )
               {
                  dtemp[2] = TMath::Abs(xfoot[j] - shfootmean);
                  break;
               }
            }

            shfootmin = dtemp[2];
            shfootmax = dtemp[1];
	 }
      }
   }

   delete slantDepth;
   delete profiledEdX;
   delete profiledEdXerr;

   delete x;
   delete xerr;
   delete itemp;
   delete[] dtemp;
}

float AdstMva::GetShowerFoot(int eye, int type)
{
   // Mean value
   if(type == 0)
      return shfootmean;
   // Mean - neg value
   else if(type == 1)
      return shfootmin;
   // Mean + pos value
   else if(type == 2)
      return shfootmax;
   else
      return -1.0;
}

float AdstMva::GetNrStations(int type)
{
   // No actual error values
   // Mean value
   if(type == 0)
      return nractstations;
   // Mean - neg value
   else if(type == 1)
      return 0.;
   // Mean + pos value
   else if(type == 2)
      return 0.;
   else
      return -1.0;
}

float AdstMva::GetShowerSize(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetShowerSize();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetShowerSizeError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetShowerSizeError();
   else
      return -1.0;
}

float AdstMva::GetSdEnergy(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetEnergy();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetEnergyError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetEnergyError();
   else
      return -1.0;
}

float AdstMva::GetBeta(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetBeta();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetBetaError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetBetaError();
   else
      return -1.0;
}

float AdstMva::GetCurvature(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetCurvature();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetCurvatureError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetCurvatureError();
   else
      return -1.0;
}

void AdstMva::InitRisetimeVariables()
{
   limitTankDistance[0] = 300.;
   limitTankDistance[1] = 1400.;
   minSignal = 5.0;
   includeSaturated = false;
   minPoints = 3;
   evalDistance = 1000.;
}

void AdstMva::RisetimeFunction(double zenith, double energy, TF1 *risetimeFit)
{
  //Karen Mora parameterisation:
  const double secTheta= 1./cos(zenith); 
  const double a_par1= -0.141152; 
  const double a_par2=0.0141074;   
  const double a_par3=1.25107;  
  const double a_par4=-0.405333;   
  
  const double b_par1=0.000904323;   
  const double b_par2=6.4291e-06; 
  const double b_par3=-1.09992; 
  const double b_par4=0.30987;  
  
  double alpha=(a_par1+a_par2*log10(energy))*exp(-0.5*pow((secTheta-a_par3)/a_par4, 2));
  double beta=(b_par1+b_par2*log10(energy))*(1+b_par3* secTheta+b_par4*pow(secTheta, 2));
  
  risetimeFit->SetParameter(0,alpha);
  risetimeFit->SetParameter(1,beta);
}

void AdstMva::CalculateRisetime()
{
//   cout << "# Entering function AdstMva::CalculateRisetime()" << endl;
   if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
   {
      SetBinary(1, 0, &rewritecode);
      risemean = -1;
      risemin = -1;
      risemax = -1;
      return;
   }

   risemean = 0;
   risemin = 0;
   risemax = 0;

   double *x;
   double *y;
   double *maxval;
   int *nrpoints;
   double *dtemp;
   int *itemp;
   bool *btemp;

   x = new double;
   y = new double;
   maxval = new double;
   nrpoints = new int;
   dtemp = new double[4];
   itemp = new int[2];
   btemp = new bool;

   int *start_bin, *stop_bin;

   start_bin = new int;
   stop_bin = new int;

   double *byrange = new double[2];
   double *bzrange = new double[2];
   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;

   vector<SdRecStation> stationVector = fRecEvent->GetSDEvent().GetStationVector();

//cout << "Number of triggered stations: " << stationVector.size() << ", Zenith angle = " << fRecEvent->GetSDEvent().GetSdRecShower().GetZenith() << ", Cos Zenith angle = " << fRecEvent->GetSDEvent().GetSdRecShower().GetCosZenith() << ", Energy = " << fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy() << endl;

   vector<float> *time = new vector<float>;
   vector<float> *vemtrace = new vector<float>;

   vector<float> *yvalue = new vector<float>;

   vector<double> *riseVect = new vector<double>;
   vector<double> *riseVectErr = new vector<double>;
   vector<double> *distVect = new vector<double>;
   vector<double> *distVectErr = new vector<double>;

   dtemp[0] = 0;
   itemp[0] = 0;

   if(TMath::Log10((fRecEvent->GetSDEvent().GetSdRecShower()).GetEnergy()) > 19.6)
      limitTankDistance[1] = 2000.;
   else
      limitTankDistance[1] = 1400.;

   // Check all stations
   for(int i = 0; i < stationVector.size(); i++)
   {
      // Only use stations that are valid candidates
      if( (stationVector[i].IsCandidate()) /*&& (stationVector[i].GetSPDistance() < 1500.)*/ )
      {
         *start_bin = stationVector[i].GetSignalStartSlot() - 4;
         *stop_bin = stationVector[i].GetSignalEndSlot();

         if( (*start_bin >= *stop_bin) || (*start_bin < 0) || (*start_bin > 5000) )
            *start_bin = 0;

//cout << "Tank " << i << " is a candidate (" << *start_bin << "," << *stop_bin << ")." << endl;

         dtemp[1] = 0;
         itemp[1] = 0;

         // Check all PMTs
         for(int j = 1; j <= 3; j++)
         {
            if(time->size() != 0)
               time->erase(time->begin(), time->end());
            if(yvalue->size() != 0)
               yvalue->erase(yvalue->begin(), yvalue->end());
            *y = 0;
            *maxval = -1.e40;

            *vemtrace = stationVector[i].GetVEMTrace(j);
            *nrpoints = vemtrace->size();

//cout << "PMT " << j << ": Number of points in the VEM trace: " << *nrpoints << " --------------------------------------------------------" << endl;

            // Continue if there is a VEM trace
            if( *nrpoints > 0 )
            {
               itemp[0]++;
               itemp[1]++;

	       dtemp[2] = 0;

               // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
               for(int k = 0; k < *nrpoints; k++)
               {
                  if( (k >= *start_bin) && (k <= *stop_bin) )
                  {
                     time->push_back((float)k*25.);

                     *y += vemtrace->at(k);
                     if(*y > *maxval) *maxval = *y;
                  
                     yvalue->push_back(*y);
		     dtemp[2] += *y;
                  }
               }

//cout << "Number of points in the signal slot: " << yvalue->size() << ", Maxval: " << *maxval << endl;

	       if(dtemp[2] < 0)
               {
                  if(DBGSIG > 0)
                     cout << "Rejected PMT " << j << " in tank " << stationVector[i].GetId() << ": Negative signal integral value = " << dtemp[2] << endl;
		  itemp[0]--;
		  itemp[1]--;
               }
	       else
               {
                  for(int k = 0; k < yvalue->size(); k++)
                  {
//cout << time->at(k)/25. << "\t" << yvalue->at(k)/(*maxval) << endl;

                     if(yvalue->at(k)/(*maxval) > 0.95)
                        break;

                     if(yvalue->at(k)/(*maxval) <= 0.10)
                     {
                        byrange[0] = yvalue->at(k)/(*maxval);
                        byrange[1] = yvalue->at(k+1)/(*maxval);

                        *y = 0.1;
                        // Find the x value of point with y value = *y = 0.1, that lies on a line between two points
                        // y = k*x + a
                        //    k = (y2 - y1)/(x2 - x1)
                        //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                        // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                        *x = ((time->at(k+1) - time->at(k))*((*y) - byrange[1]))/(byrange[1] - byrange[0]) + time->at(k+1);

                        byrange[0] = *x;
                        byrange[1] = *y;
                     }

                     if(yvalue->at(k)/(*maxval) <= 0.50)
                     {
                        bzrange[0] = yvalue->at(k)/(*maxval);
                        bzrange[1] = yvalue->at(k+1)/(*maxval);

                        *y = 0.5;
                        // Find the x value of point with y value = *y = 0.5, that lies on a line between two points
                        // y = k*x + a
                        //    k = (y2 - y1)/(x2 - x1)
                        //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                        // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                        *x = ((time->at(k+1) - time->at(k))*((*y) - bzrange[1]))/(bzrange[1] - bzrange[0]) + time->at(k+1);

                        bzrange[0] = *x;
                        bzrange[1] = *y;
                     }
                  }

//cout << /*"Reconstructed risetime = " << inRisetime <<*/ ", calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;

                  dtemp[0] += bzrange[0] - byrange[0];
                  dtemp[1] += bzrange[0] - byrange[0];
	       }
            }
         }

         dtemp[1] = dtemp[1]/itemp[1];

	 // Some ADST files don't calculate the VEM for an event
         dtemp[3] = stationVector[i].GetTotalSignal();
	 if(dtemp[3] == 0)
            dtemp[3] = TotalSignalFromPMT(&stationVector[i], true);

         // Asymmetry correction
         double eventThetaRec = fRecEvent->GetSDEvent().GetSdRecShower().GetZenith();
         double secZenith = 1/cos(eventThetaRec);
         const double alpha = 96.73 + secZenith*(-282.40 + secZenith*(241.80 - 62.61*secZenith));
         const double gamma = -0.0009572 + secZenith*(0.002068 + secZenith*(-0.001362 + 0.0002861*secZenith));
         const double g = alpha + gamma * stationVector[i].GetSPDistance()*stationVector[i].GetSPDistance();
         const double zeta = stationVector[i].GetAzimuthSP();

         risemean = dtemp[1] - g*cos(zeta);
	 risemin = fRTWeights->Eval(stationVector[i].GetSPDistance(), secZenith, dtemp[3]);
//cout << "Asymmetry corrected risetime = " << risemean << " (+- " << risemin << "), from actual data = " << stationVector[i].GetAssymCorrRiseTime(eventThetaRec) << " (+- " << stationVector[i].GetAssymCorrRiseTimeError(eventThetaRec) << ")" << endl;

	 *btemp = true;

	 // Check if station signal is below the minimum VEM signal
         if( (dtemp[3] < minSignal) )
	 {
            *btemp = false;
//            cout << "Rejected: Station signal " << dtemp[3] << " is below the minimum accepted (" << minSignal << " VEM)." << endl;
	 }

	 // Check if station is low gain saturated (making it impossible to use)
         if( (stationVector[i].IsLowGainSaturated()) || (includeSaturated) )
	 {
            *btemp = false;
//            cout << "Rejected: Station signal is low gain saturated." << endl;
	 }

	 // Check if station is outside the distance limits: (300m,1400m) for low energy and (300m,2000m) for high energy
         if( (stationVector[i].GetSPDistance() < limitTankDistance[0]) || (stationVector[i].GetSPDistance() > limitTankDistance[1]) )
	 {
            *btemp = false;
//            cout << "Rejected: Station distance " << stationVector[i].GetSPDistance() << " is outside the limits (" << limitTankDistance[0] << "m, " << limitTankDistance[1] << "m)." << endl;
	 }

	 // Check if asymmetry corrected station risetime is negative
         if(risemean < 0.)
	 {
            *btemp = false;
//            cout << "Rejected: Asymmetry corrected station risetime is negative (" << risemean << ")." << endl;
	 }

         if(*btemp)
	 {
            riseVect->push_back(risemean);
//            riseVect->push_back(stationVector[i].GetAssymCorrRiseTime(eventThetaRec));
            distVect->push_back(stationVector[i].GetSPDistance()/1000.);
            riseVectErr->push_back(risemin);
//            riseVectErr->push_back(stationVector[i].GetAssymCorrRiseTimeError(eventThetaRec));
            distVectErr->push_back(stationVector[i].GetSPDistanceError()/1000.);
	 }
      }
   }

   if(riseVect->size() >= minPoints)
   {
      TGraphErrors riseGraph(riseVect->size(), &distVect->front(), &riseVect->front(), 0, &riseVectErr->front());
      TF1 risetimeFit("RisetimeFit", "40+[0]*x+[1]*x*x", limitTankDistance[0]/1000., limitTankDistance[1]/1000.);
      risetimeFit.SetParLimits(0, 0, 10000);
      risetimeFit.SetParLimits(1, 0, 10000);
      RisetimeFunction(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(), fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy(), &risetimeFit);
      int ret = riseGraph.Fit(&risetimeFit, "Q", "", limitTankDistance[0]/1000., limitTankDistance[1]/1000.);
      if(!ret)
      {
         risemean = risetimeFit.Eval(evalDistance/1000.);
         TVirtualFitter* const fitter = TVirtualFitter::GetFitter();
         const int nPar = fitter->GetNumberTotalParameters();
         const TMatrixD covar(nPar, nPar, fitter->GetCovarianceMatrix());

         risemin = sqrt(covar[1][1] + covar[0][0] + 2*covar[0][1]);
         risemax = sqrt(covar[1][1] + covar[0][0] + 2*covar[0][1]);
      }
      else
      {
         if(DBGSIG > 0)
            cout << "Rejected: Risetime fit unsuccessful." << endl;
	 risemean = -1;
	 risemin = -1;
	 risemax = -1;
      }
   }
   else
   {
      if(DBGSIG > 0)
         cout << "Rejected: Only " << riseVect->size() << " valid tanks." << endl;
      risemean = -1;
      risemin = -1;
      risemax = -1;
   }

   delete x;
   delete y;
   delete nrpoints;
   delete maxval;
   delete[] dtemp;
   delete[] itemp;
   delete btemp;

   delete start_bin;
   delete stop_bin;

   delete[] byrange;

   delete time;
   delete vemtrace;
   delete yvalue;
   delete riseVect;
   delete riseVectErr;
   delete distVect;
   delete distVectErr;
}

float AdstMva::GetRisetime(int type, bool recalc)
{
   // Recalculated risetime
   if(recalc)
   {
      // Mean value
      if(type == 0)
         return risemean;
      // Mean - neg value
      else if(type == 1)
         return risemin;
      // Mean + pos value
      else if(type == 2)
         return risemax;
      else
         return -1.0;
   }
   else
   {
      // Mean value
      if(type == 0)
         return sdrecshw->GetRiseTimeResults().GetRiseTime1000();
      // Mean - neg value
      else if(type == 1)
         return sdrecshw->GetRiseTimeResults().GetRiseTime1000Error();
      // Mean + pos value
      else if(type == 2)
         return sdrecshw->GetRiseTimeResults().GetRiseTime1000Error();
      else
         return -1.0;
   }
}

float AdstMva::GetSdZenith(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetZenith();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetZenithError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetZenithError();
   else
      return -1.0;
}

float AdstMva::GetSdAzimuth(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetAzimuth();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetAzimuthError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetAzimuthError();
   else
      return -1.0;
}

float AdstMva::GetSdLatitude(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetGalacticLatitude();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetGalacticLatitudeError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetGalacticLatitudeError();
   else
      return -1.0;
}

float AdstMva::GetSdLongitude(int type)
{
   // Mean value
   if(type == 0)
      return sdrecshw->GetGalacticLongitude();
   // Mean - neg value
   else if(type == 1)
      return sdrecshw->GetGalacticLongitudeError();
   // Mean + pos value
   else if(type == 2)
      return sdrecshw->GetGalacticLongitudeError();
   else
      return -1.0;
}

void AdstMva::CalculateAoP()
{
   double *chpeak, *charge, *peak;
   chpeak = new double[2];
   charge = new double[2];
   peak = new double[2];

   double *midAop;
   midAop = new double[2];

   int nrpmt, nrstat = 0;

   aopmean = 0;
   aopmin = 0;
   aopmax = 0;

   vector<SdRecStation> stationVector = fRecEvent->GetSDEvent().GetStationVector();

   // Check all stations
   for(int i = 0; i < stationVector.size(); i++)
   {
      charge[0] = 0;
      charge[1] = 0;
      peak[0] = 0;
      peak[1] = 0;
      chpeak[0] = 0;
      chpeak[1] = 0;
      nrpmt = 0;

      // Only use stations that are valid candidates
      if( stationVector[i].IsCandidate() )
      {
         for(int j = 1; j <= 3; j++)
	 {
	    charge[0] = stationVector[i].GetCharge(j);
	    charge[1] = stationVector[i].GetChargeError(j);
	    peak[0] = stationVector[i].GetPeak(j);
	    peak[1] = stationVector[i].GetPeakError(j);
	    
            if( (charge[0] != 0) && (peak[0] != 0) )
	    {
	       chpeak[0] += charge[0]/peak[0];
	       chpeak[1] += (charge[0]/peak[0])*TMath::Sqrt( (charge[1]*charge[1])/(charge[0]*charge[0]) + (peak[1]*peak[1])/(peak[0]*peak[0]) );
	       nrpmt++;
	    }
	 }

         if(nrpmt > 0)
	 {
	    midAop[0] = chpeak[0]/nrpmt;
	    midAop[1] = chpeak[1]/nrpmt;

	    aopmean += midAop[0];
	    aopmin += midAop[1];
	    aopmax += midAop[1];

	    nrstat++;
	 }
      }
   }

   if(nrstat > 0)
   {
      aopmean = aopmean/nrstat;
      aopmin = aopmin/nrstat;
      aopmax = aopmax/nrstat;
   }

   delete[] chpeak;
   delete[] charge;
   delete[] peak;
   delete[] midAop;
}

float AdstMva::GetAoP(int type)
{
   // Mean value
   if(type == 0)
      return aopmean;
   // Mean - neg value
   else if(type == 1)
      return aopmin;
   // Mean + pos value
   else if(type == 2)
      return aopmax;
   else
      return -1.0;
}

float AdstMva::GetNrMuons(int type)
{
   // Mean value
   if(type == 0)
      return genshw->GetMuonNumber();
   // Mean - neg value
   else if(type == 1)
      return genshw->GetMuonNumber();
   // Mean + pos value
   else if(type == 2)
      return genshw->GetMuonNumber();
   else
      return -1.0;
}
// User defined functions for getting observable values ------------------------
