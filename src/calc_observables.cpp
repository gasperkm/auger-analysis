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

   // Go over a loop of mean, neg and pos values - also save values in all four observables, used for eyes
   for(int j = 0; j < 3; j++)
   {
      for(int i = 0; i < ALLEYES; i++)
      {
         cursig[j]->SetValue("nrmu", GetNrMuons(i, j), i);
      }
   }
   return 0;
}

int AdstMva::SetSdObservables(Observables **cursig)
{
   cout << "# SetSdObservables      #: " << "Entering function AdstMva::SetSdObservables()" << endl;
   SetBinary(1, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values - also save values in all four observables, used for eyes
   for(int j = 0; j < 3; j++)
   {
      for(int i = 0; i < ALLEYES; i++)
      {
         cursig[j]->SetValue("shwsize", GetShowerSize(i, j), i);
         cursig[j]->SetValue("energySD", GetSdEnergy(i, j), i);
         cursig[j]->SetValue("ldfbeta", GetBeta(i, j), i);
         cursig[j]->SetValue("curvature", GetCurvature(i, j), i);
         cursig[j]->SetValue("risetime", GetRisetime(i, j, false), i);
	 InitRisetimeVariables();
	 CalculateRisetime();
         cursig[j]->SetValue("risetimerecalc", GetRisetime(i, j, true), i);
         cursig[j]->SetValue("zenithSD", GetSdZenith(i, j), i);
         cursig[j]->SetValue("azimuthSD", GetSdAzimuth(i, j), i);
         cursig[j]->SetValue("latitudeSD", GetSdLatitude(i, j), i);
         cursig[j]->SetValue("longitudeSD", GetSdLongitude(i, j), i);
	 CalculateAoP();
         cursig[j]->SetValue("aop", GetAoP(i, j), i);
      }
   }
   return 0;
}

int AdstMva::SetFdObservables(Observables **cursig)
{
   cout << "# SetFdObservables      #: " << "# Entering function AdstMva::SetFdObservables()" << endl;
   int count;

   SetBinary(2, 1, &rewritecode);

   // Go over a loop of mean, neg and pos values and all eyes
   for(int j = 0; j < 3; j++)
   {
      count = 0;
      for(int i = 0; i < ALLEYES; i++)
      {
 	 if(DBGSIG > 1)
            cout << "Rewritting eye " << i << " (" << ALLEYES << ")" << endl;
         // For some reason sometimes extra eyes remain in vector, but number of eyes is correct && if eye is not active, set value to -1
	 if( (nreyes == count) || (acteyes[count].GetEyeId()-1 != i) )
	 {
 	    if(DBGSIG > 1)
	       cout << "Extra eye!" << endl;

            cursig[j]->SetValue("xmax", GetXmax(count, -1), i);
            cursig[j]->SetValue("x0", GetX0(count, -1), i);
            cursig[j]->SetValue("lambda", GetLambda(count, -1), i);
            cursig[j]->SetValue("energyFD", GetFdEnergy(count, -1), i);
            cursig[j]->SetValue("zenithFD", GetFdZenith(count, -1), i);
            cursig[j]->SetValue("azimuthFD", GetFdAzimuth(count, -1), i);
            cursig[j]->SetValue("latitudeFD", GetFdLatitude(count, -1), i);
            cursig[j]->SetValue("longitudeFD", GetFdLongitude(count, -1), i);
            cursig[j]->SetValue("shfoot", GetShowerFoot(count, -1), i);
	 }
         // If we get active eye, set the values
         else if(acteyes[count].GetEyeId()-1 == i)
         {
 	    if(DBGSIG > 1)
               cout << "Active eye!" << endl;

            cursig[j]->SetValue("xmax", GetXmax(count, j), i);
            cursig[j]->SetValue("x0", GetX0(count, j), i);
            cursig[j]->SetValue("lambda", GetLambda(count, j), i);
            cursig[j]->SetValue("energyFD", GetFdEnergy(count, j), i);
            cursig[j]->SetValue("zenithFD", GetFdZenith(count, j), i);
            cursig[j]->SetValue("azimuthFD", GetFdAzimuth(count, j), i);
            cursig[j]->SetValue("latitudeFD", GetFdLatitude(count, j), i);
            cursig[j]->SetValue("longitudeFD", GetFdLongitude(count, j), i);
            CalculateShowerFoot(count);
            cursig[j]->SetValue("shfoot", GetShowerFoot(count, j), i);
            count++;
         }
      }

      if(DBGSIG > 1)
         cout << "Finished value type" << endl;
   }

   return 0;
}
// Functions to hold rules for saving observables ------------------------------

// User defined functions for getting observable values ------------------------
double AdstMva::GetXmax(int eye, int type)
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

double AdstMva::GetX0(int eye, int type)
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

double AdstMva::GetLambda(int eye, int type)
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

double AdstMva::GetFdEnergy(int eye, int type)
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

double AdstMva::GetFdZenith(int eye, int type)
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

double AdstMva::GetFdAzimuth(int eye, int type)
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

double AdstMva::GetFdLatitude(int eye, int type)
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

double AdstMva::GetFdLongitude(int eye, int type)
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
   double *dtemp1, *dtemp2, *dtemp3;

   x = new double;
   xerr = new double;
   itemp = new int;
   dtemp1 = new double;
   dtemp2 = new double;
   dtemp3 = new double;

   *x = 0;
   *xerr = 0;

   vector<double> slantDepth;
   vector<double> profiledEdX;
   vector<double> profiledEdXerr;

   slantDepth = (acteyes[eye].GetFdRecShower()).GetDepth();
   profiledEdX = (acteyes[eye].GetFdRecShower()).GetEnergyDeposit();
   profiledEdXerr = (acteyes[eye].GetFdRecShower()).GetEnergyDepositError();

   vector<double> xfoot;
   vector<double> yfoot;
   vector<double> yerrfoot;

   if(!xfoot.empty()) xfoot.erase(xfoot.begin(),xfoot.end());
   if(!yfoot.empty()) yfoot.erase(yfoot.begin(),yfoot.end());
   if(!yerrfoot.empty()) yerrfoot.erase(yerrfoot.begin(),yerrfoot.end());

   for(int i = 0; i < slantDepth.size(); i++)
   {
      *x += profiledEdX[i];
      *xerr += profiledEdXerr[i];
   
      xfoot.push_back(slantDepth[i]);
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
            *dtemp1 = (xfoot[i] - xfoot[i-1])/(yfoot[i] - yfoot[i-1]); // 1/k = (x2 - x1)/(y2 - y1)
            *dtemp2 = shfootlimit*(yfoot[yfoot.size()-1]) - yfoot[i]; // y - y2
            *x = (*dtemp1)*(*dtemp2) + xfoot[i]; // x = (1/k)*(y - y2) + x2

            *dtemp1 = (xfoot[i] - xfoot[i-1])/(yfoot[i]+yerrfoot[i] - (yfoot[i-1]+yerrfoot[i-1])); // 1/kerr = (x2 - x1)/(y2err - y1err)
            *dtemp2 = (yfoot[i]+yerrfoot[i]) - (1/(*dtemp1))*(xfoot[i]); // aerr = y2err - (y2err - y1err)/(x2 - x1)*x2
            *dtemp3 = (1/(*dtemp1))*(*x) + (*dtemp2); // yerr = kerr*x + aerr
            *xerr = (*dtemp3) - (((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i]))); // Dy = yerr - y

            *dtemp1 = ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(*x) + (yfoot[i] - ((yfoot[i] - yfoot[i-1])/(xfoot[i] - xfoot[i-1]))*(xfoot[i])); // y = k*x + a

            shfootmean = *x;

            for(int j = i; ; j++)	// Calculate positive error
            {
               if( yfoot[j] >= (*dtemp1)+(*xerr) )
               {
                  *dtemp2 = TMath::Abs(xfoot[j] - shfootmean);
                  break;
               }
            }

            for(int j = i; ; j--)	// Calculate negative error
            {
               if(j == 0)
               {
                  *dtemp3 = TMath::Abs(xfoot[j] - shfootmean);
                  break;
               }

               if( yfoot[j] <= (*dtemp1)-(*xerr) )
               {
                  *dtemp3 = TMath::Abs(xfoot[j] - shfootmean);
                  break;
               }
            }

            shfootmin = *dtemp3;
            shfootmax = *dtemp2;
	 }
      }
   }

   delete x;
   delete xerr;
   delete itemp;
   delete dtemp1;
   delete dtemp2;
   delete dtemp3;
}

double AdstMva::GetShowerFoot(int eye, int type)
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

double AdstMva::GetShowerSize(int eye, int type)
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

double AdstMva::GetSdEnergy(int eye, int type)
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

double AdstMva::GetBeta(int eye, int type)
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

double AdstMva::GetCurvature(int eye, int type)
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
   limitTankDistance[0] = 0.;
   limitTankDistance[1] = 1800.;
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

   x = new double;
   y = new double;
   maxval = new double;
   nrpoints = new int;
   dtemp = new double[2];
   itemp = new int[2];

   int start_bin, stop_bin;
   string stemp;

   double byrange[2];
   double bzrange[2];
   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;

   vector<SdRecStation> stationVector = fRecEvent->GetSDEvent().GetStationVector();
//   cout << "Number of triggered stations: " << stationVector.size() << ", Zenith angle = " << fRecEvent->GetSDEvent().GetSdRecShower().GetZenith() << ", Cos Zenith angle = " << fRecEvent->GetSDEvent().GetSdRecShower().GetCosZenith() << ", Energy = " << fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy() << endl;

   vector<float> time;
   vector<float> vemtrace;

   vector<float> yvalue;

   vector<double> riseVect;
   vector<double> riseVectErr;
   vector<double> distVect;
   vector<double> distVectErr;

   dtemp[0] = 0;
   itemp[0] = 0;

   // Check all stations
   for(int i = 0; i < stationVector.size(); i++)
   {
      // Only use stations that are valid candidates
      if( (stationVector[i].IsCandidate()) /*&& (stationVector[i].GetSPDistance() < 1500.)*/ )
      {
         start_bin = stationVector[i].GetSignalStartSlot() - 4;
         stop_bin = stationVector[i].GetSignalEndSlot();

         if( (start_bin >= stop_bin) || (start_bin < 0) || (start_bin > 5000) )
            start_bin = 0;

//         cout << "Tank " << i << " is a candidate (" << start_bin << "," << stop_bin << ")." << endl;

         dtemp[1] = 0;
         itemp[1] = 0;

         // Check all PMTs
         for(int j = 1; j <= 3; j++)
         {
            if(time.size() != 0)
               time.erase(time.begin(), time.end());
            if(yvalue.size() != 0)
               yvalue.erase(yvalue.begin(), yvalue.end());
            *y = 0;
            *maxval = -1.e40;

            vemtrace = stationVector[i].GetVEMTrace(j);
            *nrpoints = vemtrace.size();
//          cout << "PMT " << j << ": Number of points in the VEM trace: " << *nrpoints << " --------------------------------------------------------" << endl;

            // Continue if there is a VEM trace
            if( *nrpoints > 0 )
            {
               itemp[0]++;
               itemp[1]++;

               // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
               for(int k = 0; k < *nrpoints; k++)
               {
                  if( (k >= start_bin) && (k <= stop_bin) )
                  {
                     time.push_back((float)k*25.);

                     *y += vemtrace[k];
                     if(*y > *maxval) *maxval = *y;
                  
                     yvalue.push_back(*y);
                  }
               }

//             cout << "Number of points in the signal slot: " << yvalue.size() << endl;

               for(int k = 0; k < yvalue.size(); k++)
               {
//                cout << time[k]/25. << "\t" << yvalue[k]/(*maxval) << endl;

                  if(yvalue[k]/(*maxval) <= 0.10)
                  {
                     byrange[0] = yvalue[k]/(*maxval);
                     byrange[1] = yvalue[k+1]/(*maxval);

                     *y = 0.1;
                     // Find the x value of point with y value = *y = 0.1, that lies on a line between two points
                     // y = k*x + a
                     //    k = (y2 - y1)/(x2 - x1)
                     //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                     // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                     *x = ((time[k+1] - time[k])*((*y) - byrange[1]))/(byrange[1] - byrange[0]) + time[k+1];

                     byrange[0] = *x;
                     byrange[1] = *y;
                  }

                  if(yvalue[k]/(*maxval) <= 0.50)
                  {
                     bzrange[0] = yvalue[k]/(*maxval);
                     bzrange[1] = yvalue[k+1]/(*maxval);

                     *y = 0.5;
                     // Find the x value of point with y value = *y = 0.5, that lies on a line between two points
                     // y = k*x + a
                     //    k = (y2 - y1)/(x2 - x1)
                     //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                     // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                     *x = ((time[k+1] - time[k])*((*y) - bzrange[1]))/(bzrange[1] - bzrange[0]) + time[k+1];

                     bzrange[0] = *x;
                     bzrange[1] = *y;
                  }
               }

//             cout << "Reconstructed risetime = " << inRisetime << ", calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;
               dtemp[0] += bzrange[0] - byrange[0];
               dtemp[1] += bzrange[0] - byrange[0];
            }
         }

         dtemp[1] = dtemp[1]/itemp[1];
//         cout << "Tank " << stationVector[i].GetId() << ", " << stationVector[i].GetSPDistance() << " m: Calculated average risetime (for " << itemp[1] << " PMTs in the tank) = " << dtemp[1] << ", Total signal = " << stationVector[i].GetTotalSignal() << endl;

         // Asymmetry correction
         double eventThetaRec = fRecEvent->GetSDEvent().GetSdRecShower().GetZenith();
         double secZenith = 1/cos(eventThetaRec);
         const double alpha = 96.73 + secZenith*(-282.40 + secZenith*(241.80 - 62.61*secZenith));
         const double gamma = -0.0009572 + secZenith*(0.002068 + secZenith*(-0.001362 + 0.0002861*secZenith));
         const double g = alpha + gamma * stationVector[i].GetSPDistance()*stationVector[i].GetSPDistance();
         const double zeta = stationVector[i].GetAzimuthSP();

         risemean = dtemp[1] - g*cos(zeta);
	 risemin = fRTWeights->Eval(stationVector[i].GetSPDistance(), secZenith, stationVector[i].GetTotalSignal());
//         cout << "Asymmetry corrected risetime = " << risemean << " (+- " << risemin << "), from actual data = " << stationVector[i].GetAssymCorrRiseTime(eventThetaRec) << " (+- " << stationVector[i].GetAssymCorrRiseTimeError(eventThetaRec) << ")" << endl;

         if( (stationVector[i].GetTotalSignal() > minSignal) )
	 {
            if( (!stationVector[i].IsLowGainSaturated()) && (!includeSaturated) )
	    {
               if( (stationVector[i].GetSPDistance() >= limitTankDistance[0]) && (stationVector[i].GetSPDistance() <= limitTankDistance[1]) )
	       {
                  riseVect.push_back(risemean);
//                  riseVect.push_back(stationVector[i].GetAssymCorrRiseTime(eventThetaRec));
                  distVect.push_back(stationVector[i].GetSPDistance()/1000.);
                  riseVectErr.push_back(risemin);
//                  riseVectErr.push_back(stationVector[i].GetAssymCorrRiseTimeError(eventThetaRec));
                  distVectErr.push_back(stationVector[i].GetSPDistanceError()/1000.);
	       }
/*	       else
	       {
                  cout << "Rejected: Station distance " << stationVector[i].GetSPDistance() << " is outside the limits (" << limitTankDistance[0] << ", " << limitTankDistance[1] << ")." << endl;
	       }*/
	    }
/*	    else
	    {
               cout << "Rejected: Station signal is low gain saturated." << endl;
	    }*/
	 }
/*	 else
	 {
            cout << "Rejected: Station signal " << stationVector[i].GetTotalSignal() << " is below the minimum accepted (" << minSignal << ")." << endl;
	 }*/
      }
   }

   if(riseVect.size() >= minPoints)
   {
      TGraphErrors riseGraph(riseVect.size(), &distVect.front(), &riseVect.front(), 0, &riseVectErr.front());
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
	 risemean = -1;
	 risemin = -1;
	 risemax = -1;
      }

//      cout << /*"Reconstructed risetime = " << inRisetime <<*/ "Calculated average risetime (for " << itemp[0] << " PMTs in all tanks, " << riseVect.size() << " fitting points) = " << risemean << " (+- " << risemin << ")" << endl;
   }
   else
   {
      if(DBGSIG > 0)
         cout << "Rejected: Only " << riseVect.size() << " valid tanks." << endl;
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
}

double AdstMva::GetRisetime(int eye, int type, bool recalc)
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

double AdstMva::GetSdZenith(int eye, int type)
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

double AdstMva::GetSdAzimuth(int eye, int type)
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

double AdstMva::GetSdLatitude(int eye, int type)
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

double AdstMva::GetSdLongitude(int eye, int type)
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

double AdstMva::GetAoP(int eye, int type)
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

double AdstMva::GetNrMuons(int eye, int type)
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
