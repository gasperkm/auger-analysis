#define _STANDALONE_ 1
#include "workstation.h"
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "separate_functions.h"
#include "mva_methods.h"
#include "mva_result_read.h"
#include "root_style.h"
#include "primary_type.h"
#if OFFVER == 0
   #include "OfflineIncludeOld.h"
#elif OFFVER == 1
   #include "OfflineIncludeNew.h"
#endif

using namespace std;

class AdstFile
{
private:
   string *stemp;
   int *itemp;
   double *dtemp;

   RecEventFile *fFile;
   RecEvent *fRecEvent;
   DetectorGeometry *fDetGeo;
   SdRecShower *sdrecshw;
   vector<SdRecStation> stationVector;

   // Temporary variables and holders for number of stations and events
   int nrstations, nevents;
   bool goodrec;

   // Limits for risetime calculations
   double limitTankDistance[2];
   double minSignal;
   bool includeSaturated;
   int minPoints;

   // Vectors and variables needed for calculation of risetime
   vector<double> tempVect;
   vector<float> time;
   vector<float> vemtrace;
   vector<float> yvalue;
   int start_bin, stop_bin;
   double maxval;
   double xp, yp;
   double risemean, riseerr;
   double *byrange;
   double *bzrange;
   TFormula *fRTWeights;

   double eventThetaRec;
   double secZenith;
   double alpha;
   double gamma;
   double g;
   double zeta;

public:
   AdstFile();
   virtual ~AdstFile();

   void ReadAdstFile(string inname, double *energylimitFull, double *zenithlimitFull, vector<double> *sdidVect, vector<bool> *HGsat, vector<double> *distVect, vector<double> *riseVect, vector<double> *energy, vector<double> *zenith, vector<double> *shwsize, vector<int> *eventVect, int *allcount, int *allevts);

   double GetDistanceLimit(int type);
};

AdstFile::AdstFile()
{
   fRecEvent = new RecEvent();
   fDetGeo = new DetectorGeometry();
   sdrecshw = new SdRecShower();

   stemp = new string[3];
   itemp = new int[4];
   dtemp = new double[3];

   limitTankDistance[0] = 300.;
   limitTankDistance[1] = 1400.;
   minSignal = 5.0;
   includeSaturated = false;
   minPoints = 3;

   byrange = new double[2];
   bzrange = new double[2];
   fRTWeights = new TFormula("RiseTimeWeights", "(80.0+(5.071e-7+6.48e-4*y-3.051e-4*y*y)*x*x)/z-16.46*y+36.16");

   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;
}

AdstFile::~AdstFile()
{
   delete fRTWeights;
   delete[] byrange;
   delete[] bzrange;

   delete[] dtemp;
   delete[] itemp;
   delete[] stemp;

   delete sdrecshw;
   delete fDetGeo;
   delete fRecEvent;
}

double AdstFile::GetDistanceLimit(int type)
{
   return limitTankDistance[type];
}

void AdstFile::ReadAdstFile(string inname, double *energylimitFull, double *zenithlimitFull, vector<double> *sdidVect, vector<bool> *HGsat, vector<double> *distVect, vector<double> *riseVect, vector<double> *energy, vector<double> *zenith, vector<double> *shwsize, vector<int> *eventVect, int *allcount, int *allevts)
{
   cout << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;
   cerr << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;

   // Open and prepare the ADST files for reading
   fFile = new RecEventFile(inname.c_str(), RecEventFile::eRead);
   nevents = fFile->GetNEvents();
   cout << "Number of events: " << nevents << endl;

   fFile->SetBuffers(&fRecEvent);
   fFile->ReadDetectorGeometry(*fDetGeo);

   for(int j = 0; j < nevents; j++)
   {
      if(nevents < 20)
         cerr << "Currently at " << j << "/" << nevents << endl;
      else if(j%((int)(nevents*0.05)) == 0)
         cerr << "Currently at " << j << "/" << nevents << endl;
      else
      {
         if(j%((int)(nevents*0.1)) == 0)
            cerr << "Currently at " << j << "/" << nevents << endl;
      }
   
      fFile->ReadEvent(j);
      cout << "New event (" << j+1 << ", ID = " << fRecEvent->GetEventId() << ", Time = [" << fRecEvent->GetYYMMDD() << "," << fRecEvent->GetHHMMSS() << "], E = " << TMath::Log10(fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy()) << ", sec(theta) = " <<SecTheta(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(),false)  << ") -------" << endl;
      goodrec = true;
            
      // Prepare SD station events -------------------------------------------
      *sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();
  
      // Check if there are triggered SD stations
      if(!(fRecEvent->GetSDEvent().HasTriggeredStations()))
         goodrec = false;
      // Check if there are any SD stations in the event
      if(!(fRecEvent->GetSDEvent().HasStations()))
         goodrec = false;
      // Check if SD stations have a VEM trace
      if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
         goodrec = false;

      // Limit in energy
      if(energylimitFull[0] != -1)
      {
         if(TMath::Log10(sdrecshw->GetEnergy()) < energylimitFull[0])
            goodrec = false;
      }
      if(energylimitFull[1] != -1)
      {
         if(TMath::Log10(sdrecshw->GetEnergy()) > energylimitFull[1])
            goodrec = false;
      }

      // Full limit in zenith angle (will be further binned later)
      if(zenithlimitFull[0] != -1)
      {
         if(SecTheta(sdrecshw->GetZenith(),false) < zenithlimitFull[0])
            goodrec = false;
      }
      if(zenithlimitFull[1] != -1)
      {
         if(SecTheta(sdrecshw->GetZenith(),false) > zenithlimitFull[1])
            goodrec = false;
      }

      if(TMath::Log10(sdrecshw->GetEnergy()) > 19.6)
         limitTankDistance[1] = 2000.;
      else
         limitTankDistance[1] = 1400.;
  
      if(goodrec)
      {
         dtemp[0] = 0;
         itemp[0] = 0;

         // Loop over all triggered SD stations
         stationVector = fRecEvent->GetSDEvent().GetStationVector();
         itemp[3] = 0;
         tempVect.clear();
         for(int i = 0; i < stationVector.size(); i++)
         {
            // Only use stations that are valid candidates
            if(stationVector[i].IsCandidate())
            {
               cout << "New station (" << stationVector[i].GetId() << ")" << endl;	// DEBUG
               start_bin = stationVector[i].GetSignalStartSlot() - 4;
               stop_bin = stationVector[i].GetSignalEndSlot();
   
               if( (start_bin >= stop_bin) || (start_bin < 0) || (start_bin > 5000) )
                  start_bin = 0;
   
               dtemp[1] = 0;
               itemp[1] = 0;
   
               // Check all PMTs
               for(int iPMT = 1; iPMT <= 3; iPMT++)
               {
                  time.clear();
                  yvalue.clear();
                  vemtrace.clear();
   
                  yp = 0;
                  maxval = -1.e40;
   
                  vemtrace = stationVector[i].GetVEMTrace(iPMT);
   
                  //cout << "PMT " << iPMT << ": Number of points in the VEM trace: " << vemtrace.size() << " --------------------------------------------------------" << endl;	// DEBUG
   
                  // Continue if there is a VEM trace
                  if(vemtrace.size() > 0)
                  {
                     itemp[0]++;
                     itemp[1]++;
   
                     dtemp[2] = 0;
   
                     // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
                     for(int iVEM = 0; iVEM < vemtrace.size(); iVEM++)
                     {
                        if( (iVEM >= start_bin) && (iVEM <= stop_bin) )
                        {
                           time.push_back((float)iVEM*25.);
   
                           yp += vemtrace[iVEM];
                           if(yp > maxval)
                              maxval = yp;
                        
                           yvalue.push_back(yp);
                           dtemp[2] += yp;
                        }
                     }
   
                     //cout << "Number of points in the signal slot: " << yvalue.size() << ", Maxval: " << maxval << endl;	// DEBUG
   
                     if(dtemp[2] < 0)
                     {
                        cout << "Rejected PMT " << iPMT << " in tank " << stationVector[i].GetId() << ": Negative signal integral value = " << dtemp[2] << endl;
                        itemp[0]--;
                        itemp[1]--;
                     }
                     else
                     {
                        for(int iy = 0; iy < yvalue.size(); iy++)
                        {
                           //cout << time[iy]/25. << "\t" << yvalue[iy]/maxval << endl;	// DEBUG
   
                           if(yvalue[iy]/maxval > 0.95)
                              break;
   
                           if(yvalue[iy]/maxval <= 0.10)
                           {
                              byrange[0] = yvalue[iy]/maxval;
                              byrange[1] = yvalue[iy+1]/maxval;
   
                              yp = 0.1;
                              //cout << "yp = " << yp << ", byrange = " << byrange[0] << ", " << byrange[1] << ", time = " << time[iy] << ", " << time[iy+1] << endl;	// DEBUG
                              // Find the x value of point with y value = yp = 0.1, that lies on a line between two points
                              // y = k*x + a
                              //    k = (y2 - y1)/(x2 - x1)
                              //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                              // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                              xp = ((time[iy+1] - time[iy])*(yp - byrange[1]))/(byrange[1] - byrange[0]) + time[iy+1];
   
                              byrange[0] = xp;
                              byrange[1] = yp;
                           }
   
                           if(yvalue[iy]/maxval <= 0.50)
                           {
                              bzrange[0] = yvalue[iy]/maxval;
                              bzrange[1] = yvalue[iy+1]/maxval;
   
                              yp = 0.5;
                              // Find the x value of point with y value = yp = 0.5, that lies on a line between two points
                              // y = k*x + a
                              //    k = (y2 - y1)/(x2 - x1)
                              //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                              // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                              xp = ((time[iy+1] - time[iy])*(yp - bzrange[1]))/(bzrange[1] - bzrange[0]) + time[iy+1];
   
                              bzrange[0] = xp;
                              bzrange[1] = yp;
                           }
                        }
   
                        //cout << "Calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;	// DEBUG
   
                        dtemp[0] += bzrange[0] - byrange[0];
                        dtemp[1] += bzrange[0] - byrange[0];
             	     }
                  }
               }
   
               dtemp[1] = dtemp[1]/itemp[1];
   
               //cout << "Station " << stationVector[i].GetId() << ", " << stationVector[i].GetSPDistance() << " m: Calculated average risetime (for " << itemp[1] << " PMTs in the tank) = " << dtemp[1] << ", Total signal = " << stationVector[i].GetTotalSignal() << endl;	// DEBUG
   
               // Asymmetry correction
               eventThetaRec = fRecEvent->GetSDEvent().GetSdRecShower().GetZenith();
               secZenith = 1./cos(eventThetaRec);
               alpha = 96.73 + secZenith*(-282.40 + secZenith*(241.80 - 62.61*secZenith));
               gamma = -0.0009572 + secZenith*(0.002068 + secZenith*(-0.001362 + 0.0002861*secZenith));
               g = alpha + gamma * stationVector[i].GetSPDistance()*stationVector[i].GetSPDistance();
               zeta = stationVector[i].GetAzimuthSP();
   
               risemean = dtemp[1] - g*cos(zeta);
               riseerr = fRTWeights->Eval(stationVector[i].GetSPDistance(), secZenith, stationVector[i].GetTotalSignal());
   
               if( (stationVector[i].GetTotalSignal() > minSignal) )
               {
                  if( (!stationVector[i].IsLowGainSaturated()) && (!includeSaturated) )
                  {
                     if( (stationVector[i].GetSPDistance() >= limitTankDistance[0]) && (stationVector[i].GetSPDistance() <= limitTankDistance[1]) )
                     {
                        tempVect.push_back(stationVector[i].GetId());
                        tempVect.push_back((double)stationVector[i].IsHighGainSaturated());
                        tempVect.push_back(stationVector[i].GetSPDistance());
                        tempVect.push_back(stationVector[i].GetSPDistanceError());
                        tempVect.push_back(risemean);
                        tempVect.push_back(riseerr);
                        tempVect.push_back(sdrecshw->GetEnergy());
                        tempVect.push_back(sdrecshw->GetEnergyError());
                        tempVect.push_back(sdrecshw->GetZenith());
                        tempVect.push_back(sdrecshw->GetZenithError());
                        tempVect.push_back(sdrecshw->GetShowerSize());
                        tempVect.push_back(sdrecshw->GetShowerSizeError());
   
          	        itemp[3]++;
          	     }
                     else
                        cout << "Rejected: Station distance " << stationVector[i].GetSPDistance() << " is outside the limits (" << limitTankDistance[0] << "m, " << limitTankDistance[1] << "m)." << endl;
                  }
          	  else
                     cout << "Rejected: Station signal is low gain saturated." << endl;
               }
               else
                  cout << "Rejected: Station signal " << stationVector[i].GetTotalSignal() << " is below the minimum accepted (" << minSignal << " VEM)." << endl;
            }
         }
   
         if(itemp[3] < minPoints)
         {
            goodrec = false;
            cout << "Rejected: Only " << itemp[3] << " valid tanks." << endl;
            cout << "Event reconstruction has failed." << endl;
         }
         else
         {
            for(int i = 0; i < itemp[3]; i++)
            {
               // Save SD station ID values
               sdidVect->push_back(tempVect[12*i]);
               // Save if station is high gain saturated
               HGsat->push_back((bool)tempVect[12*i+1]);
               // Save distance of station to shower axis
               distVect->push_back(tempVect[12*i+2]);
               distVect->push_back(tempVect[12*i+3]);
               // Save risetime for each station
               riseVect->push_back(tempVect[12*i+4]);
               riseVect->push_back(tempVect[12*i+5]);
               // Save event energy
               energy->push_back(tempVect[12*i+6]);
               energy->push_back(tempVect[12*i+7]);
               // Save event zenith angle
               zenith->push_back(tempVect[12*i+8]);
               zenith->push_back(tempVect[12*i+9]);
               // Save event S1000 value
               shwsize->push_back(tempVect[12*i+10]);
               shwsize->push_back(tempVect[12*i+11]);
               // Save the current event number
               eventVect->push_back(*allevts);

               (*allcount)++;

	       if(i == itemp[3]-1)
                  cout << "Risetime value from " << itemp[3] << " SD stations for event " << *allevts << " = " << tempVect[12*i+4] << " ± " << tempVect[12*i+5] << endl;
            }
   
            (*allevts)++;
         }
      }
      else
         cout << "Event reconstruction has failed." << endl;
   }
   
   fFile->Close();
   delete fFile;
}

void BinNaming(string *instring, double *energy, double *zenith)
{
   if( (energy[0] != -1) || (energy[1] != -1) )
      *instring += "_en_";
   if(energy[0] != -1)
      *instring += ToString(energy[0],2);
   if(energy[1] != -1)
      *instring += "-" + ToString(energy[1],2);
   if( (zenith[0] != -1) || (zenith[1] != -1) )
      *instring += "_zen_";
   if(zenith[0] != -1)
      *instring += ToString(zenith[0],2);
   if(zenith[1] != -1)
      *instring += "-" + ToString(zenith[1],2);
}

int ReadCustomBinning(string infile, vector<double> *outBins, double *maxRange)
{
   outBins->clear();

   ifstream *ifs = new ifstream;
   double *dtemp = new double[4];
   int itemp;

   ifs->open(infile.c_str(), ifstream::in);

   dtemp[2] = 1e+20;
   dtemp[3] = -1;

   if(ifs->is_open())
   {
      itemp = 0;
      cout << "Writing out custom binning:" << endl;
      while(1)
      {
         *ifs >> dtemp[0] >> dtemp[1];
   
         if(dtemp[0] < dtemp[2])
            dtemp[2] = dtemp[0];
         if(dtemp[1] > dtemp[3])
            dtemp[3] = dtemp[1];
   
         outBins->push_back(dtemp[0]);
         outBins->push_back(dtemp[1]);
   
         itemp++;
         ifs->ignore(1,' ');
         if(ifs->eof()) break;
      }
   
      cout << "Number of bins = " << itemp << ", Minimum bin value = " << dtemp[2] << ", Maximum bin value = " << dtemp[3] << endl;
      maxRange[0] = dtemp[2];
      maxRange[1] = dtemp[3];
   }
   
   ifs->close();

   delete ifs;
   delete[] dtemp;

   return itemp;
}

// Calculate SD fit parameters A, B and N from zenith angle
void CalculateSdFitParams(double *inpar, double *outpar, double zenith, double zenithErr)
{
   double *dtemp = new double[6];

   dtemp[0] = SecTheta(zenith,false);
   dtemp[1] = TMath::Sin(zenith);
   dtemp[2] = TMath::Cos(zenith);
   dtemp[3] = TMath::Tan(zenith);

   // A = a0 + a1*(sec(theta))^-4 = (inpar[0]+inpar[2]*TMath::Power(SecTheta(zenith,false),-4))
   outpar[0] = inpar[0]+inpar[2]*TMath::Power(dtemp[0],-4);
   // B = b0 + b1*(sec(theta))^-4 = (inpar[4]+inpar[6]*TMath::Power(SecTheta(zenith,false),-4))
   outpar[2] = inpar[4]+inpar[6]*TMath::Power(dtemp[0],-4);
   // N = n0 + n1*(sec(theta))^2 + n2*exp(sec(theta)) = (inpar[8]+inpar[10]*TMath::Power(SecTheta(zenith,false),2) + inpar[12]*TMath::Exp(SecTheta(zenith,false)))
   outpar[4] = inpar[8]+inpar[10]*TMath::Power(dtemp[0],2) + inpar[12]*TMath::Exp(dtemp[0]);

   // dA = sqrt( da0^2 + (da1/sec(theta)^4)^2 + (-4*dtheta*a1*sin(theta)*cos(theta)^3)^2 )
   //   a0 = inpar[0], da0 = inpar[1], a1 = inpar[2], da1 = inpar[3]
   dtemp[4] = inpar[1];
   dtemp[5] = TMath::Power(dtemp[4],2);
   cout << "  da0 = " << dtemp[4]; 
   dtemp[4] = inpar[3]*TMath::Power(dtemp[0],-4);
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", da1 = " << dtemp[4]; 
   dtemp[4] = -4*zenithErr*inpar[2]*dtemp[1]*TMath::Power(dtemp[2],3);
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", dtheta = " << dtemp[4] << endl; 
   outpar[1] = TMath::Sqrt(dtemp[5]);
   // dB = sqrt( db0^2 + (db1/sec(theta)^4)^2 + (-4*dtheta*b1*sin(theta)*cos(theta)^3)^2 )
   //   b0 = inpar[4], db0 = inpar[5], b1 = inpar[6], db1 = inpar[7]
   dtemp[4] = inpar[5];
   dtemp[5] = TMath::Power(dtemp[4],2);
   cout << "  db0 = " << dtemp[4]; 
   dtemp[4] = inpar[7]*TMath::Power(dtemp[0],-4);
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", db1 = " << dtemp[4]; 
   dtemp[4] = -4*zenithErr*inpar[6]*dtemp[1]*TMath::Power(dtemp[2],3);
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", dtheta = " << dtemp[4] << endl; 
   outpar[3] = TMath::Sqrt(dtemp[5]);
   // dN = sqrt( dn0^2 + (dn1*sec(theta)^2)^2 + (dn2*exp(sec(theta)))^2 + (dtheta*tan(theta)*sec(theta)*(2*n1*sec(theta)+n2*exp(sec(theta))))^2 )
   //   n0 = inpar[8], dn0 = inpar[9], n1 = inpar[10], dn1 = inpar[11], n2 = inpar[12], dn2 = inpar[13]
   dtemp[4] = inpar[9];
   dtemp[5] = TMath::Power(dtemp[4],2);
   cout << "  dn0 = " << dtemp[4]; 
   dtemp[4] = inpar[11]*TMath::Power(dtemp[0],2);
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", dn1 = " << dtemp[4]; 
   dtemp[4] = inpar[13]*TMath::Exp(dtemp[0]);
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", dn2 = " << dtemp[4]; 
   dtemp[4] = zenithErr*dtemp[3]*dtemp[0]*(2*inpar[10]*dtemp[0]+inpar[12]*TMath::Exp(dtemp[0]));
   dtemp[5] += TMath::Power(dtemp[4],2);
   cout << ", dtheta = " << dtemp[4] << endl; 
   outpar[5] = TMath::Sqrt(dtemp[5]);

   delete[] dtemp;
}

// Calculate benchmark functions
void CalculateBenchmark(double *tbench, double *fitparam, double dist, double distErr, bool sat)
{
   double *dtemp = new double[2];

   // The SD station is high-gain saturated
   //   tbench = 40 + sqrt(A^2 + B*r^2) - A
   if(sat)
   {
      tbench[0] = 40 + TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)) - fitparam[0];
      cout << "  tbench[0] = " << tbench[0];
      // (dA)^2 = (dA*(A/sqrt(A^2 + B*r^2) - 1))^2
      dtemp[0] = fitparam[1]*(fitparam[0]/TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)) - 1);
      dtemp[1] = TMath::Power(dtemp[0],2);
      cout << ", dA = " << dtemp[0];
      // (dB)^2 = (dB*r^2/(2*sqrt(A^2 + B*r^2)))^2
      dtemp[0] = fitparam[3]*TMath::Power(dist,2)/(2*TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)));
      dtemp[1] += TMath::Power(dtemp[0],2);
      cout << ", dB = " << dtemp[0];
      // (dr)^2 = (dr*B*r/(sqrt(A^2 + B*r^2)))^2
      dtemp[0] = distErr*fitparam[2]*dist/(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)));
      dtemp[1] += TMath::Power(dtemp[0],2);
      cout << ", dr = " << dtemp[0];
   }
   // The SD station has a valid high-gain trace
   //   tbench = 40 + N*(sqrt(A^2 + B*r^2) - A)
   else
   {
      tbench[0] = 40 + fitparam[4]*(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)) - fitparam[0]);
      cout << "  tbench[0] = " << tbench[0];
      // (dA)^2 = (N*dA*(A/sqrt(A^2 + B*r^2) - 1))^2
      dtemp[0] = fitparam[4]*fitparam[1]*(fitparam[0]/TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)) - 1);
      dtemp[1] = TMath::Power(dtemp[0],2);
      cout << ", dA = " << dtemp[0];
      // (dB)^2 = (N*dB*r^2/(2*sqrt(A^2 + B*r^2)))^2
      dtemp[0] = fitparam[4]*fitparam[3]*TMath::Power(dist,2)/(2*TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)));
      dtemp[1] += TMath::Power(dtemp[0],2);
      cout << ", dB = " << dtemp[0];
      // (dN)^2 = (dN*(sqrt(A^2 + B*r^2) - A))^2
      dtemp[0] = fitparam[5]*(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)) - fitparam[0]);
      dtemp[1] += TMath::Power(dtemp[0],2);
      cout << ", dN = " << dtemp[0];
      // (dr)^2 = (dr*B*N*r/(sqrt(A^2 + B*r^2)))^2
      dtemp[0] = distErr*fitparam[2]*fitparam[4]*dist/(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(dist,2)));
      dtemp[1] += TMath::Power(dtemp[0],2);
      cout << ", dr = " << dtemp[0];
   }

   tbench[1] = TMath::Sqrt(dtemp[1]);
   cout << ", tbench[1] = " << tbench[1] << endl;

   delete[] dtemp;
}

int main(int argc, char **argv)
{
   gSystem->Load("libTree.so");

   // Reading ADST files
/*   RecEventFile *fFile;
   RecEvent *fRecEvent = new RecEvent();
   DetectorGeometry *fDetGeo = new DetectorGeometry();
   SdRecShower *sdrecshw = new SdRecShower();
   vector<SdRecStation> stationVector;*/

   // Temporary variables and holders for number of stations and events
//   int nrstations, nevents;
   string *stemp = new string[3];
   int *itemp = new int[4];
   double *dtemp = new double[7];
   bool goodrec;

   // Limits for risetime calculations
/*   double limitTankDistance[2] = {300., 1400.};
   double minSignal = 5.0;
   bool includeSaturated = false;
   int minPoints = 3;*/

   // Vectors for distance, risetime and station ID
   vector<double> *distVect = new vector<double>;
   vector<double> *riseVect = new vector<double>;
   vector<double> *sdidVect = new vector<double>;
   vector<double> *energy = new vector<double>;
   vector<double> *zenith = new vector<double>;
   vector<double> *shwsize = new vector<double>;
   vector<int> *eventVect = new vector<int>;
   vector<bool> *HGsat = new vector<bool>;

   vector<double> distVectLow;
   vector<double> riseVectLow;
   vector<double> distVectHigh;
   vector<double> riseVectHigh;
   vector<double> energyVect;
   vector<double> deltaVect;
//   vector<double> tempVect;

   // Vectors and variables needed for calculation of risetime
/*   vector<float> time;
   vector<float> vemtrace;
   vector<float> yvalue;
   int start_bin, stop_bin;
   double maxval;
   double xp, yp;
   double risemean, riseerr;
   double *byrange = new double[2];
   double *bzrange = new double[2];
   byrange[0] = 1.e+40;
   byrange[1] = -1.e+40;
   bzrange[0] = 1.e+40;
   bzrange[1] = -1.e+40;
   TFormula *fRTWeights = new TFormula("RiseTimeWeights", "(80.0+(5.071e-7+6.48e-4*y-3.051e-4*y*y)*x*x)/z-16.46*y+36.16");*/

   // Graphing variables for risetimes
   RootStyle *mystyle;
   TCanvas *c1;
   TGraphErrors *allRise, *allRiseLowGain, *allRiseHighGain;
   TGraphErrors *allDelta;
   TH1D *riseErrDist, *riseErrDistLow, *riseErrDistHigh;
   TLine *line = new TLine();
   line->SetLineWidth(2);
   line->SetLineStyle(7);
   line->SetLineColor(28);
   TLine *line2 = new TLine();
   line2->SetLineWidth(2);
   line2->SetLineStyle(7);
   line2->SetLineColor(28);
   int allcount = 0;
   int allcountLow = 0;
   int allcountHigh = 0;
   int allevts = 0;
   double fitparamLow[4];
   double fitparamLowErr[4];
   double fitparamHigh[4];
   double fitparamHighErr[4];
   double fitparam[6];
   double sdfitparam[14];

   // Fitting parameters of benchmark functions from the SD analysis
   sdfitparam[0] = -72.;
   sdfitparam[1] = 10.;
   sdfitparam[2] = 410.;
   sdfitparam[3] = 30.;
   sdfitparam[4] = -0.049;
   sdfitparam[5] = 0.007;
   sdfitparam[6] = 0.36;
   sdfitparam[7] = 0.02;
   sdfitparam[8] = -0.07;
   sdfitparam[9] = 0.02;
   sdfitparam[10] = -1.14;
   sdfitparam[11] = 0.02;
   sdfitparam[12] = 0.84;
   sdfitparam[13] = 0.02;

   double tbench[2];

   // Graphing variables for A, B and N parameters
   TGraphErrors *parA, *parB, *parN;
   TF1 *parfit[3];
   TF1 *tempfunc;

   // Variables for setting energy and zenith angle limits
   double energylimit[2] = {-1,-1};
   double energylimitFull[2] = {-1,-1};
   int energynr = 12;
   vector<double> *energybins = new vector<double>;
   int energyref = -1;
   double zenithlimit[2] = {-1,-1};
   double zenithlimitFull[2] = {-1,-1};
   int zenithnr = 9;
   vector<double> *zenithbins = new vector<double>;
   int sdBenchmark = 1;

   stemp[0] = string(rootdir) + "/input/custom_energy_bins.txt";
   energynr = ReadCustomBinning(stemp[0], energybins, energylimitFull);
   stemp[0] = string(rootdir) + "/input/custom_zenith_bins.txt";
   zenithnr = ReadCustomBinning(stemp[0], zenithbins, zenithlimitFull);

   ofstream *ofs = new ofstream;
   ifstream *ifs = new ifstream;
   AdstFile *adfile = new AdstFile();

   string inname;

   if(argc > 1)
   {
      cerr << "Use fitted benchmark functions (0) or those defined in SD analysis (1): ";
      cin >> sdBenchmark;

      if(sdBenchmark == 1)
      {
         cerr << "Using benchmark functions defined in the SD analysis." << endl;
         cout << "Using benchmark functions defined in the SD analysis." << endl;
         for(int i = 0; i < energynr; i++)
	 {
            if(energybins->at(2*i) == 19.1)
               energyref = i;
	 }
      }
      else
      {
         cerr << "Using fitted benchmark functions from this analysis." << endl;
         cout << "Using fitted benchmark functions from this analysis." << endl;

         cerr << "Energy bins:" << endl;
         for(int i = 0; i < energynr; i++)
            cerr << "[" << i << "] = " << energybins->at(2*i) << ", " << energybins->at(2*i+1) << endl;
         cerr << "Select the energy bin to be used as reference (use number in square brackets): ";
         cin >> energyref;
      }

      cerr << "Energy bin used for reference = (" << energybins->at(2*energyref) << ", " << energybins->at(2*energyref+1) << ")" << endl;
      cout << "Energy bin used for reference = (" << energybins->at(2*energyref) << ", " << energybins->at(2*energyref+1) << ")" << endl;

      mystyle = new RootStyle();
      mystyle->SetBaseStyle();
      c1 = new TCanvas("c1","",1200,900);
      gStyle->SetEndErrorSize(3);

      stemp[0] = "mkdir -p ./plots";
      system(stemp[0].c_str());

      distVect->clear();
      riseVect->clear();
      sdidVect->clear();
      energy->clear();
      zenith->clear();
      shwsize->clear();
      eventVect->clear();
      HGsat->clear();

      for(int k = 0; k < argc-1; k++)
      {
         inname = string(argv[k+1]);
/*----------------------------------------------*/
	 adfile->ReadAdstFile(inname, energylimitFull, zenithlimitFull, sdidVect, HGsat, distVect, riseVect, energy, zenith, shwsize, eventVect, &allcount, &allevts);
/*----------------------------------------------*/
/*         cout << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;
         cerr << endl << "Opening file: " << inname << " -------------------------------------------------------" << endl;
         
	 // Open and prepare the ADST files for reading
         fFile = new RecEventFile(inname.c_str(), RecEventFile::eRead);
         nevents = fFile->GetNEvents();
         cout << "Number of events: " << nevents << endl;

         fFile->SetBuffers(&fRecEvent);
         fFile->ReadDetectorGeometry(*fDetGeo);

	 for(int j = 0; j < nevents; j++)
	 {
            if(nevents < 20)
	       cerr << "Currently at " << j << "/" << nevents << endl;
            else if(j%((int)(nevents*0.05)) == 0)
	       cerr << "Currently at " << j << "/" << nevents << endl;
            else
	    {
               if(j%((int)(nevents*0.1)) == 0)
	          cerr << "Currently at " << j << "/" << nevents << endl;
	    }

            fFile->ReadEvent(j);
            cout << "New event (" << j+1 << ", ID = " << fRecEvent->GetEventId() << ", Time = [" << fRecEvent->GetYYMMDD() << "," << fRecEvent->GetHHMMSS() << "], E = " << TMath::Log10(fRecEvent->GetSDEvent().GetSdRecShower().GetEnergy()) << ", sec(theta) = " <<SecTheta(fRecEvent->GetSDEvent().GetSdRecShower().GetZenith(),false)  << ") -------" << endl;
	    goodrec = true;
            
	    // Prepare SD station events -------------------------------------------
            *sdrecshw = fRecEvent->GetSDEvent().GetSdRecShower();

            // Check if there are triggered SD stations
            if(!(fRecEvent->GetSDEvent().HasTriggeredStations()))
               goodrec = false;
            // Check if there are any SD stations in the event
            if(!(fRecEvent->GetSDEvent().HasStations()))
               goodrec = false;
            // Check if SD stations have a VEM trace
            if(!(fRecEvent->GetSDEvent().HasVEMTraces()))
               goodrec = false;

	    // Limit in energy
            if(energylimit[0] != -1)
	    {
               if(TMath::Log10(sdrecshw->GetEnergy()) < energylimitFull[0])
                  goodrec = false;
	    }
            if(energylimit[1] != -1)
	    {
               if(TMath::Log10(sdrecshw->GetEnergy()) > energylimitFull[1])
                  goodrec = false;
	    }

	    // Full limit in zenith angle (will be further binned later)
            if(zenithlimitFull[0] != -1)
	    {
               if(SecTheta(sdrecshw->GetZenith(),false) < zenithlimitFull[0])
                  goodrec = false;
	    }
            if(zenithlimitFull[1] != -1)
	    {
               if(SecTheta(sdrecshw->GetZenith(),false) > zenithlimitFull[1])
                  goodrec = false;
	    }

	    if(TMath::Log10(sdrecshw->GetEnergy()) > 19.6)
               limitTankDistance[1] = 2000.;
	    else
               limitTankDistance[1] = 1400.;

	    if(goodrec)
	    {
	       dtemp[0] = 0;
               itemp[0] = 0;

	       // Loop over all triggered SD stations
	       stationVector = fRecEvent->GetSDEvent().GetStationVector();
	       itemp[3] = 0;
	       tempVect.clear();
               for(int i = 0; i < stationVector.size(); i++)
	       {
                  // Only use stations that are valid candidates
                  if(stationVector[i].IsCandidate())
	          {
cout << "New station (" << stationVector[i].GetId() << ")" << endl;
                     start_bin = stationVector[i].GetSignalStartSlot() - 4;
                     stop_bin = stationVector[i].GetSignalEndSlot();

                     if( (start_bin >= stop_bin) || (start_bin < 0) || (start_bin > 5000) )
                        start_bin = 0;

	             dtemp[1] = 0;
	             itemp[1] = 0;

                     // Check all PMTs
                     for(int iPMT = 1; iPMT <= 3; iPMT++)
                     {
                        time.clear();
	                yvalue.clear();
			vemtrace.clear();

                        yp = 0;
                        maxval = -1.e40;

                        vemtrace = stationVector[i].GetVEMTrace(iPMT);

//cout << "PMT " << iPMT << ": Number of points in the VEM trace: " << vemtrace.size() << " --------------------------------------------------------" << endl;

                        // Continue if there is a VEM trace
                        if( vemtrace.size() > 0 )
                        {
                           itemp[0]++;
                           itemp[1]++;

	                   dtemp[2] = 0;

                           // Prepare the time vector (each point is multiplied by 25 to get nanoseconds)
                           for(int iVEM = 0; iVEM < vemtrace.size(); iVEM++)
                           {
                              if( (iVEM >= start_bin) && (iVEM <= stop_bin) )
                              {
                                 time.push_back((float)iVEM*25.);

                                 yp += vemtrace[iVEM];
                                 if(yp > maxval) maxval = yp;
                              
                                 yvalue.push_back(yp);
	                         dtemp[2] += yp;
                              }
                           }

//cout << "Number of points in the signal slot: " << yvalue.size() << ", Maxval: " << maxval << endl;

	                   if(dtemp[2] < 0)
                           {
                              cout << "Rejected PMT " << iPMT << " in tank " << stationVector[i].GetId() << ": Negative signal integral value = " << dtemp[2] << endl;
	                      itemp[0]--;
	                      itemp[1]--;
                           }
	                   else
                           {
                              for(int iy = 0; iy < yvalue.size(); iy++)
                              {
//cout << time[iy]/25. << "\t" << yvalue[iy]/maxval << endl;

                                 if(yvalue[iy]/maxval > 0.95)
                                    break;

                                 if(yvalue[iy]/maxval <= 0.10)
                                 {
                                    byrange[0] = yvalue[iy]/maxval;
                                    byrange[1] = yvalue[iy+1]/maxval;

                                    yp = 0.1;
//cout << "yp = " << yp << ", byrange = " << byrange[0] << ", " << byrange[1] << ", time = " << time[iy] << ", " << time[iy+1] << endl;
                                    // Find the x value of point with y value = yp = 0.1, that lies on a line between two points
                                    // y = k*x + a
                                    //    k = (y2 - y1)/(x2 - x1)
                                    //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                                    // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                                    xp = ((time[iy+1] - time[iy])*(yp - byrange[1]))/(byrange[1] - byrange[0]) + time[iy+1];

                                    byrange[0] = xp;
                                    byrange[1] = yp;
                                 }

                                 if(yvalue[iy]/maxval <= 0.50)
                                 {
                                    bzrange[0] = yvalue[iy]/maxval;
                                    bzrange[1] = yvalue[iy+1]/maxval;

                                    yp = 0.5;
                                    // Find the x value of point with y value = yp = 0.5, that lies on a line between two points
                                    // y = k*x + a
                                    //    k = (y2 - y1)/(x2 - x1)
                                    //    a = y2 - (y2 - y1)/(x2 - x1)*x2
                                    // x = ((x2 - x1)/(y2 - y1))*(y - y2) + x2
                                    xp = ((time[iy+1] - time[iy])*(yp - bzrange[1]))/(bzrange[1] - bzrange[0]) + time[iy+1];

                                    bzrange[0] = xp;
                                    bzrange[1] = yp;
                                 }
                              }

//cout << "Calculated risetime (" << byrange[0]/25. << "," << bzrange[0]/25. << ") = " << bzrange[0] - byrange[0] << endl;

                              dtemp[0] += bzrange[0] - byrange[0];
                              dtemp[1] += bzrange[0] - byrange[0];
	           	   }
	                }
                     }

                     dtemp[1] = dtemp[1]/itemp[1];

//cout << "Station " << stationVector[i].GetId() << ", " << stationVector[i].GetSPDistance() << " m: Calculated average risetime (for " << itemp[1] << " PMTs in the tank) = " << dtemp[1] << ", Total signal = " << stationVector[i].GetTotalSignal() << endl;

                     // Asymmetry correction
                     double eventThetaRec = fRecEvent->GetSDEvent().GetSdRecShower().GetZenith();
                     double secZenith = 1/cos(eventThetaRec);
                     const double alpha = 96.73 + secZenith*(-282.40 + secZenith*(241.80 - 62.61*secZenith));
                     const double gamma = -0.0009572 + secZenith*(0.002068 + secZenith*(-0.001362 + 0.0002861*secZenith));
                     const double g = alpha + gamma * stationVector[i].GetSPDistance()*stationVector[i].GetSPDistance();
                     const double zeta = stationVector[i].GetAzimuthSP();

                     risemean = dtemp[1] - g*cos(zeta);
	             riseerr = fRTWeights->Eval(stationVector[i].GetSPDistance(), secZenith, stationVector[i].GetTotalSignal());

                     if( (stationVector[i].GetTotalSignal() > minSignal) )
	             {
                        if( (!stationVector[i].IsLowGainSaturated()) && (!includeSaturated) )
	                {
                           if( (stationVector[i].GetSPDistance() >= limitTankDistance[0]) && (stationVector[i].GetSPDistance() <= limitTankDistance[1]) )
	                   {
		              tempVect.push_back(stationVector[i].GetId());
		              tempVect.push_back((double)stationVector[i].IsHighGainSaturated());
                              tempVect.push_back(stationVector[i].GetSPDistance());
                              tempVect.push_back(stationVector[i].GetSPDistanceError());
                              tempVect.push_back(risemean);
                              tempVect.push_back(riseerr);
		              tempVect.push_back(sdrecshw->GetEnergy());
		              tempVect.push_back(sdrecshw->GetEnergyError());
		              tempVect.push_back(sdrecshw->GetZenith());
		              tempVect.push_back(sdrecshw->GetZenithError());
		              tempVect.push_back(sdrecshw->GetShowerSize());
		              tempVect.push_back(sdrecshw->GetShowerSizeError());

			      itemp[3]++;
			   }
	                   else
                              cout << "Rejected: Station distance " << stationVector[i].GetSPDistance() << " is outside the limits (" << limitTankDistance[0] << "m, " << limitTankDistance[1] << "m)." << endl;
			}
			else
                           cout << "Rejected: Station signal is low gain saturated." << endl;
		     }
		     else
                        cout << "Rejected: Station signal " << stationVector[i].GetTotalSignal() << " is below the minimum accepted (" << minSignal << " VEM)." << endl;
	          }
	       }

	       if(itemp[3] < minPoints)
	       {
                  goodrec = false;
                  cout << "Rejected: Only " << itemp[3] << " valid tanks." << endl;
                  cout << "Event reconstruction has failed." << endl;
	       }
	       else
               {
		  for(int i = 0; i < itemp[3]; i++)
		  {
		     // Save SD station ID values
		     sdidVect.push_back(tempVect[12*i]);
		     // Save if station is high gain saturated
		     HGsat.push_back((bool)tempVect[12*i+1]);
	             // Save distance of station to shower axis
                     distVect.push_back(tempVect[12*i+2]);
                     distVect.push_back(tempVect[12*i+3]);
	             // Save risetime for each station
                     riseVect.push_back(tempVect[12*i+4]);
                     riseVect.push_back(tempVect[12*i+5]);
		     // Save event energy
		     energy.push_back(tempVect[12*i+6]);
		     energy.push_back(tempVect[12*i+7]);
		     // Save event zenith angle
		     zenith.push_back(tempVect[12*i+8]);
		     zenith.push_back(tempVect[12*i+9]);
		     // Save event S1000 value
		     shwsize.push_back(tempVect[12*i+10]);
		     shwsize.push_back(tempVect[12*i+11]);
	             allcount++;
		  }

	          allevts++;
	       }
	    }
	    else
               cout << "Event reconstruction has failed." << endl;
	 }

         fFile->Close();
         delete fFile;*/
/*----------------------------------------------*/
      }

      cout << "Surviving a total of " << allcount << " risetimes from " << allevts << " valid events." << endl;
      cerr << "Surviving a total of " << allcount << " risetimes from " << allevts << " valid events." << endl;

      TF1 *fitfuncLow;
      TF1 *fitfuncHigh;
      itemp[0] = 0;

      if(sdBenchmark != 1)
      {
         parA = new TGraphErrors();
         parB = new TGraphErrors();
         parN = new TGraphErrors();

         stemp[0] = "./benchmark_functions_en_" + ToString(energybins->at(2*energyref),2) + "-" + ToString(energybins->at(2*energyref+1),2) + ".txt";
         ofs->open(stemp[0].c_str(), ofstream::out);
      }

      for(int iZen = 0; iZen < zenithnr; iZen++)
      {
         zenithlimit[0] = zenithbins->at(2*iZen);
         zenithlimit[1] = zenithbins->at(2*iZen+1);

	 energylimit[0] = energybins->at(2*energyref);
	 energylimit[1] = energybins->at(2*energyref+1);

         cout << "Reference energy limit = " << energylimit[0] << ", " << energylimit[1] << endl;
         cout << "Chosen zenith angle limit = " << zenithlimit[0] << ", " << zenithlimit[1] << endl;

//	 cout << "distVect size = " << distVect.size() << ", riseVect  size = " << riseVect.size() << ", energy size = " << energy.size() << ", zenith size = " << zenith.size() << ", sdidVect size = " << sdidVect.size() << ", HGsat size = " << HGsat.size() << endl;

         distVectLow.clear();
         riseVectLow.clear();
         distVectHigh.clear();
         riseVectHigh.clear();

	 for(int i = 0; i < sdidVect->size(); i++)
	 {
	    goodrec = true;

            if(SecTheta(zenith->at(2*i),false) < zenithlimit[0])
               goodrec = false;
            if(SecTheta(zenith->at(2*i),false) > zenithlimit[1])
               goodrec = false;
            if(TMath::Log10(energy->at(2*i)) < energylimit[0])
               goodrec = false;
            if(TMath::Log10(energy->at(2*i)) > energylimit[1])
               goodrec = false;

	    if(goodrec)
	    {
               cout << "Event " << eventVect->at(i) << ", station " << sdidVect->at(i) << ", energy = " << energy->at(2*i) << " ± " << energy->at(2*i+1)  << ", zenith = " << SecTheta(zenith->at(2*i),false) << " ± " << (TMath::Sin(zenith->at(2*i))*zenith->at(2*i+1))/TMath::Power(TMath::Cos(zenith->at(2*i)),2) << endl;
               if(HGsat->at(i))
               {
                  // Save distance of station to shower axis
                  distVectLow.push_back(distVect->at(2*i));
                  distVectLow.push_back(distVect->at(2*i+1));
                  // Save risetime for each station
                  riseVectLow.push_back(riseVect->at(2*i));
                  riseVectLow.push_back(riseVect->at(2*i+1));
                  cout << "Important: Station signal (" << sdidVect->at(i) << ") is high gain saturated." << endl;
               }
               else
               {
                  // Save distance of station to shower axis
                  distVectHigh.push_back(distVect->at(2*i));
                  distVectHigh.push_back(distVect->at(2*i+1));
                  // Save risetime for each station
                  riseVectHigh.push_back(riseVect->at(2*i));
                  riseVectHigh.push_back(riseVect->at(2*i+1));
               }
	    }
	 }

	 allcount = 0;
	 allcountLow = 0;
	 allcountHigh = 0;

         allRise = new TGraphErrors();
         allRiseLowGain = new TGraphErrors();
         allRiseHighGain = new TGraphErrors();

/*         riseErrDist = new TH1D();
         riseErrDistLow = new TH1D();
         riseErrDistHigh = new TH1D();*/

         for(int i = 0; i < distVectLow.size()/2; i++)
         {
            cout << "(" << allcount << "," << allcountLow << ") Low: " << distVectLow[2*i] << " (" << distVectLow[2*i+1] << ")\t" << riseVectLow[2*i] << " (" << riseVectLow[2*i+1] << ")" << endl;
            allRise->SetPoint(allcount, distVectLow[2*i], riseVectLow[2*i]);
            allRise->SetPointError(allcount, distVectLow[2*i+1], riseVectLow[2*i+1]);
            allRiseLowGain->SetPoint(allcountLow, distVectLow[2*i], riseVectLow[2*i]);
            allRiseLowGain->SetPointError(allcountLow, distVectLow[2*i+1], riseVectLow[2*i+1]);
            allcount++;
            allcountLow++;
         }
         
         for(int i = 0; i < distVectHigh.size()/2; i++)
         {
            cout << "(" << allcount << "," << allcountHigh << ") High: " << distVectHigh[2*i] << " (" << distVectHigh[2*i+1] << ")\t" << riseVectHigh[2*i] << " (" << riseVectHigh[2*i+1] << ")" << endl;
            allRise->SetPoint(allcount, distVectHigh[2*i], riseVectHigh[2*i]);
            allRise->SetPointError(allcount, distVectHigh[2*i+1], riseVectHigh[2*i+1]);
            allRiseHighGain->SetPoint(allcountHigh, distVectHigh[2*i], riseVectHigh[2*i]);
            allRiseHighGain->SetPointError(allcountHigh, distVectHigh[2*i+1], riseVectHigh[2*i+1]);
            allcount++;
            allcountHigh++;
         }

	 if(allcount == 0)
	 {
            cout << "No risetime values found in this zenith angle bin. Skipping plotting procedures." << endl;
            cerr << "No risetime values found in this zenith angle bin. Skipping plotting procedures." << endl;
	 }
	 else
	 {
            cout << "Plotting " << allcount << " risetimes, " << allcountLow << " are high gain saturated and " << allcountHigh << " are not." << endl;
            cerr << "Plotting " << allcount << " risetimes, " << allcountLow << " are high gain saturated and " << allcountHigh << " are not." << endl;

            mystyle->SetGraphColor(allRise, 2);
            allRise->SetMarkerStyle(20);
            allRise->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(allRise, "Distance from shower axis [m]", "t_{1/2} [ns]");
            c1->SetLogx(kFALSE);
            c1->SetLogy(kFALSE);
            allRise->Draw("AP");

            mystyle->SetGraphColor(allRiseLowGain, 0);
            allRiseLowGain->SetMarkerStyle(20);
            allRiseLowGain->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(allRiseLowGain, "Distance from shower axis [m]", "t_{1/2} [ns]");
            allRiseLowGain->Draw("P;SAME");

            mystyle->SetGraphColor(allRiseHighGain, 1);
            allRiseHighGain->SetMarkerStyle(21);
            allRiseHighGain->SetMarkerSize(0.8);
            mystyle->SetAxisTitles(allRiseHighGain, "Distance from shower axis [m]", "t_{1/2} [ns]");
            allRiseHighGain->Draw("P;SAME");

	    if(sdBenchmark == 1)
	    {
               fitfuncLow = new TF1("fitfuncLow", "40+TMath::Sqrt(TMath::Power(([0]+[1]*TMath::Power([4],-4)),2)+([2]+[3]*TMath::Power([4],-4))*TMath::Power(x,2))-([0]+[1]*TMath::Power([4],-4))", adfile->GetDistanceLimit(0), adfile->GetDistanceLimit(1));
               fitfuncLow->SetParameter(0, sdfitparam[0]);
               fitfuncLow->SetParameter(1, sdfitparam[2]);
               fitfuncLow->SetParameter(2, sdfitparam[4]);
               fitfuncLow->SetParameter(3, sdfitparam[6]);
               fitfuncLow->SetParameter(4, (zenithlimit[0]+zenithlimit[1])/2.);
               fitfuncLow->SetParError(0, sdfitparam[1]);
               fitfuncLow->SetParError(1, sdfitparam[3]);
               fitfuncLow->SetParError(2, sdfitparam[5]);
               fitfuncLow->SetParError(3, sdfitparam[7]);

	       cout << endl << "Fitting parameters from low gain fit (SD analysis):" << endl;
	       cout << "- a0 = " << fitfuncLow->GetParameter(0) << " ± " << fitfuncLow->GetParError(0) << endl;
	       cout << "- a1 = " << fitfuncLow->GetParameter(1) << " ± " << fitfuncLow->GetParError(1) << endl;
	       cout << "- b0 = " << fitfuncLow->GetParameter(2) << " ± " << fitfuncLow->GetParError(2) << endl;
	       cout << "- b1 = " << fitfuncLow->GetParameter(3) << " ± " << fitfuncLow->GetParError(3) << endl;
	       cout << endl;

	       fitfuncLow->SetLineColor(1);
	       fitfuncLow->SetLineWidth(2);
	       fitfuncLow->SetLineStyle(1);
	       fitfuncLow->Draw("SAME");

               fitfuncHigh = new TF1("fitfuncHigh", "40+([4]+[5]*TMath::Power([7],2)+[6]*TMath::Exp([7]))*(TMath::Sqrt(TMath::Power(([0]+[1]*TMath::Power([7],-4)),2)+([2]+[3]*TMath::Power([7],-4))*TMath::Power(x,2))-([0]+[1]*TMath::Power([7],-4)))", adfile->GetDistanceLimit(0), adfile->GetDistanceLimit(1));
               fitfuncHigh->SetParameter(0, sdfitparam[0]);
               fitfuncHigh->SetParameter(1, sdfitparam[2]);
               fitfuncHigh->SetParameter(2, sdfitparam[4]);
               fitfuncHigh->SetParameter(3, sdfitparam[6]);
               fitfuncHigh->SetParameter(4, sdfitparam[8]);
               fitfuncHigh->SetParameter(5, sdfitparam[10]);
               fitfuncHigh->SetParameter(6, sdfitparam[12]);
               fitfuncHigh->SetParameter(7, (zenithlimit[0]+zenithlimit[1])/2.);
               fitfuncHigh->SetParError(0, sdfitparam[1]);
               fitfuncHigh->SetParError(1, sdfitparam[3]);
               fitfuncHigh->SetParError(2, sdfitparam[5]);
               fitfuncHigh->SetParError(3, sdfitparam[7]);
               fitfuncHigh->SetParError(4, sdfitparam[9]);
               fitfuncHigh->SetParError(5, sdfitparam[11]);
               fitfuncHigh->SetParError(6, sdfitparam[13]);

	       cout << endl << "Fitting parameters from high gain fit (SD analysis):" << endl;
	       cout << "- n0 = " << fitfuncHigh->GetParameter(4) << " ± " << fitfuncHigh->GetParError(4) << endl;
	       cout << "- n1 = " << fitfuncHigh->GetParameter(5) << " ± " << fitfuncHigh->GetParError(5) << endl;
	       cout << "- n2 = " << fitfuncHigh->GetParameter(6) << " ± " << fitfuncHigh->GetParError(6) << endl;
	       cout << endl;

	       fitfuncHigh->SetLineColor(1);
	       fitfuncHigh->SetLineWidth(2);
	       fitfuncHigh->SetLineStyle(9);
	       fitfuncHigh->Draw("SAME");
	    }
	    else
	    {
               fitfuncLow = new TF1("fitfuncLow", "40+TMath::Sqrt(TMath::Power([0],2)+[1]*TMath::Power(x,2))-[0]", adfile->GetDistanceLimit(0), adfile->GetDistanceLimit(1));
               fitfuncLow->SetParameters(100.,0.1);
               fitfuncLow->SetParLimits(0,0.,1000.);
               fitfuncLow->SetParLimits(1,0.,1.);
               allRiseLowGain->Fit("fitfuncLow","0");

	       tempfunc = (TF1*)allRiseLowGain->GetFunction("fitfuncLow");
	       tempfunc->SetLineColor(1);
	       tempfunc->SetLineWidth(2);
	       tempfunc->SetLineStyle(1);
	       tempfunc->Draw("SAME");

	       for(int j = 0; j < 2; j++)
	       {
                  fitparamLow[j] = fitfuncLow->GetParameter(j);
                  fitparamLowErr[j] = fitfuncLow->GetParError(j);
	       }

	       cout << endl << "Fitting parameters from low gain fit (chi2/ndf = " << fitfuncLow->GetChisquare() << "/" << fitfuncLow->GetNDF() << " = " << (fitfuncLow->GetChisquare())/(fitfuncLow->GetNDF()) << "):" << endl;
	       cout << "- A = " << fitparamLow[0] << " ± " << fitparamLowErr[0] << endl;
	       cout << "- B = " << fitparamLow[1] << " ± " << fitparamLowErr[1] << endl;
	       cout << endl;

	       parA->SetPoint(itemp[0], (zenithlimit[0]+zenithlimit[1])/2., fitparamLow[0]);
	       parA->SetPointError(itemp[0], 0., fitparamLowErr[0]);
	       parB->SetPoint(itemp[0], (zenithlimit[0]+zenithlimit[1])/2., fitparamLow[1]);
	       parB->SetPointError(itemp[0], 0., fitparamLowErr[1]);

	       *ofs << energylimit[0] << "\t" << energylimit[1] << "\t" << zenithlimit[0] << "\t" << zenithlimit[1] << "\t" << fitparamLow[0] << "\t" << fitparamLowErr[0] << "\t" << fitparamLow[1] << "\t" << fitparamLowErr[1] << "\t";

               fitfuncHigh = new TF1("fitfuncHigh", "40+[2]*(TMath::Sqrt(TMath::Power([0],2)+[1]*TMath::Power(x,2))-[0])", adfile->GetDistanceLimit(0), adfile->GetDistanceLimit(1));
               fitfuncHigh->FixParameter(0, fitfuncLow->GetParameter(0));
               fitfuncHigh->FixParameter(1, fitfuncLow->GetParameter(1));
               fitfuncHigh->SetParameter(2, 1.1);
               allRiseHighGain->Fit("fitfuncHigh","0");

	       tempfunc = (TF1*)allRiseHighGain->GetFunction("fitfuncHigh");
	       tempfunc->SetLineColor(1);
	       tempfunc->SetLineWidth(2);
	       tempfunc->SetLineStyle(9);
	       tempfunc->Draw("SAME");

	       for(int j = 0; j < 3; j++)
	       {
                  fitparamHigh[j] = fitfuncHigh->GetParameter(j);
                  fitparamHighErr[j] = fitfuncHigh->GetParError(j);
	       }

	       cout << endl << "Fitting parameters from high gain fit (chi2/ndf = " << fitfuncHigh->GetChisquare() << "/" << fitfuncHigh->GetNDF() << " = " << (fitfuncHigh->GetChisquare())/(fitfuncHigh->GetNDF()) << "):" << endl;
	       cout << "- A = " << fitparamHigh[0] << " ± " << fitparamHighErr[0] << endl;
	       cout << "- B = " << fitparamHigh[1] << " ± " << fitparamHighErr[1] << endl;
	       cout << "- N = " << fitparamHigh[2] << " ± " << fitparamHighErr[2] << endl;
	       cout << endl;

	       parN->SetPoint(itemp[0], (zenithlimit[0]+zenithlimit[1])/2., fitparamHigh[2]);
	       parN->SetPointError(itemp[0], 0., fitparamHighErr[2]);

	       *ofs << fitparamHigh[2] << "\t" << fitparamHighErr[2] << endl;
	    }

	    itemp[0]++;

            allRise->GetXaxis()->SetRange((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRise->GetXaxis()->SetRangeUser((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRise->GetXaxis()->SetLimits((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRiseLowGain->GetXaxis()->SetRange((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRiseLowGain->GetXaxis()->SetRangeUser((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRiseLowGain->GetXaxis()->SetLimits((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRiseHighGain->GetXaxis()->SetRange((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRiseHighGain->GetXaxis()->SetRangeUser((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRiseHighGain->GetXaxis()->SetLimits((adfile->GetDistanceLimit(0)-50.), (adfile->GetDistanceLimit(1)+150.));
            allRise->GetYaxis()->SetRange(0., 1100.);
            allRise->GetYaxis()->SetRangeUser(0., 1100.);
            allRise->GetYaxis()->SetLimits(0., 1100.);

            stemp[0] = "rm ./plots/risetime_vs_distance";
            BinNaming(&stemp[0], energylimit, zenithlimit);
            stemp[0] += "*";
            system(stemp[0].c_str());

            stemp[0] = "./plots/risetime_vs_distance";
            BinNaming(&stemp[0], energylimit, zenithlimit);

            stemp[1] = stemp[0] + ".pdf";
            cout << stemp[1] << endl;
            c1->SaveAs(stemp[1].c_str());
            stemp[1] = stemp[0] + ".C";
            cout << stemp[1] << endl;
            c1->SaveAs(stemp[1].c_str());

	    delete fitfuncLow;
	    delete fitfuncHigh;
	 }

	 delete allRise;
	 delete allRiseLowGain;
	 delete allRiseHighGain;

/*	 delete riseErrDist;
	 delete riseErrDistLow;
	 delete riseErrDistHigh;*/
      }

      if(sdBenchmark != 1)
      {
         mystyle->SetGraphColor(parA, 0);
         parA->SetMarkerStyle(20);
         parA->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(parA, "SD zenith angle [sec(#theta)]", "Fitting parameter A [ns]");
         c1->SetLogx(kFALSE);
         c1->SetLogy(kFALSE);
         parA->Draw("AP");

         parfit[0] = new TF1("parfit0", "[0]+[1]*TMath::Power(x,-4)", 1.0, 2.0);
         parfit[0]->SetParameters(-50.,100.);
         parA->Fit("parfit0");

         cout << endl << "Fitting parameters from sec(theta) vs A (chi2/ndf = " << parfit[0]->GetChisquare() << "/" << parfit[0]->GetNDF() << " = " << (parfit[0]->GetChisquare())/(parfit[0]->GetNDF()) << "):" << endl;
         cout << "- a0 = " << parfit[0]->GetParameter(0) << " ± " << parfit[0]->GetParError(0) << endl;
         cout << "- a1 = " << parfit[0]->GetParameter(1) << " ± " << parfit[0]->GetParError(1) << endl;
         cout << endl;

         stemp[0] = "rm ./plots/fitting_parA.*";
         system(stemp[0].c_str());
         stemp[0] = "./plots/fitting_parA.pdf";
         c1->SaveAs(stemp[0].c_str());
//         stemp[0] = "./plots/fitting_parA.C";
//         c1->SaveAs(stemp[0].c_str());

         mystyle->SetGraphColor(parB, 0);
         parB->SetMarkerStyle(20);
         parB->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(parB, "SD zenith angle [sec(#theta)]", "Fitting parameter B [ns^{2}/m^{2}]");
         c1->SetLogx(kFALSE);
         c1->SetLogy(kFALSE);
         parB->Draw("AP");

         parfit[1] = new TF1("parfit1", "[0]+[1]*TMath::Power(x,-4)", 1.0, 2.0);
         parfit[1]->SetParameters(-0.1,0.1);
         parB->Fit("parfit1");

         cout << endl << "Fitting parameters from sec(theta) vs B (chi2/ndf = " << parfit[1]->GetChisquare() << "/" << parfit[1]->GetNDF() << " = " << (parfit[1]->GetChisquare())/(parfit[1]->GetNDF()) << "):" << endl;
         cout << "- b0 = " << parfit[1]->GetParameter(0) << " ± " << parfit[1]->GetParError(0) << endl;
         cout << "- b1 = " << parfit[1]->GetParameter(1) << " ± " << parfit[1]->GetParError(1) << endl;
         cout << endl;

         stemp[0] = "rm ./plots/fitting_parB.*";
         system(stemp[0].c_str());
         stemp[0] = "./plots/fitting_parB.pdf";
         c1->SaveAs(stemp[0].c_str());
//         stemp[0] = "./plots/fitting_parB.C";
//         c1->SaveAs(stemp[0].c_str());

         mystyle->SetGraphColor(parN, 0);
         parN->SetMarkerStyle(20);
         parN->SetMarkerSize(0.8);
         mystyle->SetAxisTitles(parN, "SD zenith angle [sec(#theta)]", "Fitting parameter N");
         c1->SetLogx(kFALSE);
         c1->SetLogy(kFALSE);
         parN->Draw("AP");

         parfit[2] = new TF1("parfit2", "[0]+[1]*TMath::Power(x,2)+[2]*TMath::Exp(x)", 1.0, 2.0);
         parfit[2]->SetParameters(-0.1,-1.1,1.);
         parN->Fit("parfit2");

         cout << endl << "Fitting parameters from sec(theta) vs N (chi2/ndf = " << parfit[2]->GetChisquare() << "/" << parfit[2]->GetNDF() << " = " << (parfit[2]->GetChisquare())/(parfit[2]->GetNDF()) << "):" << endl;
         cout << "- n0 = " << parfit[2]->GetParameter(0) << " ± " << parfit[2]->GetParError(0) << endl;
         cout << "- n1 = " << parfit[2]->GetParameter(1) << " ± " << parfit[2]->GetParError(1) << endl;
         cout << "- n2 = " << parfit[2]->GetParameter(2) << " ± " << parfit[2]->GetParError(2) << endl;
         cout << endl;

         stemp[0] = "rm ./plots/fitting_parN.*";
         system(stemp[0].c_str());
         stemp[0] = "./plots/fitting_parN.pdf";
         c1->SaveAs(stemp[0].c_str());
//         stemp[0] = "./plots/fitting_parN.C";
//         c1->SaveAs(stemp[0].c_str());

         ofs->close();
      }
      delete ofs;

      if(sdBenchmark != 1)
      {
         // Convert risetimes into delta, using benchmark functions from before
         stemp[0] = "./benchmark_functions_en_" + ToString(energybins->at(2*energyref),2) + "-" + ToString(energybins->at(2*energyref+1),2) + ".txt";
         ifs->open(stemp[0].c_str(), ifstream::in);
      }

      allcount = 0;
      allevts = 0;
      cout << "Converting risetimes into delta values" << endl;

      allDelta = new TGraphErrors();

      energyVect.clear();
      deltaVect.clear();

      dtemp[6] = 0.;

      for(int iZen = 0; iZen < zenithnr; iZen++)
      {
         zenithlimit[0] = zenithbins->at(2*iZen);
         zenithlimit[1] = zenithbins->at(2*iZen+1);

         cout << "Chosen zenith angle limit = " << zenithlimit[0] << ", " << zenithlimit[1] << endl;

	 if(sdBenchmark != 1)
	 {
            for(int i = 0; i < 4; i++)
               *ifs >> dtemp[0];
            for(int i = 0; i < 6; i++)
               *ifs >> fitparam[i];

            cout << "Fitting parameters for benchmark function:" << endl;
            for(int i = 0; i < 3; i++)
               cout << "- fitparam[" << i << "] = " << fitparam[2*i] << " ± " << fitparam[2*i+1] << endl;
            cout << endl;
	 }

         dtemp[0] = 0.;
         itemp[1] = 0;
         itemp[2] = 0;

         for(int i = 0; i < sdidVect->size(); i++)
         {
            goodrec = true;
        
            if(SecTheta(zenith->at(2*i),false) < zenithlimit[0])
               goodrec = false;
            if(SecTheta(zenith->at(2*i),false) > zenithlimit[1])
               goodrec = false;

            if(i == sdidVect->size()-1)
            {
	       if(itemp[1] == 0)
                  cout << "Error! No valid events found!" << endl;
	       else
	       {
                  dtemp[1] = dtemp[1]/itemp[0];
                  dtemp[2] = TMath::Sqrt(dtemp[2])/itemp[0];

	          if(deltaVect[2*(allevts-1)] == dtemp[1])
                     cout << "Error! Last value already inserted!" << endl;
		  else
		  {
                     cout << "Last: Delta value from " << itemp[0] << " SD stations for event " << allevts << " (evt " << itemp[2] << ") = " << dtemp[1] << " ± " << dtemp[2] << endl;
 
                     // Save energy of the event
                     energyVect.push_back(energy->at(2*i));
                     energyVect.push_back(energy->at(2*i+1));
      	             // Save Delta of the event
                     deltaVect.push_back(dtemp[1]);
                     deltaVect.push_back(dtemp[2]);
                     allevts++;

                     if(TMath::Abs(dtemp[1]) > dtemp[6])
                        dtemp[6] = TMath::Abs(dtemp[1]);
		  }
	       }
            }
        
            if(goodrec)
            {
               // Calculate parameters A, B and N for the event
	       if(sdBenchmark == 1)
               {
		  CalculateSdFitParams(sdfitparam, fitparam, zenith->at(2*i), zenith->at(2*i+1));

                  cout << "Fitting parameters for benchmark function (from SD analysis):" << endl;
                  for(int i = 0; i < 3; i++)
                     cout << "- fitparam[" << i << "] = " << fitparam[2*i] << " ± " << fitparam[2*i+1] << endl;
                  cout << endl;
	       }

               // Check if this is already a new event or not
//               if(dtemp[0] != energy->at(2*i))
               if(itemp[2] != eventVect->at(i))
               {
                  if(allcount > 0)
                  {
                     dtemp[1] = dtemp[1]/itemp[0];
                     dtemp[2] = TMath::Sqrt(dtemp[2])/itemp[0];
                     cout << "Delta value from " << itemp[0] << " SD stations for event " << allevts << " (evt " << itemp[2] << ") = " << dtemp[1] << " ± " << dtemp[2] << endl;

                     // Save energy of the event
                     energyVect.push_back(energy->at(2*i));
                     energyVect.push_back(energy->at(2*i+1));
             	     // Save Delta of the event
                     deltaVect.push_back(dtemp[1]);
                     deltaVect.push_back(dtemp[2]);
                     allevts++;

             	     if(TMath::Abs(dtemp[1]) > dtemp[6])
                        dtemp[6] = TMath::Abs(dtemp[1]);
                  }

                  itemp[0] = 0;
                  dtemp[1] = 0.;
                  dtemp[2] = 0.;
                  cout << "New event (" << allevts << ")" << endl;
               }

               cout << "Event " << eventVect->at(i) << ", station " << sdidVect->at(i) << ", energy = " << energy->at(2*i)/1.e+18 << " ± " << energy->at(2*i+1)/1.e+18 << ", zenith = " << SecTheta(zenith->at(2*i),false) << " ± " << (TMath::Sin(zenith->at(2*i))*zenith->at(2*i+1))/TMath::Power(TMath::Cos(zenith->at(2*i)),2) << ", A = " << fitparam[0] << " ± " << fitparam[1] << ", B = " << fitparam[2] << " ± " << fitparam[3] << ", N = " << fitparam[4] << " ± " << fitparam[5] << endl;

               CalculateBenchmark(tbench, fitparam, distVect->at(2*i), distVect->at(2*i+1), HGsat->at(i));
               // Benchmark function (HGsat) and its uncertainty
               //   tbench = 40 + sqrt(A^2 + B*r^2) - A
/*               if(HGsat->at(i))
               {
                  tbench[0] = 40 + TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)) - fitparam[0];
                  cout << "  tbench[0] = " << tbench[0];
                  // (dA)^2 = (dA*(A/sqrt(A^2 + B*r^2) - 1))^2
                  dtemp[3] = fitparam[1]*(fitparam[0]/TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)) - 1);
                  dtemp[4] = TMath::Power(dtemp[3],2);
                  cout << ", dA = " << dtemp[3];
                  // (dB)^2 = (dB*r^2/(2*sqrt(A^2 + B*r^2)))^2
                  dtemp[3] = fitparam[3]*TMath::Power(distVect->at(2*i),2)/(2*TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)));
                  dtemp[4] += TMath::Power(dtemp[3],2);
                  cout << ", dB = " << dtemp[3];
                  // (dr)^2 = (dr*B*r/(sqrt(A^2 + B*r^2)))^2
                  dtemp[3] = (distVect->at(2*i+1))*fitparam[2]*(distVect->at(2*i))/(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)));
                  dtemp[4] += TMath::Power(dtemp[3],2);
                  cout << ", dr = " << dtemp[3];
               }
               // Benchmark function and its uncertainty
               //   tbench = 40 + N*(sqrt(A^2 + B*r^2) - A)
               else
               {
                  tbench[0] = 40 + fitparam[4]*(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)) - fitparam[0]);
                  cout << "  tbench[0] = " << tbench[0];
                  // (dA)^2 = (N*dA*(A/sqrt(A^2 + B*r^2) - 1))^2
                  dtemp[3] = fitparam[4]*fitparam[1]*(fitparam[0]/TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)) - 1);
                  dtemp[4] = TMath::Power(dtemp[3],2);
                  cout << ", dA = " << dtemp[3];
                  // (dB)^2 = (N*dB*r^2/(2*sqrt(A^2 + B*r^2)))^2
                  dtemp[3] = fitparam[4]*fitparam[3]*TMath::Power(distVect->at(2*i),2)/(2*TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)));
                  dtemp[4] += TMath::Power(dtemp[3],2);
                  cout << ", dB = " << dtemp[3];
                  // (dN)^2 = (dN*(sqrt(A^2 + B*r^2) - A))^2
                  dtemp[3] = fitparam[5]*(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)) - fitparam[0]);
                  dtemp[4] += TMath::Power(dtemp[3],2);
                  cout << ", dN = " << dtemp[3];
                  // (dr)^2 = (dr*B*N*r/(sqrt(A^2 + B*r^2)))^2
                  dtemp[3] = (distVect->at(2*i+1))*fitparam[2]*fitparam[4]*(distVect->at(2*i))/(TMath::Sqrt(TMath::Power(fitparam[0],2) + fitparam[2]*TMath::Power(distVect->at(2*i),2)));
                  dtemp[4] += TMath::Power(dtemp[3],2);
                  cout << ", dr = " << dtemp[3];
               }

               tbench[1] = TMath::Sqrt(dtemp[4]);
               cout << ", tbench[1] = " << tbench[1] << endl;*/

               dtemp[3] = riseVect->at(2*i);
               dtemp[4] = riseVect->at(2*i+1);

               cout << itemp[0] << ": distance = " << distVect->at(2*i) << " ± " << distVect->at(2*i+1) << ", t_bench = " << tbench[0] << " ± " << tbench[1] << ", risetime = " << dtemp[3] << " ± " << dtemp[4] << ", delta = " << /*(*/dtemp[3] - tbench[0]/*)/dtemp[4]*/ << " ± " << TMath::Sqrt(TMath::Power(dtemp[4],2) + TMath::Power(tbench[1],2)) << endl;

//             dtemp[1] += (dtemp[3] - tbench[0])/dtemp[4];
               dtemp[1] += dtemp[3] - tbench[0];
               dtemp[2] += TMath::Power(dtemp[4],2) + TMath::Power(tbench[1],2);

               itemp[0]++;
	       itemp[1]++;
               allcount++;
//               dtemp[0] = energy->at(2*i);
	       itemp[2] = eventVect->at(i);
            }
         }
      }

      cout << "Surviving a total of " << allcount << " deltas from " << allevts << " valid events." << endl;
      cerr << "Surviving a total of " << allcount << " deltas from " << allevts << " valid events." << endl;

      for(int i = 0; i < energyVect.size()/2; i++)
      {
         allDelta->SetPoint(i, energyVect[2*i]/1.e+18, deltaVect[2*i]);
         allDelta->SetPointError(i, energyVect[2*i+1]/1.e+18, deltaVect[2*i+1]);

//	 cout << i << ": " << energyVect[2*i]/1.e+18 << " ± " << energyVect[2*i+1]/1.e+18 << ", " << deltaVect[2*i] << " ± " << deltaVect[2*i+1] << endl;
      }

      mystyle->SetGraphColor(allDelta, 0);
      allDelta->SetMarkerStyle(21);
      allDelta->SetMarkerSize(0.8);
//      mystyle->SetAxisTitles(allDelta, "SD energy (EeV)", "#Delta_{s}");
      mystyle->SetAxisTitles(allDelta, "SD energy (EeV)", "t_{1/2} - t_{1/2}^{bench} (ns)");
      c1->SetLogx(kTRUE);
      c1->SetLogy(kFALSE);
      allDelta->GetXaxis()->SetMoreLogLabels(kTRUE);
/*      dtemp[6] = 1.1*dtemp[6];
      allDelta->GetYaxis()->SetRange(-dtemp[6], dtemp[6]);
      allDelta->GetYaxis()->SetRangeUser(-dtemp[6], dtemp[6]);*/
      allDelta->Draw("AP");
      c1->Update();
      cout << TMath::Power(10.,energylimit[0])/1.e+18 << "\t" << TMath::Power(10.,energylimit[1])/1.e+18 << "\t" << gPad->GetUymin() << "\t" << gPad->GetUymax() << endl;
      line->DrawLine(TMath::Power(10.,energylimit[0])/1.e+18, gPad->GetUymin(), TMath::Power(10.,energylimit[0])/1.e+18, gPad->GetUymax());
      line2->DrawLine(TMath::Power(10.,energylimit[1])/1.e+18, gPad->GetUymin(), TMath::Power(10.,energylimit[1])/1.e+18, gPad->GetUymax());

      stemp[0] = "rm ./plots/delta_vs_energySD.pdf";
//      BinNaming(&stemp[0], energylimit, zenithlimit);
//      stemp[0] += "*";
      system(stemp[0].c_str());

      stemp[1] = "./plots/delta_vs_energySD.pdf";
//      BinNaming(&stemp[0], energylimit, zenithlimit);

//      stemp[1] = stemp[0] + ".pdf";
      cout << stemp[1] << endl;
      c1->SaveAs(stemp[1].c_str());
      stemp[1] = "./plots/delta_vs_energySD.C";
      c1->SaveAs(stemp[1].c_str());

      delete allDelta;

      delete adfile;

      if(sdBenchmark != 1)
      {
         ifs->close();

         delete parfit[0];
         delete parfit[1];
         delete parfit[2];

         delete parA;
         delete parB;
         delete parN;
      }

      delete ifs;

      delete mystyle;
      delete c1;
   }
   else
   {
      cerr << "Error! No input files supplied. Rerun program and add input files as arguments (ADST files)." << endl;
      return 1;
   }
/*
   delete[] byrange;
   delete[] bzrange;

   delete fRecEvent;
   delete fDetGeo;
   delete sdrecshw;*/

   delete distVect;
   delete riseVect;
   delete sdidVect;
   delete energy;
   delete zenith;
   delete shwsize;
   delete HGsat;

   delete[] stemp;
   delete[] itemp;
   delete[] dtemp;

   delete line;
   delete line2;

   return 0;
}
