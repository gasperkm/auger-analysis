#include "observables.h"
#include "separate_functions.h"
#include "smearing_selection.h"
#include "smearing_resolution.h"

Observables::Observables(vector<string> obs)
{
   nrobs = obs.size();

   obsstruct = new ObserStruct[nrobs];
   for(int i = 0; i < nrobs; i++)
      obsstruct[i].name = obs[i];
   Zero();

   minimum = new float[nrobs];
   maximum = new float[nrobs];
   desc = new string[nrobs];

   timer = new time_t;
   *timer = time(NULL);
   rand = new TRandom(*timer);
}

Observables::~Observables()
{
   delete timer;
   delete rand;
   delete[] obsstruct;
   delete[] minimum;
   delete[] maximum;
   delete[] desc;
}

string Observables::GetName(int obs)
{
   return obsstruct[obs].name;
}

string Observables::GetName(string obsname)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         return GetName(i);
   }

   return "";
}

float Observables::GetValue(int obs)
{
   return obsstruct[obs].value;
}

float Observables::GetValue(string obsname)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         return GetValue(i);
   }

   return -1;
}

int Observables::GetInt(string obsname)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         return i;
   }

   return -1;
}

void Observables::SetValue(int obs, float val)
{
   obsstruct[obs].value = val;
}

void Observables::SetValue(string obsname, float val)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         SetValue(i, val);
   }
}

void Observables::Zero(int obs)
{
   obsstruct[obs].value = 0;
}

void Observables::Zero(string obsname)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         Zero(i);
   }
}

void Observables::Zero()
{
   for(int i = 0; i < nrobs; i++)
      Zero(i);
}

int Observables::Count()
{
   return nrobs;
}

void Observables::SetMin(int obs, float value)
{
   minimum[obs] = value;
}

void Observables::SetMin(string obsname, float value)
{
   minimum[GetInt(obsname)] = value;
}

void Observables::SetMax(int obs, float value)
{
   maximum[obs] = value;
}

void Observables::SetMax(string obsname, float value)
{
   maximum[GetInt(obsname)] = value;
}

float Observables::GetMin(int obs)
{
   return minimum[obs];
}

float Observables::GetMin(string obsname)
{
   return minimum[GetInt(obsname)];
}

float Observables::GetMax(int obs)
{
   return maximum[obs];
}

float Observables::GetMax(string obsname)
{
   return maximum[GetInt(obsname)];
}

void Observables::SetLabel(int obs, string value)
{
   desc[obs] = value;
}

void Observables::SetLabel(string obsname, string value)
{
   desc[GetInt(obsname)] = value;
}

string Observables::GetLabel(int obs)
{
   return desc[obs];
}

string Observables::GetLabel(string obsname)
{
   return desc[GetInt(obsname)];
}

// Applying Xmax correction to Auger FD standard data
void Observables::ApplyCorrectionFD()
{
   float *logE = new float;
   float *ftemp = new float[4];

   ftemp[0] = GetValue("energyFD");
   ftemp[1] = GetValue("xmax");
   // Only apply correction, if we have a valid value
   if( (ftemp[0] != -1) && (ftemp[1] != -1) )
   {
      *logE = TMath::Log10(ftemp[0]); 
      ftemp[2] = 6.5/(TMath::Exp((*logE - 18.23)/0.41) + 1.);
      ftemp[3] = 0.93*(*logE - 18.);
      ftemp[1] = ftemp[1] - ftemp[2] + 3.4 - ftemp[3];

      SetValue("xmax", ftemp[1]);
   }

   delete logE;
   delete[] ftemp;
}

// Applying Xmax and energy corrections to Auger HECO data
void Observables::ApplyCorrectionHECO()
{
   float *logE = new float;
   float *ftemp = new float[5];

   ftemp[0] = GetValue("energyFD");
   ftemp[1] = GetValue("xmax");
   // Only apply correction, if we have a valid value
   if( (ftemp[0] != -1) && (ftemp[1] != -1) )
   {
      *logE = TMath::Log10(ftemp[0]); 
      ftemp[2] = TMath::Log(*logE - 16.5);
      ftemp[0] = (ftemp[0])*(1 + 0.04536243 - 0.08317193*ftemp[2]);
      SetValue("energyFD", ftemp[0]);

      ftemp[0] = GetValue("energyFD");
      *logE = TMath::Log10(ftemp[0]); 
      ftemp[2] = -2.96976097 - 0.99199218*(*logE - 17.80937342);
      ftemp[3] = 6.5/(TMath::Exp((*logE - 18.23)/0.41) + 1.);
      if(*logE < 17.55)
         ftemp[4] = 0.12 - 6.43*(*logE - 17.55);
      else
         ftemp[4] = 0.12 - 0.27*(*logE - 17.55);
      ftemp[1] = ftemp[1] - (ftemp[2] + ftemp[3] + 0.5*ftemp[4]);
      SetValue("xmax", ftemp[1]);
   }

   delete logE;
   delete[] ftemp;
}

// Applying Xmax and energy corrections to Auger HECO data (Errors only)
void Observables::ApplyCorrectionHECOErrors(Observables *mean, int type)
{
   float *logE = new float;
   float *ftemp = new float[3];

   ftemp[0] = mean->GetValue("energyFD");
   ftemp[1] = GetValue("energyFD");
   // Only apply correction, if we have a valid value
   if( (ftemp[0] != -1) && (ftemp[1] != -1) )
   {
      *logE = TMath::Log10(ftemp[0]); 
      ftemp[2] = TMath::Log(*logE - 16.5);
      ftemp[1] = (ftemp[1])*(1 + 0.04536243 - 0.08317193*ftemp[2]);
      SetValue("energyFD", ftemp[1]);
   }

   delete logE;
   delete[] ftemp;
}

// Applying systematic uncertainty to the selected observable
void Observables::ApplyUncertainty(Observables *errneg, Observables *errpos, string uncerObs, int type)
{
   if((type == 0) || (type == 1))
   {
      float *ftemp = new float[3];

      for(int i = 0; i < nrobs; i++)
      {
	 if(GetName(i) == uncerObs)
	 {
            ftemp[0] = GetValue(i);
            ftemp[1] = errneg->GetValue(i);
            ftemp[2] = errpos->GetValue(i);

	    if(ftemp[0] != -1)
	    {
               if(uncerObs == "xmax")
	       {
	          if(type == 0)
                     ftemp[0] -= 10.;//ftemp[1];
	          if(type == 1)
                     ftemp[0] += 10.;//ftemp[2];
	          SetValue(i, ftemp[0]);
	       }
               else if((uncerObs == "risetimerecalc") || (uncerObs == "risetime"))
	       {
	          if(type == 0)
                     ftemp[0] -= 25.0;//ftemp[1];
	          if(type == 1)
                     ftemp[0] += 25.0;//ftemp[2];
	          SetValue(i, ftemp[0]);
	       }
               /*------ Define systematic uncertainty shifts for other observables here ------*/
               /*------ Define systematic uncertainty shifts for other observables here ------*/
	    }
	 }
      }

      delete[] ftemp;
   }
}

// Applying smearing to MC simulations
void Observables::ApplySmearing()
{
   float *logE = new float;
   float *ftemp = new float[4];

   ftemp[0] = GetValue("energyFD");
   ftemp[1] = GetValue("xmax");
   // Only apply correction, if we have a valid value
   if( (ftemp[0] != -1) && (ftemp[1] != -1) )
   {
      *logE = TMath::Log10(ftemp[0]);
//      cout << "------" << endl;
//      cout << "Energy = " << *logE << ", xmax = " << ftemp[1] << endl;

      // Sigma value for resolution
      erAna::reso::ESigma *sys = new erAna::reso::ESigma;
      *sys = erAna::reso::eZero;
      erAna::ESelection *selType = new erAna::ESelection;
      *selType = erAna::eDefault;

      // Smearing from multiple scattering
      ftemp[2] = XmaxVariance(*logE, erAna::reso::eMultipleScattering, erAna::reso::eMultipleScattering, *selType, *sys);
      ftemp[3] = ftemp[2];
      // Smearing from correlated VAOD (Vertical Aerosol Optical Depth)
      ftemp[2] = XmaxVariance(*logE, erAna::reso::eCorrVAOD, erAna::reso::eCorrVAOD, *selType, *sys);
      ftemp[3] += ftemp[2];
      // Smearing from molecular atmosphere
      ftemp[2] = XmaxVariance(*logE, erAna::reso::eMolAtm, erAna::reso::eMolAtm, *selType, *sys);
      ftemp[3] += ftemp[2];
      // Smearing from VAOD uniformity
      ftemp[2] = XmaxVariance(*logE, erAna::reso::eVAODUniformity, erAna::reso::eVAODUniformity, *selType, *sys);
      ftemp[3] += ftemp[2];
      // Smearing from alignment
      ftemp[2] = XmaxVariance(*logE, erAna::reso::eAlignment, erAna::reso::eAlignment, *selType, *sys);
      ftemp[3] += ftemp[2];
      // Smearing from statistical VAOD
      ftemp[2] = XmaxVariance(*logE, erAna::reso::eStatVAOD, erAna::reso::eStatVAOD, *selType, *sys);
      ftemp[3] += ftemp[2];
//      cout << "Additional Variance = " << ftemp[3] << endl;

      // Bias correction
      ftemp[2] = 3.4 - 0.93*(*logE - 18.);	// Should this be applied?
//      cout << "Bias correction = " << ftemp[2] << endl;

      // Smeared value of Xmax
      ftemp[0] = rand->Gaus(0, TMath::Sqrt(ftemp[3]));
//      cout << "Random gaussian = " << ftemp[0] << endl;
      ftemp[1] += ftemp[0] + ftemp[2];
//      cout << "Smeared Xmax = " << ftemp[1] << endl;
      SetValue("xmax", ftemp[1]);

      delete selType;
      delete sys;
   }

/*  double additionalVariance =
      XmaxVariance(profLgE, eMultipleScattering, eMultipleScattering, selType, sys) +
      XmaxVariance(profLgE, eCorrVAOD, eCorrVAOD, selType, sys) +
     XmaxVariance(profLgE, eMolAtm, eMolAtm, selType, sys) +
     XmaxVariance(profLgE, eVAODUniformity, eVAODUniformity, selType, sys) +
      XmaxVariance(profLgE, eAlignment, eAlignment, selType, sys) +
      XmaxVariance(profLgE, eStatVAOD, eStatVAOD, selType, sys);
     cout << "additionalVariance = " << additionalVariance << endl;
   const double biasCorr = 3.4 - 0.93  * (profLgE-18);
   cout << "biasCorr = " << biasCorr << endl;
  const double xmaxSmeared =
     xmax + rand.Gaus(0, sqrt(additionalVariance)) + biasCorr;
    cout << "xmaxSmeared = " << xmaxSmeared << endl;*/

   delete logE;
   delete[] ftemp;
}

// Setting up zenith as sec(theta) for analysis
void Observables::SetupZenith(int type, Observables *errneg, Observables *errpos)
{
   float *ftemp = new float[4];

   ftemp[0] = GetValue(type);
   ftemp[1] = errneg->GetValue(type);
   ftemp[2] = errpos->GetValue(type);

   // Only apply calculation, if we have a valid value
   if( (ftemp[0] != -1) && (ftemp[1] != -1) && (ftemp[2] != -1) )
   {
//      cout << "theta = " << RadToDeg(ftemp[0]) << ", " << ftemp[0] << " (" << ftemp[1] << ", " << ftemp[2] << ")" << ", sec(theta) = " << SecTheta(ftemp[0], false) << endl;

      ftemp[3] = tan(ftemp[0])/cos(ftemp[0]);

//      cout << "sin(theta)/cos2(theta) = " << ftemp[3] << ", dsec(theta) = " << ftemp[1]*ftemp[3] << endl;

      ftemp[0] = SecTheta(ftemp[0], false);
      ftemp[1] = ftemp[1]*ftemp[3];
      ftemp[2] = ftemp[2]*ftemp[3];

      SetValue(type, ftemp[0]);
      errneg->SetValue(type, ftemp[1]);
      errpos->SetValue(type, ftemp[2]);

//      cout << "final = " << ftemp[0] << " (" << ftemp[1] << ", " << ftemp[2] << ")" << endl;

   }

   delete[] ftemp;
}

// Converting S1000 to S38
void Observables::ConvertToS38(int type, Observables *errneg, Observables *errpos, vector<float> *fitresults)
{
//   cout << "# Convert to S38  #: " << "Converting S1000 to S38 for " << GetName(type) << endl;

   int *itemp = new int;
   *itemp = 0;

   float *ftemp = new float[4];
   double *seczenith = new double[3];
   float *s1000temp = new float[3];
   float *entemp = new float[3];
   float *fitpar = new float[4];
   float *fitparErr = new float[4];

   ftemp[0] = GetValue("zenithFD");
   ftemp[1] = errneg->GetValue("zenithFD");
   ftemp[2] = errpos->GetValue("zenithFD");

   s1000temp[0] = GetValue("shwsize");
   s1000temp[1] = errneg->GetValue("shwsize");
   s1000temp[2] = errpos->GetValue("shwsize");

   entemp[0] = GetValue("energyFD");
   entemp[1] = errneg->GetValue("energyFD");
   entemp[2] = errpos->GetValue("energyFD");

//   cout << "zenith = " << ftemp[0] << ", " << ftemp[1] << ", " << ftemp[2] << ", S1000 = " << s1000temp[0] << ", " << s1000temp[1] << ", " << s1000temp[2] << ", energy = " << TMath::Log10(entemp[0]) << ", " << TMath::Log10(entemp[1]) << ", " << TMath::Log10(entemp[2]) << endl;

   for(int i = 0; i < fitresults->size()/10; i++)
   {
      if( (TMath::Log10(entemp[0]) > fitresults->at(10*i)) && (TMath::Log10(entemp[0]) <= fitresults->at(10*i+1)) )
      {
         for(int j = 2; j < 10; j++)
         {
            fitpar[(j-2)/2] = fitresults->at(10*i+j);
            fitparErr[(j-2)/2] = fitresults->at(10*i+j+1);
	    j++;
	 }
	 *itemp = i;

         break;
      }
   }

   // Calculate sec(theta) from zenith angle
   seczenith[0] = SecTheta(ftemp[0],false);
   seczenith[1] = (TMath::Sin(ftemp[0])*ftemp[1])/TMath::Power(TMath::Cos(ftemp[0]),2);
   seczenith[2] = (TMath::Sin(ftemp[0])*ftemp[2])/TMath::Power(TMath::Cos(ftemp[0]),2);

//   cout << "sec(theta) = " << seczenith[0] << ", " << seczenith[1] << ", " << seczenith[2];

   if(s1000temp[0] == -1.)
   {
      SetValue(type, -1.);
      errneg->SetValue(type, -1.);
      errpos->SetValue(type, -1.);
//      cout << "S1000 value is negative for this event, setting S38 values to -1." << endl;
   }
   else
   {
      // x = (1/sec(theta))^2 - (cos(thetaref))^2
      ftemp[0] = TMath::Power(1./seczenith[0],2) - TMath::Power(TMath::Cos(DegToRad(38.)),2);
//      cout << ", x = " << ftemp[0];
      // fCIC = 1 + a*x + b*x^2 + c*x^3
      ftemp[1] = 1.+fitpar[1]*ftemp[0]+fitpar[2]*TMath::Power(ftemp[0],2)+fitpar[3]*TMath::Power(ftemp[0],3);
//      cout << ", fCIC = " << ftemp[1];
      // S38 = S1000/fCIC
      ftemp[2] = s1000temp[0]/ftemp[1];
//      cout << ", S38 = " << ftemp[2] << endl;
      
      SetValue(type, ftemp[2]);
      
      // (dS1000/fCIC)^2
      ftemp[2] = s1000temp[1]/ftemp[1];
      ftemp[3] = TMath::Power(ftemp[2],2);
//      cout << "  S1000 err = " << ftemp[2];
      // (-S1000*da*x/fCIC^2)^2
      ftemp[2] = -(s1000temp[0]*fitparErr[1]*ftemp[0])/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", a err = " << ftemp[2];
      // (-S1000*db*x^2/fCIC^2)^2
      ftemp[2] = -(s1000temp[0]*fitparErr[2]*TMath::Power(ftemp[0],2))/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", b err = " << ftemp[2];
      // (-S1000*dc*x^3/fCIC^2)^2
      ftemp[2] = -(s1000temp[0]*fitparErr[3]*TMath::Power(ftemp[0],3))/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", c err = " << ftemp[2];
      // (-S1000*dx*(a+2b*x+3c*x^2)/fCIC^2)^2
      // dx = 2*cos(theta)*sin(theta)*dtheta
      ftemp[2] = 2*seczenith[1]/TMath::Power(seczenith[0],3);
//      cout << ", dx = " << ftemp[2];
      ftemp[2] = -(s1000temp[0]*ftemp[2]*(fitpar[1]+2*fitpar[2]*ftemp[0]+3*fitpar[3]*TMath::Power(ftemp[0],2)))/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", x err = " << ftemp[2];
      
      ftemp[3] = TMath::Sqrt(ftemp[3]);
//      cout << ", dS38 neg = " << ftemp[3] << endl;
      errneg->SetValue(type, ftemp[3]);
      
      // (dS1000/fCIC)^2
      ftemp[2] = s1000temp[2]/ftemp[1];
      ftemp[3] = TMath::Power(ftemp[2],2);
//      cout << "  S1000 err = " << ftemp[2];
      // (-S1000*da*x/fCIC^2)^2
      ftemp[2] = -(s1000temp[0]*fitparErr[1]*ftemp[0])/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", a err = " << ftemp[2];
      // (-S1000*db*x^2/fCIC^2)^2
      ftemp[2] = -(s1000temp[0]*fitparErr[2]*TMath::Power(ftemp[0],2))/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", b err = " << ftemp[2];
      // (-S1000*dc*x^3/fCIC^2)^2
      ftemp[2] = -(s1000temp[0]*fitparErr[3]*TMath::Power(ftemp[0],3))/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", c err = " << ftemp[2];
      // (-S1000*dx*(a+2b*x+3c*x^2)/fCIC^2)^2
      // dx = 2*cos(theta)*sin(theta)*dtheta
      ftemp[2] = 2*seczenith[2]/TMath::Power(seczenith[0],3);
//      cout << ", dx = " << ftemp[2];
      ftemp[2] = -(s1000temp[0]*ftemp[2]*(fitpar[1]+2*fitpar[2]*ftemp[0]+3*fitpar[3]*TMath::Power(ftemp[0],2)))/TMath::Power(ftemp[1],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", x err = " << ftemp[2];
      
      ftemp[3] = TMath::Sqrt(ftemp[3]);
//      cout << ", dS38 pos = " << ftemp[3] << endl;
      errpos->SetValue(type, ftemp[3]);
   }

   delete itemp;
   delete[] fitpar;
   delete[] fitparErr;
   delete[] s1000temp;
   delete[] entemp;
   delete[] seczenith;
   delete[] ftemp;
}

// Converting S38 to DeltaS38
void Observables::ConvertToDeltaS38(int type, Observables *errneg, Observables *errpos, float *fitresults, int syst, bool *isdata)
{
//   cout << "# Convert to Delta S38  #: " << "Converting S38 to DeltaS38 for " << GetName(type) << endl;

   float *ftemp = new float[4];
   float *s38temp = new float[3];
   float *entemp = new float[3];

   // E, dE, S38, dS38
   s38temp[0] = GetValue(type);
   s38temp[1] = errneg->GetValue(type);
   s38temp[2] = errpos->GetValue(type);

   entemp[0] = GetValue("energyFD")/1.e+18;
   entemp[1] = errneg->GetValue("energyFD")/1.e+18;
   entemp[2] = errpos->GetValue("energyFD")/1.e+18;

//   cout << "E = " << entemp[0] << ", " << entemp[1] << ", " << entemp[2] << ", S38 = " << s38temp[0] << ", " << s38temp[1] << ", " << s38temp[2]; 

   // A, B
//   cout << ", A = " << fitresults[0] << " ± " << fitresults[1] << ", B = " << fitresults[2] << " ± " << fitresults[3] << endl; 

   // (E/A)^(1/B)
   ftemp[0] = TMath::Power(entemp[0]/fitresults[0], 1./fitresults[2]);
//   cout << "  (E/A)^(1/B) = " << ftemp[0]; 

   // DeltaS38
   if(s38temp[0] == -1.)
   {
      SetValue(type, -1.);
      errneg->SetValue(type, -1.);
      errpos->SetValue(type, -1.);
//      cout << "S38 value is negative for this event, setting DeltaS38 values to -1." << endl;
   }
   else
   {
      if(*isdata)
      {
         if(syst == 0)
            ftemp[1] = s38temp[0] - ftemp[0];
         else if(syst == 1)
            ftemp[1] = (s38temp[0]-0.24) - ftemp[0];
         else if(syst == 2)
            ftemp[1] = (s38temp[0]+0.24) - ftemp[0];
      }
      else
         ftemp[1] = s38temp[0] - ftemp[0];

//      ftemp[1] = s38temp[0] - ftemp[0];
      SetValue(type, ftemp[1]);
//      cout << ", DeltaS38 = " << ftemp[1]; 

      // (dS38)^2
      ftemp[2] = s38temp[1];
      ftemp[3] = TMath::Power(ftemp[2],2);
//      cout << ", S38 err = " << ftemp[2];
      // (-((E/A)^(1/B))*dE/(B*E))^2
      ftemp[2] = ftemp[0]*entemp[1]/(fitresults[2]*entemp[0]);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", E err = " << ftemp[2];
      // (((E/A)^(1/B))*dA/(A*B))^2
      ftemp[2] = ftemp[0]*fitresults[1]/(fitresults[0]*fitresults[2]);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", A err = " << ftemp[2];
      // (((E/A)^(1/B))*ln(E/A)*dB/B^2)^2
      ftemp[2] = ftemp[0]*TMath::Log(entemp[0]/fitresults[0])*fitresults[3]/TMath::Power(fitresults[2],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", B err = " << ftemp[2];
      
      ftemp[3] = TMath::Sqrt(ftemp[3]);
//      cout << ", dDeltaS38 neg = " << ftemp[3] << endl;
      errneg->SetValue(type, ftemp[3]);

      // (dS38)^2
      ftemp[2] = s38temp[2];
      ftemp[3] = TMath::Power(ftemp[2],2);
//      cout << ", S38 err = " << ftemp[2];
      // (-((E/A)^(1/B))*dE/(B*E))^2
      ftemp[2] = ftemp[0]*entemp[2]/(fitresults[2]*entemp[0]);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", E err = " << ftemp[2];
      // (((E/A)^(1/B))*dA/(A*B))^2
      ftemp[2] = ftemp[0]*fitresults[1]/(fitresults[0]*fitresults[2]);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", A err = " << ftemp[2];
      // (((E/A)^(1/B))*ln(E/A)*dB/B^2)^2
      ftemp[2] = ftemp[0]*TMath::Log(entemp[0]/fitresults[0])*fitresults[3]/TMath::Power(fitresults[2],2);
      ftemp[3] += TMath::Power(ftemp[2],2);
//      cout << ", B err = " << ftemp[2];
      
      ftemp[3] = TMath::Sqrt(ftemp[3]);
//      cout << ", dDeltaS38 pos = " << ftemp[3] << endl;
      errpos->SetValue(type, ftemp[3]);
   }

   delete[] s38temp;
   delete[] entemp;
   delete[] ftemp;
}

// Converting risetime to Delta Risetime
void Observables::ConvertToDelta(int type, Observables *errneg, Observables *errpos, vector<float> **SDdist, vector<float> **SDrise, vector<bool> *SDsat, vector<float> *fitresults, int syst, bool *isdata)
{
//   cout << "# Convert to Delta      #: " << "Converting risetime to Risetime Delta for " << GetName(type) << endl;

   float *ftemp = new float[4];
   float *eventInfo = new float[2];
   int *itemp = new int[2];
   float *tbench = new float[3];
   float *deltaVal = new float[3];

   eventInfo[0] = GetValue("energyFD")/1.e+18;
   eventInfo[1] = SecTheta(GetValue("zenithFD"),false);
//   cout << "Energy = " << eventInfo[0] << ", Sec(Zenith) = " << eventInfo[1]  << ", Zenith = " << GetValue("zenithFD") << endl;

   itemp[0] = -1;
   // Find fitting parameters for selected zenith angle bin
   for(int i = 0; i < fitresults->size()/10; i++)
   {
      ftemp[0] = fitresults->at(10*i+2);
      ftemp[1] = fitresults->at(10*i+3);
      
      if((eventInfo[1] > ftemp[0]) && (eventInfo[1] <= ftemp[1]))
      {
         itemp[0] = i;
         break;
      }
   }

//   cout << "Selected zenith angle bin = " << itemp[0] << endl;

   if(itemp[0] != -1)
   {
      deltaVal[0] = 0;
      deltaVal[1] = 0;
      deltaVal[2] = 0;

      itemp[1] = GetValue("nrstations");

      if(itemp[1] == 0)
      {
         SetValue(type, -1.);
         errneg->SetValue(type, -1.);
         errpos->SetValue(type, -1.);
//	 cout << "No stations in this event, setting risetime values to -1." << endl;
      }
      else
      {
         for(int i = 0; i < itemp[1]; i++)
         {
//            cout << "  Converting to delta for station " << i+1 << "/" << itemp[1] << endl;
//            cout << "    distance = " << SDdist[0]->at(i) << ", " << SDdist[1]->at(i) << ", " << SDdist[2]->at(i) << endl;
//            cout << "    risetime = " << SDrise[0]->at(i) << ", " << SDrise[1]->at(i) << ", " << SDrise[2]->at(i) << endl;
//            cout << "    saturation = " << SDsat->at(i) << endl;

            // Calculate the benchmarking functions (high-gain saturated and non-saturated) for event energy
            //   tbench = 40 + sqrt(A^2 + B*r^2) - A
            if(SDsat->at(i))
            {
               tbench[0] = 40. + TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)) - fitresults->at(10*itemp[0]+4);
//               cout << "  tbench[0] = " << tbench[0];
               // (dA)^2 = (dA*(A/sqrt(A^2 + B*r^2) - 1))^2
               ftemp[0] = fitresults->at(10*itemp[0]+5)*(fitresults->at(10*itemp[0]+4)/TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)) - 1);
               ftemp[1] = TMath::Power(ftemp[0],2);
               ftemp[3] = TMath::Power(ftemp[0],2);
//               cout << ", dA = " << ftemp[0];
               // (dB)^2 = (dB*r^2/(2*sqrt(A^2 + B*r^2)))^2
               ftemp[0] = fitresults->at(10*itemp[0]+7)*TMath::Power(SDdist[0]->at(i),2)/(2*TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)));
               ftemp[1] += TMath::Power(ftemp[0],2);
               ftemp[3] += TMath::Power(ftemp[0],2);
//               cout << ", dB = " << ftemp[0];
               // (dr)^2 = (dr*B*r/(sqrt(A^2 + B*r^2)))^2
               ftemp[0] = SDdist[1]->at(i)*fitresults->at(10*itemp[0]+6)*SDdist[0]->at(i)/(TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)));
               ftemp[1] += TMath::Power(ftemp[0],2);
               ftemp[2] = SDdist[2]->at(i)*fitresults->at(10*itemp[0]+6)*SDdist[0]->at(i)/(TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)));
               ftemp[3] += TMath::Power(ftemp[2],2);
//               cout << ", dr = (" << ftemp[0] << ", " << ftemp[2] << ")";
            }
            //   tbench = 40 + N*(sqrt(A^2 + B*r^2) - A)
            else
            {
               tbench[0] = 40. + fitresults->at(10*itemp[0]+8)*(TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)) - fitresults->at(10*itemp[0]+4));
//               cout << "  tbench[0] = " << tbench[0];
               // (dA)^2 = (N*dA*(A/sqrt(A^2 + B*r^2) - 1))^2
               ftemp[0] = fitresults->at(10*itemp[0]+8)*fitresults->at(10*itemp[0]+5)*(fitresults->at(10*itemp[0]+4)/TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)) - 1);
               ftemp[1] = TMath::Power(ftemp[0],2);
               ftemp[3] = TMath::Power(ftemp[0],2);
//               cout << ", dA = " << ftemp[0];
               // (dB)^2 = (N*dB*r^2/(2*sqrt(A^2 + B*r^2)))^2
               ftemp[0] = fitresults->at(10*itemp[0]+8)*fitresults->at(10*itemp[0]+7)*TMath::Power(SDdist[0]->at(i),2)/(2*TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)));
               ftemp[1] += TMath::Power(ftemp[0],2);
               ftemp[3] += TMath::Power(ftemp[0],2);
//               cout << ", dB = " << ftemp[0];
               // (dN)^2 = (dN*(sqrt(A^2 + B*r^2) - A))^2
               ftemp[0] = fitresults->at(10*itemp[0]+9)*(TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)) - fitresults->at(10*itemp[0]+4));
               ftemp[1] += TMath::Power(ftemp[0],2);
               ftemp[3] += TMath::Power(ftemp[0],2);
//               cout << ", dN = " << ftemp[0];
               // (dr)^2 = (dr*B*N*r/(sqrt(A^2 + B*r^2)))^2
               ftemp[0] = SDdist[1]->at(i)*fitresults->at(10*itemp[0]+8)*fitresults->at(10*itemp[0]+6)*SDdist[0]->at(i)/(TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)));
               ftemp[1] += TMath::Power(ftemp[0],2);
               ftemp[2] = SDdist[2]->at(i)*fitresults->at(10*itemp[0]+8)*fitresults->at(10*itemp[0]+6)*SDdist[0]->at(i)/(TMath::Sqrt(TMath::Power(fitresults->at(10*itemp[0]+4),2) + fitresults->at(10*itemp[0]+6)*TMath::Power(SDdist[0]->at(i),2)));
               ftemp[3] += TMath::Power(ftemp[2],2);
//               cout << ", dr = (" << ftemp[0] << ", " << ftemp[2] << ")";
            }

            tbench[1] = TMath::Sqrt(ftemp[1]);
            tbench[2] = TMath::Sqrt(ftemp[3]);
//            cout << ", tbench[1-2] = (" << tbench[1] << ", " << tbench[2]  << ")" << endl;

	    if(*isdata)
	    {
	       if(syst == 0)
                  ftemp[0] = SDrise[0]->at(i) - tbench[0];
	       else if(syst == 1)
                  ftemp[0] = (SDrise[0]->at(i)-25.0) - tbench[0];
	       else if(syst == 2)
                  ftemp[0] = (SDrise[0]->at(i)+25.0) - tbench[0];
	    }
	    else
               ftemp[0] = SDrise[0]->at(i) - tbench[0];
//          cout << "deltaVal_i = " << ftemp[0];
            deltaVal[0] += ftemp[0];

            ftemp[0] = TMath::Sqrt(TMath::Power(SDrise[1]->at(i),2) + TMath::Power(tbench[1],2));
//          cout << ", " << ftemp[0];
            deltaVal[1] += TMath::Power(ftemp[0],2);

            ftemp[0] = TMath::Sqrt(TMath::Power(SDrise[2]->at(i),2) + TMath::Power(tbench[2],2));
//          cout << ", " << ftemp[0] << endl;
            deltaVal[2] += TMath::Power(ftemp[0],2);
         }

         deltaVal[0] = deltaVal[0]/itemp[1];
         deltaVal[1] = TMath::Sqrt(deltaVal[1])/itemp[1];
         deltaVal[2] = TMath::Sqrt(deltaVal[2])/itemp[1];

         SetValue(type, deltaVal[0]);
         errneg->SetValue(type, deltaVal[1]);
         errpos->SetValue(type, deltaVal[2]);
//         cout << "Final Delta = " << deltaVal[0] << ", " << deltaVal[1] << ", " << deltaVal[2] << endl;
      }
   }
   else
   {
      SetValue(type, -1.);
      errneg->SetValue(type, -1.);
      errpos->SetValue(type, -1.);
//      cout << "Event is not valid for reconstruction (risetime Delta)" << endl;
      delete[] deltaVal;
      delete[] tbench;
      delete[] itemp;
      delete[] eventInfo;
      delete[] ftemp;
      return;
   }

   delete[] deltaVal;
   delete[] tbench;
   delete[] itemp;
   delete[] eventInfo;
   delete[] ftemp;
}
