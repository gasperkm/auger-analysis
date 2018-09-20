#include "observables.h"
#include "separate_functions.h"
#include "root_include.h"

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
}

Observables::~Observables()
{
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

float Observables::GetValue(int obs, int eye)
{
   return obsstruct[obs].value[eye];
}

float Observables::GetValue(string obsname, int eye)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         return GetValue(i, eye);
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

void Observables::SetValue(int obs, float val, int eye)
{
   obsstruct[obs].value[eye] = val;
}

void Observables::SetValue(string obsname, float val, int eye)
{
   for(int i = 0; i < nrobs; i++)
   {
      if(obsname == obsstruct[i].name)
         SetValue(i, val, eye);
   }
}

void Observables::Zero(int obs)
{
   for(int i = 0; i < ALLEYES; i++)
      obsstruct[i].value[i] = 0;
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

   for(int i = 0; i < ALLEYES; i++)
   {
      ftemp[0] = GetValue("energyFD", i);
      ftemp[1] = GetValue("xmax", i);
      // Only apply correction, if we have a valid value
      if( (ftemp[0] != -1) && (ftemp[1] != -1) )
      {
         *logE = TMath::Log10(ftemp[0]); 
	 ftemp[2] = 6.5/(TMath::Exp((*logE - 18.23)/0.41) + 1.);
	 ftemp[3] = 0.93*(*logE - 18.);
	 ftemp[1] = ftemp[1] - ftemp[2] + 3.4 - ftemp[3];

         SetValue("xmax", ftemp[1], i);
      }
   }

   delete logE;
   delete[] ftemp;
}

// Applying Xmax and energy corrections to Auger HECO data
void Observables::ApplyCorrectionHECO()
{
   float *logE = new float;
   float *ftemp = new float[5];

   for(int i = 0; i < ALLEYES; i++)
   {
      ftemp[0] = GetValue("energyFD", i);
      ftemp[1] = GetValue("xmax", i);
      // Only apply correction, if we have a valid value
      if( (ftemp[0] != -1) && (ftemp[1] != -1) )
      {
         *logE = TMath::Log10(ftemp[0]); 
	 ftemp[2] = TMath::Log(*logE - 16.5);
	 ftemp[0] = (ftemp[0])*(1 + 0.04536243 - 0.08317193*ftemp[2]);
         SetValue("energyFD", ftemp[0], i);

	 ftemp[0] = GetValue("energyFD", i);
         *logE = TMath::Log10(ftemp[0]); 
	 ftemp[2] = -2.96976097 - 0.99199218*(*logE - 17.80937342);
	 ftemp[3] = 6.5/(TMath::Exp((*logE - 18.23)/0.41) + 1.);
	 if(*logE < 17.55)
            ftemp[4] = 0.12 - 6.43*(*logE - 17.55);
	 else
            ftemp[4] = 0.12 - 0.27*(*logE - 17.55);
         ftemp[1] = ftemp[1] - (ftemp[2] + ftemp[3] + 0.5*ftemp[4]);
         SetValue("xmax", ftemp[1], i);
      }
   }

   delete logE;
   delete[] ftemp;
}

// Applying Xmax and energy corrections to Auger HECO data (Errors only)
void Observables::ApplyCorrectionHECOErrors(Observables *mean, int type)
{
   float *logE = new float;
   float *ftemp = new float[3];

   for(int i = 0; i < ALLEYES; i++)
   {
      ftemp[0] = mean->GetValue("energyFD", i);
      ftemp[1] = GetValue("energyFD", i);
      // Only apply correction, if we have a valid value
      if( (ftemp[0] != -1) && (ftemp[1] != -1) )
      {
         *logE = TMath::Log10(ftemp[0]); 
	 ftemp[2] = TMath::Log(*logE - 16.5);
	 ftemp[1] = (ftemp[1])*(1 + 0.04536243 - 0.08317193*ftemp[2]);
         SetValue("energyFD", ftemp[1], i);
      }
   }

   delete logE;
   delete[] ftemp;
}

// Applying Xmax and energy corrections to Auger HECO data
void Observables::ApplyUncertainty(Observables *errneg, Observables *errpos, int type)
{
   if((type == 0) || (type == 1))
   {
//      cout << "Changing value (" << type << "):" << endl;
      float *ftemp = new float[3];

      for(int i = 0; i < nrobs; i++)
      {
         // TESTING!
	 if(GetName(i) == "xmax")
	 {
            for(int j = 0; j < ALLEYES; j++)
	    {
               ftemp[0] = GetValue(i, j);
               ftemp[1] = errneg->GetValue(i, j);
               ftemp[2] = errpos->GetValue(i, j);
//               ftemp[1] = 10.;
//               ftemp[2] = 10.;
//	       cout << i << "\teye" << j << ", " << GetName(i) << "\t" << ftemp[0] << "\t" << ftemp[1] << "\t" << ftemp[2] << "\tnewVal = ";

	       if(ftemp[0] != -1)
	       {
	          if(type == 0)
                     ftemp[0] -= ftemp[1];

	          if(type == 1)
                     ftemp[0] += ftemp[2];

	          SetValue(i, ftemp[0], j);
	       }

//	       cout << GetValue(i, j) << endl;
	    }
	 }
      }
//      cout << endl;

      delete[] ftemp;
   }
}

// Setting up zenith as sec(theta) for analysis
void Observables::SetupZenith(int type, Observables *errneg, Observables *errpos)
{
   float *ftemp = new float[4];

   for(int i = 0; i < ALLEYES; i++)
   {
      ftemp[0] = GetValue(type, i);
      ftemp[1] = errneg->GetValue(type, i);
      ftemp[2] = errpos->GetValue(type, i);

      // Only apply calculation, if we have a valid value
      if( (ftemp[0] != -1) && (ftemp[1] != -1) && (ftemp[2] != -1) )
      {
//         cout << "theta = " << RadToDeg(ftemp[0]) << ", " << ftemp[0] << " (" << ftemp[1] << ", " << ftemp[2] << ")" << ", sec(theta) = " << SecTheta(ftemp[0], false) << endl;

	 ftemp[3] = tan(ftemp[0])/cos(ftemp[0]);

//         cout << "sin(theta)/cos2(theta) = " << ftemp[3] << ", dsec(theta) = " << ftemp[1]*ftemp[3] << endl;

	 ftemp[0] = SecTheta(ftemp[0], false);
	 ftemp[1] = ftemp[1]*ftemp[3];
	 ftemp[2] = ftemp[2]*ftemp[3];

	 SetValue(type, ftemp[0], i);
	 errneg->SetValue(type, ftemp[1], i);
	 errpos->SetValue(type, ftemp[2], i);

//         cout << "final = " << ftemp[0] << " (" << ftemp[1] << ", " << ftemp[2] << ")" << endl;

      }
   }

   delete[] ftemp;

/*   float *logE = new float;
   float *ftemp = new float[3];

   for(int i = 0; i < ALLEYES; i++)
   {
      ftemp[0] = mean->GetValue("energyFD", i);
      ftemp[1] = GetValue("energyFD", i);
      // Only apply correction, if we have a valid value
      if( (ftemp[0] != -1) && (ftemp[1] != -1) )
      {
         *logE = TMath::Log10(ftemp[0]); 
	 ftemp[2] = TMath::Log(*logE - 16.5);
	 ftemp[1] = (ftemp[1])*(1 + 0.04536243 - 0.08317193*ftemp[2]);
         SetValue("energyFD", ftemp[1], i);
      }
   }

   delete logE;
   delete[] ftemp;*/
}
