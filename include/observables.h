#ifndef _OBSERVABLES_H_
#define _OBSERVABLES_H_

#include "workstation.h"
#include "root_include.h"
#include <vector>
#include <string>
#include <ctime>

using namespace std;

class Observables
{
private:
   float *minimum;
   float *maximum;
   string *desc;
   time_t *timer;
   TRandom *rand;
public:
   Observables(vector<string> obs);
   virtual ~Observables();

   struct ObserStruct
   {
      string name;
      float value;
   };
   int nrobs;
   ObserStruct *obsstruct;

   // Function for getting an observable through name or number
   string GetName(int obs);
   string GetName(string obsname);
   float GetValue(int obs);
   float GetValue(string obsname);
   int GetInt(string obsname);

   // Function for setting an observable through name or number
   void SetValue(int obs, float val);
   void SetValue(string obsname, float val);

   // Zero all observables or just one
   void Zero();
   void Zero(int obs);
   void Zero(string obsname);

   // Get number of observables
   int Count();

   // Functions for setting and retrieving bin limits for observables
   void SetMin(int obs, float value);
   void SetMin(string obsname, float value);
   void SetMax(int obs, float value);
   void SetMax(string obsname, float value);

   float GetMin(int obs);
   float GetMin(string obsname);
   float GetMax(int obs);
   float GetMax(string obsname);

   // Functions for setting and retrieving axis labels for observables
   void SetLabel(int obs, string value);
   void SetLabel(string obsname, string value);

   string GetLabel(int obs);
   string GetLabel(string obsname);

   // Functions for applying specific corrections to observables
   void ApplyCorrectionFD();
   void ApplyCorrectionHECO();
   void ApplyCorrectionHECOErrors(Observables *mean, int type);
   void ApplyUncertainty(Observables *errneg, Observables *errpos, int type);
   void ApplySmearing();
   void SetupZenith(int type, Observables *errneg, Observables *errpos);

   void ConvertToS38(int type, Observables *errneg, Observables *errpos, vector<float> *fitresults);
   void ConvertToDeltaS38(int type, Observables *errneg, Observables *errpos, float *fitresults);
   void ConvertToDelta(int type, Observables *errneg, Observables *errpos, vector<float> **SDdist, vector<float> **SDrise, vector<bool> *SDsat, vector<float> *fitresults);
};

#endif
