#ifndef _MVA_RESULT_READ_H_
#define _MVA_RESULT_READ_H_

#include "workstation.h"
#include "root_include.h"
#include <vector>
#include <string>

using namespace std;

class ResultRead
{
private:
   float *ebin;
   int *valtype;
   float *mvacut;
   int *nrtrees;

   vector<int> treeType;
   vector<float> allEvents;
   vector<float> siglikeEvents;
   vector<float> bgdlikeEvents;
   vector<string> treeName;

   int *itemp;
   string *stemp;
   float *ftemp;

   string filename;
   float *fraction;
public:
   ResultRead();
   virtual ~ResultRead();

   int ReadFile(string inname);
   void ZeroVectors();
   
   // 0 = background, 1 = signal, 2 = data
   float GetFraction(int sigbackdata, float norm);
   void GetFractionError(float *err);
   float GetEnergy();
   void GetEnergyError(float *err);
   float GetLowEnergy();
   float GetHighEnergy();
   float GetMvaCut(int type);
   int GetNrTrees(int type);

   void PrintVectors();
   void PrintVectors(int type);
   int FindPos(int sigbackdata, int type);
   void FindPos(int sigbackdata, int type, vector<int> *out);

   string GetTreeName(int nr);
   int GetFileType();
   string GetObservableType();
};

#endif
