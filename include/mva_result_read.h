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
   /* NEWREMOVE - TODO
   int *valtype;
   */
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
   /* NEWREMOVE - TODO
   float GetFraction(int tree);
   */
   void GetFractionError(float *err);
   float GetEnergy();
   void GetEnergyError(float *err);
   float GetLowEnergy();
   float GetHighEnergy();
   /* NEWREMOVE - TODO
   float GetMvaCut(int type);
   */
   float GetMvaCut();
   /* NEWREMOVE - TODO
   int GetNrTrees(int type);
   */
   int GetNrTrees();

   void PrintVectors(int output);
   /* NEWREMOVE - TODO
   void PrintVectors(int type, int output);
   int FindPos(int sigbackdata, int type);
   void FindPos(int sigbackdata, int type, vector<int> *out);
   */
   int FindPos(int sigbackdata);
   void FindPos(int sigbackdata, vector<int> *out);

   string GetTreeName(int nr);
   int GetTreeType(int nr);
   /* NEWREMOVE - TODO
   int GetFileType();
   string GetObservableType();
   */
};

#endif
