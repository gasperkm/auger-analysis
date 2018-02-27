//#include "frame.h"
#include "separate_functions.h"
#include <fstream>
#include "mva_result_read.h"
//#include <algorithm>

ResultRead::ResultRead()
{
   ebin = new float[2];
   valtype = new int[3];
   mvacut = new float[3];
   nrtrees = new int[3];

   itemp = new int[6];
   stemp = new string;
   ftemp = new float[5];
   fraction = new float[3];
}

ResultRead::~ResultRead()
{
   delete[] ebin;
   delete[] valtype;
   delete[] mvacut;
   delete[] nrtrees;
   delete[] itemp;
   delete stemp;
   delete[] ftemp;
   delete[] fraction;
}

void ResultRead::ZeroVectors()
{
   if(!treeType.empty())
      treeType.erase(treeType.begin(), treeType.end());
   if(!allEvents.empty())
      allEvents.erase(allEvents.begin(), allEvents.end());
   if(!siglikeEvents.empty())
      siglikeEvents.erase(siglikeEvents.begin(), siglikeEvents.end());
   if(!bgdlikeEvents.empty())
      bgdlikeEvents.erase(bgdlikeEvents.begin(), bgdlikeEvents.end());
   if(!treeName.empty())
      treeName.erase(treeName.begin(), treeName.end());
}

int ResultRead::ReadFile(string inname)
{
   ZeroVectors();
   fraction[0] = -1;
   fraction[1] = -1;
   fraction[2] = -1;

   ifstream ifile;
   ifile.open(inname.c_str(), ifstream::in);

   filename = inname;

   if(ifile.is_open())
   {
      ifile >> ebin[0] >> ebin[1];

      for(int i = 0; i < 3; i++)
      {
         ifile >> valtype[i] >> mvacut[i] >> nrtrees[i];

         for(int j = 0; j < *nrtrees; j++)
         {
            ifile >> itemp[0] >> itemp[1] >> itemp[2] >> itemp[3] >> *stemp;

/*            if(itemp[0] > 0)
            {*/
               treeType.push_back(itemp[0]);
               allEvents.push_back((float)itemp[1]);
               siglikeEvents.push_back((float)itemp[2]);
               bgdlikeEvents.push_back((float)itemp[3]);
               treeName.push_back(*stemp);
//            }
         }
      }

      ifile.close();
      return 0;
   }
   else
      return 1;
}

float ResultRead::GetEnergy()
{
   return TMath::Log10(ebin[0]) + (TMath::Log10(ebin[1]) - TMath::Log10(ebin[0]))/2.;
}

void ResultRead::GetEnergyError(float *err)
{
   err[0] = (TMath::Log10(ebin[1]) - TMath::Log10(ebin[0]))/2.1;
   err[1] = (TMath::Log10(ebin[1]) - TMath::Log10(ebin[0]))/2.1;
}

float ResultRead::GetFraction(int sigbackdata, float rebal)
{
   // User selects background
   if(sigbackdata == 0)
      itemp[0] = 2;
   // User selects signal
   else if(sigbackdata == 1)
      itemp[0] = 1;
   // User selects data
   else if(sigbackdata == 2)
      itemp[0] = 3;

   cout << "Selected value = " << itemp[0] << endl;

   // Find the number of all signal, background and data trees per type (mean, negative, positive)
   itemp[1] = treeType.size()/3;
   cout << "Number of trees = " << itemp[1] << endl;

   // Some trees are missing
   if(treeType.size() < 9)
   {
      cout << "Error! Some trees seem to be missing. Please check file " << filename << " for errors." << endl;
      return -1;
   }

   // Loop over signal, background and data values (only select values that are selected by sigbackdata values
   for(int i = 0; i < treeType.size(); i++)
   {
      // Zero values at beginning
      if(i == 0)
      {
         ftemp[0] = 0;
         ftemp[1] = 0;
         ftemp[2] = 0;
      }

      if(treeType[i] == itemp[0])
      {
	 // Background events (save background fraction)
         if(itemp[0] == 2)
	 {
            ftemp[0] += 100.*bgdlikeEvents[i]/allEvents[i];
            ftemp[1] += 100.*bgdlikeEvents[i+itemp[1]]/allEvents[i+itemp[1]];
            ftemp[2] += 100.*bgdlikeEvents[i+2*itemp[1]]/allEvents[i+2*itemp[1]];
            cout << i << ": Background mean =\t" << ftemp[0] << endl;
            cout << i << ": Background neg =\t" << ftemp[1] << endl;
            cout << i << ": Background pos =\t" << ftemp[2] << endl;
	 }
	 // Signal events (save signal fraction)
	 else if(itemp[0] == 1)
	 {
            ftemp[0] += 100.*siglikeEvents[i]/allEvents[i];
            ftemp[1] += 100.*siglikeEvents[i+itemp[1]]/allEvents[i+itemp[1]];
            ftemp[2] += 100.*siglikeEvents[i+2*itemp[1]]/allEvents[i+2*itemp[1]];
            cout << i << ": Sig mean =\t" << ftemp[0] << endl;
            cout << i << ": Sig neg =\t" << ftemp[1] << endl;
            cout << i << ": Sig pos =\t" << ftemp[2] << endl;
	    break;
	 }
	 // Data events (save signal fraction)
	 else if(itemp[0] == 3)
	 {
            // No rebalancing of data
            if(rebal == -1)
	    {
               ftemp[0] += 100.*siglikeEvents[i]/allEvents[i];
               ftemp[1] += 100.*siglikeEvents[i+itemp[1]]/allEvents[i+itemp[1]];
               ftemp[2] += 100.*siglikeEvents[i+2*itemp[1]]/allEvents[i+2*itemp[1]];
               cout << i << ": Data mean =\t" << ftemp[0] << endl;
               cout << i << ": Data neg =\t" << ftemp[1] << endl;
               cout << i << ": Data pos =\t" << ftemp[2] << endl;
	       break;
	    }
            // Rebalancing of data
	    else if(rebal >= 0)
	    {
	       // Mean, negative error, positive error
	       for(int k = 0; k < 3; k++)
	       {
                  // Sum background trees together
	          itemp[2] = 0;
	          itemp[3] = 0;
	          itemp[4] = 0;
                  for(int j = k*itemp[1]; j < (k+1)*itemp[1]; j++)
	          {
                     if(treeType[j] == 2)
	             {
                        itemp[2] += allEvents[j];
                        itemp[3] += siglikeEvents[j];
                        itemp[4] += bgdlikeEvents[j];
	             }
	          }
		  cout << k << ": " << itemp[2] << ", " << itemp[3] << ", " << itemp[4] << endl;
	          // Fraction of signal-like data events - Fraction of wrongly classified signal events
		  itemp[5] = FindPos(3, k);
		  cout << "Data tree is at: " << itemp[5] << endl;
	          ftemp[4] = (siglikeEvents[itemp[5]]/allEvents[itemp[5]]);
		  itemp[5] = FindPos(1, k);
		  cout << "Signal tree is at: " << itemp[5] << endl;
		  ftemp[4] -= (bgdlikeEvents[itemp[5]]/allEvents[itemp[5]]);
//	          ftemp[4] = (siglikeEvents[(k+1)*itemp[1]-1]/allEvents[(k+1)*itemp[1]-1]) - (bgdlikeEvents[k*itemp[1]]/allEvents[k*itemp[1]]);
		  cout << k << ": f4 = " << ftemp[4] << endl;
	          // Fraction of signal-like data events - Fraction of wrongly classified background events
		  itemp[5] = FindPos(3, k);
		  cout << "Data tree is at: " << itemp[5] << endl;
	          ftemp[5] = (siglikeEvents[itemp[5]]/allEvents[itemp[5]]) - ((float)itemp[3]/(float)itemp[2]);
		  cout << k << ": f5 = " << ftemp[5] << endl;
	          // Rebalanced value
		  itemp[5] = FindPos(1, k);
		  cout << "Signal tree is at: " << itemp[5] << endl;
	          ftemp[k] = 100.*( (rebal*allEvents[itemp[5]]*ftemp[4]/(siglikeEvents[itemp[5]] - bgdlikeEvents[itemp[5]])) + ((1.-rebal)*(float)itemp[2]*ftemp[5]/((float)itemp[4] - (float)itemp[3])) );
	       }

               cout << i << ": Data mean =\t" << ftemp[0] << " (rebalanced)" << endl;
               cout << i << ": Data neg =\t" << ftemp[1] << " (rebalanced)" << endl;
               cout << i << ": Data pos =\t" << ftemp[2] << " (rebalanced)" << endl;
	    }
	 }
      }

      if(i == itemp[1]-1)
         break;
   }

   fraction[0] = ftemp[0];
   fraction[1] = TMath::Abs(ftemp[0] - ftemp[2]);
   fraction[2] = TMath::Abs(ftemp[0] - ftemp[1]);

   cout << "Fraction mean =\t" << fraction[0] << endl;
   cout << "Fraction neg =\t" << fraction[1] << endl;
   cout << "Fraction pos =\t" << fraction[2] << endl;

   return fraction[0];
}

void ResultRead::GetFractionError(float *err)
{
   if( (fraction[0] == -1) || (fraction[1] == -1) || (fraction[2] == -1) )
      cout << "Values for fractions have not been set yet. Please run function ResultRead::GetFraction first." << endl;

   err[0] = fraction[1];
   err[1] = fraction[2];
}

float ResultRead::GetMvaCut(int type)
{
   if(type == 0)
      return mvacut[0];
   else if(type == -1)
      return mvacut[1];
   else if(type == 1)
      return mvacut[2];
}

int ResultRead::GetNrTrees(int type)
{
   if(type == 0)
      return nrtrees[0];
   else if(type == -1)
      return nrtrees[1];
   else if(type == 1)
      return nrtrees[2];
}

void ResultRead::PrintVectors()
{
   cout << "Reading information:" << endl;
   PrintVectors(0);
   PrintVectors(1);
   PrintVectors(2);
}

void ResultRead::PrintVectors(int type)
{
   for(int j = 0; j < 3; j++)
   {
      if(j == type)
      {
         if(j == 0)
            cout << "- Mean cut (" << GetMvaCut(0) << ")" << endl;
         else if(j == 1)
            cout << "- Negative cut (" << GetMvaCut(-1) << ")" << endl;
         else if(j == 2)
            cout << "- Positive cut (" << GetMvaCut(1) << ")" << endl;

         itemp[0] = treeType.size()/3;
         for(int i = j*itemp[0]; i < (j+1)*itemp[0]; i++)
         {
            cout << "  [" << i << "]: " << treeType[i] << "\t";
            cout << allEvents[i] << "\t";
            cout << siglikeEvents[i] << "\t";
            cout << bgdlikeEvents[i] << "\t";
            cout << treeName[i] << endl;
         }
      }
   }
}

// For background it only works when we have one background tree
//   sigbackdata: 1 = signal, 2 = background, 3 = data
//   type: 0 = mean, 1 = negerror, 2 = poserror
int ResultRead::FindPos(int sigbackdata, int type)
{
   for(int i = type*(treeType.size()/3); i < (type+1)*(treeType.size()/3); i++)
   {
      if(sigbackdata == treeType[i])
         return i;
   }

   return -1;
}

// Find the tree position and save it to output vector
//   sigbackdata: 1 = signal, 2 = background, 3 = data
//   type: 0 = mean, 1 = negerror, 2 = poserror
void ResultRead::FindPos(int sigbackdata, int type, vector<int> *out)
{
   for(int i = type*(treeType.size()/3); i < (type+1)*(treeType.size()/3); i++)
   {
      if(sigbackdata == treeType[i])
         out->push_back(i);
   }
}

// Get the name of the tree using the identifier for the tree
string ResultRead::GetTreeName(int nr)
{
   return treeName[nr];
}
