//#include "frame.h"
#include "separate_functions.h"
#include <fstream>
#include "mva_result_read.h"
//#include <algorithm>

ResultRead::ResultRead()
{
   ebin = new float[2];
   /* NEWREMOVE - TODO
   valtype = new int[3];
   mvacut = new float[3];
   nrtrees = new int[3];
   */
   mvacut = new float;
   nrtrees = new int;

   itemp = new int[6];
   stemp = new string[2];
   ftemp = new float[6];
   fraction = new float[3];
}

ResultRead::~ResultRead()
{
   delete[] ebin;
   /* NEWREMOVE - TODO
   delete[] valtype;
   delete[] mvacut;
   delete[] nrtrees;
   */
   delete mvacut;
   delete nrtrees;
   delete[] itemp;
   delete[] stemp;
   delete[] ftemp;
   delete[] fraction;
}

void ResultRead::ZeroVectors()
{
   treeType.clear();
   allEvents.clear();
   siglikeEvents.clear();
   bgdlikeEvents.clear();
   treeName.clear();
}

int ResultRead::ReadFile(string inname)
{
   ZeroVectors();
   fraction[0] = -1;
   fraction[1] = -1;
   fraction[2] = -1;

   char *ctemp = new char[1024];
   ifstream *ifile = new ifstream;
   ifile->open(inname.c_str(), ifstream::in);

   filename = inname;

   if(ifile->is_open())
   {
      *ifile >> ebin[0] >> ebin[1];

      /* NEWREMOVE - TODO
      for(int i = 0; i < 3; i++)
      {
      */
         /* NEWREMOVE - TODO
         *ifile >> valtype[i] >> mvacut[i] >> nrtrees[i];
	 */
         *ifile >> *mvacut >> *nrtrees;

         for(int j = 0; j < *nrtrees; j++)
         {
            *ifile >> itemp[0] >> itemp[1] >> itemp[2] >> itemp[3];
	    ifile->getline(ctemp, 1024, '\n');
	    stemp[0] = string(ctemp);
	    RemoveLeadingSpaces(&stemp[0]);

            treeType.push_back(itemp[0]);
            allEvents.push_back((float)itemp[1]);
            siglikeEvents.push_back((float)itemp[2]);
            bgdlikeEvents.push_back((float)itemp[3]);
            treeName.push_back(stemp[0]);
         }
      /* NEWREMOVE - TODO
      }
      */

      ifile->close();
      delete ifile;
      delete[] ctemp;
      return 0;
   }
   else
   {
      delete ifile;
      delete[] ctemp;
      return 1;
   }
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

float ResultRead::GetLowEnergy()
{
   return TMath::Log10(ebin[0]);
}

float ResultRead::GetHighEnergy()
{
   return TMath::Log10(ebin[1]);
}

float ResultRead::GetFraction(int sigbackdata, float norm)
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

//   cout << "Selected value = " << itemp[0] << endl;

   /* NEWREMOVE - TODO
   // Find the number of all signal, background and data trees per type (mean, negative, positive)
   itemp[1] = treeType.size()/3;
   */
   itemp[1] = *nrtrees;
//   cout << "Number of trees = " << itemp[1] << endl;

   // Some trees are missing
   /* NEWREMOVE - TODO
   if(treeType.size() < 3)
   */
   if(*nrtrees < 1)
   {
      cout << "Error! Some trees seem to be missing. Please check file " << filename << " for errors." << endl;
      return -1;
   }

   int nrbacktrees = 0;

   // Loop over signal, background and data values (only select values that are selected by sigbackdata values
   /* NEWREMOVE - TODO
   for(int i = 0; i < treeType.size(); i++)
   */
   for(int i = 0; i < *nrtrees; i++)
   {
      // Zero values at beginning
      if(i == 0)
      {
         ftemp[0] = 0;
         ftemp[1] = 0;
         ftemp[2] = 0;
         ftemp[3] = 0;
         ftemp[4] = 0;
         ftemp[5] = 0;
      }

      if(treeType[i] == itemp[0])
      {
         // No normalizing of data
         if(norm == -1)
	 {
	    // Background events (save background fraction)
            if(itemp[0] == 2)
	    {
               ftemp[0] += bgdlikeEvents[i];
               ftemp[1] += bgdlikeEvents[i+itemp[1]];
               ftemp[2] += bgdlikeEvents[i+2*itemp[1]];
               ftemp[3] += allEvents[i];
               ftemp[4] += allEvents[i+itemp[1]];
               ftemp[5] += allEvents[i+2*itemp[1]];
	       nrbacktrees++;
//               cout << i << ": Background mean =\t" << ftemp[0] << endl;
//               cout << i << ": Background neg =\t" << ftemp[1] << endl;
//               cout << i << ": Background pos =\t" << ftemp[2] << endl;
	    }
	    // Signal events (save signal fraction)
	    else if( (itemp[0] == 1) || (itemp[0] == 3) )
	    {
               ftemp[0] += siglikeEvents[i]/allEvents[i];
               ftemp[1] += siglikeEvents[i+itemp[1]]/allEvents[i+itemp[1]];
               ftemp[2] += siglikeEvents[i+2*itemp[1]]/allEvents[i+2*itemp[1]];
//               cout << i << ": Sig mean =\t" << ftemp[0] << endl;
//               cout << i << ": Sig neg =\t" << ftemp[1] << endl;
//               cout << i << ": Sig pos =\t" << ftemp[2] << endl;
	       break;
	    }
	 }
         // Normalizing of data
         else if(norm > -1)
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
//	       cout << k << ": " << itemp[2] << ", " << itemp[3] << ", " << itemp[4] << endl;
	       // Fraction of signal-like data events - Fraction of wrongly classified signal events
	       if(itemp[0] == 1)	// signal-like signal events
	       {
		  /* NEWREMOVE - TODO
	          itemp[5] = FindPos(1, k);
		  */
	          itemp[5] = FindPos(1);
	          ftemp[4] = (siglikeEvents[itemp[5]]/allEvents[itemp[5]]);
	       }
	       else if(itemp[0] == 2)	// signal-like background events
	       {
	          ftemp[4] = itemp[3]/itemp[2];
	       }
	       else if(itemp[0] == 3)	// signal-like data events
	       {
		  /* NEWREMOVE - TODO
	          itemp[5] = FindPos(3, k);
		  */
	          itemp[5] = FindPos(3);
	          ftemp[4] = (siglikeEvents[itemp[5]]/allEvents[itemp[5]]);
	       }
//	       cout << "Data tree is at: " << itemp[5] << endl;
	       /* NEWREMOVE - TODO
	       itemp[5] = FindPos(1, k);
	       */
	       itemp[5] = FindPos(1);
//	       cout << "Signal tree is at: " << itemp[5] << endl;
	       ftemp[4] -= (bgdlikeEvents[itemp[5]]/allEvents[itemp[5]]);
//	       cout << k << ": f4 = " << ftemp[4] << endl;
	       // Fraction of signal-like data events - Fraction of wrongly classified background events
	       if(itemp[0] == 1)	// signal-like signal events
	       {
		  /* NEWREMOVE - TODO
	          itemp[5] = FindPos(1, k);
		  */
	          itemp[5] = FindPos(1);
	          ftemp[5] = (siglikeEvents[itemp[5]]/allEvents[itemp[5]]) - ((float)itemp[3]/(float)itemp[2]);
	       }
	       else if(itemp[0] == 2)	// signal-like background events
	       {
	          ftemp[5] = itemp[3]/itemp[2];
	       }
	       else if(itemp[0] == 3)	// signal-like data events
	       {
		  /* NEWREMOVE - TODO
	          itemp[5] = FindPos(3, k);
		  */
	          itemp[5] = FindPos(3);
	          ftemp[5] = (siglikeEvents[itemp[5]]/allEvents[itemp[5]]) - ((float)itemp[3]/(float)itemp[2]);
	       }
//	       cout << "Data tree is at: " << itemp[5] << endl;
//	       cout << k << ": f5 = " << ftemp[5] << endl;
	       // Normalized value
	       /* NEWREMOVE - TODO
	       itemp[5] = FindPos(1, k);
	       */
	       itemp[5] = FindPos(1);
//	       cout << "Signal tree is at: " << itemp[5] << endl;

	       if(itemp[0] == 1)	// signal-like signal events
	          ftemp[k] = (allEvents[itemp[5]]*ftemp[4]/(siglikeEvents[itemp[5]] - bgdlikeEvents[itemp[5]]));
	       else if(itemp[0] == 2)	// background-like background events
	          ftemp[k] = (1. - (float)itemp[2]*ftemp[5]/((float)itemp[4] - (float)itemp[3]));
	       else if(itemp[0] == 3)	// signal-like data events
	          ftemp[k] = ( (norm*allEvents[itemp[5]]*ftemp[4]/(siglikeEvents[itemp[5]] - bgdlikeEvents[itemp[5]])) + ((1.-norm)*(float)itemp[2]*ftemp[5]/((float)itemp[4] - (float)itemp[3])) );
	    }
	 }
      }

      if(i == itemp[1]-1)
         break;
   }

   if(sigbackdata == 0)
   {
      if(norm == -1)
      {
         fraction[0] = ftemp[0]/ftemp[3];
         fraction[1] = TMath::Abs(ftemp[0]/ftemp[3] - ftemp[1]/ftemp[4]);
         fraction[2] = TMath::Abs(ftemp[0]/ftemp[3] - ftemp[2]/ftemp[5]);
      }
   }
   else
   {
      fraction[0] = ftemp[0];
      fraction[1] = TMath::Abs(ftemp[0] - ftemp[2]);
      fraction[2] = TMath::Abs(ftemp[0] - ftemp[1]);
   }

//   cout << "Fraction mean =\t" << fraction[0] << endl;
//   cout << "Fraction neg =\t" << fraction[1] << endl;
//   cout << "Fraction pos =\t" << fraction[2] << endl;

   return fraction[0];
}

/* NEWREMOVE -TODO
// Same as above, but with no normalization and to make it possible to select any tree
float ResultRead::GetFraction(int tree)
{
   for(int i = 0; i < treeType.size(); i++)
   {
	cout << i << ", treeType = " << treeType.at(i) << ", allEvents = " << allEvents.at(i) << ", siglikeEvents = " << siglikeEvents.at(i) << ", bgdlikeEvents = " << bgdlikeEvents.at(i) << ", treeName = " << treeName.at(i) << endl;
   }

   return 0.;
}
*/

void ResultRead::GetFractionError(float *err)
{
   if( (fraction[0] == -1) || (fraction[1] == -1) || (fraction[2] == -1) )
      cout << "Values for fractions have not been set yet. Please run function ResultRead::GetFraction first." << endl;

   err[0] = fraction[1];
   err[1] = fraction[2];
}

/* NEWREMOVE - TODO
float ResultRead::GetMvaCut(int type)
{
   if(type == 0)
      return mvacut[0];
   else if(type == -1)
      return mvacut[1];
   else if(type == 1)
      return mvacut[2];
}
*/
float ResultRead::GetMvaCut()
{
   return *mvacut;
}

/* NEWREMOVE - TODO
int ResultRead::GetNrTrees(int type)
{
   if(type == 0)
      return nrtrees[0];
   else if(type == -1)
      return nrtrees[1];
   else if(type == 1)
      return nrtrees[2];
}
*/
int ResultRead::GetNrTrees()
{
   return *nrtrees;
}

void ResultRead::PrintVectors(int output)
{
   cout << "Reading information:" << endl;
   /* NEWREMOVE - TODO
   PrintVectors(0, output);
   PrintVectors(1, output);
   PrintVectors(2, output);
   */

   if(output == 1) // printout to stdout
   {
      cout << "- MVA cut (" << GetMvaCut() << ")" << endl;
      for(int i = 0; i < *nrtrees; i++)
      {
         cout << "  [" << i << "]: " << treeType[i] << "\t";
         cout << allEvents[i] << "\t";
         cout << siglikeEvents[i] << "\t";
         cout << bgdlikeEvents[i] << "\t";
         cout << treeName[i] << endl;
      }
   }
   if(output == 0) // printout to stderr
   {
      cerr << "- MVA cut (" << GetMvaCut() << ")" << endl;
      for(int i = 0; i < *nrtrees; i++)
      {
         cerr << "  [" << i << "]: " << treeType[i] << "\t";
         cerr << allEvents[i] << "\t";
         cerr << siglikeEvents[i] << "\t";
         cerr << bgdlikeEvents[i] << "\t";
         cerr << treeName[i] << endl;
      }
   }
}

/* NEWREMOVE - TODO
void ResultRead::PrintVectors(int type, int output)
{
   if(output == 1) // printout to stdout
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
   if(output == 0) // printout to stderr
   {
      for(int j = 0; j < 3; j++)
      {
         if(j == type)
         {
            if(j == 0)
               cerr << "- Mean cut (" << GetMvaCut(0) << ")" << endl;
            else if(j == 1)
               cerr << "- Negative cut (" << GetMvaCut(-1) << ")" << endl;
            else if(j == 2)
               cerr << "- Positive cut (" << GetMvaCut(1) << ")" << endl;

            itemp[0] = treeType.size()/3;
            for(int i = j*itemp[0]; i < (j+1)*itemp[0]; i++)
            {
               cerr << "  [" << i << "]: " << treeType[i] << "\t";
               cerr << allEvents[i] << "\t";
               cerr << siglikeEvents[i] << "\t";
               cerr << bgdlikeEvents[i] << "\t";
               cerr << treeName[i] << endl;
            }
         }
      }
   }
}
*/

/* NEWREMOVE - TODO
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
*/

// For background it only works when we have one background tree
//   sigbackdata: 1 = signal, 2 = background, 3 = data
//   type: 0 = mean, 1 = negerror, 2 = poserror
int ResultRead::FindPos(int sigbackdata)
{
   for(int i = 0; i < *nrtrees; i++)
   {
      if(sigbackdata == treeType[i])
         return i;
   }

   return -1;
}

// Find the tree position and save it to output vector
//   sigbackdata: 1 = signal, 2 = background, 3 = data
//   type: 0 = mean, 1 = negerror, 2 = poserror
void ResultRead::FindPos(int sigbackdata, vector<int> *out)
{
   for(int i = 0; i < *nrtrees; i++)
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

// Get the tree type (1 = signal, 2 = background, 3 = data, 0 = not used in analysis)
int ResultRead::GetTreeType(int nr)
{
   return treeType[nr];
}

/* NEWREMOVE - TODO
// Get the type of the file (0 = individual observable analysis, 1 = mva analysis)
int ResultRead::GetFileType()
{
   itemp[0] = filename.find("individual_results");
   if(itemp[0] != string::npos)
      return 0;

   itemp[0] = filename.find("application_results");
   if(itemp[0] != string::npos)
      return 1;

   return -1;
}

// Get observable type from individual results file
string ResultRead::GetObservableType()
{
   stemp[0] = RemovePath(&filename);
   stemp[1] = RemoveExtension(&stemp[0]);
   stemp[0] = "individual_results_";
   stemp[1].erase(0, (size_t)stemp[0].length()); 
   if(stemp[1].length() > 0)
      return stemp[1];
   else
      return "mva-analysis";
}
*/
