#include "separate_functions.h"
#include "workstation.h"

// Convert Double to string
string ToString(double value, int prec)
{
   ostringstream ss;
   ss << fixed << setprecision(prec) << value;
   return ss.str();
}

// Convert Integer to string
string ToString(int value)
{
   ostringstream ss;
   ss << value;
   return ss.str();
}

// Convert Double in scientific format to string
string ToSciString(double value, int prec)
{
   ostringstream ss;
   ss << scientific << setprecision(prec) << value;
   return ss.str();
}

// Conversions between degrees and radians
double RadToDeg(double rad)
{
   return rad*180./TMath::Pi();
}

double DegToRad(double deg)
{
   return deg*TMath::Pi()/180.;
}

// Calculate the sin square of an angle (if degree = true, use degrees)
double SinSquare(double input, bool degree)
{
   if(degree)
      return sin(DegToRad(input))*sin(DegToRad(input));
   else
      return sin(input)*sin(input);
}

// Calculate the arcsin square of a value (if degree = true, result is in degrees)
double AsinSqrt(double input, bool degree)
{
   if(degree)
      return RadToDeg(asin(sqrt(input)));
   else
      return asin(sqrt(input));
}

// Calculate secant of an angle (if degree = true, use degrees)
double SecTheta(double input, bool degree)
{
   if(degree)
      return 1./cos(DegToRad(input));
   else
      return 1./cos(input);
}

// Calculate secant of a value (if degree = true, result is in degrees)
double InvSecTheta(double input, bool degree)
{
   if(degree)
      return RadToDeg(acos(1./input));
   else
      return acos(1./input);
}

// Remove extension from a string
string RemoveExtension(string *inname)
{
   size_t lastchar = inname->find_last_of(".");
   if(lastchar == string::npos)
      return *inname;
   return inname->substr(0, lastchar);
}

// Remove path from a string
string RemovePath(string *inname)
{
   return RemoveBeforeLast(inname, "/");
}

// Remove filename from a string
string RemoveFilename(string *inname)
{
   return RemoveFromLast(inname, "/");
}

// Replace extension in a string
string ReplaceExtension(string *inname, string ext)
{
   return (RemoveExtension(inname) + ext);
}

// Check if file has extension - if not, add it
string CheckExtension(string *inname, string ext)
{
   if(DBGSIG > 1)
      cout << "# CheckExtension        #: " << "Remove before last = " << RemoveBeforeLast(inname, ".") << endl;

   if(RemoveBeforeLast(inname, ".") == ext)
   {
      if(DBGSIG > 1)
         cout << "# CheckExtension        #: " << "Extension correct: " << *inname << endl;
      return *inname;
   }
   else
   {
      if(DBGSIG > 1)
         cout << "# CheckExtension        #: " << "Extension wrong: " << *inname << endl;
      return (*inname + "." + ext);
   }
}

// Compare the extension extension
int CompareExtension(string *inname, string ext)
{
   if(DBGSIG > 1)
      cout << "# CompareExtension        #: " << "Remove before last = " << RemoveBeforeLast(inname, ".") << endl;

   if(RemoveBeforeLast(inname, ".") == ext)
   {
      if(DBGSIG > 1)
         cout << "# CompareExtension        #: " << "Extensions equal: " << *inname << endl;
      return 1;
   }
   else
   {
      if(DBGSIG > 1)
         cout << "# CompareExtension        #: " << "Extensions not equal: " << *inname << endl;
      return 0;
   }
}

// Remove everything from the last occurence of a character to the end
string RemoveFromLast(string *inname, string search)
{
   size_t lastchar = inname->find_last_of(search);
   if(lastchar == string::npos)
      return *inname;
   return inname->substr(0, lastchar);
}

// Remove everything from the beginning to the last occurence of a character
string RemoveBeforeLast(string *inname, string search)
{
   size_t lastchar = inname->find_last_of(search);
   if(lastchar == string::npos)
      return *inname;
   return inname->substr(lastchar+1, inname->length());
}

// Remove everything from the beginning to the first occurence of a character
string RemoveBeforeFirst(string *inname, string search)
{
   size_t firstchar = inname->find_first_of(search);
   if(firstchar == string::npos)
      return *inname;
   return inname->substr(firstchar+1, inname->length());
}

// Replace a part of a string with something else
void ReplacePart(string *inname, string toreplace, string replace)
{
   size_t *searchchar;
   searchchar = new size_t;
  
   *searchchar = inname->find(toreplace);
   if(*searchchar != string::npos)
      inname->replace(*searchchar, toreplace.length(), replace);

   delete searchchar;
}

// Replace all occurrences of a character inside another string
void ReplaceAllCharacters(string *inname, string toreplace, string replace)
{
   size_t *start_pos;
   start_pos = new size_t;
   *start_pos = 0;

   while((*start_pos = inname->find(toreplace, *start_pos)) != string::npos)
   {
      inname->replace(*start_pos, toreplace.length(), replace);
      *start_pos += replace.length();
   }

   delete start_pos;
}

// Remove all leading spaces before the string
void RemoveLeadingSpaces(string *inname)
{
   size_t *pos = new size_t;
   *pos = inname->find_first_not_of(" \t");

   inname->erase(0, *pos);

   delete pos;
}

// Set the bin at set place for a binary number and return this number
void SetBinary(int place, int val, int *number)
{
   string *stemp;
   int *itemp;

   stemp = new string;
   *stemp = DecimalToBinary(*number);

   itemp = new int;
   *itemp = stemp->length();

   // Make sure we have leading zeroes, if place is above the current string length
   if(place >= *itemp)
   {
      for(int i = 0; i < (place-(*itemp)+1); i++)
         *stemp = "0" + *stemp;
   }
   *itemp = stemp->length();

   // Place the new value at the designated position (bins go from right to left, string positions go from left to right)
   if(stemp->at(*itemp-place-1) == '0')
   {
      if(val == 1)
         *number += pow(2, place);
   }
   else if(stemp->at(*itemp-place-1) == '1')
   {
      if(val == 0)
         *number -= pow(2, place);
   }

   delete stemp;
   delete itemp;
}

// Convert binary to decimal and read it
int BinaryToDecimal(string bin)
{
   int dec = 0;

   for(int i = 0; i <= bin.length(); i++)
   {
      if(bin[bin.length()-i-1] == '1')
      {
         dec += pow(2, i);
         if(DBGSIG > 1)
            cout << "# BinaryToDecimal       #: " << "Bin " << i << " is 1: " << pow(2, i) << " -> dec = " << dec << endl;
      }
   }

   return dec;
}

// Convert decimal to binary and read it
string DecimalToBinary(int dec)
{
   if(dec == 0)
      return "0";

   int *itemp;
   itemp = new int[3];
   itemp[0] = 0;
   itemp[1] = 0;
   itemp[2] = 0;

   string bin = "";

   while( dec > itemp[1] )
   {
      itemp[1] += pow(2, itemp[0]);
      itemp[0]++;
   }

   itemp[2] = dec;

   for(int i = itemp[0]-1; i >= 0; i--)
   {
      itemp[1] = itemp[2]/pow(2, i);
      if(itemp[1] == 1)
      {
         bin = bin + "1";
	 itemp[2] -= pow(2, i);
         if(DBGSIG > 1)
            cout << "# DecimalToBinary       #: " << "Bin " << i << " is 1 -> bin = " << bin << ", leftover = " << itemp[2] << endl;
      }
      else
         bin = bin + "0";
   }

   delete[] itemp;

   return bin;
}

// Searching for parts in a string vector
int Find(vector<string> strvec, string search)
{
   int retval = -1;

   for(int i = 0; i < strvec.size(); i++)
   {
      if( strvec[i].find(search) != string::npos )
      {
         retval = i;
	 break;
      }
   }

   return retval;
}

int FindVecPos(vector<int> intvec, int search)
{
   int retval = -1;

   for(int i = 0; i < intvec.size(); i++)
   {
      if( intvec[i] == search )
      {
         retval = i;
	 break;
      }
   }

   return retval;
}

int FindMinElement(vector<double> *invec)
{
   int retval = -1;
   double *dtemp = new double;
   *dtemp = 1.e+10;

   for(int i = 0; i < invec->size(); i++)
   {
      if( invec->at(i) < (*dtemp) )
      {
         *dtemp = invec->at(i);
	 retval = i;
      }
   }

   delete dtemp;

   return retval;
}

int SelectionPass(vector<double> *invec, double lowVal, double highVal)
{
   int retval = 0;

   if(invec->size() == 0)
      return -1;

   for(int i = 0; i < invec->size(); i++)
   {
      if( (invec->at(i) > lowVal) && (invec->at(i) <= highVal) )
         retval++;
   }

   return retval;
}

int FindStringPart(string instr, string tofind)
{
   if( instr.find(tofind) != string::npos )
      return 1;
   else
      return 0;
}

void SplitDelimitedList(string inlist, vector<int> *outvec)
{
   stringstream ss(inlist);
   int *itemp = new int;
   while (ss >> *itemp)
   {
      outvec->push_back(*itemp);
      if (ss.peek() == ',')
         ss.ignore();
   }

   for (int i = 0; i < outvec->size(); i++)
      cout << outvec->at(i) << endl;

   delete itemp;
}

void SplitDelimitedList(string inlist, vector<double> *outvec)
{
   stringstream ss(inlist);
   double *dtemp = new double;
   while (ss >> *dtemp)
   {
      outvec->push_back(*dtemp);
      if (ss.peek() == ',')
         ss.ignore();
   }

   for (int i = 0; i < outvec->size(); i++)
      cout << outvec->at(i) << endl;

   delete dtemp;
}

void SplitDelimitedList(string inlist, vector<string> *outvec)
{
   size_t *pos = new size_t;
   string *stemp = new string;
   string *delimit = new string;
   *delimit = ",";
   while((*pos = inlist.find(*delimit)) != std::string::npos)
   {
      *stemp = inlist.substr(0, *pos);
      outvec->push_back(*stemp);
      inlist.erase(0, *pos + delimit->length());
   }
   outvec->push_back(inlist);

   for (int i = 0; i < outvec->size(); i++)
      cout << outvec->at(i) << endl;

   delete delimit;
   delete stemp;
   delete pos;
}

// Get the correct number of keys, even if file was updated
int GetRootKeys(TFile *ifile, string keytemp)
{
   int ret;
   string *stemp = new string;
//   TList *keys;

//   ret = ifile->GetNkeys();
//   cout << "Number of keys reported by ROOT = " << ret << endl;
//   keys = (TList*) ifile->GetListOfKeys();

   ret = 0;
   for(int i = 1; i <= ifile->GetNkeys(); i++)
   {
      *stemp = keytemp + ToString(i);
      if(ifile->GetListOfKeys()->Contains(stemp->c_str()))
         ret++;
   }

   cout << "# GetRootKeys           #: " << "Number of " << keytemp << " keys in file = " << ret << endl;

/*   if(ifile->GetListOfKeys()->Contains("TreeA"))
   {
      cout << "# GetRootKeys           #: " << "Number of TreeA keys in file = " << 1 << endl;
      ret++;
   }*/

   delete stemp;

   return ret;
}

string GetSavedMethod(string *location)
{
   string nmeth;
   ifstream *instr = new ifstream;

   nmeth = *location + "/mva_method.dat";
   instr->open(nmeth.c_str(), ifstream::in);

   if(instr->is_open())
      *instr >> nmeth;
   else
   {
      cout << "No information on method name found." << endl;
      nmeth = "Fisher";
   }

   instr->close();
   delete instr;

   return nmeth;
}
