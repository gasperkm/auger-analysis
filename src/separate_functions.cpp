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
      return sin(input*TMath::Pi()/180.)*sin(input*TMath::Pi()/180.);
   else
      return sin(input)*sin(input);
}

// Calculate the arcsin square of a value (if degree = true, result is in degrees)
double AsinSqrt(double input, bool degree)
{
   if(degree)
      return asin(sqrt(input))*180./TMath::Pi();
   else
      return asin(sqrt(input));
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
   size_t searchchar = inname->find(toreplace);
   if(searchchar != string::npos)
      inname->replace(searchchar+1, (searchchar + toreplace.length()), replace);
}

// Replace all occurrences of a character inside another string
void ReplaceAllCharacters(string *inname, string toreplace, string replace)
{
   size_t start_pos = 0;

   while((start_pos = inname->find(toreplace, start_pos)) != string::npos)
   {
      inname->replace(start_pos, toreplace.length(), replace);
      start_pos += replace.length();
   }
}

// Set the bin at set place for a binary number and return this number
void SetBinary(int place, int val, int *number)
{
   string *stemp;
   stemp = new string;
   *stemp = DecimalToBinary(*number);

   int *itemp;
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
