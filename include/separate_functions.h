#ifndef _SEPARATE_FUNCTIONS_H_
#define _SEPARATE_FUNCTIONS_H_

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "root_include.h"

using namespace std;

// Convert numbers to string
string ToString(double value, int prec);
string ToString(int value);
string ToSciString(double value, int prec);
double RadToDeg(double rad);
double DegToRad(double deg);

// Mathematical functions
double SinSquare(double input, bool degree);
double AsinSqrt(double input, bool degree);
double SecTheta(double input, bool degree);
double InvSecTheta(double input, bool degree);

// String manipulation
string RemoveExtension(string *inname);
string RemovePath(string *inname);
string RemoveFilename(string *inname);
string ReplaceExtension(string *inname, string ext);
string CheckExtension(string *inname, string ext);
int CompareExtension(string *inname, string ext);
string RemoveFromLast(string *inname, string search);
string RemoveBeforeLast(string *inname, string search);
string RemoveBeforeFirst(string *inname, string search);
void ReplacePart(string *inname, string toreplace, string replace);
void ReplaceAllCharacters(string *inname, string toreplace, string replace);
void RemoveLeadingSpaces(string *inname);

// Conversion binary <> decimal
void SetBinary(int place, int val, int *number);
int BinaryToDecimal(string bin);
string DecimalToBinary(int dec);

// Searching for parts in a string vector
int Find(vector<string> strvec, string search);
// Searching for location of integer in vector
int FindVecPos(vector<int> intvec, int search);

// Search vector for a specific value
int FindMinElement(vector<double> *invec);

// Search vector for the number of elements passing a selection
int SelectionPass(vector<double> *invec, double lowVal, double highVal);

// Searching for a part of a string
int FindStringPart(string instr, string tofind);

// Split comma delimited list
void SplitDelimitedList(string inlist, vector<int> *outvec);
void SplitDelimitedList(string inlist, vector<double> *outvec);
void SplitDelimitedList(string inlist, vector<string> *outvec);

// Get the correct number of keys, even if file was updated
int GetRootKeys(TFile *ifile, string keytemp);

// Read the method name
string GetSavedMethod(string *location);

#endif
