#ifndef _PRIMARY_TYPE_H_
#define _PRIMARY_TYPE_H_

#include "workstation.h"
//#include "root_include.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class PrimPart
{
private:
   vector<string> *names;
   vector<string> *shortnames;
   vector<double> *masses;

   int nrelem;

   string stemp;
   double dtemp;
public:
   PrimPart();
   virtual ~PrimPart();

   string GetName(string *inname);
   string GetName(int inz);
   string GetShortName(string *inname);
   string GetShortName(int inz);
   int GetZ(string *inname);
   double GetA(string *inname);
   double GetA(int inz);
   int Nr();

/*   int Z(string name);
   string Name(string shortname);
   string ShortName(string name);
   string Name(int z);
   string ShortName(int z);
   double A(string name);
   double A(int z);
   int N();*/
};

#endif
