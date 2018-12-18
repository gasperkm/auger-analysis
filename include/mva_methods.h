#ifndef _MVA_METHODS_H_
#define _MVA_METHODS_H_

#include <string>
#include <iostream>
#include "root_include.h"

using namespace std;

void PrintMethods();
double GetMethodMin(TString name);
double GetMethodMax(TString name);
bool GetMethodDist(TString name);

#endif
