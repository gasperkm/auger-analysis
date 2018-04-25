#include "frame.h"
#include <fstream>
#include <algorithm>

// Compare method names from values on the drop-down list
string MyFrame::GetMethodName(string name)
{
   vector<string>::iterator it;
   it = find(methodsDesc.begin(), methodsDesc.end(), name);
//   cout << "Finding 1 = " << *it << ", " << methodsDesc[distance(methodsDesc.begin(), it)] << endl;
   if(it != methodsDesc.end())
      return methods[distance(methodsDesc.begin(), it)];
   it = find(methods.begin(), methods.end(), name);
//   cout << "Finding 2 = " << *it << ", " << methods[distance(methods.begin(), it)] << endl;
   if(it != methods.end())
      return methodsDesc[distance(methods.begin(), it)];
}
