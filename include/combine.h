#ifndef _COMBINE_H_
#define _COMBINE_H_

#include "frame.h"
#include "workstation.h"
#include "separate_functions.h"

void CombineRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames, vector<string> &titles );
void MergeRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames, int opt );
void CheckKeys( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts );
void hmerge(int nrfiles, string *files, string *outname);
void hadd(int nrfiles, string *files, string *outname, string *titlename, bool settitles);

#endif
