#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cstring>
#include <cstdio>
#include <iostream>

#include "parse_directory_files.h"
#include "separate_functions.h"
#include "workstation.h"

int ParseDirectory(string dir, vector<string> *outFiles)
{
   DIR *tDir;
   tDir = opendir(dir.c_str());
   if(tDir == nullptr)
   {
      cout << "Error opening directory " << dir << endl;
      return -1;
   }

   struct dirent *dirP;
   struct stat filestat;
   string *stemp = new string;

   while(dirP = readdir(tDir))
   {
      // Skip current object if it is this directory or parent directory
      if( !strncmp(dirP->d_name, ".", 1) || !strncmp(dirP->d_name, "..", 2) )
         continue;

      if(dir == ".")
         *stemp = dirP->d_name;
      else
         *stemp = dir + "/" + dirP->d_name;

      // Skip current file/directory if it is invalid in some way
      if(stat(stemp->c_str(), &filestat))
         continue;

      // Recursively call this function if current object is a directory
      if(S_ISDIR(filestat.st_mode))
      {
         ParseDirectory(*stemp, outFiles);
	 continue;
      }

      // Check if the current file is the file you are searching
      if(FindStringPart(*stemp, "/mvatree_file.root"))
      {
//         cout << "Found file: " << *stemp << endl;
         outFiles->push_back(*stemp);
      }
   }

   closedir(tDir);

   delete stemp;

   return 0;
}
