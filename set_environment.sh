#!/bin/bash

   function cleanAdstEnvironment
   {
     cleandir=$(echo $1 | sed 's/\//\\\//g')
     newpath=$(echo ${PATH} | awk -v RS=: -v ORS=: "/${cleandir}/ {next} {print}" | sed 's/:*$//')
     newldpath=$(echo ${LD_LIBRARY_PATH} | awk -v RS=: -v ORS=: "/${cleandir}/ {next} {print}" | sed 's/:*$//')

     export PATH=$newpath
     export LD_LIBRARY_PATH=$newldpath
     unset ADSTROOT
   }

   function cleanLibEnvironment
   {
     cleandir=$(echo $1 | sed 's/\//\\\//g')
     newldpath=$(echo ${LD_LIBRARY_PATH} | awk -v RS=: -v ORS=: "/${cleandir}/ {next} {print}" | sed 's/:*$//')

     export LD_LIBRARY_PATH=$newldpath
   }

   function cleanTmvaEnvironment
   {
     cleandir=$(echo $1 | sed 's/\//\\\//g')
     cleanrootdir=$(echo $ROOTSYS | sed 's/\//\\\//g')
     newldpath=$(echo ${LD_LIBRARY_PATH} | awk -v RS=: -v ORS=: "/${cleandir}/ {next} {print}" | sed 's/:*$//')
     newpypath=$(echo ${PYTHONPATH} | awk -v RS=: -v ORS=: "/${cleandir}/ {next} {print}" | sed 's/:*$//')
     newpypath=$(echo ${newpypath} | awk -v RS=: -v ORS=: "/${cleanrootdir}/ {next} {print}" | sed 's/:*$//')

     export LD_LIBRARY_PATH=$newldpath
     export PYTHONPATH=$newpypath
     unset TMVASYS
   }

   function helptext
   {
     echo ""
     echo "Run script with one of the two options below (depending on version of ADST):"
     echo "  source set_environment.sh old"
     echo "  source set_environment.sh new"
     echo ""
     echo "Clean all environment variables with:"
     echo "  source set_environment.sh clean"
   }

   function finishtext
   {
     echo "Done!"
     echo ""
     echo "Can now use Makefile to compile:"
     echo "- Main program (with old ADST): make auger-analysis-gui-old"
     echo "- Main program (with new ADST): make auger-analysis-gui-new"
     echo "- Support programs (plotting,...): make support"
     echo "- Substructure library: make library"
   }

   function checkEnvironment
   {
     case ":$1:" in
       *:$2:*) val=0;;
       *) val=1;;
     esac

     return $val
   }

   # Check for arguments
   if [ "$1" == "" ]; then
     echo "Error! No arguments supplied (old/new)."
     helptext
   else
     startdir=$PWD

     # Cleaning old environment variables
     echo "Cleaning old TMVA environment variables"
#     checkEnvironment $PATH $startdir/setup/TMVA
#     if [ $? == 1 ]; then
     cleanTmvaEnvironment $startdir/setup/TMVA
#     fi
     echo "Cleaning old custom library path"
#     checkEnvironment $LD_LIBRARY_PATH $startdir/lib
#     echo "$startdir/lib"
#     if [ $? == 1 ]; then
     cleanLibEnvironment $startdir/lib
#     fi
     echo "Cleaning old ADST environment variables"
#     checkEnvironment $PATH $startdir/setup/ADST
#     if [ $? == 1 ]; then
     cleanAdstEnvironment $startdir/setup/ADST
#     fi

     if [ "$1" != "clean" ]; then
       cd $startdir
       # Set TMVA environment variables
       echo "Setting TMVA environment variables"
       source $startdir/setup/setup_TMVA.sh $startdir
       cd $TMVASYS
       source $TMVASYS/test/setup.sh $TMVASYS

       cd $startdir
       # Set custom library path
       echo "Setting custom library path"
       export LD_LIBRARY_PATH=$startdir/lib:$LD_LIBRARY_PATH
  
       cd $startdir
       # Set ADST environment variables
       echo "Setting ADST environment variables"
       if [ "$1" == "old" ]; then
         source $startdir/setup/setup_ADST.sh $startdir old
	 finishtext
       elif [ "$1" == "new" ]; then
         source $startdir/setup/setup_ADST.sh $startdir new
	 finishtext
       else
         echo "Error! No arguments supplied (old/new)."
         helptext
       fi
     fi
   fi
