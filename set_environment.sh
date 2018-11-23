#!/bin/bash

   startdir=$PWD

   # Check for available ADST versions
   adstcount=0
   for adst in $startdir/setup/ADST_*.tar.gz
   do
     baseadst=$(basename $adst)
     baseadst=$(echo ${baseadst%%.*})
     baseadst=$(echo ${baseadst##*_})
     
     adstver[$adstcount]=$baseadst
     adstcount=$(( $adstcount + 1 ))
   done

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
     for verHelp in ${adstver[@]}
     do
       echo "  source set_environment.sh $verHelp"
     done
     echo ""
     echo "Clean all environment variables with:"
     echo "  source set_environment.sh clean"
   }

   function finishtext
   {
     echo "Done!"
     echo ""
     echo "Can now use Makefile to compile:"
     for verFin in ${adstver[@]}
     do
       echo "- Main program (with ADST version $verFin): make auger-analysis-gui-$verFin"
     done
     for verFin in ${adstver[@]}
     do
       echo "- Batch rewrite program (with ADST version $verFin): make bin/batch-rewrite-$verFin"
     done
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
     echo "Error! No arguments supplied (ADST version)."
     helptext
   else
     # Cleaning old environment variables
     echo "Cleaning old TMVA environment variables"
     cleanTmvaEnvironment $startdir/setup/TMVA
     echo "Cleaning old custom library path"
     cleanLibEnvironment $startdir/lib
     echo "Cleaning old ADST environment variables"
     cleanAdstEnvironment $startdir/setup/ADST

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
       err=1
       for ver in ${adstver[@]}
       do
	 if [ "$1" == "$ver" ]; then
           source $startdir/setup/setup_ADST.sh $startdir $ver
	   finishtext
	 fi
	 err=0
       done
        
       if [ $err == 1 ]; then
         echo "Error! No arguments supplied (ADST version)."
         helptext
       fi
     fi
   fi
