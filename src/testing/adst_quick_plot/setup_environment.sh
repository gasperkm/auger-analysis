#!/bin/bash

basedir=$PWD
cd ../../../
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

if [ "$1" == "" ]; then
  echo "ADST setup error! No ADST version supplied. Possible options:"
  for ver in ${adstver[@]}
  do
    echo "- $ver"
  done
else
  count=0
  for ver in ${adstver[@]}
  do
    if [ "$1" == "$ver" ]; then
      export ADSTROOT=$startdir/setup/ADST_${ver}
      export LD_LIBRARY_PATH=${ADSTROOT}/lib:${LD_LIBRARY_PATH}
      export PATH=${ADSTROOT}/bin:${PATH}
    fi

    count=$(( $count + 1 ))
  done
  echo "Using ADST in folder: $ADSTROOT"
fi

cd $basedir
