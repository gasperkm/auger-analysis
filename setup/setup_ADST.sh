#!/bin/bash

if [ "$1" == "" ]; then
  echo "ADST setup error! No starting directory supplied as first argument (same folder where configure is at)."
  exit 1
else
  startdir=$1
fi

if [ "$2" == "" ]; then
  echo "ADST setup error! No ADST version supplied as second argument."
  exit 1
fi

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

count=0
for ver in ${adstver[@]}
do
  if [ "$2" == "$ver" ]; then
    export ADSTROOT=$startdir/setup/ADST_${ver}
    export LD_LIBRARY_PATH=${ADSTROOT}/lib:${LD_LIBRARY_PATH}
    export PATH=${ADSTROOT}/bin:${PATH}
  fi

  count=$(( $count + 1 ))
done
