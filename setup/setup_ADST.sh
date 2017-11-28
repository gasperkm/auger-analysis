#!/bin/bash

if [ "$1" == "" ]; then
  echo "ADST setup error! No starting directory supplied as first argument (same folder where configure is at)."
  exit 1
else
  startdir=$1
fi

adstold=$startdir/setup/ADST.r29701
adstnew=$startdir/setup/ADST_trunk

if [ "$2" == "old" ]; then
  export ADSTROOT=$adstold
  export LD_LIBRARY_PATH=${ADSTROOT}/lib:${LD_LIBRARY_PATH}
  export PATH=${ADSTROOT}/bin:${PATH}
elif [ "$2" == "new" ]; then
  export ADSTROOT=$adstnew
  export LD_LIBRARY_PATH=${ADSTROOT}/lib:${LD_LIBRARY_PATH}
  export PATH=${ADSTROOT}/bin:${PATH}
else
  echo "ADST setup error! No ADST version supplied as second argument (old or new)."
  exit 1
fi

#exit 0
