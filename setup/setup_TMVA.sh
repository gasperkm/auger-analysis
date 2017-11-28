#!/bin/bash

if [ "$1" == "" ]; then
  echo "TMVA setup error! No starting directory supplied as first argument (same folder where configure is at)."
  exit 1
else
  startdir=$1
fi

tmvadir=$startdir/setup/TMVA-v4.2.0

export TMVASYS=$tmvadir
