#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "Please specify a file name"
  exit 1
fi

INFILENAME=$1
OUTFILENAME=`echo $1 | sed s/"unboundedrv"/"boundedrv"/`

echo "Generating $OUTFILENAME from $INFILENAME"
cat $INFILENAME | sed s/"boundMutation = False"/"boundMutation = True"/ > $OUTFILENAME
