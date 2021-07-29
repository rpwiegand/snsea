#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "Please specify sigma and rho-min:"
  echo "  $0  01 02"
  exit 1
fi

OLDBASENAME="unboundedrv-s$1-r$2"
NEWBASENAME="unboundedrv-conv-s$1-r$2"

if [ ! -f ${OLDBASENAME}.raw ] ; then
  echo "${OLDBASENAME}.raw doesn't exist"
  exit 1
fi

echo "Moving ${OLDBASENAME}.raw to ${NEWBASENAME}.raw"
mv ${OLDBASENAME}.raw  ${NEWBASENAME}.raw

