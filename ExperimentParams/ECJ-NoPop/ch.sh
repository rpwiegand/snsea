#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "Please specify sigma and rho-min:"
  echo "  $0  01 02"
  exit 1
fi

OLDBASENAME="boundedrv-500-s$1-r$2.ini"
NEWBASENAME="unboundedrv-500-s$1-r$2.ini"

if [ ! -f ${OLDBASENAME} ] ; then
  echo "${OLDBASENAME} doesn't exist"
  exit 1
fi

#echo "Moving ${OLDBASENAME}.ini to ${NEWBASENAME}.ini"
#git mv ${OLDBASENAME}.ini  ${NEWBASENAME}.ini

cat ${OLDBASENAME} | sed s/"boundMutation = True"/"boundMutation = False"/ > ${NEWBASENAME}
 
