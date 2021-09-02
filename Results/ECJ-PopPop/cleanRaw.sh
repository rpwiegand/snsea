#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "Please specify a file (omit the extension)"
    exit 1
fi

INFILE=${1}.raw

if [ ! -f $INFILE ] ; then
    echo "$INFILE does not exist."
    exit 1
fi

echo "Reading $INFILE, writing ${1}.XX and ${1}.YY"
grep XX $INFILE | cut -d':' -f2 > ${1}.XX
#grep YY $INFILE | cut -d':' -f2 > ${1}.YY
