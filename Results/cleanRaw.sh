#!/bin/bash

PREFIX=XX

if [ $# -lt 1 ] ; then
    echo "Please specify a file (omit the extension)"
    exit 1
fi


INFILE=${1}.raw

if [ ! -f $INFILE ] ; then
    echo "$INFILE does not exist."
    exit 1
fi

if [ $# -gt 1 ] ; then
    PREFIX=$2
fi

OUTFILE=${1}.${PREFIX}

echo "Reading $INFILE, writing $OUTFILE"
grep $PREFIX $INFILE | cut -d':' -f2 > $OUTFILE
