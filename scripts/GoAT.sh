#!/bin/bash

PROG=$(readlink -f ${0})
GOATDIR=$(dirname "${PROG}")

CONFIG=$GOATDIR/../configfiles/$1
INFILE=$(readlink -f "${2}")
OUTFILE=$(readlink -f "${3}")

cd $GOATDIR
cd ..

build/bin/goat $CONFIG -f "${INFILE}" -F "${OUTFILE}"
