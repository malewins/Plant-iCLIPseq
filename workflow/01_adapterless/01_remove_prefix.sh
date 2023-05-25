#!/bin/bash

for FQGZ in *.gz;
do
  TARGET="${FQGZ#_barcode_}"

  mv $FQGZ $TARGET

  echo " -moving $TARGET"


done

