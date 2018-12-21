#!/bin/bash

# args

tempfile=$(mktemp)
./config_gen.py > $tempfile

pcsim -c $tempfile -o OUTDIR -r1

rm "$tempfile"


