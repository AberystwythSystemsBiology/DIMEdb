#!/bin/sh

FD=~/.data/dimedb/jsons/*

for f in $FD
do
	mongoimport --db dimedb --collection metabolites --type json --file $f --jsonArray
done
