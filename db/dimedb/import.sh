#!/bin/sh

FD=~/Data/dimedb/output/dimedb_jsons/*

for f in $FD
do
	mongoimport --db dimedb --collection metabolites --type json --file $f --jsonArray
done
