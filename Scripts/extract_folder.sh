#!/bin/bash

for i in "$1/"*.cif
do
	echo "Beginning $i" 1>&2
	echo "----------------------------" 1>&2
	python Python/extract_moffles.py $i
done
