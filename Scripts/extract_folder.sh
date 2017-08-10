#!/bin/bash

OUTPUT_COPY="results_part.txt"

echo "Analyzing $1 with MOFid commit:" 1>&2
git rev-parse --verify HEAD 1>&2
echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2

for i in "$1/"*.cif
do
	echo "Beginning $i" 1>&2
	echo "----------------------------" 1>&2
	python Python/extract_moffles.py $i | tee -a $OUTPUT_COPY
done

echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2
