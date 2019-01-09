#!/bin/bash

OUTPUT_COPY="$2/results_part.txt"
mkdir -p "$2"

echo "Analyzing $1 with MOFid commit:" 1>&2
git rev-parse --verify HEAD 1>&2
echo "Writing results to intermediate directory: $2" 1>&2
echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2

for i in "$1/"*.cif
do
	echo "Beginning $i" 1>&2
	echo "----------------------------" 1>&2
	python Python/extract_moffles.py $i "$2" | tee -a $OUTPUT_COPY
done

echo "----------------------------" 1>&2
echo "Time:" 1>&2
date 1>&2
echo "----------------------------" 1>&2
