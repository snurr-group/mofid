#!/bin/bash
# Combines TSV tables written by Python/convert_smi_to_tables.py
# into a single combined_table.tsv.
# Requires installation of [sqlite3](https://www.sqlite.org/index.html)
# and specifying the output directory to simplify.

TABLE_DIR="$1"

if [ $# -ne 1 ]
then
	if [ -f "smiles.tsv" ]
	then
		TABLE_DIR="."
	elif [ -d "TableOutput" ]
	then
		TABLE_DIR="TableOutput"
	else
		echo "ERROR: Cannot find relevant files.  Need to specify the directory to analyze" 1>&2 && exit
	fi
fi

echo "Combining TSV files in $TABLE_DIR" 1>&2
cd "$TABLE_DIR"

sqlite3 <<EMBEDDED_SQL_HEREDOC
.mode tabs
.headers off

CREATE TABLE smiles (refcode TEXT, smiles TEXT);
CREATE TABLE topology (refcode TEXT, topology TEXT);
CREATE TABLE cat (refcode TEXT, cat INTEGER);
CREATE TABLE parts (refcode TEXT, part TEXT);

.import smiles.tsv smiles
.import topology.tsv topology
.import cat.tsv cat
.import smiles_part.tsv parts

CREATE TABLE mofid AS
	SELECT s.refcode refcode, smiles, topology, cat
	FROM smiles s, topology t, cat c
	WHERE s.refcode = t.refcode AND s.refcode=c.refcode;


.headers on
.output combined_table.tsv
SELECT * FROM mofid;
EMBEDDED_SQL_HEREDOC
