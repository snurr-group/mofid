#!/bin/bash
# Calculates duplicates within a database based on the MOFid and/or MOFkeys
# TODO: consider between databases as a new DB_OP argument

# TODO: implement a new function for overlap based on this script,
# and document it as such.  Or, include it as a separate command line argument to this script.
# Note, for the common MOFs, can I just run an overlap followed by duplicates search?
# Perhaps the duplicates code should still export two outputs.
# But in that case, both the left and right would be concatenated in case of duplicates (with additional fields counting the number of duplicates)
# So in that case, the overlap script would actually form two duplicates databases separately, then join them. (so building out two duplicate detection strings with special L/R field names, then slightly adapting the end step

# And then the find_duplicates.sh script would just be a specialization of find_overlap.sh, with fewer args???

# TODO: make a new DB_OP as the first argument as duplicates OR overlap, and list both usages
# The initial part will probably be similar.  Just some minor differences in the column aggregation, etc.
# (just importing two databases, then filtering by an overlap qty???, and reporting qty_left and qty_right?)
# Also then a probably a new DB_NAME column in the overlap output file?  Or not exporting an explicit file like that?

if [ $# -ne 3 ]
then
	echo "ERROR: Usage: Analysis/find_duplicates.sh in_mofkey.tsv_OR_in_mofid.smi out_with_names.tsv out_summary.tsv" 1>&2 && exit
fi
INPUT_FILE="$1"
OUTPUT_NAMES_FILE="$2"
OUTPUT_SUMMARY_FAMILIES="$3"

# Define sqlite3 commands per filetype using multiline strings
# https://stackoverflow.com/questions/23929235/multi-line-string-with-extra-space-preserved-indentation
if [[ $INPUT_FILE == *.smi ]]  # MOFid input data
then
	read -r -d '' RUN_IMPORT <<IMPORT_HEREDOC
.mode tabs
.separator ";"
.headers off
CREATE TABLE mofs (identifier TEXT, filename TEXT);
.import ${INPUT_FILE} mofs
IMPORT_HEREDOC
elif [[ $INPUT_FILE == *.tsv ]]  # MOFkey input data
then
	read -r -d '' RUN_IMPORT <<IMPORT_HEREDOC
.mode tabs
.headers on
.import ${INPUT_FILE} mofs
-- Rename TSV column, keeping "filename" the same
ALTER TABLE mofs RENAME COLUMN mofkey TO identifier;
IMPORT_HEREDOC
else
	echo "Unknown input data type.  Should be either raw MOFid or MOFkey input data" 1>&2
	exit 2
fi


sqlite3 <<EMBEDDED_SQL_HEREDOC
${RUN_IMPORT}

-- reset tabs and other I/O params
.mode tabs
.headers on

--SELECT * FROM mofs LIMIT 10;

-- Note that this query is rather similar to polymorphs, actually
CREATE TABLE duplicates AS
	SELECT identifier, COUNT(*) AS qty, GROUP_CONCAT(DISTINCT filename) AS duplicates
	FROM mofs
	-- could add a WHERE topology filter here
	GROUP BY identifier
	HAVING qty > 1
	ORDER BY qty DESC;

.headers off
SELECT "Number of duplicates families:", COUNT(*) FROM duplicates;
.headers on

.output ${OUTPUT_SUMMARY_FAMILIES}
SELECT identifier, qty, duplicates
	FROM duplicates
	ORDER BY qty DESC, identifier;

.output ${OUTPUT_NAMES_FILE}
SELECT mofs.filename, duplicates.qty, mofs.identifier
	FROM duplicates
	LEFT JOIN mofs on duplicates.identifier = mofs.identifier
	ORDER BY mofs.filename;

/*
-- Old overlap code here:
.mode tabs
.headers on

.import CoreTables/combined_table.tsv core
.import GaTables/combined_table.tsv ga

.headers off
SELECT "Number of original GA MOFs:", COUNT(*) FROM ga;
DELETE FROM ga WHERE smiles = '*' OR topology = 'NA';
SELECT "Number of GA MOFs after removing NA's/*:", COUNT(*) FROM ga;
SELECT "";


CREATE TABLE overlap AS
	SELECT core.smiles, core.topology FROM
	(SELECT DISTINCT smiles, topology FROM core) AS core
	INNER JOIN
	(SELECT DISTINCT smiles, topology FROM ga) AS ga
	ON core.smiles = ga.smiles AND core.topology = ga.topology;


.output overlap_results.smi
SELECT * FROM overlap ORDER BY topology, smiles;
*/

EMBEDDED_SQL_HEREDOC
