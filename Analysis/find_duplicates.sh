#!/bin/bash
# Runs duplicate-related operations based on the MOFid and/or MOFkeys,
# either duplicates within a database or the overlap between databases


function validate_input_file() {
	# Avoid sqlite3 errors by verifying that input files exist
	if [ ! -e "$1" ]
	then
		echo "Input file $1 does not exist!" 1>&2 && exit
	fi
	if [[ "$1" =~ ' ' ]]
	then
		echo "Input file \'$1\' cannot contain spaces" 1>&2 && exit
	fi
}

function build_sql_import() {
	# Build the $RUN_IMPORT command on a specified input file
	# WARNING: sqlite3.import will not accept paths containing spaces
	# Usage: build_sql_import input_file.tsv_or_smi sql_table_name
	
	# Using multiline strings: https://stackoverflow.com/questions/23929235/multi-line-string-with-extra-space-preserved-indentation
	if [[ "$1" == *.smi ]]  # MOFid input data
	then
		read -r -d '' temp_import <<IMPORT_HEREDOC
.mode tabs
.separator ";"
.headers off
CREATE TABLE $2 (identifier TEXT, filename TEXT);
.import $1 $2
IMPORT_HEREDOC
	elif [[ "$1" == *.tsv ]]  # MOFkey input data
	then
		read -r -d '' temp_import <<IMPORT_HEREDOC
.mode tabs
.headers on
.import $1 $2
-- Rename TSV column, keeping "filename" the same
ALTER TABLE $2 RENAME COLUMN mofkey TO identifier;
IMPORT_HEREDOC
	else
		echo "Unknown input data type.  Should be either raw MOFid or MOFkey input data" 1>&2
		exit 2
	fi
	
	# Build $RUN_IMPORT from the global environment
	RUN_IMPORT+=$'\n'  # https://stackoverflow.com/questions/3005963/how-can-i-have-a-newline-in-a-string-in-sh
	RUN_IMPORT+="$temp_import"
}


function build_sql_duplicates() {
	# Build the $RUN_DUPLICATES command similar to add_to_import()
	# Args: build_sql_duplicates orig_table duplicates_table
	# NOTE: this iteration of the duplicates analysis only deduplicates structures, not filtering out singletons
	
	read -r -d '' temp_duplicates <<DUPLICATES_HEREDOC
-- Note that this query is rather similar to polymorphs, actually
CREATE TABLE $2 AS
	SELECT identifier, COUNT(*) AS qty, GROUP_CONCAT(DISTINCT filename) AS duplicates
	FROM $1
	-- could add a WHERE topology filter here
	GROUP BY identifier
	-- HAVING qty > 1  -- doing this post-processing
	ORDER BY qty DESC;
DUPLICATES_HEREDOC
	
	# Build $RUN_DUPLICATES from the global environment
	RUN_DUPLICATES+=$'\n'
	RUN_DUPLICATES+="$temp_duplicates"
}


# Define sqlite3 variables
RUN_IMPORT=""  # build via add_to_import
RUN_DUPLICATES=""  # build via add_to_sql_duplicates
RUN_OUTPUT=""  # define below for overlap vs. duplicates


# Parse command line arguments
if [ $# -eq 4 ]
then
	if [[ "$1" != "duplicates" ]]
	then
		echo "Usage: \"duplicates\" must be specified as the first argument" 1>&2 && exit
	fi
	INPUT_FILE="$2"
	OUTPUT_NAMES_FILE="$3"
	OUTPUT_SUMMARY_FAMILIES="$4"
	validate_input_file "$INPUT_FILE"
	build_sql_import "$INPUT_FILE" mofs
	build_sql_duplicates mofs duplicates
	read -r -d '' RUN_OUTPUT <<OUTPUT_HEREDOC
.headers off
SELECT "Number of duplicates families:", COUNT(*) FROM duplicates;
.headers on

.output ${OUTPUT_SUMMARY_FAMILIES}
SELECT identifier, qty, duplicates
	FROM duplicates
	WHERE qty > 1
	ORDER BY qty DESC, identifier;

.output ${OUTPUT_NAMES_FILE}
SELECT mofs.filename, duplicates.qty, mofs.identifier
	FROM duplicates
	LEFT JOIN mofs on duplicates.identifier = mofs.identifier
	WHERE duplicates.qty > 1
	ORDER BY mofs.filename;
OUTPUT_HEREDOC

elif [ $# -eq 5 ]
then
	if [[ "$1" != "overlap" ]]
	then
		echo "Usage: \"overlap\" must be specified as the first argument" 1>&2 && exit
	fi
	INPUT_LEFT="$2"
	INPUT_RIGHT="$3"
	OUTPUT_NAMES_FILE="$4"
	OUTPUT_SUMMARY_FAMILIES="$5"
	validate_input_file "$INPUT_LEFT"
	validate_input_file "$INPUT_RIGHT"
	build_sql_import "$INPUT_LEFT" mofs_left
	build_sql_import "$INPUT_RIGHT" mofs_right
	build_sql_duplicates mofs_left dup_left
	build_sql_duplicates mofs_right dup_right
	read -r -d '' RUN_OUTPUT <<OUTPUT_HEREDOC
CREATE TABLE overlap AS
	SELECT
		dup_left.identifier AS identifier,
		dup_left.qty AS left_qty,
		dup_right.qty AS right_qty,
		dup_left.duplicates AS left_filenames,
		dup_right.duplicates AS right_filenames
	FROM dup_left
	INNER JOIN dup_right
		ON dup_left.identifier = dup_right.identifier;
ALTER TABLE overlap ADD COLUMN total_qty INTEGER;
UPDATE overlap SET total_qty = left_qty + right_qty;

-- SELECT identifier, total_qty, left_qty, right_qty FROM overlap LIMIT 10;
.headers off
SELECT "Number of overlap families:", COUNT(*) FROM overlap;
.headers on

.output ${OUTPUT_SUMMARY_FAMILIES}
SELECT * FROM overlap ORDER BY total_qty DESC, identifier;

.output ${OUTPUT_NAMES_FILE}
-- From what I recall, sqlite3 doesn't allow us to actually use CTE's with joins (WITH statements), so just create temporary tables
CREATE TABLE left_cte AS
	SELECT
		"left" AS db,
		mofs_left.filename AS filename,
		overlap.identifier AS identifier,
		left_qty AS db_qty, total_qty
	FROM overlap
	LEFT JOIN mofs_left ON overlap.identifier = mofs_left.identifier;
CREATE TABLE right_cte AS
	SELECT
		"right" AS db,
		mofs_right.filename AS filename,
		overlap.identifier AS identifier,
		right_qty AS db_qty, total_qty
	FROM overlap
	LEFT JOIN mofs_right ON overlap.identifier = mofs_right.identifier;
SELECT * FROM left_cte UNION ALL SELECT * FROM right_cte;
OUTPUT_HEREDOC

else
	echo "ERROR: incorrect usage." 1>&2
	echo "\"duplicates\" or \"overlap\" must be specified as the first argument, e.g." 1>&2
	echo "Analysis/find_duplicates.sh duplicates in_mofkey.tsv_OR_in_mofid.smi out_with_names.tsv out_summary.tsv" 1>&2
	echo "Analysis/find_duplicates.sh overlap left_mofkey.tsv_OR_in_mofid.smi right_input out_with_names.tsv out_summary.tsv" 1>&2
	exit
fi


# Run the built SQL query
sqlite3 <<EMBEDDED_SQL_HEREDOC
${RUN_IMPORT}

-- reset tabs and other I/O params
.mode tabs
.headers on

${RUN_DUPLICATES}

${RUN_OUTPUT}
EMBEDDED_SQL_HEREDOC
