#!/bin/bash
# Extracts polymorphs from a specified combined_table.tsv

if [ $# -ne 3 ]
then
	echo "ERROR: Usage: Analysis/find_polymorphs.sh in_combined_table.tsv out_with_names.tsv out_summary.tsv" 1>&2 && exit
fi
TSV_FILE="$1"
OUTPUT_NAMES_FILE="$2"
OUTPUT_SUMMARY_FAMILIES="$3"

sqlite3 <<EMBEDDED_SQL_HEREDOC
.mode tabs
.headers on

.import ${TSV_FILE} mofs

.headers off
CREATE TABLE polymorphic AS
	SELECT smiles, COUNT(DISTINCT topology) AS qty, GROUP_CONCAT(DISTINCT topology) AS topologies
	FROM mofs
	WHERE topology NOT IN ('ERROR', 'NEW', 'NA', 'MISMATCH', 'TIMEOUT')
	GROUP BY smiles
	HAVING qty > 1
	ORDER BY qty DESC;

-- SELECT "Number of polymorphic families:", COUNT(*) FROM polymorphic;

.headers on
.output ${OUTPUT_NAMES_FILE}
SELECT mofs.refcode, mofs.smiles, mofs.topology, polymorphic.qty
	FROM polymorphic
	LEFT JOIN mofs ON mofs.smiles = polymorphic.smiles
	ORDER BY polymorphic.qty DESC, mofs.smiles, mofs.topology, refcode;

.output ${OUTPUT_SUMMARY_FAMILIES}
SELECT * FROM polymorphic;

.output stdout
EMBEDDED_SQL_HEREDOC
