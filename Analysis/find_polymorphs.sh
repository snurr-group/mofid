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

ALTER TABLE mofs
	ADD COLUMN sn_topology TEXT;
UPDATE mofs
	SET sn_topology = SUBSTR(topology, 1, COALESCE(
		NULLIF((INSTR(topology, ',')-1), -1),  -- returns the comma position minus 1 if exists, else 0-1=-1
		LENGTH(topology)  -- else, coalesce and return the full string, because no comma
		));
ALTER TABLE mofs RENAME COLUMN topology TO all_topology;
ALTER TABLE mofs RENAME COLUMN refcode TO filename;

/*
-- Spot checking the block above, for debugging
SELECT topology, sn_topology FROM mofs WHERE INSTR(topology, ',') != 0 LIMIT 10;
SELECT topology, sn_topology FROM mofs WHERE INSTR(topology, ',') == 0 LIMIT 10;
*/

.headers off
CREATE TABLE polymorphic AS
	SELECT smiles, COUNT(DISTINCT sn_topology) AS qty, GROUP_CONCAT(DISTINCT sn_topology) AS topologies
	FROM mofs
	WHERE sn_topology NOT IN ('ERROR', 'NEW', 'NA', 'MISMATCH', 'TIMEOUT', 'UNKNOWN')
	GROUP BY smiles
	HAVING qty > 1
	ORDER BY qty DESC;

-- SELECT "Number of polymorphic families:", COUNT(*) FROM polymorphic;

.headers on
.output ${OUTPUT_NAMES_FILE}
SELECT mofs.filename, mofs.smiles, mofs.sn_topology, polymorphic.qty
	FROM polymorphic
	LEFT JOIN mofs ON mofs.smiles = polymorphic.smiles
	ORDER BY polymorphic.qty DESC, mofs.smiles, mofs.sn_topology, mofs.filename;

.output ${OUTPUT_SUMMARY_FAMILIES}
SELECT * FROM polymorphic;

.output stdout
EMBEDDED_SQL_HEREDOC
