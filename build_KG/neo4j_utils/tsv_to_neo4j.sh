#!/usr/bin/env bash
# tsv_to_neo4j.sh: Import TSV files into Neo4j

if [[ "${1:-}" == "" || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo Usage: "$0 <database-name>=MetagenomicsKG <path_TSV_folder>='../../neo4j/input_files'"
    exit 2
fi

# Usage: tsv_to_neo4j.sh <database-name> <path_TSV_folder>

echo "================= starting tsv-to-neo4j.sh =================="
date

database=${1:-"MetagenomicsKG"}
tsv_dir=${2:-"../../neo4j/input_files"}

# import TSV files into Neo4j as Neo4j
neo4j-admin database import full --nodes="${tsv_dir}/nodes.tsv" \
--relationships="${tsv_dir}/edges.tsv" \
--delimiter="\t" \
--array-delimiter="ǂ" "${database}.db"