#!/usr/bin/env bash
# tsv_to_neo4j.sh: Import TSV files into Neo4j

if [[ "${1:-}" == "" || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo Usage: "$0 <database-name>=graph.db"
    exit 2
fi

# Usage: tsv_to_neo4j.sh <database-name> <path_neo4j.conf> <path_TSV_folder>

echo "================= starting tsv-to-neo4j.sh =================="
date

database=${1:-"graph.db"}
neo4j_config=${2:-"../neo4j/neo4j-community-4.4.21/conf/neo4j.conf"}
tsv_dir=${3:-"../neo4j/input_files"}



# make a backup for the original configuration file
random_name=`date | sed 's/ /_/g' | sed 's/__/_/'`
neo4j_original_config=`echo ${neo4j_config} | sed 's/neo4j.conf.*/neo4j.conf.original/'`
cp ${neo4j_config} ${neo4j_original_config}_${random_name}

# change the default database to a given database name
sed -i "s/#dbms.active_database=.*/dbms.active_database=${database}.db/" ${neo4j_config}

# stop neo4j
neo4j_command=`echo ${neo4j_config} | sed 's/conf\/neo4j.conf/bin\/neo4j/'`
${neo4j_command} stop

# delete the old log file and create a new one
rm -rf ${tsv_dir}/import.report
touch ${tsv_dir}/import.report

# increae dbms memory max size
sed -i "s/#dbms.memory.heap.max_size=512m/dbms.memory.heap.max_size=10000m/" ${neo4j_config}

# setup http connection
# sed -i "s/#dbms.connector.http.listen_address=:7474/dbms.connector.http.listen_address=:7474/" ${neo4j_config}
# sed -i "s/#dbms.connector.https.listen_address=:7473/dbms.connector.https.listen_address=:7473/" ${neo4j_config}

# import TSV files into Neo4j as Neo4j
neo4j_admin_command=`echo ${neo4j_config} | sed 's/conf\/neo4j.conf/bin\/neo4j-admin/'`
${neo4j_admin_command} import --nodes "${tsv_dir}/nodes_header.tsv,${tsv_dir}/nodes.tsv" \
                              --relationships "${tsv_dir}/edges_header.tsv,${tsv_dir}/edges.tsv" \
                              --max-memory=20G \
                              --multiline-fields=true --delimiter "\t" \
                              --array-delimiter="ǂ" --report-file="${tsv_dir}/import.report" \
                              --database=${database}.db --ignore-missing-nodes=true


