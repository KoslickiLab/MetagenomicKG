#!/usr/bin/env bash
# tsv_to_neo4j.sh: Import TSV files into Neo4j

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo Usage: "$0 <database-name>=graph.db [test]"
    exit 2
fi

echo "================= starting tsv-to-neo4j.sh =================="
date

neo4j_config=${2:-"/data/shared_data/neo4j-community-3.5.26/conf/neo4j.conf"}
database_path=`grep dbms.directories.data ${neo4j_config} | cut -f2 -d=`
database=${1:-"graph.db"}
tsv_dir=${3:-"/home/grads/cqm5886/work/pathon_detection/data/KEGG/KG_info/"}

# change database and database paths to current database and database path in config file
# random_name=`date | sed 's/ /_/g' | sed 's/__/_/'`
# sed -i.${random_name}_bk "s/dbms.active_database=.*/dbms.active_database=${database}/" ${neo4j_config}
sed -i.bk "s/dbms.active_database=.*/dbms.active_database=${database}/" ${neo4j_config}
rm -rf ${neo4j_config}.bk

# stop neo4j
neo4j_command=`echo ${neo4j_config} | sed 's/conf\/neo4j.conf/bin\/neo4j/'`
${neo4j_command} stop

# change the database to write mode
sed -i.temp "s/dbms.read_only=true/dbms.read_only=false/" ${neo4j_config}
rm -rf ${neo4j_config}.temp

# delete the old log file and create a new one
rm -rf ${tsv_dir}/import.report
touch ${tsv_dir}/import.report

# import TSV files into Neo4j as Neo4j
neo4j_admin_command=`echo ${neo4j_config} | sed 's/conf\/neo4j.conf/bin\/neo4j-admin/'`
${neo4j_admin_command} import --nodes "${tsv_dir}/nodes_header.tsv,${tsv_dir}/nodes.tsv" \
    --relationships "${tsv_dir}/edges_header.tsv,${tsv_dir}/edges.tsv" \
    --max-memory=20G --multiline-fields=true --delimiter "\t" \
    --array-delimiter="ǂ" --report-file="${tsv_dir}/import.report" \
    --database=${database} --ignore-missing-nodes=true

# wait while neo4j boots up
echo "Sleeping for 1 minute, please do not SIGINT...."
sleep 1m
${neo4j_command} stop
${neo4j_command} start

# add indexes and constraints to the graph database
# python /home/grads/cqm5886/work/pathon_detection/scripts/process_data/KEGG/create_indexes_contrains.py
