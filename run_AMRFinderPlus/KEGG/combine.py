import pandas as pd

# Read mapping file
mapping_file = 'KEGG_microbe_assembly_mapping.tsv'
mapping = pd.read_csv(mapping_file, sep='\t', header=0)

# Read accession file with sccessful runs
accession_file = 'finished_accession_list.txt'
accession_file = pd.read_csv(accession_file, sep='\t', header=None)

# Filter some genomes without AMR results
filtered_accessions = mapping.loc[mapping['assembly_id'].isin(accession_file[0].to_list()),:].reset_index(drop=True)

# combine AMR results
combined_df = []
for row in filtered_accessions.to_numpy():
    gn_id, assembly_id = row
    temp_df = pd.read_csv(f'seqs/{assembly_id}/amrfinder_results.txt', sep='\t', header=0)
    temp_df['genome_id'] = gn_id
    temp_df['source'] = 'KEGG'
    combined_df += [temp_df]
combined_df = pd.concat(combined_df, axis=0).reset_index(drop=True)
combined_df.to_csv('amrfinderplus_all_results.txt', sep='\t', index=None)
    