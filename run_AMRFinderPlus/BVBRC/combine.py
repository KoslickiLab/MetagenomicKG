import pandas as pd
from tqdm import tqdm

# Read accession file with sccessful runs
accession_file = 'finished_accession_list.txt'
accession_file = pd.read_csv(accession_file, sep='\t', header=None)

# combine AMR results
combined_df = []
for genome_id in tqdm(accession_file[0].to_list()):
    temp_df = pd.read_csv(f'genomes/{genome_id}/amrfinder_results.txt', sep='\t', header=0)
    temp_df['genome_id'] = genome_id
    temp_df['source'] = 'PATRIC'
    combined_df += [temp_df]
combined_df = pd.concat(combined_df, axis=0).reset_index(drop=True)
combined_df.to_csv('amrfinderplus_all_results.txt', sep='\t', index=None)
    