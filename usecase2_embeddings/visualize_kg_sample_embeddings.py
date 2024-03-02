"""
This script is used to visualize the sample-specific embeddings using t-SNE.
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import argparse
import logging
from functools import reduce
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

## Import custom libraries
from utils import get_logger, read_tsv_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize the sample-specific embeddings using t-SNE')
    parser.add_argument('--metadata', type=str, help='path of the metadata file')
    parser.add_argument('--raw_abundance', type=str, help='path of the raw abundance file')
    parser.add_argument('--top_values', type=int, help='number of top values to extract', default=1000)
    parser.add_argument('--embeddings', type=str, help='path of the sample-specific embeddings')
    parser.add_argument('--random_state', type=int, help='random state for t-SNE', default=42)
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()

    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.DEBUG)
    
    # Read metadata
    logger.info('Reading metadata file...')
    metadata = read_tsv_file(args.metadata)
    metadata = pd.DataFrame(metadata[1:], columns=metadata[0])
    sample_to_label = {row[0]:row[3] for row in metadata.to_numpy()}
    
    # Read raw abundance
    logger.info('Reading raw abundance file...')
    raw_abundance = pd.read_csv(args.raw_abundance, sep=',', header=0)
    X1 = raw_abundance.iloc[:,1:].to_numpy().T
    
    # Read embeddings
    logger.info('Reading sample-specific embeddings...')
    embeddings = pd.read_csv('data/usecase2/results/sample_embeddings.csv', sep=',', header=0)
    
    # Extract top values and merge table
    logger.info('Extracting top values and merging tables...')
    top_dfs = []

    # Loop through each samples
    for sample in tqdm(embeddings.columns[1:]):
        # Sort the DataFrame based on the current column in descending order and take the top values
        top_values = embeddings[['node_id', sample]].sort_values(by=sample, ascending=False).head(1000)
        # Append to the list
        top_dfs.append(top_values)

    # Merge all DataFrames in the list on 'node_id'
    merged_df = reduce(lambda left, right: pd.merge(left, right, on='node_id', how='outer'), top_dfs)
    
    # Fill NaN values with its original value
    merged_df = embeddings.loc[embeddings['node_id'].isin(merged_df['node_id'].to_list()),:]
    X2 = merged_df.iloc[:,1:].to_numpy().T
    
    # Perform t-SNE
    logger.info('Performing t-SNE...')
    tsne = TSNE(n_components=2, random_state=args.random_state)
    X1_tsne = tsne.fit_transform(X1)
    X2_tsne = tsne.fit_transform(X2)
    tsne_df_X1 = pd.DataFrame(X1_tsne, columns=['tsne_dim1', 'tsne_dim2'])
    tsne_df_X2 = pd.DataFrame(X2_tsne, columns=['tsne_dim1', 'tsne_dim2'])
    
    # Generate plot data
    body_site_label = pd.DataFrame([sample_to_label[sample] for sample in merged_df.columns[1:]], columns=['body_site'])
    plot_data_X1 = pd.concat([body_site_label, tsne_df_X1], axis=1)
    plot_data_X2 = pd.concat([body_site_label, tsne_df_X2], axis=1)
    
    # Generate plot
    logger.info('Generating plot...')
    # Create a figure and a set of subplots
    fig, axs = plt.subplots(1, 2, figsize=(20,10))

    # Combine the datasets for easier iteration
    plot_data_combined = [plot_data_X1, plot_data_X2]

    # Titles for each subplot
    titles = ['tSNE cluster based on raw aundance', 'tSNE cluster based on sample-specific embeddings']

    for ax, plot_data, title in zip(axs, plot_data_combined, titles):
        # Loop through each body site in the current plot_data
        for body_site in plot_data['body_site'].unique():
            # Filter data for the current body site
            subset = plot_data[plot_data['body_site'] == body_site]
            # Plot the filtered data
            ax.scatter(subset['tsne_dim1'], subset['tsne_dim2'], label=body_site)
        
        # Set title, labels and legend for each subplot
        ax.set_title(title, fontsize=20)
        ax.set_xlabel('t-SNE Dimension 1', fontsize=15)
        ax.set_ylabel('t-SNE Dimension 2', fontsize=15)
        ax.legend(title='Body Site')

    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(args.output_dir, 'tSNE_results.png'))