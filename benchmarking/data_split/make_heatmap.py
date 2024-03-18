from typing import Any
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys, csv, os, argparse, re, glob
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import fcluster
sys.setrecursionlimit(1000000)  # to avoid max recur for k=7 matrix

if len(sys.argv) < 2:
    print("Usage: python script.py workdir")
    sys.exit(1)

### read df
df_ani = pd.read_csv("ani_matrix_for_all_genomes.csv")
df_ani.columns = [re.sub(r".*/", "", col) for col in df_ani.columns]

### read metadata and prepare color map
# only has 1 label col now, extendable
df_meta = pd.read_csv("file_label.txt", header=None, sep="\t", names=['file', 'label'])
df_meta['file'] = df_meta['file'].apply(lambda x: re.sub(r".*/", "", x))
color_map_dict = dict()
legend_list = []
for col_name in ['label']:
    temp_distinct_category = df_meta[col_name].unique()
    lut = dict(zip(temp_distinct_category, sns.color_palette("hls", len(temp_distinct_category))))
    legend_patches = [mpatches.Patch(color=color, label=label) for label, color in lut.items()]
    legend_list.append(legend_patches)
    # need to get a group information in the same order as the matrix go
    records_in_sim_matrix = []
    for item in list(df_ani.columns):   # it's symmatric so ok to use cols
        records_in_sim_matrix.append(df_meta.loc[df_meta['file'] == item, col_name].values[0])
    color_map_dict[col_name] = pd.Series(records_in_sim_matrix, index=list(df_meta.index), name=col_name).map(lut)
color_map_df = pd.DataFrame(color_map_dict)


# generate clustermap
fig = sns.clustermap(df_ani, figsize=(20, 16), cmap="Greens", row_colors=color_map_df, metric='euclidean')
plt.legend(handles=legend_list[0], bbox_to_anchor=(1.05, 1), loc=2, borderaxespa
d=0.)
fig.fig.suptitle("Corr_ani_matrix_euclidean")
fig.savefig("ANI_clustermap_corrlation.png", dpi=300)

# extract cluster label
row_linkage = fig.dendrogram_row.linkage
n_clusters = 10
cluster_labels = fcluster(row_linkage, n_clusters, criterion='maxclust')
df_ani.index = df_ani.columns
# Map the cluster labels to the original data (use original data index)
df_ani['cluster'] = pd.Series(cluster_labels, index=df_ani.index)
df_ani.to_csv("ani_matric_with_cluster_label_of_10.csv")
del fig


# add group label into figure, code similar to above
df_ani_w_cluster = pd.read_csv("ani_matric_with_cluster_label_of_10.csv", index_col=0)
df_ani_w_cluster.columns = [re.sub(r".*/", "", col) for col in df_ani_w_cluster.columns]
# add the cluster label to metadata
df_meta = pd.read_csv("file_label.txt", header=None, sep="\t", names=['file', 'label'])
df_meta['file'] = df_meta['file'].apply(lambda x: re.sub(r".*/", "", x))
# cluster label
temp_df_cluster = df_ani_w_cluster.reset_index()
temp_df_cluster = temp_df_cluster[['index', 'cluster']].copy()
# merge
new_meta = pd.merge(df_meta, temp_df_cluster, left_on="file", right_on="index", how="inner")
new_meta.drop("index", axis=1, inplace=True)
# drop cluster col
df_ani_w_cluster.drop("cluster", axis=1, inplace=True)
df_ani_w_cluster.reset_index(drop=True, inplace=True)


new_color_map_dict = dict()
new_legend_list = []
for col_name in ['label', 'cluster']:
    temp_distinct_category = new_meta[col_name].unique()
    lut = dict(zip(temp_distinct_category, sns.color_palette("hls", len(temp_distinct_category))))
    legend_patches = [mpatches.Patch(color=color, label=label) for label, color in lut.items()]
    new_legend_list.append(legend_patches)
    # need to get a group information in the same order as the matrix go
    records_in_sim_matrix = []
    for item in list(df_ani_w_cluster.columns):   # it's symmatric so ok to use cols
        records_in_sim_matrix.append(new_meta.loc[new_meta['file'] == item, col_name].values[0])
    new_color_map_dict[col_name] = pd.Series(records_in_sim_matrix, index=list(df_ani_w_cluster.index), name=col_name).map(lut)
new_color_map_df = pd.DataFrame(new_color_map_dict)
new_color_map_df.reset_index(drop=True, inplace=True)

# generate clustermap
fig = sns.clustermap(df_ani_w_cluster, figsize=(20, 16), cmap="Greens", row_colors=new_color_ma
p_df, metric='euclidean')
# Get the axes for the heatmap
ax = fig.ax_heatmap
# Create the first legend and add it to the plot
first_legend = ax.legend(handles=new_legend_list[0], bbox_to_anchor=(1.5, 1), loc=2, borderaxes
pad=0.)
ax.add_artist(first_legend)
# Create the second legend and add it to the plot
ax.legend(handles=new_legend_list[1], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()

# plt.legend(handles=new_legend_list[1], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
fig.fig.suptitle("Corr_ani_matrix_euclidean_w_cluster")
fig.savefig("ANI_clustermap_corrlation_w_cluster.png", dpi=300)
