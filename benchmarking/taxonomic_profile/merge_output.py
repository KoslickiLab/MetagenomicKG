# temp usage: merge Sourmash gather results
import os, glob
import pandas as pd
import re

gather_profile_dict = {}

for sub_file in glob.glob('sourmash*.csv'):
	filename = re.sub("(.csv)|(sourmash_gather_out_s1000_k31_)", "", os.path.basename(sub_file))
	# load data
	temp_df = pd.read_csv(sub_file, header=0, index_col='name', usecols=['name', 'f_unique_weighted'])
	temp_df.columns = [filename]
	# clean index col and keep GCx id only (exclude the string name)
	temp_df.index = temp_df.index.str.split(' ').str[0]
	gather_profile_dict[filename] = temp_df

# merge dataframes by "KO" column
out_df = pd.concat(gather_profile_dict.values(), axis=1, join='outer', sort=True)
out_df.fillna(0, inplace=True)
out_df.to_csv("merged_funiqweighted_fillna0.csv")
