import pandas as pd
import glob
file_list = glob.glob("result_*.xlsx")
sheet_list = ["min_coverage0.1", "min_coverage0.05", "min_coverage0.01"]

for sheet_name in sheet_list:
	print("Merging "+ sheet_name)
	out_list = []
	for file in file_list:
		temp_df = pd.read_excel(file, sheet_name=sheet_name)
		temp_list = temp_df['organism_name'].to_list()
		temp_out = [x.split(" ")[0] for x in temp_list]
		out_list.extend(temp_out)
	out_list = list(set(out_list))
	out_df = pd.DataFrame(out_list)
	out_df.to_csv("merged_genomes_"+sheet_name+".txt", header=False, index=False)
