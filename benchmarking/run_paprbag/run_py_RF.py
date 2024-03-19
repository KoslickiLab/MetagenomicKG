import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
import os, sys
from joblib import dump


training_file="training.tsv"
valid_file="validation.tsv"


# load data
df_train = pd.read_csv(training_file, index_col=['RowName'], sep="|")
df_valid = pd.read_csv(valid_file, index_col=['RowName'], sep="|")
X_train = df_train.drop('Labels', axis=1)
y_train = df_train['Labels']
X_test = df_valid


# train RF
clf = RandomForestClassifier(n_estimators=100, random_state=33)
clf.fit(X_train, y_train)

# fit validation data
y_proba = clf.predict_proba(X_test)  # col0: F, col1: T

# load and merge true
out_df = pd.DataFrame(y_proba, index=X_test.index, columns=["Pred_False_Pr", "Pred_True_Pr"])
out_df["contig_name"] = out_df.index
df_label = pd.read_csv("merged_test_data/temp_contig_label_validation.tsv", sep="\t", names=["contig_name", "True_label", "genome_id"])

final_out = pd.merge(out_df, df_label,  left_on="contig_name", right_on="contig_name", how='inner')
final_out = final_out[["True_label", "Pred_True_Pr", "Pred_False_Pr", "genome_id", "contig_name"]]
final_out.to_csv("py_contig_level_RF_prediction.tsv", sep="\t", header=True, index=False)
