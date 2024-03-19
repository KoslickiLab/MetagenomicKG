import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
import os, sys
from joblib import dump


training_file=sys.argv[1]
valid_file=sys.argv[2]
label_file=sys.argv[3]


# load data
df_train = pd.read_csv(training_file, index_col=['filename'], sep=",")
df_valid = pd.read_csv(valid_file, index_col=['filename'], sep=",")
X_train = df_train.drop('label', axis=1)
y_train = df_train['label']
X_test = df_valid.drop('label', axis=1)


# sort cols by colname
X_train.sort_index(axis=1, inplace=True)
X_test.sort_index(axis=1, inplace=True)


# train RF
clf = RandomForestClassifier(n_estimators=100, random_state=33)
clf.fit(X_train, y_train)

# fit validation data
y_proba = clf.predict_proba(X_test)  # col0: F, col1: T

# load and merge true
out_df = pd.DataFrame(y_proba, index=X_test.index, columns=["Pred_False_Pr", "Pred_True_Pr"])
df_label = pd.read_csv(label_file, header=None, sep="\t", index_col=0, names=['label'])

final_out = pd.merge(out_df, df_label, left_index=True, right_index=True, how='inner')
final_out.to_csv("py_genome_level_RF_prediction.tsv", sep="\t", header=True, index=True)
