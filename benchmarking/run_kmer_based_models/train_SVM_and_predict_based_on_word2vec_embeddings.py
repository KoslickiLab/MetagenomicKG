from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
import numpy as np
import os, sys
import pandas as pd


metric_file = sys.argv[1]
training_list=sys.argv[2]
test_list=sys.argv[3]



# load data
df_matrics = pd.read_csv(metric_file, index_col=['filename'], sep=",")
df_matrics.index = df_matrics.index.map(lambda x: '_'.join(x.split('_')[:2]))
label_train = pd.read_csv(training_list, header=None, sep="\t", index_col=0, names=['label'])
label_test = pd.read_csv(test_list, header=None, sep="\t", index_col=0, names=['label'])

# merge and split
df_train = pd.merge(df_matrics, label_train, left_index=True, right_index=True, how='inner')
df_test = pd.merge(df_matrics, label_test, left_index=True, right_index=True, how='inner')

X_train = df_train.drop('label', axis=1)
y_train = df_train['label']
X_test = df_test.drop('label', axis=1)
y_test = df_test['label']


# Initialize and train the SVM classifier
svm_classifier = SVC(kernel='linear', probability=True)  # Enable probability estimates
svm_classifier.fit(X_train, y_train)


# Predict the labels for the test set
y_proba = svm_classifier.predict_proba(X_test)

# output
out_df = pd.DataFrame(y_proba, index=X_test.index, columns=["Pred_False_Pr", "Pred_True_Pr"])
final_out = pd.merge(out_df, label_test, left_index=True, right_index=True, how='inner')
final_out.index.name="filename"
final_out.to_csv("py_genome_level_SVM_prediction.tsv", sep="\t", header=True, index=True)

