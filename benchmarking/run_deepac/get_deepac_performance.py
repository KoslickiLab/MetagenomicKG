import sys, os, torch
import pandas as pd
metric_code=sys.argv[1]

# load calculation function
script_directory = os.path.dirname(metric_code)
if script_directory not in sys.path:
    sys.path.append(script_directory)

from metrics import ModelMetrics
from torchmetrics import Accuracy


# prepare data
merged_df = pd.read_csv("prediction/merged_deepac_prediction.tsv", header=None, sep="\t", names=['genome_id', 'Label', 'Pred_True_Pr'])
merged_df['Label'] = merged_df['Label'].astype('bool')
merged_df['Pred_False_Pr'] = 1 - merged_df['Pred_True_Pr']


metrics_obj = ModelMetrics(num_classes=2)
pos_accuracy = Accuracy(task="binary")
neg_accuracy = Accuracy(task="binary")


### input list
true_label = torch.tensor(merged_df['Label'].values, dtype=torch.float32)
posi_pr = torch.tensor(merged_df['Pred_True_Pr'].values, dtype=torch.float32)
nega_pr = torch.tensor(merged_df['Pred_False_Pr'].values, dtype=torch.float32)

# calculate metrics
metrics_obj.update(torch.tensor(merged_df['Pred_True_Pr'].values, dtype=torch.float32),
                                                                        torch.tensor(merged_df['Label'].values, dtype=torch.int32))

# pos data
pos_accuracy.update(torch.tensor(merged_df[merged_df['Label'] == True]['Pred_True_Pr'].values, dtype=torch.float32),
                                                                        torch.tensor(merged_df[merged_df['Label'] == True]['Label'].values, dtype=torch.int32))

# neg data
neg_accuracy.update(torch.tensor(merged_df[merged_df['Label'] == False]['Pred_True_Pr'].values, dtype=torch.float32),
                                                                        torch
.tensor(merged_df[merged_df['Label'] == False]['Label'].values, dtype=torch.i
nt32))


metrics_results = metrics_obj.compute()
pos_acc = pos_accuracy.compute()
neg_acc = neg_accuracy.compute()
acc = metrics_results['accuracy']
auroc = metrics_results['auroc']
ap = metrics_results['average_precision']
f1score = metrics_results['f1_score']

# output results
out_df = pd.DataFrame([acc, pos_acc, neg_acc, auroc, ap, f1score]).transpose(
)
out_df.columns = ['ACC', 'Pos_acc', 'Neg_acc', 'AUROC', 'AP', 'F1score']
out_df.to_csv("performance.tsv", sep="\t", header=True, index=False)
print(out_df)


