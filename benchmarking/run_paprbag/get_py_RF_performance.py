import sys, os, torch
import pandas as pd
metric_code=sys.argv[1]
target_dir=sys.argv[2]

# load calculation function
script_directory = os.path.dirname(metric_code)
if script_directory not in sys.path:
    sys.path.append(script_directory)

from metrics import ModelMetrics
from torchmetrics import Accuracy


# prepare data
df = pd.read_csv(target_dir+"/py_contig_level_RF_prediction.tsv", header=0, sep="\t", usecols=['Pred_True_Pr', 'Pred_False_Pr', 'genome_id'])
### get group means for contig predictions
grouped_means = df.groupby('genome_id').mean()
### merge with true labels
valid_genome = pd.read_csv(target_dir+"/merged_test_data/test_genome_label.tsv", header=None, sep="\t", index_col=0, names=['Label'])
valid_genome['Label'] = valid_genome['Label'].astype('bool')
merged_df = pd.merge(grouped_means, valid_genome, left_index=True, right_index=True, how='inner')


# set up model metrics
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
                                                                        torch
.tensor(merged_df[merged_df['Label'] == True]['Label'].values, dtype=torch.in
t32))

# neg data
neg_accuracy.update(torch.tensor(merged_df[merged_df['Label'] == False]['Pred
_True_Pr'].values, dtype=torch.float32),
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
out_df.to_csv(target_dir+"/py_RF_performance.tsv", sep="\t", header=True, ind
ex=False)
print(target_dir)
print(out_df)
