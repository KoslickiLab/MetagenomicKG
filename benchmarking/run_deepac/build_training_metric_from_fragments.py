# run under the "input_npy" dir
# file: build_training_metric_from_fragments.py
import numpy as np

### merge training data
arr_train1 = np.load("merged_neg_training_fragmented_genomes.npy")
arr_label1 = np.full(arr_train1.shape[0], 0, dtype=np.int32)
arr_train2 = np.load("merged_pos_training_fragmented_genomes.npy")
arr_label2 = np.full(arr_train2.shape[0], 1, dtype=np.int32)
# merge
out_train = np.vstack((arr_train1, arr_train2))
out_t_label = np.concatenate((arr_label1, arr_label2))
# save
np.save("train.npy", out_train)
np.save("label_train.npy", out_t_label)


### merge valid data
arr_valid1 = np.load("merged_neg_valid_fragmented_genomes.npy")
arr_v_label1 = np.full(arr_valid1.shape[0], 0, dtype=np.int32)
arr_valid2 = np.load("merged_pos_valid_fragmented_genomes.npy")
arr_v_label2 = np.full(arr_valid2.shape[0], 1, dtype=np.int32)
# merge
out_valid = np.vstack((arr_valid1, arr_valid2))
out_v_label = np.concatenate((arr_v_label1, arr_v_label2))
# save
np.save("valid.npy", out_valid)
np.save("label_valid.npy", out_v_label)

print("done")
