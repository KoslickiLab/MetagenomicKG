# ref:
# Github page: https://github.com/crarlus/paprbag
# Document: https://rdrr.io/github/crarlus/paprbag/man/Predict.ReadSet.fromFiles.html

# parameters
input_training_data_dir = "merged_training_data"
input_validation_data_dir = "merged_test_data"
training_data = "temp_merged_training_contigs.fasta"
training_label = "temp_contig_label_training.tsv"
validation_data = "temp_merged_validation_contigs.fasta"
validation_label = "temp_contig_label_validation.tsv"
validation_genome_label = "test_genome_label.tsv"

# ENV setup
library(paprbag)
library(ranger) # load ranger package explicitly

# output training data
Path2FeatureFile = paste0("Features/Features_", gsub(".fasta","", training_data), ".rds")
trainingData = readRDS(Path2FeatureFile)
# read labels
Path2LabelFile = paste0(input_training_data_dir, "/", training_label)
df_label = read.csv(Path2LabelFile, sep="\t", header= FALSE, row.names = 1)
df_label$V3 = NULL
colnames(df_label) = "Labels" # this name is required for the command

# merge df and put Lebels at first
trainingData$RowName <- rownames(trainingData)
df_label$RowName <- rownames(df_label)
merged_df <- merge(df_label, trainingData, by.x = "RowName", by.y = "RowName")
write.table(merged_df, "training.tsv", quote=FALSE, sep="|", row.names = FALSE)


# output valid data
Path2valid_data = paste0("Features/Features_", gsub(".fasta","", validation_data), ".rds")
feature_validation_data = readRDS(Path2valid_data)
feature_validation_data$RowName = row.names(feature_validation_data)
write.table(feature_validation_data, "validation.tsv", quote=FALSE, sep="|", row.names = FALSE)
