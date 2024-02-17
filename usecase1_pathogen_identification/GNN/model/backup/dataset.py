import os, sys
import dgl
import torch
import logging
import pickle
import numpy as np
import random
import pandas as pd
import sklearn.model_selection as ms
from sklearn.model_selection import train_test_split 
from dgl.data.utils import save_graphs

## Import custom libraries
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from utils import read_tsv_file, check_files, check_directory

class ProcessedDataset():
    def __init__(self, processed_dir, data_path, ratio, random_state, k_fold, logger):
        self.logger = logger
        
        self.processed_dir = processed_dir
        if not check_directory(processed_dir, self.logger):
            self.logger.error(f"Directory {processed_dir} does not exist. Please check it.")
            exit(1)
        self.data_path = data_path
        if not check_directory(data_path, self.logger):
            self.logger.error(f"Directory {data_path} does not exist. Please check it.")
            exit(1)
        self.ratio = ratio
        self.k_fold = k_fold
        self.random_state = random_state
        self.process_data()
    
    def __processed_file_names(self):
        return ["mask.pt", "info.pt", "kfold_index.pkl", "holdout_index.pkl", "mkg.bin"]
    
    def __create_mask(self, info_list):
        info_index = [self.id2index[item[0]] for item in info_list]
        mask = torch.zeros(len(self.id2index))
        mask[info_index] = 1
        return mask
    
    # def __split_data(self, pathogen_info, non_pathogen_info, holdout_pathogen_info):
        
    #     # exclude holdout pathogen from pathogen info as training data
    #     temp_pathogen_info = pd.DataFrame(pathogen_info)
    #     temp_pathogen_info = temp_pathogen_info.loc[~temp_pathogen_info[0].isin(pd.DataFrame(holdout_pathogen_info)[0]),:].sample(frac=1)
    #     training_pathogen_info = [tuple(row) for row in temp_pathogen_info.itertuples(index=False, name=None)]
        
    #     # get number of positive and negative data
    #     pos_num = len(training_pathogen_info)
    #     neg_num = int(pos_num * self.ratio)
        
    #     # get training negative data
    #     if neg_num > len(non_pathogen_info):
    #         self.logger.error(f"Ratio {self.ratio} is too large. Please check it.")
    #         exit(1)
    #     random.shuffle(non_pathogen_info)
    #     training_non_pathogen_info = non_pathogen_info[:neg_num]
    #     holdout_non_pathogen_info = non_pathogen_info[neg_num:]

    #     ## create pathogen mask
    #     all_pathogen_mask = self.__create_mask(pathogen_info)
        
    #     ## create non-pathogen mask
    #     all_non_pathogen_mask = self.__create_mask(non_pathogen_info)

    #     ## create training pathogen mask
    #     training_pathogen_mask = self.__create_mask(training_pathogen_info)
        
    #     ## create training non-pathogen mask
    #     training_non_pathogen_mask = self.__create_mask(training_non_pathogen_info)

    #     ## create holdout pathogen mask
    #     holdout_pathogen_mask = self.__create_mask(holdout_pathogen_info)
        
    #     ## create holdout non-pathogen mask
    #     holdout_non_pathogen_mask = self.__create_mask(holdout_non_pathogen_info)

    #     ## split data
    #     pos_index = torch.where(training_pathogen_mask == 1)[0].numpy().tolist()
    #     neg_index = torch.where(training_non_pathogen_mask == 1)[0].numpy().tolist()
    #     X = np.array(neg_index + pos_index)
    #     y = np.array([0] * neg_num + [1] * pos_num)

    #     cv = ms.StratifiedKFold(n_splits=10, random_state=self.random_state, shuffle=True)
    #     kfold_index = []
    #     for i, (train_val_index, test_index) in enumerate(cv.split(X, y)):
    #         X_train_val, X_test = X[train_val_index], X[test_index]
    #         y_train_val, y_test = y[train_val_index], y[test_index]
            
    #         # shuffle data
    #         X_train_val, y_train_val = self.__make_shuffle_index(X_train_val, y_train_val)
    #         X_test, y_test = self.__make_shuffle_index(X_test, y_test)
            
    #         # split train data into train and validation
    #         X_train, X_val, y_train, y_val = train_test_split(X_train_val, y_train_val, test_size=int(len(X) * (1/self.k_fold)), random_state=self.random_state, stratify=y_train_val, shuffle=True)
            
    #         # shuffle data
    #         X_train, y_train = self.__make_shuffle_index(X_train, y_train)
    #         X_val, y_val = self.__make_shuffle_index(X_val, y_val)
            
    #         kfold_index += [((X_train, y_train), (X_val, y_val), (X_test, y_test))]
        
    #     ## create holdout index
    #     holdout_pos_index = torch.where(holdout_pathogen_mask == 1)[0].numpy().tolist()
    #     holdout_neg_index = torch.where(holdout_non_pathogen_mask == 1)[0].numpy().tolist()
    #     holdout_index = np.array(holdout_neg_index + holdout_pos_index)
    #     holdout_y = np.array([0] * len(holdout_neg_index) + [1] * len(holdout_pos_index))
    #     holdout_index, holdout_y = self.__make_shuffle_index(holdout_index, holdout_y)
    #     holdout_index = (holdout_index, holdout_y)
        
    #     return kfold_index, holdout_index, (all_pathogen_mask, all_non_pathogen_mask, training_pathogen_mask, training_non_pathogen_mask, holdout_pathogen_mask, holdout_non_pathogen_mask), (training_pathogen_info, training_non_pathogen_info, holdout_non_pathogen_info)

    def __split_data(self, pathogen_info, non_pathogen_info, training_pathogen_info, training_non_pathogen_info, holdout_pathogen_info, holdout_non_pathogen_info):
        
        # get number of positive and negative data
        pos_num = len(training_pathogen_info)
        neg_num = len(training_non_pathogen_info)

        ## create pathogen mask
        all_pathogen_mask = self.__create_mask(pathogen_info)
        
        ## create non-pathogen mask
        all_non_pathogen_mask = self.__create_mask(non_pathogen_info)

        ## create training pathogen mask
        training_pathogen_mask = self.__create_mask(training_pathogen_info)
        
        ## create training non-pathogen mask
        training_non_pathogen_mask = self.__create_mask(training_non_pathogen_info)

        ## create holdout pathogen mask
        holdout_pathogen_mask = self.__create_mask(holdout_pathogen_info)
        
        ## create holdout non-pathogen mask
        holdout_non_pathogen_mask = self.__create_mask(holdout_non_pathogen_info)

        ## split data
        pos_index = torch.where(training_pathogen_mask == 1)[0].numpy().tolist()
        neg_index = torch.where(training_non_pathogen_mask == 1)[0].numpy().tolist()
        X = np.array(neg_index + pos_index)
        y = np.array([0] * neg_num + [1] * pos_num)

        cv = ms.StratifiedKFold(n_splits=10, random_state=self.random_state, shuffle=True)
        kfold_index = []
        for i, (train_val_index, test_index) in enumerate(cv.split(X, y)):
            X_train_val, X_test = X[train_val_index], X[test_index]
            y_train_val, y_test = y[train_val_index], y[test_index]
            
            # shuffle data
            X_train_val, y_train_val = self.__make_shuffle_index(X_train_val, y_train_val)
            X_test, y_test = self.__make_shuffle_index(X_test, y_test)
            
            # split train data into train and validation
            X_train, X_val, y_train, y_val = train_test_split(X_train_val, y_train_val, test_size=int(len(X) * (1/self.k_fold)), random_state=self.random_state, stratify=y_train_val, shuffle=True)
            
            # shuffle data
            X_train, y_train = self.__make_shuffle_index(X_train, y_train)
            X_val, y_val = self.__make_shuffle_index(X_val, y_val)
            
            kfold_index += [((X_train, y_train), (X_val, y_val), (X_test, y_test))]
        
        ## create holdout index
        holdout_pos_index = torch.where(holdout_pathogen_mask == 1)[0].numpy().tolist()
        holdout_neg_index = torch.where(holdout_non_pathogen_mask == 1)[0].numpy().tolist()
        holdout_index = np.array(holdout_neg_index + holdout_pos_index)
        holdout_y = np.array([0] * len(holdout_neg_index) + [1] * len(holdout_pos_index))
        holdout_index, holdout_y = self.__make_shuffle_index(holdout_index, holdout_y)
        holdout_index = (holdout_index, holdout_y)
        
        return kfold_index, holdout_index, (all_pathogen_mask, all_non_pathogen_mask, training_pathogen_mask, training_non_pathogen_mask, holdout_pathogen_mask, holdout_non_pathogen_mask)
    

    @staticmethod
    def __make_shuffle_index(X, y):
        xy_pairs = list(zip(X, y))
        random.shuffle(xy_pairs)
        X, y = zip(*xy_pairs)
        X = np.array(X)
        y = np.array(y)
        
        return X, y
    
    
    def load_data(self):        

        ## load pathogen information
        self.logger.info("Loading pathogen information...")
        pathogen_info = read_tsv_file(os.path.join(self.data_path, "pathogen_info.tsv"))
        pathogen_info_header = pathogen_info[0]
        pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in pathogen_info[1:] if 'species' not in item[5]]
        # temp_pathogen_info = pd.DataFrame(pathogen_info)[0].to_list()

        ## load non-pathogen information
        self.logger.info("Loading non-pathogen information...")
        non_pathogen_info = read_tsv_file(os.path.join(self.data_path, "non_pathogen_info.tsv"))
        non_pathogen_info_header = non_pathogen_info[0]
        non_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in non_pathogen_info[1:]]

        ## load training pathogen information
        self.logger.info("Loading training pathogen information...")
        training_pathogen_info = read_tsv_file(os.path.join(self.data_path, "pathogen_training_info.tsv"))
        training_pathogen_info_header = training_pathogen_info[0]
        training_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in training_pathogen_info[1:]]
        
        ## load training non-pathogen information
        self.logger.info("Loading training non-pathogen information...")
        training_non_pathogen_info = read_tsv_file(os.path.join(self.data_path, "non_pathogen_training_info.tsv"))
        training_non_pathogen_info_header = training_non_pathogen_info[0]
        training_non_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in training_non_pathogen_info[1:]]
        
        ## load holdout pathogen information
        self.logger.info("Loading holdout pathogen information...")
        holdout_pathogen_info = read_tsv_file(os.path.join(self.data_path, "pathogen_holdout_info.tsv"))
        holdout_pathogen_info_header = holdout_pathogen_info[0]
        holdout_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in holdout_pathogen_info[1:]]
        
        ## load holdout non-pathogen information
        self.logger.info("Loading holdout non-pathogen information...")
        holdout_non_pathogen_info = read_tsv_file(os.path.join(self.data_path, "non_pathogen_holdout_info.tsv"))
        holdout_non_pathogen_info_header = holdout_non_pathogen_info[0]
        holdout_non_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in holdout_non_pathogen_info[1:]]

        # ## load non-pathogen information
        # self.logger.info("Loading non-pathogen information...")
        # non_pathogen_info = read_tsv_file(os.path.join(self.data_path, "non_pathogen_info.tsv"))
        # non_pathogen_info_header = non_pathogen_info[0]
        # non_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in non_pathogen_info[1:]]
        # # remove non-pathogen that is actually pathogen identified by KEGG
        # temp_non_pathogen_info = pd.DataFrame(non_pathogen_info)
        # temp_non_pathogen_info = temp_non_pathogen_info.loc[~temp_non_pathogen_info[0].isin(temp_pathogen_info),:].sample(frac=1)
        # non_pathogen_info = [tuple(row) for row in temp_non_pathogen_info.itertuples(index=False, name=None)]

        # ## load holdout pathogen information
        # self.logger.info("Loading holdout pathogen information...")
        # holdout_pathogen_info = read_tsv_file(os.path.join(self.data_path, "hold_out_pathogen_info.tsv"))
        # holdout_pathogen_info_header = holdout_pathogen_info[0]
        # holdout_pathogen_info = [(item[0], item[1], eval(item[2]), item[3], eval(item[4]), eval(item[5]), eval(item[6]), eval(item[7]), eval(item[8]), eval(item[9])) for item in holdout_pathogen_info[1:]]

        ## Load edge information
        self.logger.info("Loading edge information...")
        edge_info = read_tsv_file(os.path.join(self.data_path, "edge_info.tsv"))
        edge_info_header = edge_info[0]
        edge_info = edge_info[1:]
        
        ## Load node embedding
        self.logger.info("Loading node embedding...")
        with open(os.path.join(self.data_path, "text_embedding", "embedding_biobert_namecat.pkl"), "rb") as f:
            embeddings_pl = pickle.load(f)
        
        node_embeddings = torch.tensor(np.array([embeddings_pl[k] for k, v in sorted(self.id2index.items(), key=lambda item: item[1])]))
        
        return (pathogen_info_header, pathogen_info), (non_pathogen_info_header, non_pathogen_info), \
               (training_pathogen_info_header, training_pathogen_info), (training_non_pathogen_info_header, training_non_pathogen_info), \
               (holdout_pathogen_info_header, holdout_pathogen_info), (holdout_non_pathogen_info_header, holdout_non_pathogen_info), \
               (edge_info_header, edge_info), node_embeddings
    
    def process_data(self):
        
        # check if processed data exists
        processed_file_path = [os.path.join(self.processed_dir, name) for name in self.__processed_file_names()]
        flag = False
        for file_path in processed_file_path:
            if not check_files(file_path, self.logger):
                self.logger.info("Re-run data processing step.")
                flag = True
                break

        self.logger.info("Loading node to index mapping...")
        node_to_index = read_tsv_file(os.path.join(self.data_path, "node_to_index.tsv"))
        node_to_index_header = node_to_index[0]
        node_to_index = {item[0]: int(item[1]) for item in node_to_index[1:]}
        self.id2index = node_to_index
        self.index2id = {value:key for key, value in self.id2index.items()}

        if flag:
            # load data
            (_, pathogen_info), (_, non_pathogen_info), (_, training_pathogen_info), (_, training_non_pathogen_info), (_, holdout_pathogen_info), (_, holdout_non_pathogen_info), (edge_info_header, edge_info), node_embeddings = self.load_data()
            
            ## create microbe mask
            microbe_info_index = [v for k, v in self.id2index.items() if k.split(":")[0] == "Microbe"]
            microbe_mask = torch.zeros(len(self.id2index))
            microbe_mask[microbe_info_index] = 1
            
            source_node_index, target_node_index = [], []
            for item in edge_info:
                source_node_index.append(self.id2index[item[edge_info_header.index('source_node')]])
                target_node_index.append(self.id2index[item[edge_info_header.index('target_node')]])
            # create graph
            mkg = dgl.graph((source_node_index, target_node_index), num_nodes=len(self.id2index), idtype=torch.int32)
            # make graph bidirectional
            mkg = dgl.to_bidirected(mkg)
            # set up node features
            mkg.ndata['feat'] = node_embeddings.float()
            
            # split training data into train, validation, and test by k-fold
            kfold_index, holdout_index, mask_list = self.__split_data(pathogen_info, non_pathogen_info, training_pathogen_info, training_non_pathogen_info, holdout_pathogen_info, holdout_non_pathogen_info)
            
            ## save processed data
            self.logger.info("Saving processed data...")
            
            # save mask
            mask_list = [microbe_mask] + list(mask_list)
            mask_file_names = ["microbe_mask", "pathogen_mask", "non_pathogen_mask", "training_pathogen_mask", "training_non_pathogen_mask", "holdout_pathogen_mask", "holdout_non_pathogen_mask"]
            torch.save([mask_file_names, mask_list], os.path.join(self.processed_dir, "mask.pt"))
            
            # save info
            info_list = [pathogen_info, non_pathogen_info, training_pathogen_info, training_non_pathogen_info, holdout_pathogen_info, holdout_non_pathogen_info]
            info_file_names = ["pathogen_info", "non_pathogen_info", "training_pathogen_info", "training_non_pathogen_info", "holdout_pathogen_info", "holdout_non_pathogen_info"]
            with open(os.path.join(self.processed_dir, "info.pt"), "wb") as f:
                torch.save([info_file_names, info_list], f)
            
            # save k-fold index
            with open(os.path.join(self.processed_dir, "kfold_index.pkl"), "wb") as f:
                pickle.dump(kfold_index, f)
                
            # save holdout index
            with open(os.path.join(self.processed_dir, "holdout_index.pkl"), "wb") as f:
                pickle.dump(holdout_index, f)
                
            # save graph
            save_graphs(os.path.join(self.processed_dir, "mkg.bin"), mkg)
            
            
        else:
            self.logger.info(f"Processed data exists in {self.processed_dir}, skip processing step.")
            return
        
    def load_mask_data(self):
        # check if processed mask data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[0])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
            
        self.logger.info(f"Loading processed mask data from {self.processed_dir}...")
        _, mask_list = torch.load(processed_file_path)
        all_pathogen_mask = mask_list[1]
        all_non_pathogen_mask = mask_list[2]
        training_pathogen_mask = mask_list[3]
        training_non_pathogen_mask = mask_list[4]
        holdout_pathogen_mask = mask_list[5]
        holdout_non_pathogen_mask = mask_list[6]
        
        return all_pathogen_mask, all_non_pathogen_mask, training_pathogen_mask, training_non_pathogen_mask, holdout_pathogen_mask, holdout_non_pathogen_mask

    def load_info_data(self):
        # check if processed info data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[1])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
            
        self.logger.info(f"Loading processed info data from {self.processed_dir}...")
        with open(processed_file_path, "rb") as f:
            _, info_list = torch.load(f)
        pathogen_info = info_list[0]
        non_pathogen_info = info_list[1]
        training_pathogen_info = info_list[2]
        training_non_pathogen_info = info_list[3]
        holdout_pathogen_info = info_list[4]
        holdout_non_pathogen_info = info_list[5]
        
        return pathogen_info, non_pathogen_info, training_pathogen_info, training_non_pathogen_info, holdout_pathogen_info, holdout_non_pathogen_info

    def load_kfold_index(self):
        # check if processed k-fold index data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[2])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None

        self.logger.info(f"Loading processed k-fold index data from {self.processed_dir}...")
        with open(os.path.join(self.processed_dir, "kfold_index.pkl"), "rb") as f:
            kfold_index = pickle.load(f)

        return kfold_index
    
    def load_holdout_index(self):
        # check if processed holdout index data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[3])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
        
        self.logger.info(f"Loading processed holdout index data from {self.processed_dir}...")
        with open(os.path.join(self.processed_dir, "holdout_index.pkl"), "rb") as f:
            holdout_index = pickle.load(f)
            
        return holdout_index
        
    
    def load_graph(self):
        # check if processed graph data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[4])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
        
        self.logger.info(f"Loading processed graph data from {self.processed_dir}...")
        mkg = dgl.load_graphs(os.path.join(self.processed_dir, "mkg.bin"))[0][0]
        return mkg
    