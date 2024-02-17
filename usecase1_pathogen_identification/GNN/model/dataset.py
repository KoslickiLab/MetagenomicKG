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
    def __init__(self, processed_dir, data_path, ratio, random_state, logger):
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
        self.random_state = random_state
        self.process_data()
    
    def __processed_file_names(self):
        return ["mask.pt", "info.pt", "mkg.bin", "X_y.pt"]
    
    def __create_mask(self, info_list):
        info_index = [self.id2index[item[0]] for item in info_list]
        mask = torch.zeros(len(self.id2index))
        mask[info_index] = 1
        return mask
    
    def __split_data(self, pathogen_info, non_pathogen_info, training_info, valid_info, test_info):
        
        train_set_header, train_set = training_info
        valid_set_header, valid_set = valid_info
        test_set_header, test_set = test_info
        
        ## create pathogen mask
        all_pathogen_mask = self.__create_mask(pathogen_info)
        
        ## create non-pathogen mask
        all_non_pathogen_mask = self.__create_mask(non_pathogen_info)

        ## create training set mask
        train_mask = self.__create_mask(train_set)
        pos_index = [self.id2index[x[0]] for x in train_set if x[3] == 'True']
        neg_index = [self.id2index[x[0]] for x in train_set if x[3] == 'False']
        X_train = np.array(neg_index + pos_index)
        y_train = np.array([0] * len(neg_index) + [1] * len(pos_index))
        X_train, y_train = self.__make_shuffle_index(X_train, y_train)
        
        ## create validation set mask
        valid_mask = self.__create_mask(valid_set)
        pos_index = [self.id2index[x[0]] for x in valid_set if x[3] == 'True']
        neg_index = [self.id2index[x[0]] for x in valid_set if x[3] == 'False']
        X_valid = np.array(neg_index + pos_index)
        y_valid = np.array([0] * len(neg_index) + [1] * len(pos_index))
        X_valid, y_valid = self.__make_shuffle_index(X_valid, y_valid)
        
        ## create test set mask
        test_mask = self.__create_mask(test_set)
        pos_index = [self.id2index[x[0]] for x in test_set if x[3] == 'True']
        neg_index = [self.id2index[x[0]] for x in test_set if x[3] == 'False']
        X_test = np.array(neg_index + pos_index)
        y_test = np.array([0] * len(neg_index) + [1] * len(pos_index))
        X_test, y_test = self.__make_shuffle_index(X_test, y_test)
        
        return (all_pathogen_mask, all_non_pathogen_mask, train_mask, valid_mask, test_mask), ((X_train, y_train), (X_valid, y_valid), (X_test, y_test))
    

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

        ## load training set
        self.logger.info("Loading training set...")
        train_set = read_tsv_file(os.path.join(self.data_path, "train_set.tsv"))
        train_set_header = train_set[0][3:]
        train_set = [(item[3], item[4], eval(item[5]), item[6], eval(item[7]), eval(item[8]), eval(item[9]), eval(item[10]), eval(item[11]), eval(item[12])) for item in train_set[1:]]
        
        ## load validation set
        self.logger.info("Loading validation set...")
        valid_set = read_tsv_file(os.path.join(self.data_path, "valid_set.tsv"))
        valid_set_header = valid_set[0][3:]
        valid_set = [(item[3], item[4], eval(item[5]), item[6], eval(item[7]), eval(item[8]), eval(item[9]), eval(item[10]), eval(item[11]), eval(item[12])) for item in valid_set[1:]]
        
        ## load test set
        self.logger.info("Loading test set...")
        test_set = read_tsv_file(os.path.join(self.data_path, "test_set.tsv"))
        test_set_header = test_set[0][3:]
        test_set = [(item[3], item[4], eval(item[5]), item[6], eval(item[7]), eval(item[8]), eval(item[9]), eval(item[10]), eval(item[11]), eval(item[12])) for item in test_set[1:]]
        
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
        
        return (pathogen_info_header, pathogen_info), (non_pathogen_info_header, non_pathogen_info), (train_set_header, train_set), (valid_set_header, valid_set), (test_set_header, test_set), (edge_info_header, edge_info), node_embeddings
    
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
            (pathogen_info_header, pathogen_info), (non_pathogen_info_header, non_pathogen_info), (train_set_header, train_set), (valid_set_header, valid_set), (test_set_header, test_set), (edge_info_header, edge_info), node_embeddings = self.load_data()
            
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
            
            # split training data into train, validation, and test sets
            mask_list, X_y_list = self.__split_data(pathogen_info, non_pathogen_info, (train_set_header, train_set), (valid_set_header, valid_set), (test_set_header, test_set))
            
            ## save processed data
            self.logger.info("Saving processed data...")
            
            # X and y pairs
            with open(os.path.join(self.processed_dir, "X_y.pt"), "wb") as f:
                torch.save(X_y_list, f)
            
            # save mask
            mask_list = [microbe_mask] + list(mask_list)
            mask_file_names = ["microbe_mask", "pathogen_mask", "non_pathogen_mask", "train_mask", "valid_mask", "test_mask"]
            torch.save([mask_file_names, mask_list], os.path.join(self.processed_dir, "mask.pt"))
            
            # save info
            info_list = [(pathogen_info_header, pathogen_info), (non_pathogen_info_header, non_pathogen_info), (train_set_header, train_set), (valid_set_header, valid_set), (test_set_header, test_set)]
            info_file_names = ["pathogen_info", "non_pathogen_info", "train_info", "valid_info", "test_info"]
            with open(os.path.join(self.processed_dir, "info.pt"), "wb") as f:
                torch.save([info_file_names, info_list], f)
            
            # save graph
            save_graphs(os.path.join(self.processed_dir, "mkg.bin"), mkg)
            
            
        else:
            self.logger.info(f"Processed data exists in {self.processed_dir}, skip processing step.")
            return
        
    def load_X_y_data(self):
        # check if processed X_y data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[3])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
            
        self.logger.info(f"Loading processed X_y data from {self.processed_dir}...")
        X_y_list = torch.load(processed_file_path)
        train_set = X_y_list[0]
        valid_set = X_y_list[1]
        test_set = X_y_list[2]
        
        return train_set, valid_set, test_set
        
        
    def load_mask_data(self):
        # check if processed mask data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[0])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
            
        self.logger.info(f"Loading processed mask data from {self.processed_dir}...")
        with open(processed_file_path, "rb") as f:
            _, mask_list = torch.load(processed_file_path)
        microbe_mask = mask_list[0]
        all_pathogen_mask = mask_list[1]
        all_non_pathogen_mask = mask_list[2]
        train_mask = mask_list[3]
        valid_mask = mask_list[4]
        test_mask = mask_list[5]
        
        return microbe_mask, all_pathogen_mask, all_non_pathogen_mask, train_mask, valid_mask, test_mask

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
        train_info = info_list[2]
        valid_info = info_list[3]
        test_info = info_list[4]
        
        return pathogen_info, non_pathogen_info, train_info, valid_info, test_info
    
    def load_graph(self):
        # check if processed graph data exists
        processed_file_path = os.path.join(self.processed_dir, self.__processed_file_names()[2])
        if not check_files(processed_file_path, self.logger):
            self.logger.warning(f"File {processed_file_path} does not exist. Please run data processing step first.")
            return None
        
        self.logger.info(f"Loading processed graph data from {self.processed_dir}...")
        mkg = dgl.load_graphs(os.path.join(self.processed_dir, "mkg.bin"))[0][0]
        return mkg
    