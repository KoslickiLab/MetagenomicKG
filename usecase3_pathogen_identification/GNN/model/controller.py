import os, sys
import logging
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

# import model associated classes
from .dataset import ProcessedDataset
from .model import MyModel
from .trainer import Trainer
from .tester import Tester
from .metrics import ModelMetrics

## Import PyTorch functions
import torch
import torch.optim as optim
from torch.utils.tensorboard import SummaryWriter

## Import DGL functions
import dgl
import dgl.nn as dglnn
import dgl.function as fn
from dgl.dataloading import (
    DataLoader,
    MultiLayerFullNeighborSampler,
    NeighborSampler,
)

# import custom libraries
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils import comments

class ModelController(object):

    def __init__(self, args):

        # Create a logger object
        self.logger = args.logger

        # Set up general parameters
        self.args = args
        self.data_path = args.data_path
        self.experiment_name = args.experiment_name
        self.device = args.device
        self.random_state = args.random_seed
        
        # set up training parameters
        self.ratio = args.ratio
        self.lr = args.lr
        self.num_neighbors = args.num_neighbors
        self.num_workers = args.num_workers
        self.num_epochs = args.num_epochs
        self.train_batch_size = args.train_batch_size
        self.eval_batch_size = args.eval_batch_size
        self.opt_method = args.opt_method
        self.patience = args.patience
        self.factor = args.factor
        self.early_stop_n = args.early_stop_n
        self.dropout_p = args.dropout_p
        
        # set up model parameters
        self.emb_size = args.emb_size
        self.num_layers = args.num_layers
        
        # setup output directory
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        self.output_dir = args.output_dir
        
        # initialize result
        self.summary_cv_test_dict = dict()
        self.summary_table_test_cv = os.path.join(self.output_dir, f"performance_summary_test.tsv")

    def run(self):
        """
        Run model training
        """
        comments(f"Train Model", logger=self.logger)
        self.run_cv_training()
        comments(f"Finished all CV iterations, saving model training summary", logger=self.logger)
        self.save_model_stats()

    def run_cv_training(self):
        """
        Perform cross validation training runs.
        Returns: None
        """
        # Load processed data
        processdata_dir = os.path.join(self.data_path, 'processed_model_data', self.experiment_name)
        # create a processed dataset directory if it doesn't exist
        if not os.path.exists(processdata_dir):
            os.makedirs(processdata_dir)
            
        # create a ProcessedDataset object
        data_obj = ProcessedDataset(processdata_dir, self.data_path, self.ratio, self.random_state, self.logger)
        
        # load processed data
        X_y_data = data_obj.load_X_y_data()
        mkg = data_obj.load_graph()
        mkg = mkg.to(self.device)
        
        # set up negihbor sampler
        train_step_sampler = NeighborSampler(self.num_neighbors, prefetch_node_feats=['feat'])
        eval_step_sampler = MultiLayerFullNeighborSampler(self.num_layers, prefetch_node_feats=['feat'])

        # set up model metrics
        metrics_obj = ModelMetrics(num_classes=2)

        ## set up dataloader
        train_set, valid_set, test_set = X_y_data
        
        # train set
        X_train, y_train = train_set
        X_train = torch.tensor(X_train, dtype=torch.int32)
        y_train = torch.tensor(y_train, dtype=torch.int32)
        sorted_indices = torch.argsort(X_train)
        train_sorted_X = X_train[sorted_indices]
        train_sorted_y = y_train[sorted_indices]

        if self.args.use_gpu:
            X_train = X_train.to(self.device)
            train_dataloader = DataLoader(mkg, X_train, eval_step_sampler, device=self.device, batch_size=self.eval_batch_size, \
                shuffle=True, drop_last=False, num_workers=0)
        else:
            train_dataloader = DataLoader(mkg, X_train, eval_step_sampler, batch_size=self.eval_batch_size, \
                shuffle=True, drop_last=False, num_workers=self.num_workers)

        # validation set
        X_val, y_val = valid_set
        X_val = torch.tensor(X_val, dtype=torch.int32)
        y_val = torch.tensor(y_val, dtype=torch.int32)
        sorted_indices = torch.argsort(X_val)
        val_sorted_X = X_val[sorted_indices]
        val_sorted_y = y_val[sorted_indices]

        if self.args.use_gpu:
            X_val = X_val.to(self.device)
            val_dataloader = DataLoader(mkg, X_val, eval_step_sampler, device=self.device, batch_size=self.eval_batch_size, \
                shuffle=True, drop_last=False, num_workers=0)
        else:
            val_dataloader = DataLoader(mkg, X_val, eval_step_sampler, batch_size=self.eval_batch_size, \
                shuffle=True, drop_last=False, num_workers=self.num_workers)

        # test set
        X_test, y_test = test_set
        X_test = torch.tensor(X_test, dtype=torch.int32)
        y_test = torch.tensor(y_test, dtype=torch.int32)
        sorted_indices = torch.argsort(X_test)
        test_sorted_X = X_test[sorted_indices]
        test_sorted_y = y_test[sorted_indices]
        
        if self.args.use_gpu:
            X_test = X_test.to(self.device)
            test_dataloader = DataLoader(mkg, X_test, eval_step_sampler, device=self.device, batch_size=self.eval_batch_size, \
                shuffle=True, drop_last=False, num_workers=0)
        else:
            test_dataloader = DataLoader(mkg, X_test, eval_step_sampler, batch_size=self.eval_batch_size, \
                shuffle=True, drop_last=False, num_workers=self.num_workers)

        # create a model object
        model = MyModel(self.num_layers, mkg.ndata['feat'].shape[1], self.emb_size, 1, self.dropout_p)
        model = model.to(self.device)

        ## set up optimizer and scheduler
        if self.opt_method == "adam":
            optimizer = optim.Adam(
                model.parameters(),
                lr = self.lr
            )
        else:
            optimizer = optim.SGD(
                model.parameters(),
                lr = self.lr
            )
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=self.factor, patience=self.patience, threshold=0.00001, threshold_mode='rel')

        ## set up tensorboard
        writer = SummaryWriter(log_dir=os.path.join(self.output_dir, 'tensorboard_runs'))

        data = {"train_dataloader": train_dataloader,
                "train_sorted_X": train_sorted_X,
                "train_sorted_y": train_sorted_y,
                "val_dataloader": val_dataloader,
                "val_sorted_X": val_sorted_X,
                "val_sorted_y": val_sorted_y,
                "test_dataloader": test_dataloader,
                "test_sorted_X": test_sorted_X,
                "test_sorted_y": test_sorted_y
                }


        ## train model
        self.logger.info("Initialize traner object")
        train_obj = Trainer(self.args, model, data, optimizer, scheduler, metrics_obj, writer)
        self.logger.info("Model controller starts training model")
        train_obj.run()
        ## save best model and data
        self.logger.info(f"Saving best model and data")
        self.save_model(train_obj.best_model_state, train_obj.best_model_name)
        self.save_data(data, 'data.pkl')

        ## evaluate model
        ## load best model
        model_state = self.load_best_model(train_obj.best_model_name)
        model = MyModel(self.num_layers, mkg.ndata['feat'].shape[1], self.emb_size, 1, self.dropout_p)
        model.load_state_dict(model_state)
        model = model.to(self.device)
        
        ## evaluate model with test data
        self.logger.info("Evaluate model with test data")
        test_obj = Tester(self.args, model, data, metrics_obj, mode='test')
        acc, auroc, ap, f1score, precision, recall, fpr, tpr, pos_acc, neg_acc = test_obj.run()
        self.update_model_stats(acc, pos_acc, neg_acc, auroc, ap, f1score, 'test')
        self.plot_performance_metrics(precision, recall, fpr, tpr, 'test')
            

    def save_model(self, model_state, model_name):
        """
        Save the model parameters
        Args:
            model_state (model.state_dict): model parameters
            model_name (str): model filename name
        Returns: None
        """
        save_dir = self.output_dir
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        # save model
        torch.save(model_state, os.path.join(save_dir, model_name))
        
    def save_data(self, data, data_name):
        """
        Save the data
        Args:
            data (dict): data dictionary
            data_name (str): data filename name
        Returns: None
        """
        save_dir = self.output_dir
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        # save data
        # remove dataloader objects
        data = {key:value for key, value in data.items() if 'dataloader' not in key}
        with open(os.path.join(save_dir, data_name), 'wb') as f:
            pickle.dump(data, f)

    def load_best_model(self, model_name):
        """
        Load the best model parameters
        Args:
            model_name (str): model filename name
        Returns:
            model_state (model.state_dict): model parameters
        """
        save_dir = self.output_dir
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        # load model
        model_state = torch.load(os.path.join(self.output_dir, model_name))
        return model_state

    def update_model_stats(self, acc, pos_acc, neg_acc, auroc, ap, f1score, type):
        """
        Save the model performance metrics
        Args:
            acc (float): Accuracy
            pos_acc (float): Accuracy for positive samples
            neg_acc (float): Accuracy for negative samples
            auroc (float): AUROC
            ap (float): Average Precision
            f1score (float): F1 score
            type (str): type of data (train, val, test, or holdout)
        Returns: None
        """
        
        if f'{type}_acc' in self.summary_cv_test_dict:
            self.summary_cv_test_dict[f'{type}_acc'].append(acc)
        else:
            self.summary_cv_test_dict[f'{type}_acc'] = [acc]
            
        if f'{type}_pos_acc' in self.summary_cv_test_dict:
            self.summary_cv_test_dict[f'{type}_pos_acc'].append(pos_acc)
        else:
            self.summary_cv_test_dict[f'{type}_pos_acc'] = [pos_acc]
            
        if f'{type}_neg_acc' in self.summary_cv_test_dict:
            self.summary_cv_test_dict[f'{type}_neg_acc'].append(neg_acc)
        else:
            self.summary_cv_test_dict[f'{type}_neg_acc'] = [neg_acc]
        
        if f'{type}_auroc' in self.summary_cv_test_dict:
            self.summary_cv_test_dict[f'{type}_auroc'].append(auroc)
        else:
            self.summary_cv_test_dict[f'{type}_auroc'] = [auroc]
        
        if f'{type}_ap' in self.summary_cv_test_dict:
            self.summary_cv_test_dict[f'{type}_ap'].append(ap)
        else:
            self.summary_cv_test_dict[f'{type}_ap'] = [ap]
            
        if f'{type}_f1score' in self.summary_cv_test_dict:
            self.summary_cv_test_dict[f'{type}_f1score'].append(f1score)
        else:
            self.summary_cv_test_dict[f'{type}_f1score'] = [f1score]

    def save_model_stats(self):
        """
        Save the model performance metrics
        Args: None
        Returns: None
        """
        # save model performance metrics
        df = pd.DataFrame(self.summary_cv_test_dict)
        df.to_csv(self.summary_table_test_cv, sep='\t', index=None, header=True)


    def plot_performance_metrics(self, precision, recall, fpr, tpr, type):
        """
        Plot the metric curves
        Args:
            precision (numpy.ndarray): Precision
            recall (numpy.ndarray): Recall
            fpr (numpy.ndarray): False Positive Rate
            tpr (numpy.ndarray): True Positive Rate
            type (str): type of data (train, val, test, or holdout)
        Returns: None
        """
        pr_path = os.path.join(self.output_dir, f"{type}_precision-recall_curve.png")
        p = sns.lineplot(x=recall, y=precision)
        p.set(xlabel="Recall", ylabel="Precision")
        plt.savefig(pr_path, bbox_inches='tight')
        plt.clf()

        roc_path = os.path.join(self.output_dir, f"{type}_roc_curve.png")
        p = sns.lineplot(x=fpr, y=tpr)
        p.set(xlabel="FPR", ylabel="TPR")
        plt.savefig(roc_path, bbox_inches='tight')
        plt.clf()

