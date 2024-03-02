import os, sys
import numpy as np
from tqdm import tqdm
import logging

## Import PyTorch functions
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchmetrics import Accuracy

## Import custom libraries
from .metrics import ModelMetrics
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from utils import empty_gpu_cache, EVAL_STEP_MESSAGE


class Tester(object):

    def __init__(self, args, 
                model, 
                data,
                metrics_obj,
                mode='test'):
        """
        This class controls model test step
        Attributes:
            args (obj): an argparse object contains user-defined parameters
            model (obj): a model object
            data (obj): model-specific data preprocessing class object
            metrics_obj (obj): a TorchMetrics object
            mode (str): mode of the model (train, val, test, or holdout)
        Returns: None
        """
        # Create a logger object
        self.logger = args.logger

        # set up parameter
        self.args = args
        self.model = model
        self.data = data
        self.metrics_obj = metrics_obj
        self.mode = mode

    def run(self):
        """
        run Tester class to predict label
        Returns:
            acc (numpy.ndarray): Accuracy
            auroc (numpy.ndarray): AUROC
            ap (numpy.ndarray): Average Precision
            f1score (numpy.ndarray): F1 score
            recall (numpy.ndarray): Recalls
            precision (numpy.ndarray): Precisions
            fpr (numpy.ndarray): False Positive Rates
            tpr (numpy.ndarray): True Positive Rates
            pos_acc (numpy.ndarray): Accuracy for positive samples
            neg_acc (numpy.ndarray): Accuracy for negative samples
        """

        if self.mode == 'train':
            data_loader = self.data['train_dataloader']
            sorted_X = self.data['train_sorted_X']
            sorted_y = self.data['train_sorted_y']
        elif self.mode == 'val':
            data_loader = self.data['val_dataloader']
            sorted_X = self.data['val_sorted_X']
            sorted_y = self.data['val_sorted_y']
        elif self.mode == 'test':
            data_loader = self.data['test_dataloader']
            sorted_X = self.data['test_sorted_X']
            sorted_y = self.data['test_sorted_y']
        elif self.mode == 'holdout':
            data_loader = self.data['holdout_dataloader']
            sorted_X = self.data['holdout_sorted_X']
            sorted_y = self.data['holdout_sorted_y']
        else:
            raise ValueError("Only train, val, and test modes are supported")

        self.logger.info(f"Evaluate best model using {self.mode} dataset")
        
        acc, auroc, ap, f1score, recall, precision, fpr, tpr, pos_acc, neg_acc = self.predict(loader=data_loader, X=sorted_X, y=sorted_y)
        self.logger.info(EVAL_STEP_MESSAGE.format(
                                            eval_acc = acc,
                                            eval_auroc = auroc,
                                            eval_ap = ap,
                                            eval_f1score = f1score
                                        ))

        return acc, auroc, ap, f1score, recall, precision, fpr, tpr, pos_acc, neg_acc

    def predict(self, loader, X, y):
        """
        predict based on the given batch data
        Attributes:
            loader (obj): a PyTorch DataLoader object
            X (torch.tensor): a PyTorch tensor containing the sorted microbe node indices
            y (torch.tensor): a PyTorch tensor containing the sorted node label
        Returns:
            acc (numpy.ndarray): Accuracy
            auroc (numpy.ndarray): AUROC
            ap (numpy.ndarray): Average Precision
            f1score (numpy.ndarray): F1 score
            recall (numpy.ndarray): Recalls
            precision (numpy.ndarray): Precisions
            fpr (numpy.ndarray): False Positive Rates
            tpr (numpy.ndarray): True Positive Rates
            pos_acc (numpy.ndarray): Accuracy for positive samples
            neg_acc (numpy.ndarray): Accuracy for negative samples
        """ 

        ## switch model to evaluation mode
        self.model.eval()
        ## reset metrics
        self.metrics_obj.reset()

        X = X.to(self.args.device)
        y = y.to(self.args.device)

        # set up accuracy class for positive and negative samples
        if self.metrics_obj.num_classes == 2:
            pos_accuracy = Accuracy(task="binary")
            neg_accuracy = Accuracy(task="binary")
        else:
            pos_accuracy = Accuracy(task="multiclass", num_classes=self.metrics_obj.num_classes)
            neg_accuracy = Accuracy(task="multiclass", num_classes=self.metrics_obj.num_classes)

        with torch.no_grad():
            for _, output_nodes, blocks in tqdm(loader):

                input_feats = blocks[0].srcdata['feat']
                output_labels = y[torch.searchsorted(X, output_nodes)].float()
                pos_output_labels = output_labels[output_labels == 1]
                neg_output_labels = output_labels[output_labels == 0]
                output_preds= self.model(blocks, input_feats)
                pos_output_preds = output_preds[output_labels == 1]
                neg_output_preds = output_preds[output_labels == 0]
                empty_gpu_cache(self.args.device)
                self.metrics_obj.update(torch.sigmoid(output_preds).cpu().detach(), output_labels.int().cpu())
                if len(pos_output_preds) > 0:
                    pos_accuracy.update(torch.sigmoid(pos_output_preds).cpu().detach(), pos_output_labels.int().cpu())
                if len(neg_output_preds) > 0:
                    neg_accuracy.update(torch.sigmoid(neg_output_preds).cpu().detach(), neg_output_labels.int().cpu())

            metrics_results = self.metrics_obj.compute()
            pos_acc = pos_accuracy.compute()
            neg_acc = neg_accuracy.compute()
            acc = metrics_results['accuracy']
            auroc = metrics_results['auroc']
            ap = metrics_results['average_precision']
            f1score = metrics_results['f1_score']
            recall = metrics_results['recall']
            precision = metrics_results['precision']
            fpr = metrics_results['fpr']
            tpr = metrics_results['tpr']

        return acc, auroc, ap, f1score, recall, precision, fpr, tpr, pos_acc.numpy(), neg_acc.numpy()
