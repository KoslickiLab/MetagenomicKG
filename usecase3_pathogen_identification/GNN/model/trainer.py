import os, sys
from copy import deepcopy
import numpy as np
from tqdm import tqdm
import logging

## Import PyTorch functions
import torch
import torch.nn as nn
import torch.nn.functional as F

## Import custom libraries
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from utils import empty_gpu_cache, TRAINING_STEP_MESSAGE, VAL_STEP_MESSAGE

class Trainer(object):

    def __init__(self, args,
                model,
                data,
                optimizer,
                scheduler,
                metrics_obj,
                writer):
        """
        This class controls model training step
        Attributes:
            args (obj): an argparse object contains user-defined parameters
            model (obj): a model object
            data (obj): model-specific data preprocessing class object
            optimizer (obj): a PyTorch optimizer object
            scheduler (obj): a PyTorch scheduler object
            metrics_obj (obj): a TorchMetrics object
            writer (obj): a PyTorch SummaryWriter object
        Returns: None
        """
        # Create a logger object
        self.logger = args.logger

        # set up parameter
        self.args = args
        self.model = model
        self.data = data
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.metrics_obj = metrics_obj
        self.writer = writer
        self.best_model_state = None
        self.best_model_name = None
        self.flag = True

    def run(self):
        count = 0
        current_min_val_loss = float('inf')
        training_range = tqdm(range(1, self.args.num_epochs + 1))
        self.logger.info("Start training model")
        for epoch in training_range:
            if count > self.args.early_stop_n:
                break
            self.logger.info('Training model...')
            train_loss, train_acc, train_auroc, train_ap, train_f1score  = self.__train_one_step()
            self.logger.info(TRAINING_STEP_MESSAGE.format(
                                                epoch = epoch,
                                                num_epochs = len(training_range),
                                                train_loss = train_loss,
                                                train_acc = train_acc,
                                                train_auroc = train_auroc,
                                                train_ap = train_ap,
                                                train_f1score = train_f1score
                                            ))
            self.logger.info('Validating model...')
            val_loss, val_acc , val_auroc, val_ap, val_f1score = self.__val_one_step()
            self.logger.info(VAL_STEP_MESSAGE.format(
                                                epoch = epoch,
                                                num_epochs = len(training_range),
                                                val_loss = val_loss,
                                                val_acc = val_acc,
                                                val_auroc = val_auroc,
                                                val_ap = val_ap,
                                                val_f1score = val_f1score
                                            ))
            self.scheduler.step(val_loss)
            self.writer.add_scalars('Loss', {'train': train_loss, 'val': val_loss}, epoch)
            self.writer.add_scalars('Accuracy', {'train': train_acc, 'val': val_acc}, epoch)
            self.writer.add_scalars('AUROC', {'train': train_auroc, 'val': val_auroc}, epoch)
            self.writer.add_scalars('Average Precision', {'train': train_ap, 'val': val_ap}, epoch)
            self.writer.add_scalars('F1 Score', {'train': train_f1score, 'val': val_f1score}, epoch)
            count += 1
            if val_loss < current_min_val_loss:
                count = 0
                current_min_val_loss = val_loss
                self.best_model_state = deepcopy(self.model.state_dict())
                self.best_model_name = f'{self.args.experiment_name}_best_epoch{epoch}.pt'

        self.writer.close()

    def __train_one_step(self):
        """
        This method is used for training step in one epoch
        Returns:
            train_loss (numpy.ndarray): the training loss calculated in this epoch
            train_acc (numpy.ndarray): the training accuracy calculated in this epoch
            train_auroc (numpy.ndarray): the training AUROC calculated in this epoch
            train_ap (numpy.ndarray): the training average precision calculated in this epoch
            train_f1score (numpy.ndarray): the training F1 score calculated in this epoch
        """
        ## switch model to training mode
        self.model.train()
        ## reset metrics
        self.metrics_obj.reset()

        total_loss = 0
        dataloader = self.data['train_dataloader']
        sorted_X = self.data['train_sorted_X'].to(self.args.device)
        sorted_y = self.data['train_sorted_y'].to(self.args.device)
        
        for _, output_nodes, blocks in tqdm(dataloader):

            input_feats = blocks[0].srcdata['feat']
            output_labels = sorted_y[torch.searchsorted(sorted_X, output_nodes)].float()
            output_preds= self.model(blocks, input_feats)
            empty_gpu_cache(self.args.device) # empty GPU cache to avoid memory leak
            train_loss = F.binary_cross_entropy_with_logits(output_preds, output_labels)
            self.metrics_obj.update(torch.sigmoid(output_preds).cpu().detach(), output_labels.int().cpu())
            self.optimizer.zero_grad()
            train_loss.backward()
            empty_gpu_cache(self.args.device) # empty GPU cache to avoid memory leak
            self.optimizer.step()
            total_loss += train_loss.detach().cpu().numpy()

        train_loss = total_loss / len(dataloader)
        train_metrics_results = self.metrics_obj.compute()
        train_acc = train_metrics_results['accuracy']
        train_auroc = train_metrics_results['auroc']
        train_ap = train_metrics_results['average_precision']
        train_f1score = train_metrics_results['f1_score']

        return train_loss, train_acc, train_auroc, train_ap, train_f1score

    def __val_one_step(self):
        """
        This method is model-specific and used for validation step in one epoch
        Returns:
            val_loss (numpy.ndarray): the validation loss calculated in this epoch
            val_acc (numpy.ndarray): the validation accuracy calculated in this epoch
            val_auroc (numpy.ndarray): the validation AUROC calculated in this epoch
            val_ap (numpy.ndarray): the validation average precision calculated in this epoch
            val_f1score (numpy.ndarray): the validation F1 score calculated in this epoch
        """
        ## switch model to evaluation mode
        self.model.eval()
        ## reset metrics
        self.metrics_obj.reset()

        total_loss = 0
        dataloader = self.data['val_dataloader']
        sorted_X = self.data['val_sorted_X'].to(self.args.device)
        sorted_y = self.data['val_sorted_y'].to(self.args.device)
        
        with torch.no_grad():
            for _, output_nodes, blocks in tqdm(dataloader):

                input_feats = blocks[0].srcdata['feat']
                output_labels = sorted_y[torch.searchsorted(sorted_X, output_nodes)].float()
                output_preds= self.model(blocks, input_feats)
                empty_gpu_cache(self.args.device)
                val_loss = F.binary_cross_entropy_with_logits(output_preds, output_labels)
                self.metrics_obj.update(torch.sigmoid(output_preds).cpu().detach(), output_labels.int().cpu())
                total_loss += val_loss.detach().cpu().numpy()

            val_loss = total_loss / len(dataloader)
            val_metrics_results = self.metrics_obj.compute()
            val_acc = val_metrics_results['accuracy']
            val_auroc = val_metrics_results['auroc']
            val_ap = val_metrics_results['average_precision']
            val_f1score = val_metrics_results['f1_score']

        return val_loss, val_acc, val_auroc, val_ap, val_f1score