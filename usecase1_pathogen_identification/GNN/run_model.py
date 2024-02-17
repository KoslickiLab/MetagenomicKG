import os
import argparse
import torch
import numpy as np
import logging
import random

## Import custom libraries
from model import ModelController
from utils import get_logger, set_seed, comments

def get_args():
    """
    Parse the arguments from command line
    """
    parser = argparse.ArgumentParser(description='Main script for running model for pathogen detection')
    
    ## required arguments
    parser.add_argument('--data_path', required=True, type=str, help='path to the data folder containing preprocessed data')
    parser.add_argument('--output_dir', required=True, type=str, help='path of the output directory')
    
    ## optional general arguments
    parser.add_argument('--experiment_name', required=False, type=str, help='name of the experiment', default='pathogen_detection')
    parser.add_argument('--gpu', required=False, type=int, help='gpu device (default: 0)', default=0)
    parser.add_argument("--use_gpu", action="store_true", help="Whether use GPU or not", default=False)
    parser.add_argument('--random_seed', required=False, type=int, default=100, help='random seed (default: 100)')
    
    ## optional training arguments
    parser.add_argument('--ratio', required=False, type=float, help='ratio of negative and positive samples', default=2.0)
    parser.add_argument('--cv_num', required=False, type=int, help='number of folds for cross validation', default=10)
    parser.add_argument("--lr", required=False, type=float, help="Learning ratio", default=0.001)
    parser.add_argument("--num_neighbors", required=False, type=int, nargs='+', help="List of neighbors to sample", default=[100, 100, 100])
    parser.add_argument("--num_workers", required=False, type=int, help="Number of workers for data loader", default=32)
    parser.add_argument("--num_epochs", required=False, type=int, help="Number of epochs to train model", default=100)
    parser.add_argument("--train_batch_size", required=False, type=int, help="Batch size of training step", default=256)
    parser.add_argument("--eval_batch_size", required=False, type=int, help="Batch size of evaluation step", default=128)
    parser.add_argument("--opt_method", required=False, type=str, choices=["adam", "sgd"], help="Optimization method", default="adam")
    parser.add_argument("--patience", required=False, type=int, help="Number of epochs with no improvement after which learning rate will be reduced", default=10)
    parser.add_argument("--factor", required=False, type=float, help="The factor for learning rate to be reduced", default=0.1)
    parser.add_argument("--early_stop_n", required=False, type=int, help="Early stop if validation loss doesn't further decrease after n step", default=50)
    parser.add_argument("--dropout_p", required=False, type=float, help="Drop out proportion", default=0)
    
    ## optional model arguments
    parser.add_argument("--emb_size", required=False, type=int, help="Embedding vertor dimension", default=512)
    parser.add_argument("--num_layers", required=False, type=int, help="Number of GNN layers to train model", default=3)
    
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    
    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.DEBUG)
    args.logger = logger

    # set the seed for reproducibility
    set_seed(args.random_seed)
    
    # check for any conflicting parameter settings
    if len(args.num_neighbors) != args.num_layers:
        raise ValueError(f"Number of neighbors ({len(args.num_neighbors)}) must be the same as number of layers ({args.num_layers})")
    
    # set up device
    if args.use_gpu and torch.cuda.is_available():
        use_gpu = True
        device = torch.device(f'cuda:{args.gpu}')
        torch.cuda.reset_peak_memory_stats()
        torch.cuda.set_device(args.gpu)
    elif args.use_gpu:
        logger.warning("GPU is not available. Use CPU instead.")
        use_gpu = False
        device = 'cpu'
    else:
        use_gpu = False
        device = 'cpu'
    args.use_gpu = use_gpu
    args.device = device
    logger.info(f"Arugments: {args}")
    
    comments(f"Start an experiment '{args.experiment_name}'", logger=logger)
    # create a model controller object to train and test the model
    model_controller = ModelController(args)
    model_controller.run()


if __name__ == "__main__":
    main()