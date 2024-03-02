
# Import Python libraries
import os
import sys
import logging
import pandas as pd
import csv
import torch
import numpy as np
import random
from tqdm import tqdm, trange
from typing import List, Dict, Tuple, Union, Any, Optional, Set
csv.field_size_limit(100000000)  # set the maximum field size limit to 100MB or any value you need

# Set up constants
TRAINING_STEP_MESSAGE = "Epoch {epoch} / {num_epochs} | Train Metrics -- Loss: {train_loss:.5f}; Accuracy: {train_acc:.5f}; AUROC: {train_auroc:.5f}; Average Precision: {train_ap:.5f}; F1 Score: {train_f1score:.5f}"
VAL_STEP_MESSAGE = "Epoch {epoch} / {num_epochs} | Validation Metrics -- Loss: {val_loss:.5f}; Accuracy: {val_acc:.5f}; AUROC: {val_auroc:.5f}; Average Precision: {val_ap:.5f}; F1 Score: {val_f1score:.5f}"
EVAL_STEP_MESSAGE = "Evaluation Metrics -- Accuracy: {eval_acc:.5f}; AUROC: {eval_auroc:.5f}; Average Precision: {eval_ap:.5f}; F1 Score: {eval_f1score:.5f}"
COMM_STRING = "{delimiter}{left_whitespace}{message}{right_whitespace}{delimiter}\n"

def get_logger():
    """
    Setup a logger object
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def read_tsv_file(file_path, encoding='utf-8'):
    """
    Read a tsv file and return a list of lists
    """
    
    with open(file_path, newline='', encoding=encoding) as f:
        reader = csv.reader(f, delimiter='\t')
        data = [row for row in reader]
    return data

def check_files(file_path: str, logger):
    """
    Check if file path exists.
    Args:
        file_path: path of the file
        logger: logger object
    :return: True if the file exists, False otherwise
    """
    
    if not os.path.exists(file_path):
        logger.error(f"File {file_path} does not exist.")
        return False
    else:
        return True

def check_directory(dir_path: str, logger):
    """
    Check if directory path exists.
    Args:
        dir_path: path of the directory
        logger: logger object
    :return: True if the directory exists, False otherwise
    """
    
    if not os.path.exists(dir_path):
        logger.error(f"Directory {dir_path} does not exist.")
        return False
    else:
        return True
    

def set_seed(seed):
    """
    Set the seed for both PyTorch and NumPy to ensure reproducible results.
    Args:
        seed (int): The seed value to set for random number generation.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)

def empty_gpu_cache(device):
    """
    Empty the GPU cache to avoid memory overflow
    Args:
        device (torch.device object): the device object
    """
    if isinstance(device, torch.device):
        with torch.cuda.device(device):
            torch.cuda.empty_cache()


def comments(comment, logger, delim="#", max_len=90):
    """
    Print out pretty formatted comments
    Args:
        comment (str): the comment string to print
        logger (obj): the logger object
        delim (str): character to be used to wrap the comment
        max_len (int): maximum line length before break
    Returns: None
    """
    
    out_str = '\n' + delim * 120 + "\n"
    cut = comment.split(" ")
    line = ""
    for word in cut:
        if (len(line) + 1 + len(word)) > max_len:
            edge1 = (120 - len(line)) // 2 - 5
            edge2 = 120 - edge1 - len(line) - 10
            comm_info = {"delimiter": 5 * delim,
                         "left_whitespace": edge1 * " ",
                         "right_whitespace": edge2 * " ",
                         "message": line}
            out_str += COMM_STRING.format(**comm_info)
            line = word
        else:
            line = line + " " + word
    edge1 = (120 - len(line)) // 2 - 5
    edge2 = 120 - edge1 - len(line) - 10
    comm_info = {"delimiter": 5 * delim,
                 "left_whitespace": edge1 * " ",
                 "right_whitespace": edge2 * " ",
                 "message": line}
    out_str += COMM_STRING.format(**comm_info)
    out_str += delim * 120 + "\n"

    if logger is not None:
        logger.info(out_str)

def write_datatable(df_table, dist_path, separator="\t", header=True, index=False):
    """
    Write a pandas dataframe out to a file
    Args:
        df_table (pandas datafram): the pandas dataframe to write to a file
        dist_path (str): path to the table to write the dataframe to
        separator (str): the separator to use for column delimiting
        header (bool): include header in output
        index (bool): include index in output
    """
    # set up logger
    logger = get_logger()

    # create directory if it does not exist
    dir_name = os.path.dirname(dist_path)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    
    # write dataframe to file
    logger.info(f"Wrote dataframe with shape {df_table.shape} to file {dist_path}")
    df_table.to_csv(dist_path, sep=separator, index=index, header=header)
    
