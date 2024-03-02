
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