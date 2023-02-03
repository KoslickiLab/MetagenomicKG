import os, sys
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
import random
# import torch
import numpy as np
import pandas as pd
import logging
import logging.handlers
from subprocess import Popen, PIPE, STDOUT

## constants
COMM_STRING = "{delimiter}{left_whitespace}{message}{right_whitespace}{delimiter}\n"

## utility functions

# def set_random_seed(seed=1234):
#     random.seed(seed)
#     np.random.seed(seed)
#     torch.manual_seed(seed)
#     if torch.cuda.is_available():
#         torch.cuda.manual_seed_all(seed)

def get_logger(logname):
    logger = logging.getLogger(logname)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    fh = logging.handlers.RotatingFileHandler(logname, mode='w')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def strip_extension(name, extensions=[".tar.bz2", ".tar.gz", ".fastq.gz", ".gz", ".daa", ".rma"]):
    for ext in extensions:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name

def empty_gpu_cache(gpu):
    with torch.cuda.device(f'cuda:{gpu}'):
        torch.cuda.empty_cache()

def comm(comment, delim="*", max_len=90, logger=None):

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
    return out_str

def load_edge(path, sep='\t', symmetrize=False):
    edge_df = pd.read_csv(path, sep=sep, header=0)
    edge_df.columns = ['id1','id2','weight']
    edge_df.dropna(inplace=True)

    if symmetrize:
        rev_edge_df = edge_df.copy().rename(columns={'id1' : 'id2', 'id2' : 'id1'})
        edge_df = pd.concat([edge_df, rev_edge_df])

    idx_map = {name:idx for idx, name in enumerate(sorted(set(list(edge_df['id1']) + list(edge_df['id2']))))}
    idx_df = np.array(list(map(idx_map.get, np.array(edge_df[['id1', 'id2']]).flatten())), dtype=np.int32).reshape(np.array(edge_df[['id1', 'id2']]).shape).astype('int')
    weights = edge_df.weight.values.astype('float')

    return idx_df, idx_map, weights

def execute(cmd, logger=None, exit_on_fail=True, verbose=False):
    """
    Runs a system-level command with subprocess and makes sure it completes successfully.
    Args:
        cmd (str): command to execute
        logger (logging object): logger to write info to
        exit_on_fail (bool): should the script exit if the command returns a non-0
        verbose (bool): should the script print extra detailed messages
    Returns:
        exit_code (int): exit code of command (should be 0)

    """
    if logger is None:
        logger = get_logger('/'.join(path_list[:(project_index+1)] + ['log_folder', 'temp.log']))
    else:
        logger = logger

    logger.info("Executing command:\n{0}".format(cmd))
    process = Popen(cmd, stdout=PIPE, stderr=STDOUT, shell=True, executable="/bin/bash")

    while process.poll() is None:
        line = process.stdout.readline()
        if str(line) != "" and str(line) != b"":
            if "\\r" not in str(line):
                logger.info(line.rstrip())
    exitcode = process.wait()

    if exitcode != 0 and exit_on_fail:
        logger.info("Subprocess failed with error code %s" % exitcode)
        sys.exit(exitcode)
    else:
        if verbose:
            logger.info(f"Subprocess finished with exit code {exitcode}!")
        return exitcode

