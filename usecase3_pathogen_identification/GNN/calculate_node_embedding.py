import os
import argparse
import json
import time
import pickle
from transformers import AutoModel, AutoTokenizer
import torch
import numpy as np
import logging
import random
from sklearn.decomposition import PCA

## Import custom libraries
from utils import get_logger, read_tsv_file, set_seed

def get_bert_embedding(texts, tokenizer, model, device):
    """
    Get BERT embedding for a list of texts
    """
    inputs = tokenizer(texts, padding=True, truncation=True, return_tensors="pt", max_length=512)
    inputs = inputs.to(device)
    with torch.no_grad():
        embeddings = model(**inputs, output_hidden_states=True, return_dict=True).pooler_output
    return embeddings.detach().to("cpu").numpy()

def get_args():
    """
    Parse the arguments from command line
    """
    parser = argparse.ArgumentParser(description='Convert node information to node embedding')
    parser.add_argument('--node_info', type=str, help='path to a file containing node information')
    parser.add_argument('--node_to_index', type=str, help='path to a file containing node to index mapping')
    parser.add_argument('--gpu', type=int, help='gpu device (default: 0)', default=0)
    parser.add_argument("--use_gpu", action="store_true", help="Whether use GPU or not", default=False)
    parser.add_argument('--random_seed', type=int, default=100, help='random seed (default: 100)')
    parser.add_argument('--final_embedding_dim', type=int, default=100, help='final embedding dimension (default: 100)')
    parser.add_argument("--batch_size", type=int, help="Batch size of bert embedding calculation", default=50)
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    
    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.DEBUG)
    set_seed(args.random_seed)

    # Load node information
    logger.info("Loading node information...")
    # load node information
    node_info = read_tsv_file(args.node_info)
    node_info_header = node_info[0]
    node_to_texts = {item[0]: item[1] + " " + " ".join(eval(item[2])) for item in node_info[1:]}
    # load node to index mapping
    node_to_index = read_tsv_file(args.node_to_index)
    node_to_index_header = node_to_index[0]
    id2index = {item[0]: int(item[1]) for item in node_to_index[1:]}
    index2id = {value:key for key, value in id2index.items()}
    # set up text descriptions for each node
    texts = [node_to_texts[node_id] for node_id, _ in id2index.items()]

    # create directory for saving embeddings
    if not os.path.exists(os.path.join(args.output_dir, "text_embedding")):
        os.makedirs(os.path.join(args.output_dir, "text_embedding"))

    # writing files
    with open(os.path.join(args.output_dir, "text_embedding", "id2index.json"), "w") as f:
        json.dump(id2index, f)
    with open(os.path.join(args.output_dir, "text_embedding", "index2id.json"), "w") as f:
        json.dump(index2id, f)

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

    # set up tokenizer and model
    tokenizer = AutoTokenizer.from_pretrained("microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext")
    model = AutoModel.from_pretrained("microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext")
    model.to(args.device)
    # set up original embedding
    ori_embedding = np.zeros([len(texts), 768])

    logger.info(f"Calculating BERT embedding on {args.device} with batch size: {args.batch_size}")

    start_time = time.time()

    for i in range(len(texts) // args.batch_size):
        if (i * args.batch_size) % 1000 == 0:
            logger.info(f"Finished: {i * args.batch_size} in {time.time() - start_time}")
            start_time = time.time()
        batch_text = texts[i*args.batch_size:(i+1)*args.batch_size]
        batch_embeddings = get_bert_embedding(batch_text, tokenizer, model, args.device)
        ori_embedding[i*args.batch_size:(i+1)*args.batch_size] = batch_embeddings
    
    if (i+1)*args.batch_size < len(texts):
        batch_text = texts[(i+1)*args.batch_size:]
        batch_embeddings = get_bert_embedding(batch_text, tokenizer, model, args.device)
        ori_embedding[(i+1)*args.batch_size:] = batch_embeddings
            
    # use PCA to reduce dimension
    logger.info("Fitting new embedding with PCA")
    pca = PCA(n_components=args.final_embedding_dim)
    pca_embedding = pca.fit_transform(ori_embedding)

    # save embeddings
    logger.info("Generating and saving data")
    id2embedding = {}
    for n_id in id2index.keys():
        id2embedding[n_id] = pca_embedding[id2index[n_id]]

    with open(os.path.join(args.output_dir, "text_embedding", "embedding_biobert_namecat.pkl"), 'wb') as f:
        pickle.dump(id2embedding, f, protocol=pickle.HIGHEST_PROTOCOL)
        
    logger.info("Finished")

if __name__ == "__main__":
    main()
