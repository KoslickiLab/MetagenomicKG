import os
import numpy as np

## Import PyTorch functions
import torch
import torch.nn as nn
import torch.nn.functional as F

## Import DGL functions
import dgl
import dgl.nn as dglnn
import dgl.function as fn


class MyModel(nn.Module):
    def __init__(self, num_layers, in_size, hid_size, out_size, dropout, aggregator_type="mean"):
        super().__init__()
        self.num_layers = num_layers
        self.convs = nn.ModuleList()
        self.dropout = dropout
        self.hid_size = hid_size
        self.out_size = out_size
        
        self.fina_lin = nn.Linear(hid_size, out_size, bias=False)
        self.transform_lin = torch.nn.Linear(in_size, hid_size)
        # set up GraphSAGE layers
        for _ in range(num_layers):
            self.convs.append(dglnn.SAGEConv(hid_size, hid_size, aggregator_type))
        
    def reset_parameters(self):
        self.transform_lin.reset_parameters()
        for conv in self.convs:
            conv.reset_parameters()
        self.fina_lin.reset_parameters()

    def forward(self, blocks, x):
        h = x
        ## transform input features
        h = self.transform_lin(h)
        for _, (layer, block) in enumerate(zip(self.convs, blocks)):
            h = layer(block, h)
            h = F.relu(h)
            h = F.dropout(h, self.dropout, training=self.training)
        h = self.fina_lin(h)
            
        return h.squeeze()
