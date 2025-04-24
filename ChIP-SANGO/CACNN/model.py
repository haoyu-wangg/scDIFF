import torch
from torch import Tensor
import torch.nn as nn
from typing import Iterable, Optional
from ECA_layer import eca_layer

ONEHOT = torch.cat((
    torch.ones(1, 4) / 4, # 
    torch.eye(4), # A, C, G, T
    torch.zeros(1, 4), # padding
), dim=0).float()

class ConvTower(nn.Module):
    """
    CNN layer with ECA-Net

    in_channel: input dim
    out_channel: output dim
    kernel_size: kernel size of CNN

    return: embedding after CNN layer
    """
    def __init__(self, in_channel, out_channel: int, kernel_size) -> None:
        super().__init__()
        self.conv1 = nn.Conv1d(in_channel, out_channel, kernel_size, padding=kernel_size//2, stride=1, bias=False)
        self.bn1 = nn.BatchNorm1d(out_channel)
        self.relu = nn.ReLU()

        self.conv2 = nn.Conv1d(out_channel, out_channel, kernel_size, padding=kernel_size//2, stride=1, bias=False)
        self.bn2 = nn.BatchNorm1d(out_channel)
        self.eca = eca_layer(out_channel)

        self.maxpool = nn.MaxPool1d(kernel_size=kernel_size//2)

        self.downsample = nn.Sequential(
            nn.Conv1d(in_channel, out_channel, kernel_size=1, stride=1, bias=False),
            nn.BatchNorm1d(out_channel),
        )
    
    def forward(self, x: Tensor) -> Tensor:
        residual = x

        # conv1
        y = self.conv1(x)
        y = self.bn1(y)
        y = self.relu(y)

        # conv2 + eca
        y = self.conv2(y)
        y = self.bn2(y)
        y = self.eca(y)

        residual = self.downsample(residual)

        y += residual
        y = self.maxpool(y)
        y = self.relu(y)

        return y

class CACNN(nn.Module):
    """
    CACNN model implementation
    """
    def __init__(self, n_cells: int, batch_ids: Optional[Iterable[int]] = None, use_reg_cell=False, hidden_size=32, seq_len: int = 1344, epifeature_dim: int = 3):
        super().__init__()
        self.config = {
            "n_cells": n_cells,
            "hidden_size": hidden_size,
            "seq_len": seq_len,
            "epifeature_dim": epifeature_dim,
        }
        if batch_ids is None:
            self.batch_ids = None
        else:
            self.batch_embedding = nn.Embedding(max(batch_ids) + 1, hidden_size)
            self.batch_ids = nn.Parameter(torch.as_tensor(batch_ids), requires_grad=False)
            assert self.batch_ids.ndim == 1
        self.onehot = nn.Parameter(ONEHOT, requires_grad=False)
        self.seq_len = seq_len
        self.use_reg_cell = use_reg_cell
        self.epifeature_dim = epifeature_dim  # 动态调整通道数

        # 1
        input_channels = 4 + epifeature_dim  # 输入通道数为4 (one-hot) + epifeature通道数
        current_len = seq_len
        self.pre_conv = nn.Sequential(                  # input: (batch_size, 4 + epifeature_dim, seq_len)
            nn.Conv1d(input_channels, out_channels=288, kernel_size=17, padding=8),
            nn.BatchNorm1d(288),
            nn.MaxPool1d(kernel_size=3),                # output: (batch_size, 288, 448)
            nn.ReLU(),
        )
        current_len = current_len // 3

        # 2
        self.conv_towers = []
        self.conv_towers.append(ConvTower(288, 64, 5))  # output: (batch_size, 64, 224)
        current_len = current_len // 2
        self.conv_towers.append(ConvTower(64, 128, 5))  # output: (batch_size, 128, 112)
        current_len = current_len // 2
        self.conv_towers.append(ConvTower(128, 256, 5)) # output: (batch_size, 256, 56)
        current_len = current_len // 2
        self.conv_towers.append(ConvTower(256, 512, 5)) # output: (batch_size, 512, 28)
        current_len = current_len // 2
        self.conv_towers = nn.Sequential(*self.conv_towers)

        # 3
        self.post_conv = nn.Sequential(
            nn.Conv1d(512, 256, kernel_size=1),
            nn.BatchNorm1d(256),
            nn.ReLU(),                                  # output: (batch_size, 256, 28)
        )
        current_len = current_len // 1

        # 4
        self.flatten = nn.Flatten()                     # output: (batch_size, 7168)

        current_len = current_len * 256

        # 5
        self.dense = nn.Sequential(
            nn.Linear(current_len, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.ReLU(),
            nn.Dropout(0.2),                            # output: (batch_size, 64)
        )

        # 6 
        self.cell_embedding = nn.Linear(hidden_size, n_cells)   # output: (batch_size, 10537)
    
    def get_embedding(self):
        return self.cell_embedding.state_dict()["weight"]
    
    
    def forward(self, sequence: Tensor, epifeature: Tensor) -> Tensor:
        """
        sequence: (batch_size, seq_len) 输入DNA序列
        epifeature: (batch_size, epifeature_dim, seq_len ) 输入附加特征
        """
        # transform sequence to one-hot vector
        sequence = self.onehot[sequence.long()].transpose(1, 2)  # (batch_size, 4, seq_len)

        # 将sequence和epifeature拼接
        x = torch.cat((sequence, epifeature), dim=1)  # 拼接后的shape: (batch_size, 4 + epifeature_dim, seq_len)

        x = self.pre_conv(x)
        x = self.conv_towers(x)
        x = self.post_conv(x)
        x = self.flatten(x)
        x = self.dense(x)
        logits = self.cell_embedding(x)
        
        if self.use_reg_cell:
            lr_reg_cell = torch.norm(self.cell_embedding.weight, p=2) + torch.norm(self.cell_embedding.bias, p=2)
        else:
            lr_reg_cell = None
        return logits, lr_reg_cell

