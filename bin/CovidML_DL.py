import os
from os.path import join
import sys

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
import torch.optim as optim
from torchvision import datasets, transforms
from torch.autograd import Variable
from sklearn.metrics import roc_auc_score
from sklearn import metrics


# Here we store the ML model for the deep learning task.
class C_FCClassifier(nn.Module):
    def __init__(self):
        super(C_FCClassifier, self).__init__()
        self.fc1 = nn.Linear(251, 256)
        self.bn1 = nn.BatchNorm1d(256)
        self.fc2 = nn.Linear(256, 256)
        self.bn2 = nn.BatchNorm1d(256)
        self.fc3 = nn.Linear(256, 256)
        self.bn3 = nn.BatchNorm1d(256)
        self.fc4 = nn.Linear(256, 128)
        self.bn4 = nn.BatchNorm1d(128)
        self.fc5 = nn.Linear(128, 64)
        self.bn5 = nn.BatchNorm1d(64)
        self.fc6 = nn.Linear(64, 32)
        self.bn6 = nn.BatchNorm1d(32)
        self.fc7 = nn.Linear(32, 1)

    def forward(self, x):
        x = self.fc1(x)
        self.bn1(x)
        x = F.leaky_relu(x)

        x = self.fc2(x)
        self.bn2(x)
        x = F.leaky_relu(x)

        x = self.fc3(x)
        self.bn3(x)
        x = F.leaky_relu(x)

        x = self.fc4(x)
        self.bn4(x)
        x = F.leaky_relu(x)

        x = self.fc5(x)
        self.bn5(x)
        x = F.leaky_relu(x)

        x = self.fc6(x)
        self.bn6(x)
        x = F.leaky_relu(x)

        x = self.fc7(x)

        return x.squeeze()


class CA_FCClassifier(nn.Module):
    def __init__(self):
        super(CA_FCClassifier, self).__init__()
        self.fc1 = nn.Linear(252, 256)
        self.bn1 = nn.BatchNorm1d(256)
        self.fc2 = nn.Linear(256, 256)
        self.bn2 = nn.BatchNorm1d(256)
        self.fc3 = nn.Linear(256, 256)
        self.bn3 = nn.BatchNorm1d(256)
        self.fc4 = nn.Linear(256, 128)
        self.bn4 = nn.BatchNorm1d(128)
        self.fc5 = nn.Linear(128, 64)
        self.bn5 = nn.BatchNorm1d(64)
        self.fc6 = nn.Linear(64, 32)
        self.bn6 = nn.BatchNorm1d(32)
        self.fc7 = nn.Linear(32, 1)

    def forward(self, x):
        x = self.fc1(x)
        self.bn1(x)
        x = F.leaky_relu(x)

        x = self.fc2(x)
        self.bn2(x)
        x = F.leaky_relu(x)

        x = self.fc3(x)
        self.bn3(x)
        x = F.leaky_relu(x)

        x = self.fc4(x)
        self.bn4(x)
        x = F.leaky_relu(x)

        x = self.fc5(x)
        self.bn5(x)
        x = F.leaky_relu(x)

        x = self.fc6(x)
        self.bn6(x)
        x = F.leaky_relu(x)

        x = self.fc7(x)

        return x.squeeze()


def map_cohort_labels(cohort_label):
    mapped_label = 0

    if cohort_label == 'Severe':
        mapped_label = 1
    elif cohort_label == 'Mild':
        mapped_label = 0
    else:
        mapped_label = None
    return mapped_label


def get_torch_data(dataframe_data, data_col_index_start=2, batch_size_in=128):
    TensorX = torch.Tensor(dataframe_data.iloc[:, data_col_index_start:].to_numpy())
    TensorY = torch.Tensor(dataframe_data.apply(lambda x: map_cohort_labels(x['Cohort']), axis=1).to_numpy())
    CovidDataSet = TensorDataset(TensorX, TensorY)  # create your datset
    CovidTrainDataLoader = DataLoader(CovidDataSet, batch_size=batch_size_in, shuffle=True, drop_last=True)

    return CovidDataSet, CovidTrainDataLoader, TensorX, TensorY

