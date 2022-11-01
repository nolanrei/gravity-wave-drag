import numpy as np
import torch
from torch import nn
import matplotlib.pyplot as plt
import xarray as xr

# Basic neural network architecture for learning gravity wave drag
class GravNet(nn.Module):
    def __init__(self, params):
        super().__init__()
        self.radius = params[0]         # number of grid cells we're allowed to see in each direction
        self.nz = params[1]             # number of grid cells in vertical column
        self.n2dvars = params[2]        # number of 2d input fields
        self.n3dvars = params[3]        # number of 3d input fields
        self.has_loc = params[4]        # 0 or 1; can this model see exact lat/lon

        # Number of features to expect as input
        # For computational ease, uses L_\infty ball to determine distance on grid WRT self.radius
        self.ninputs = (2*self.radius+1)**2*(self.nz*self.n3dvars + self.n2dvars + 2*self.has_loc)
        self.nhid1 = 10*self.ninputs
        self.nhid2 = self.nhid1
        self.noutputs = 3*self.nz       # (drag_x, drag_y, drag_z) at each level in the column

        # Network layers
        self.fc1 = nn.Linear(self.ninputs, self.nhid1)
        self.fc2 = nn.Linear(self.nhid1, self.nhid2)
        self.fc3 = nn.Linear(self.nhid2, self.noutputs)
        self.relu = nn.ReLU()

    def forward(self, x):
        return self.fc3(self.relu(self.fc2(self.relu(self.fc1(x)))))


