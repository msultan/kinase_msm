#!/bin/evn python
from kinase_msm.data_loader import load_yaml_file
import numpy as np
import os
from kinase_msm.mdl_analysis import ProteinSeries, Protein
from kinase_msm.data_loader import load_frame
from kinase_msm.data_transformer import create_assignment_matrix, create_tics_array