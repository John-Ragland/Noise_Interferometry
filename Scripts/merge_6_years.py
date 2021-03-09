import os
import sys

from matplotlib import pyplot as plt
import datetime
import numpy as np
import pickle
import scipy
import seaborn as sns
import progressbar

from Noise_Interferometry.Modules import analysis

file_name = "/Volumes/Ocean_Acoustics/NCCFs/MJ03F-MJ03E/all_years"
exp1 = analysis.NCCF_experiment(file_name, verbose=False)

NCCFs = exp1.MA_TDGF(201,1,start_hour=0,end_hour=52607)

# Write to pickle file
file_name = "/Volumes/Ocean_Acoustics/NCCFs/MJ03F-MJ03E/long_term_MA_allYears_201.pkl"

with open(file_name, 'wb') as f:
    pickle.dump(NCCFs, f)