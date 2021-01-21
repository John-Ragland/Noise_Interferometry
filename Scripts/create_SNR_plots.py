import os
import sys
cwd = os.getcwd()
ooipy_dir = os.path.dirname(os.path.dirname(cwd)) + '/ooipy'
sys.path.append(ooipy_dir)
# print(ooipy_dir)
import ooipy

from matplotlib import pyplot as plt
import datetime
import numpy as np
from obspy import read,Stream, Trace
from obspy.core import UTCDateTime
# from ooipy.hydrophone import Noise_Interferometry as NI
import pickle
import scipy
from gwpy.timeseries import TimeSeries
import seaborn as sns
import gwpy
import progressbar
import scipy

cwd = os.getcwd()
ni_dir = os.path.dirname(os.path.dirname(cwd))
sys.path.append(ni_dir)
# print(ni_dir)
from Noise_Interferometry.Modules import analysis

import time as tm


# Create Spectrogram for Time Range
file_base = "/Users/jhrag/Code/ocean_acoustics/Noise_Interferometry/NCFs/MJ03F-MJ03E/"
file_name_2016 = file_base + '2016'
file_name_2017 = file_base + '2017'
file_name_2018 = file_base + '2018'

# Create Experiment for 2017
exp1 = analysis.NCCF_experiment(file_name_2017)

peak_names = ['dA', 's1b0A', 's2b1A', 's3b2A', 'dB', 's1b0B', 's2b1B', 's3b2B']

avg_time = 584
stride = 146
overlap = 1/4

SNRs = np.ma.zeros((8,int(8760/stride - 3),avg_time+1))

# change to 8
for n in range(8):
    #change to 47
    for k in range(int(8760/stride - 3)):
        start = round(k*avg_time*overlap)
        end = start + avg_time
        print(f'\nCalculating SNR for {peak_names[n]}: hours {start} to {end} ({k+1}/{int(8760/stride - 3)})')
        tm.sleep(0.01)
        SNRs[n,k,:] = (exp1.SNR_plot(start,end, peak_id=peak_names[n], plot=False))
    
# Save SNRs to .pkl file
with open('SNRs.pkl', 'wb') as f:
    pickle.dump(SNRs, f)