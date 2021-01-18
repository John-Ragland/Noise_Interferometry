import os
import sys
cwd = os.getcwd()
ooipy_dir = os.path.dirname(os.path.dirname(cwd)) + '/ooipy'
sys.path.append(ooipy_dir)

from matplotlib import pyplot as plt
import datetime
import numpy as np
from obspy import read,Stream, Trace
from obspy.core import UTCDateTime
from ooipy.hydrophone import Noise_Interferometry as NI
import pickle
import scipy
from gwpy.timeseries import TimeSeries
import seaborn as sns
import gwpy
import progressbar
import scipy

cwd = os.getcwd()
ni_dir = os.path.dirname(cwd)
sys.path.append(ni_dir)

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
SNRs = []
for n in range(8):
    SNR = []
    for k in range(47):
        start = round(k*365/2)
        end = start + 365
        print(f'\nCalculating SNR for {peak_names[n]}: hours {start} to {end} ({k+1}/47)')
        tm.sleep(0.01)
        SNR.append(exp1.SNR_plot(start,end, peak_id=peak_names[n], plot=False))
    SNR = np.array(SNR)
    SNRs.append(SNR)

mdic = {"SNRs": SNRs}
scipy.io.savemat('2017_SNR_plot.mat',mdic)
# SNR = exp1.SNR_plot(600,900, peak_id='s1b0B', plot=False, )