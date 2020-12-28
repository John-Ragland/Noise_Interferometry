import os
import sys
from matplotlib import pyplot as plt
import datetime
import numpy as np
from obspy import read,Stream, Trace
from obspy.core import UTCDateTime
import pickle

cwd = os.getcwd()
ni_dir = os.path.dirname(cwd)
sys.path.append(ni_dir)

from Noise_Interferometry import calculate
# sys.path.append(ooipy_dir)

num_periods = 5
avg_time = 60  #minutes 
start_time = datetime.datetime(2017,1,1,0,0,0) # time of first sample
node1 = 'Central_Caldera'
node2 = 'Eastern_Caldera'
filter_cutoffs = np.array([1, 90])
W = 30
htype = 'low_frequency'
whiten= True
kstart = 0
other_notes = 'Butterworth and Hann Window for whitening. Longterm Analysis'
sp_method = 'sabra_b'

calculate.calculate_NCF_loop(num_periods, node1, node2, avg_time, start_time, W, filter_cutoffs, verbose=True, whiten=whiten, htype=htype, kstart=kstart, sp_method = sp_method, other_notes=other_notes)

