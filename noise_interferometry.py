import os
import sys
cwd = os.getcwd()
ooipy_dir = os.path.dirname(cwd) + '/ooipy'
sys.path.append(ooipy_dir)
from matplotlib import pyplot as plt
import datetime
import numpy as np
from obspy import read,Stream, Trace
from obspy.core import UTCDateTime
from ooipy.hydrophone import Noise_Interferometry as NI
import pickle

num_periods = 600
avg_time = 60  #minutes 
start_time = datetime.datetime(2017,7,1,0,0,0) # time of first sample
node1 = 'Central_Caldera'
node2 = 'Eastern_Caldera'
filter_cutoffs = np.array([1, 90])
W = 30
htype = 'low_frequency'
whiten= True
kstart= 0
other_notes = 'Butterworth and Hann Window for whitening. No Low bound on frequency (HARDCODED CHANGE!)'
sp_method = 'phase_weight'

NI.calculate_NCF_loop(num_periods, node1, node2, avg_time, start_time, W, filter_cutoffs, verbose=True, whiten=whiten, htype=htype, kstart=kstart, sp_method = sp_method, other_notes=other_notes)

