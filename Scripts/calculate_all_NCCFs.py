import os
import sys
from matplotlib import pyplot as plt
import datetime
import numpy as np
# from obspy import read,Stream, Trace
# from obspy.core import UTCDateTime
import pickle

cwd = os.getcwd()
code_dir = os.path.dirname(os.path.dirname(cwd))
print(code_dir)
sys.path.append(code_dir)
from Noise_Interferometry.Modules import calculate
# sys.path.append(ooipy_dir)


file_base = '/Volumes/Ocean_Acoustics/NCCFs/MJ03F-MJ03E/'
years = [2015,2016,2017,2018,2019,2020]

for year in years:
    file_dir = f'{file_base}{year}/'
    if year % 4 == 0:
        leap = True
    else:
        leap = False

    if leap:
        num_periods = 8784
    else:
        num_periods = 8760

    avg_time = 60  #minutes 
    start_time = datetime.datetime(year,1,1) # time of first sample
    node1 = 'Central_Caldera'
    node2 = 'Eastern_Caldera'
    filter_cutoffs = np.array([1, 90])
    W = 30
    htype = 'low_frequency'
    whiten= True
    if year == 2015:
        kstart = 4371
    else:
        kstart = 0

    other_notes = 'Looking at Airgun Experiments'
    sp_method = 'sabra_b'

    calculate.calculate_NCF_loop(num_periods, node1, node2, avg_time, start_time, W, filter_cutoffs, file_dir=file_dir, verbose=True, whiten=whiten, htype=htype, kstart=kstart, sp_method = sp_method, other_notes=other_notes)
