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
from ooipy.hydrophone import Noise_Interferometry
import time
import gc

num_periods = 168
avg_time = 10 #minutes

start_time = datetime.datetime(2017,3,10,0,0,0) # time of first sample
node1 = '/LJ01C'
node2 = '/PC01A'
filter_data = True
'''
#Initialize Hydrophone Xcorr Class
hyd_xcorr = Noise_Interferometry.Hydrophone_Xcorr(node1=node1, node2=node2, avg_time=avg_time, mp=True, W=60, filter_data=filter_data, ckpts=True)

# Average Cross Correlation Over Specified time
t, xcorr_stack, bearing_max = hyd_xcorr.avg_over_mult_periods(num_periods, start_time)

'''




'''
Computes average over num_periods of averaging periods

Inputs:
num_periods - number of periods to average over
start_time - start time for data

Outputs:
xcorr - average xcorr over num_periods of averaging
'''
count = 0

for k in range(num_periods):
    #Initialize Hydrophone Xcorr Class
    hyd_xcorr = Noise_Interferometry.Hydrophone_Xcorr(node1=node1, node2=node2, avg_time=avg_time, mp=True, W=60, filter_data=filter_data)
    hyd_xcorr.count = count
    stopwatch_start = time.time()
    print('\n\nTime Period: ',k + 1)
    
    h1, h2, flag = hyd_xcorr.get_audio(start_time)

    if flag == None:
        print(f'{k+1}th average period skipped, no data available')
        continue
    
    h1_processed, h2_processed = hyd_xcorr.preprocess_audio(h1,h2)
    
    hyd_xcorr.xcorr_over_avg_period(h1_processed, h2_processed)

    count = count+1

    del hyd_xcorr
    gc.collect()

