import os
import sys
cwd = os.getcwd()
ooipy_dir = os.path.dirname(os.path.dirname(cwd)) + '/ooipy'
sys.path.append(ooipy_dir)

import ooipy
from matplotlib import pyplot as plt
import datetime
import numpy as np
from obspy import read,Stream, Trace
from obspy.core import UTCDateTime
import pickle
import scipy
from gwpy.timeseries import TimeSeries
import seaborn as sns
import gwpy
import progressbar
import scipy.io

cwd = os.getcwd()
ni_dir = os.path.dirname(cwd)
sys.path.append(ni_dir)

from Modules import analysis

code_dir = os.path.dirname(os.path.dirname(cwd))

sys.path.append(code_dir)
from compute_spectrograms import spec_tools

def plot_SNR(time, freq, values, SNR, file_name, start_time):


    res_reduction_freq = 1


    font = {'size': 22}
    matplotlib.rc('font', **font)
            
    v = values[::res_reduction_time, ::res_reduction_freq]

    if len(time) != len(values):
        t = np.linspace(0, len(values) - 1,
                        int(len(values) / res_reduction_time))
    else:
        t = time[::res_reduction_time]
                
            
    if len(freq) != len(values[0]):
        f = np.linspace(0, len(values[0]) - 1,
                        int(len(values[0]) / res_reduction_freq))
    else:
        f = freq[::res_reduction_freq]
                
    cbarticks = np.arange(vmin, vmax + vdelta, vdelta)
    fig, ax = plt.subplots(figsize=(figsize), dpi=dpi)
    plt.xticks(rotation=xlabel_rot)
    im = ax.contourf(t, f, np.transpose(v), cbarticks,
                    norm=colors.Normalize(vmin=vmin, vmax=vmax),
                    cmap=plt.cm.jet)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.ylim([fmin, fmax])
    if time_limits is not None:
        plt.xlim(time_limits)
    
    title = f'SNR Plot starting {start_time.strftime("%m/%d/%Y, %H:%M:%S")}'
    plt.title(title)
    cbar = plt.colorbar(im, ax=ax,
                ticks=np.arange(vmin, vmax + vdelta, vdelta_cbar), pad=0.1)

    cbar.ax.set_ylabel('PSD Amplitude (dB)', rotation=90, fontsize='16')
    plt.tick_params(axis='y')

    # Add SNR Plot
    
    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'purple'
    ax2.set_ylabel('SNR (dB)', color=color)  # we already handled the x-label with ax1
    plt.xticks(rotation=xlabel_rot)
    ax2.plot(SNR, color=color, linewidth=3)
    ax2.tick_params(axis='y', labelcolor=color)
    
    fig.savefig(filename,dpi=dpi)

    plt.close('all')

    return

# Load NCCF Experiment
file_name = "/Users/jhrag/Code/Noise_Interferometry/NCFs/MJ03F-MJ03E/2016"
exp1 = analysis.NCCF_experiment(file_name, verbose=True)


# Calculate SNR and Merge Spectrograms for Time Period

for k in range(24):
    start = k*366
    end = (k+1)*366
    print(f'\nCalculating SNR for hours {start} to {end}...')
    SNR = (exp1.SNR_plot(start,end, plot=False))
    
    filename = f'SNR_Plots/SNR_only/SNR_hr{k*365:03}-{(k+1)*365}_2016_Central.png'

    fig = plt.figure(figsize=(7,5))
    plt.plot(SNR)
    plt.xlabel('Average Time (Hours)')
    plt.ylabel('SNR (dB)')
    plot_SNR(times_new, freq, values, SNR, file_name=filename, start_time = times[0])






