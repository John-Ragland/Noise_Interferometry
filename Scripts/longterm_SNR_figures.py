#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Noise_Interferometry
import pickle
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import math as m

# In[ ]:

peak_names = ['dA','s1b0A','s2b1A','dB','s1b0B','s2b1B']
avg_times = [201, 601]
years = [2015, 2016, 2017, 2018, 2019, 2020]
SNRs = []
peak_idxs = []
noise_stds = []
dates = []
num_availables = []

for avg_time in avg_times:
    print(f'Calculating for {avg_time}...')
    for peak_name in peak_names:
        print(f'Calculating for {peak_name}...')

        for year in years:
            print(f'Loading {year}...')
            file_name = f'/Volumes/Ocean_Acoustics/NCCFs/MJ03F-MJ03E/Longterm_MA_{year}_{avg_time}.pkl'
            long_term = pickle.load(open(file_name,'rb'))
            dates.append(long_term.dates)
            num_availables.append(np.array(long_term.num_available))
            SNR, peak_idx, noise_std = long_term.snr_of_peak_amp(peak_name)
            SNRs.append(SNR)
            peak_idxs.append(peak_idx)
            noise_stds.append(noise_std)


        # ### Set Number of Hours Threshold to 100 and create SNR masked arrays

        # In[ ]:


        threshold = 100
        mask = []
        SNRs_mask = []
        for k in range(6):
            mask.append(~(num_availables[k] >= threshold))
            SNRs_mask.append(np.ma.masked_array(SNRs[k],mask=mask[k]))


        # In[ ]:


        SNRs_mask[0].shape


        # In[ ]:


        fig, axes = plt.subplots(2,3,figsize=(15,5))
        k = 0
        for ax in fig.get_axes():
            ax.plot(dates[k], SNRs_mask[k])

            plt.sca(ax)
            plt.grid()
            plt.title(f'SNR for {years[k]}')
            plt.xticks(rotation=20)
            plt.ylabel('SNR [dB]')
            k += 1
            
        plt.tight_layout()

        fig.savefig(f'longterm_SNR_figures/all_SNRs_201_{peak_name}.png',dpi=300)


        # In[ ]:


        fig, axes = plt.subplots(2,3,figsize=(15,5))
        k = 0
        for ax in fig.get_axes():
            ax.plot(dates[k], num_availables[k])

            plt.sca(ax)
            plt.grid()
            plt.title(f'Available Hours {years[k]}')
            
            k += 1
            
        plt.tight_layout()


        # # Seperate into Months

        # In[ ]:


        SNR_month_all = []
        for k in range(6):
            num_hours = len(SNRs_mask[k])
            num_split = m.trunc(num_hours/12)*12
            
            x = np.split(SNRs_mask[k][:num_split],12)
            
            SNR_month = np.zeros((12,1))
            for n in range(12):
                SNR_month[n] = np.mean(x[n])
                
            SNR_month_all.append(SNR_month)


        # In[ ]:


        SNR_month_all_array = np.ndarray.flatten(np.array(SNR_month_all))
        years = [2015, 2016, 2017, 2018, 2019, 2020]
        years2 = [val for val in years for _ in range(12)]
        months = ['January','February','March','April','May','June','July','August','September','October','November','December']*6

        d = {'year':years2, 'month':months, 'SNR':SNR_month_all_array}


        # In[ ]:


        df = pd.DataFrame.from_dict(data=d)
        df


        # In[ ]:


        cp = sns.color_palette('Set2')
        fig, ax = plt.subplots(1,1,figsize=(10,7))

        ax = sns.barplot(x="month", y="SNR", hue="year", data=df, palette=cp)
        plt.sca(ax)
        plt.xticks(rotation=45)
        plt.legend(loc="lower right", ncol=len(df.columns))
        plt.title('Average SNR for Each Month')
        plt.ylabel('SNR [dB]')

        fig.savefig(f'longterm_SNR_figures/SNR_for_every_month_{peak_name}.png',dpi=200)


        # In[ ]:


        cp = sns.color_palette('Set2')
        fig, ax = plt.subplots(1,1,figsize=(10,7))

        ax = sns.barplot(x="month", y="SNR", data=df, palette=cp)
        plt.sca(ax)
        plt.xticks(rotation=45)
        plt.title('Average SNR for Each Month')
        plt.ylabel('SNR [dB]')

        fig.savefig(f'longterm_SNR_figures/SNR_for_every_month_avg_{peak_name}.png',dpi=200)


        # In[ ]:




