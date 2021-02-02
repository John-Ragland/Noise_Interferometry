'Create Moving Average NCCF_array object for each year [2015-2020]'
import sys
sys.path.append('/Users/jhrag/Code/ocean_acoustics')
from Noise_Interferometry.Modules import analysis
import pickle


years = [2015, 2016, 2017, 2018, 2019, 2020]

nccfs_array_ls = []

avg_time = 601

for year in years:
    if year % 4 == 0:
        end_hour = 366*24
    else:
        end_hour = 365*24
    file_name = "/Volumes/Ocean_Acoustics/NCCFs/MJ03F-MJ03E/" + str(year)
    exp1 = analysis.NCCF_experiment(file_name, verbose=True)
    nccfs_array = exp1.MA_TDGF(avg_time,1,start_hour=0,end_hour=end_hour, verbose=True)

    # Save NCCFs to Pickle File
    pickle_name = f"/Volumes/Ocean_Acoustics/NCCFs/MJ03F-MJ03E/Longterm_MA_{year}_{avg_time}.pkl"
    with open(pickle_name, 'wb') as f:
        pickle.dump(nccfs_array, f)
