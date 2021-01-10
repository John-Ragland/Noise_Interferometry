'''
Library of functions (maybe will be changed to a class in future) to do peak
analysis on long term ambient noise cross correlation functions
'''

def init_from_NCCF_experiment(exp):
    '''
    init_from_NCCF_experiment - create NCCF_exp object from NCCF_experiment data
    object

    personal notes - it might make sense to combine these two data types in the
        future (or at least rename them so they don't basically have the 
        same exact name)

    Parameters
    ----------
    exp : Noise_Interferometry.Modules.analysis.NCCF_experiment object
        data type containing Noise Interferometry Experiement

    Returns
    -------
    NCCF_exp : Noise_Interferometry.Modules.peak_analysis.NCCF_exp
        data type containing Noise Interferometry Experiement...
            (Kind of confusing I know......)
    '''

    NCCF = exp.



    