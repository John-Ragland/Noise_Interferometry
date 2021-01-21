import numpy as np
import os
import pickle
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.animation as animation
import matplotlib
import progressbar
import scipy
from scipy import signal
import datetime
import warnings

# Define manually determined peak slices (in delay time dimension)
peak_windows = [
    [5535, 5635],
    [5303, 5503],
    [4969, 5169],
    [4622, 4782],
    [6367, 6467],
    [6498, 6698],
    [6833, 7033],
    [7242, 7402],
]

slices = []
for k in range(len(peak_windows)):
    slices.append(np.s_[peak_windows[k][0]:peak_windows[k][1]])

peak_names = ['dA', 's1b0A', 's2b1A', 's3b2A', 'dB', 's1b0B', 's2b1B', 's3b2B']
peak_slices = dict(zip(peak_names,slices))


def snr_of_single_NCCF(NCCF, peak_id, t):
    '''
    Parameters
    ----------
    NCCF : numpy array
        single NCCF of shape [1,n]
    peak_id : str
        peak identifier string
    t : numpy array
        time array (usually shape [19999,])

    Returns
    -------
    SNR : float
        SNR of peak in dB
    '''
    
    noise_bounds = np.array([-1.5, 1.5])
    NCCF_c = scipy.signal.hilbert(NCCF)

    noise_idx = np.array([
        np.argmin(np.abs(t - noise_bounds[0])),
        np.argmin(np.abs(t - noise_bounds[1]))
    ])
    
    noise_slice = np.s_[noise_idx[0]:noise_idx[1]]
    noise_std = np.std(np.abs(NCCF_c[noise_slice]))

    peak_amp = np.max(np.abs(NCCF_c[peak_slices[peak_id]]))
    SNR = 20*np.log10(peak_amp/noise_std)

    return SNR


class NCCF_experiment:
    '''
    class with tools for analyzing NCFs. NCF experiments are stored in ./NCFs.

    Attributes
    ----------
    ckpt_dir : str
        name of specific NCF directory. NCF .pkl files will be located in ckpt_dir
        where ckpt_dir is an absolute path
    start_time : datetime.datetime
        specifies start time of Noise Interferometry experiment
    filter_cutoffs : list
        specifies low and high frequencies cutoffs used
    htype : str
        specifies whether hydrophone in 'broadband' or 'low_frequency'
    num_periods : int
        specifices number of hours that noise interferometry experiment lasted
    node1 : str
        location code for node 1
    node2 : str
        location code for node 2
    W : int
        length of short time window in seconds
    distance_between_nodes : float
        distance between the two nodes (in km)
    xcorr : numpy array
        NCCF created from pickle files

    Methods
    -------
    __read_header()
        reads header from ckpt_dir and prints pertanent information

    NCF_plot(num_hours)
        create average NCF for num_hours (starting from beginning) and plot
        the results

    average_NCF(num_hours)
        averages NCF for num_hours (starting from beginning)

    animate_NCF(num_hours, direction, animation_len)

    available_hours(max_number)
        loops through all ckpts in ckpt_dir and reports number of valid hours.
        it also creates a plot of availability for the time frame.

    edit_header()
        changes header. This method contains hard coding and should be used with extrem caution

    animate_phase()
        sweeps through 360 degree phase shift using Hilbert Transform and creates animation.
    
    SNR_plot(start, end, mode=0, t_peak=None, plot=False, savefig=False)
        creates plot of SNR verses average time

    
    '''
    def __init__(self, ckpt_dir, verbose=False):
        self.ckpt_dir = ckpt_dir

        self.__read_header(verbose=verbose)

        self.__distance_between_hydrophones()

        self.__bearing_between_hydrophones()


        if verbose:
            print(f'Distance Between Nodes: {self.distance_between_nodes/1000:4.4} km')
            print(f'Bearing from node 1 to 2: {self.bearing_node_1_2:4.4} °')
        return


    def __read_header(self, verbose=False):
        header_name = self.ckpt_dir + '/0HEADER.pkl'
        with open(header_name, 'rb') as f:
            header = pickle.load(f)

        self.start_time = header.start_time
        self.filter_cutoffs = header.filter_cutoffs
        self.htype = header.htype
        self.num_periods = header.num_periods
        self.node1 = header.node1
        self.node2 = header.node2
        self.W = header.W
        if header.htype == 'low_frequency':
            self.Fs = 200
        elif header.htype == 'broadband':
            self.Fs = 64000
        if verbose:
            print(f'Start Time: {header.start_time}')
            print(f'Frequency Range: {header.filter_cutoffs}')
            print(f'Hydrophone Type: {header.htype}')
            print(f'Number of Average Periods: {header.num_periods}')
            print(f'Node 1: {header.node1}')
            print(f'Node 2: {header.node2}')
            print(f'Window Length (s): {header.W}')
            print(f'Signal Processing Method: {header.sp_method}')
            try:
                print(f'Specific Notes: {header.other_notes}')
                return header
            except:
                return header
        else:
            return header


    def available_hours(self, savefig=False):
        NCF_available = np.zeros((self.num_periods,1))
        num_available = 0

        for k in range(self.num_periods):
            ckpt_name = self.ckpt_dir + '/ckpt_' + str(k) + '.pkl'
            try:
                with open(ckpt_name, 'rb') as f:
                    xcorr_1hr = pickle.load(f)
                NCF_available[k] = 1
                num_available = num_available + 1
            except:
                NCF_available[k] = 0

        fig = plt.figure(figsize=(20,0.5))
        ax = plt.imshow(NCF_available.T, aspect='auto')
        ax.axes.get_yaxis().set_visible(False)
        plt.xlabel('Hours')
        self.num_available = num_available

        if savefig:
            fig.savefig('available_hours.png')


        print(f'Number of Available hours: {num_available}')


    def average_NCF(self, hour_start, hour_end, plot=False, verbose=False):
        '''
        reads pickle file and averages between hour_start and hour_end and
        returns numpy array of shape [M,] where M is the number of points in
        time dimension (normally 11999 for low frequency and 30 s w)

        Parameters
        ----------
        hour_start : int
            start hour
        hour_end : int
            end hour
        plot : bool
            indicates whether to plot NCCF or not
        verbose : bool
            indicates whether to print updates or not

        Returns
        -------
        xcorr_avg : numpy array
            average NCCF between hour_start and hour_end (shape [M,])
        num_available : float
            number of hours within average period that are not invalid
        '''
        
        invalid = 0

        for k in range(hour_start,hour_end+1):
            # Read Pickle File
            ckpt_name = self.ckpt_dir + '/ckpt_'+ str(k) +'.pkl'
            try:
                with open(ckpt_name, 'rb') as f:
                    xcorr_1hr = pickle.load(f)
            # if pickle file read results in error
            except:
                invalid = invalid + 1
                if k == hour_end:
                    self.SNR_plot_flag = True
                continue
            
            # if file contains any NaN Value
            if np.isnan(np.sum(xcorr_1hr)):
                invalid = invalid + 1
                if k == hour_end:
                    self.SNR_plot_flag = True
                continue

            # No errors (either add to xcorr or create xcorr)
            else:
                try:
                    xcorr = xcorr + xcorr_1hr
                except NameError:
                    xcorr = xcorr_1hr

                self.SNR_plot_flag = False

        # Check if variable xcorr exists, if not create NaN entry
        try:
            xcorr
        except NameError:
            xcorr_len = 2*self.W*self.Fs - 1
            xcorr = np.empty(xcorr_len)
            xcorr[:] = np.NaN
            if verbose: print('\n Entire Average Period Invalid')
        
        num_available = (hour_end+1) - hour_start - invalid
        self.num_available = num_available
        # make average instead of sum
        xcorr_avg = xcorr / num_available #changed from num_available

        self.xcorr = xcorr_avg

        # create time variable in NCCF_experiment class
        dt = 1/self.Fs
        self.t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)

        if plot:
            self.NCF_plot(xcorr_avg, num_available)

        return xcorr_avg, num_available


    def NCF_plot(self, xcorr, num_valid, save_fig=False, file_name=None, frequency=False, xlimits=None, ylimits=None, symetric=False, print_time_delay=False, title=None):

        if 'self.xcorr' not in locals():
            Exception: 'Must Calculate average NCF first. use average_NCF(num_hours) method'
        t = self.t
        f = np.linspace(0,self.Fs,xcorr.shape[0])
        node1 = self.node1.replace("_"," ")
        node2 = self.node2.replace("_"," ")

        sns.set()
        f1 = plt.figure(figsize=(7,5))
        if frequency:
            plt.plot(f,np.abs(scipy.fft.fft(xcorr)))
            plt.xlim([0,self.Fs/2])
            plt.title(f'{num_valid} Hours - {node1}/{node2} - Spectrum', fontsize=18, y=1.08)
            plt.xlabel('Frequency (Hz)', fontsize=14)
            plt.ylabel('Average NCF Spectrum Amplitude', fontsize=14)

        elif symetric:
            import math as m
            xcorr1 = xcorr[:m.floor(len(xcorr)/2)]
            xcorr2 = xcorr[m.ceil(len(xcorr)/2):]

            xcorr_sym = xcorr2 + np.flip(xcorr1)
            t_sym = t[m.ceil(len(xcorr)/2):]
            plt.plot(t_sym, xcorr_sym)
            plt.title(f'{num_valid} Hours - {node1}/{node2} - Average NCF', fontsize=18, y=1.08)

            plt.xlabel('Delay τ (s)', fontsize=14)
            plt.ylabel('Average NCF Amplitude', fontsize=14)
            if xlimits is not None:
                plt.xlim(xlimits)
            if ylimits is not None:
                plt.ylim(ylimits)
        else:
            xcorr_c = signal.hilbert(xcorr)
            plt.plot(t, np.abs(xcorr_c))
            # plt.title(f'{num_valid} Hours - {node1}/{node2} - Average NCF', fontsize=18, y=1.08)

            plt.xlabel('Delay τ (s)', fontsize=14)
            plt.ylabel('Average NCF Amplitude', fontsize=14)
            if xlimits is not None:
                plt.xlim(xlimits)
            if ylimits is not None:
                plt.ylim(ylimits)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        if save_fig:
            f1.savefig(file_name + '.png', dpi=500)

        if print_time_delay:
            print(f'Time of Max Peak: {t[np.argmax(np.abs(xcorr))]}')


    def animate_NCF(self, num_hours, length, file_name, xlim=None, frequency=False, interval=1, dpi=None):
        interval=int(interval)

        # set y limit
        xcorr, num_available = self.average_NCF(1,2)
        ylim = np.max(np.abs(xcorr))

        # matplotlib.rcParams[figure.max_open_warning'] = 0
        if xlim == None:
            xlim = [-self.W,self.W]
        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line,time_text,

        def animate(i):
            i = interval*i
            dt = 1/self.Fs
            t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)
            f = np.linspace(0,self.Fs,xcorr.shape[0])
            try:
                y, num_available = self.average_NCF(1,i)
            except:
                y = np.squeeze(np.zeros((xcorr.shape[0],1)))

            line.set_data(t, y)
            if frequency:
                line.set_data(f,np.abs(scipy.fft.fft(y)))
            time_text.set_text(f'Hour: {i}')

            ax.set_xlim(xlim)
            bar.update(i)
            return line, time_text,

        xcorr = self.xcorr
        fig = plt.figure(figsize=(13.33,7.5))
        ax = plt.axes(xlim=(-self.W, self.W), ylim=(-ylim, ylim))
        line, = ax.plot([], [], lw=1)
        time_text = ax.text(0.01, 0.95, '', transform=ax.transAxes)

        bar = progressbar.ProgressBar(maxval=num_hours, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=int(num_hours/interval), interval=1, blit=False)
        if dpi == None:
            anim.save(file_name+'.mp4',fps=int((num_hours/interval)/length))
        else:
            anim.save(file_name+'.mp4',fps=int((num_hours/interval)/length), dpi=dpi)


    def animate_tdoa2(self, num_hours, length, file_name, xlim=None, frequency=False, interval=1, dpi=None):
        interval=int(interval)

        # set y limit
        xcorr, num_available = self.average_NCF(1,2)

        # matplotlib.rcParams[figure.max_open_warning'] = 0
        if xlim == None:
            xlim = [-self.W,self.W]
        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line,time_text,

        def animate(i):
            i = interval*i
            dt = 1/self.Fs
            t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)
            f = np.linspace(0,self.Fs,xcorr.shape[0])
            try:
                y, num_available = self.average_NCF(i,i+1)
            except:
                y = np.squeeze(np.zeros((xcorr.shape[0],1)))

            line.set_data(t, y)
            if frequency:
                line.set_data(f,np.abs(scipy.fft.fft(y)))
            time_text.set_text(f'Hour: {i}')

            ax.set_xlim(xlim)
            bar.update(i)
            return line, time_text,

        xcorr = self.xcorr
        fig = plt.figure(figsize=(13.33,7.5))
        ax = plt.axes(xlim=(-self.W, self.W))
        line, = ax.plot([], [], lw=1)
        time_text = ax.text(0.01, 0.95, '', transform=ax.transAxes)

        bar = progressbar.ProgressBar(maxval=num_hours, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=int(num_hours/interval), interval=1, blit=False)
        if dpi == None:
            anim.save(file_name+'.mp4',fps=int((num_hours/interval)/length))
        else:
            anim.save(file_name+'.mp4',fps=int((num_hours/interval)/length), dpi=dpi)


    def animate_tdoa(self, num_hours, length, file_name, xlim=None, frequency=False, interval=1, dpi=None):
        interval=int(interval)

        # set y limit
        xcorr, num_available = self.average_NCF(1,2)
        ylim = np.max(np.abs(xcorr))

        # matplotlib.rcParams[figure.max_open_warning'] = 0
        if xlim == None:
            xlim = [-self.W,self.W]
        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line,time_text,

        def animate(i):
            
            i = interval*i
            dt = 1/self.Fs
            t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)
            freq = np.linspace(0,self.Fs,xcorr.shape[0])
            
            ckpt_name = self.ckpt_dir + '/ckpt_'+ str(i) +'.pkl'
            try:
                with open(ckpt_name, 'rb') as f:
                    xcorr_1hr = pickle.load(f)
                    y=xcorr_1hr

            except:
                y = np.zeros(np.shape(t))

            if (np.sum(y) != 0) and (np.isnan(np.sum(xcorr_1hr))):
                y = np.zeros(np.shape(t))
            
            try:
                y, num_available = self.average_NCF(1,i)
            except:
                y = np.squeeze(np.zeros((xcorr.shape[0],1)))

            line.set_data(t, y)
            if frequency:
                line.set_data(freq,np.abs(scipy.fft.fft(y)))
            time_text.set_text(f'Hour: {i}')

            ax.set_xlim(xlim)
            bar.update(i)
            return line, time_text,

        xcorr = self.xcorr
        fig = plt.figure(figsize=(13.33,7.5))
        ax = plt.axes(xlim=(-self.W, self.W))
        line, = ax.plot([], [], lw=1)
        time_text = ax.text(0.01, 0.95, '', transform=ax.transAxes)

        bar = progressbar.ProgressBar(maxval=num_hours, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=int(num_hours/interval), interval=1, blit=False)
        if dpi == None:
            anim.save(file_name+'.mp4',fps=int((num_hours/interval)/length))
        else:
            anim.save(file_name+'.mp4',fps=int((num_hours/interval)/length), dpi=dpi)


    def edit_header(self):
        header = self.read_header()
        header.num_periods = 600

        header_name = self.ckpt_dir + '/0HEADER.pkl'

        with open(header_name, 'wb') as f:
            pickle.dump(header, f)


    def animate_phase(self, xcorr, file_name, xlim=None):
        '''
        animate full 360 degree phase shift of signal

        first method is just shifting phase in Frequecy domain
        '''
        length = 10 #length of animation in seconds
        # matplotlib.rcParams[figure.max_open_warning'] = 0
        if xlim == None:
            xlim = [-self.W,self.W]

        ylim = 1.5*np.max(np.abs(xcorr))

        def init():
            line.set_data([], [])
            phase_text.set_text('')
            return line, phase_text

        def animate(i):

            dt = 1/self.Fs
            t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)
            f = np.linspace(0,self.Fs,xcorr.shape[0])

            # calculate phase shifted signal
            xcorr_f = scipy.fft.fft(signal.hilbert(xcorr))
            xcorr_mag = np.abs(xcorr_f)
            xcorr_phase = np.angle(xcorr_f)
            xcorr_phase_new = xcorr_phase + np.deg2rad(i)

            xcorr_new = scipy.fft.ifft(xcorr_mag*np.exp(1j*xcorr_phase_new))

            y = xcorr_new

            line.set_data(t, y)

            phase_text.set_text(f'Phase Change: {i} $^\circ$')

            ax.set_xlim(xlim)
            bar.update(i)
            return line, phase_text,

        fig = plt.figure(figsize=(13.33,7.5))
        ax = plt.axes(xlim=(-self.W, self.W), ylim=(-ylim,ylim))
        line, = ax.plot([], [], lw=1)
        phase_text = ax.text(0.01, 0.95, '', transform=ax.transAxes)

        bar = progressbar.ProgressBar(maxval=360, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=360, interval=1, blit=False)

        anim.save(file_name+'.mp4',fps=360/length)

        return


    def SNR_plot(self, start, end, peak_id, plot=False, savefig=False,
        file_name=None):
        '''
        Attributes
        ----------
        start : int
            index of hour to start within directory
        end : int
            index of hour to end within directory
        peak_id : str
            specified which peak to calculate SNR for
        plot : bool
            specifices whether to plot SNR or not
        savefig : bool
            specifies whether to save plot to file
        file_name : str
            name of file to save. default None

        Returns
        -------
        t : numpy array
            time variable
        SNR_ma : numpy array
            masked array of valid SNR values
        '''
        count = 0
        signal_amp = np.zeros(((end-start)+1,1))
        noise_std = np.zeros(((end-start)+1,1))

        bar = progressbar.ProgressBar(maxval=end-start+1, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()

        t = np.linspace(-self.W,self.W,self.Fs*2*self.W-1)
        SNR = []
        last_valid = []
        for k in range(start, end+1):
            xcorr, num_available = self.average_NCF(start,k)
            last_valid.append(self.SNR_plot_flag)
            
            # Get SNR for peak_id
            SNR.append(snr_of_single_NCCF(xcorr, peak_id, t))

            count = count+1
            bar.update(count)

        if plot:
            sns.reset_orig()
            fig = plt.figure(figsize=(7,5))
            plt.plot(SNR, linewidth='2')
            #plt.title(f'Plot of SNR for {(end-start)+1} Hours', fontsize=18, y=1.08)

            plt.xlabel('Average Time (Hours)', fontsize=16)
            plt.ylabel('SNR (dB)', fontsize=16)

            if savefig:
                if file_name == None:
                    raise Exception('Invalid Filename')
                fig.savefig(file_name, dpi=300)

        # make SNR_ma - masked array
        SNR_ma = np.ma.masked_array(SNR,last_valid)
        
        return SNR_ma


    def xcorr_loop(self, start, end, plot=False, savefig=False,
        file_name=None):
        '''
        Attributes
        ----------
        start : int
            index of hour to start within directory
        end : int
            index of hour to end within directory
        plot : bool
            specifices whether to plot 2D xcorr or not
        savefig : bool
            specifies whether to save plot to file
        file_name : str
            name of file to save. default None

        Returns
        -------
        tau : numpy array
            array of time delay from cross correlation
        avg_t : numpy array
            array of average time
        xcorr2d : numpy array
            xcorr for different average times (2D array)
        '''
        count = 0
        num_available_list = []
        for k in range(start, end+1):


            xcorr, num_available = self.average_NCF(0,k)
            xcorr_single, num_available = self.average_NCF(0,k)

            count = count+1

            num_available_list.append(self.num_available)

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(xcorr2d)
            ax.set_aspect(aspect=100)
            plt.show()


        return num_available_list


    def __distance_between_hydrophones(self):
        """
        distance_between_hydrophones(coord1, coord2) - calculates the distance
        in meters between two global coordinates

        Inputs:
        coord1 - numpy array of shape [2,1] containing
            latitude and longitude of point 1
        coord2 - numpy array of shape [2,1] containing
            latitude and longitude of point 2

        Outputs:
        self.distance - distance between 2 hydrophones in meters
        self.time_delay - approximate time delay between 2 hydrophones
            (assuming speed of sound = 1480 m/s)

        """

        coord1, coord2 = self.__get_node_coords()

        from math import radians, cos, sin, asin, sqrt
        # The math module contains a function named
        # radians which converts from degrees to radians.
        lon1 = radians(coord1[1])
        lon2 = radians(coord2[1])
        lat1 = radians(coord1[0])
        lat2 = radians(coord2[0])

        # Haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2

        c = 2 * asin(sqrt(a))

        # Radius of earth in kilometers. Use 3956 for miles
        r = 6371000
        D = c * r

        self.distance_between_nodes = D
        return


    def __bearing_between_hydrophones(self):
        """
        bearing_between_hydrophones(coord1, coord2) - calculates the
            bearing in degrees (NSEW) between coord1 and coord2

        Inputs:
        coord1 - numpy array
            of shape [2,1] containing latitude and longitude of point1
        coord2 - numpy array
            of shape [2,1] containing latitude and longitude of point2

        Outputs:
        self.bearing_d_1_2 - float
            bearing in degrees between node 1 and node 2
        """
        coord1, coord2 = self.__get_node_coords()


        psi1 = np.deg2rad(coord1[0])
        lambda1 = np.deg2rad(coord1[1])
        psi2 = np.deg2rad(coord2[0])
        lambda2 = np.deg2rad(coord2[1])
        del_lambda = lambda2 - lambda1

        y = np.sin(del_lambda)*np.cos(psi2)
        x = np.cos(psi1)*np.sin(psi2) - np.sin(psi1)*np.cos(psi2)*np.cos(del_lambda)

        theta_bearing_rad = np.arctan2(y, x)
        theta_bearing_d_1_2 = (np.rad2deg(theta_bearing_rad) + 360) % 360

        self.bearing_node_1_2 = theta_bearing_d_1_2


    def __get_node_coords(self):
        #get coordinates from location codes
        hydrophone_locations = {
            '/LJ01D': [44.63714, -124.30598],
            '/LJ01C': [44.36943, -124.95357],
            '/PC01A': [44.52897, -125.38967],
            '/LJ01A': [44.51512, -125.38992],
            '/LJ03A': [45.81668, -129.75435],
            '/PC03A': [45.83049, -129.75327],
            'Slope_Base': [44.509781, -125.405296],
            'Southern_Hydrate': [44.569118,-125.147903],
            'Axial_Base': [45.820179,-129.736694],
            'Central_Caldera': [45.954681,-130.008896],
            'Eastern_Caldera': [45.939671,-129.973801]}

        coord1 = hydrophone_locations[self.node1]
        coord2 = hydrophone_locations[self.node2]

        return coord1, coord2


    def get_bearing_angle(self, delay_T):
        # speed of sound manually set from SSP
        c = 1485 #m/s

        #bearing is with respect to node1 (where node2 is at 0 deg)
        bearing_local = np.array([np.rad2deg(np.arccos(c*delay_T/self.distance_between_nodes)),
            -np.rad2deg(np.arccos(c*delay_T/self.distance_between_nodes))])
        
        #convert to global (NSEW) degrees
        bearing_global = self.bearing_node_1_2 + bearing_local
        
        #make result between 0 and 360
        bearing_global = bearing_global % 360

        return bearing_global


    def plot_map_bearing(self, bearing_angle):
        
        # get coordinates of hydrophone locations
        coord1, coord2 = self.__get_node_coords()

        thetaB1 = bearing_angle[0]
        thetaB2 = bearing_angle[1]

        midpoint, phantom_point1 = self.__find_phantom_point(coord1, coord2, thetaB1)
        midpoint, phantom_point2 = self.__find_phantom_point(coord1, coord2, thetaB2)

        import plotly.graph_objects as go

        hyd_lats = [coord1[0], coord2[0]]
        hyd_lons = [coord1[1], coord2[1]]

        antmidpoint = self.__get_antipode(midpoint)
        fig = go.Figure()

        fig.add_trace(go.Scattergeo(
            lat=[midpoint[0], phantom_point1[0], antmidpoint[0]],
            lon=[midpoint[1], phantom_point1[1], antmidpoint[1]],
            mode='lines',
            line=dict(width=1, color='blue')
        ))

        fig.add_trace(go.Scattergeo(
            lat=[midpoint[0], phantom_point2[0], antmidpoint[0]],
            lon=[midpoint[1], phantom_point2[1], antmidpoint[1]],
            mode='lines',
            line=dict(width=1, color='green')
        ))

        fig.add_trace(go.Scattergeo(
            lon=hyd_lons,
            lat=hyd_lats,
            hoverinfo='text',
            text=[self.node1,
                  self.node2],
            mode='markers',
            marker=dict(
                size=5,
                color='rgb(255, 0, 0)',
                line=dict(
                    width=3,
                    color='rgba(68, 68, 68, 0)'
                )
            )))

        fig.update_layout(
            title_text='Possible Bearings of Max Correlation Peak',
            showlegend=False,
            geo=dict(
                resolution=50,
                showland=True,
                showlakes=True,
                landcolor='rgb(204, 204, 204)',
                countrycolor='rgb(204, 204, 204)',
                lakecolor='rgb(255, 255, 255)',
                projection_type="natural earth",
                coastlinewidth=1,
                lataxis=dict(
                    # range=[20, 60],
                    showgrid=True,
                    dtick=10
                ),
                lonaxis=dict(
                    # range=[-100, 20],
                    showgrid=True,
                    dtick=20
                ),
            )
        )

        fig.show()
        fig.write_html('Map_Bearing.html')

        return


    def __find_phantom_point(self, coord1, coord2, thetaB):
        """
        find_phantom_point

        Inputs:
        coord1 - list
            coordinate of first hydrophone
        coord2 - list
            coordinate of second hydrophone

        Output:
        midpoint, phantom_point
        """
        midpoint = [coord1[0] - (coord1[0] - coord2[0]) / 2,
                    coord1[1] - (coord1[1] - coord2[1]) / 2]

        del_lat = 0.01 * np.cos(np.deg2rad(thetaB))
        del_lon = 0.01 * np.sin(np.deg2rad(thetaB))

        phantom_point = [midpoint[0] + del_lat, midpoint[1] + del_lon]

        return midpoint, phantom_point


    def __get_antipode(self, coord):
        # get antipodes
        antlon = coord[1] + 180
        if antlon > 360:
            antlon = antlon - 360
        antlat = -coord[0]
        antipode_coord = [antlat, antlon]
        return antipode_coord


    def MA_TDGF(self, avg_time, stride, start_hour = 0, end_hour = None,
                verbose=True):
        '''
        MA_TDGF - calculates the moving average of a NCCF experiment. Returns
            a 2D array with dimensions defined as [date, Tau]

        Parameters:
        -----------
        avg_time : int
            number of hours each TDGF is estimated from (must be odd)
        stride : int
            stride of moving average in hours
        start_hour : int
            hour of start for moving average in experiment. Default is zero
        end_hour : int
            hour of end for moving average in experiment. Default is last hour
            entry

        Returns
        -------

        NCCF_array : NCCF_array
            Data type for storing NCCFs array data
        '''
        if end_hour is None:
            end_hour = self.num_periods

        exp_length = self.num_periods

        # Create dates variable
        dates = []
        for k in range(end_hour):
            hours_to_add = datetime.timedelta(hours=k)
            dates.append(self.start_time + hours_to_add)

        # Check if avg time is valid
        if avg_time > (end_hour-start_hour):
            raise Exception ('average time is longer than experiment length')
        if avg_time % 2 is not 1:
            raise Exception ('average time must be odd')

        m = avg_time
        n = end_hour

        bar = progressbar.ProgressBar(maxval=int(((end_hour-start_hour)-m)/stride), \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        if verbose: bar.start()

        tdgfs = []
        num_available = []
        dates_new = []
        count = 0

        # k is now middle average
        for k in range(int(start_hour + (m-1)/2), int(end_hour-(m-1)/2) ,stride):
            ma_start = int(k - (m-1)/2)
            ma_mid = int(k)
            ma_end = int(k+(m-1)/2)
            

            # print(ma_start, ma_mid, ma_end)
            dates_new.append(dates[ma_mid])
            try:
                nccf, x = self.average_NCF(ma_start, ma_end)
                num_available.append(x)
                tdgfs.append(nccf)
            except UnboundLocalError:
                print('test')

            if verbose: bar.update(count)
            count = count + 1
        
        NCCFs = np.array(tdgfs)
        dates = np.array(dates_new)
        if self.htype == 'broadband':
            Fs = 64000
        elif self.htype == 'low_frequency':
            Fs = 200
        else:
            raise Exception ('Invalid Hydrophone Type')
        Ts = 1/Fs

        NCCFs_object = NCCF_array(
            NCCFs, stride, avg_time, num_available, dates=dates)
        return NCCFs_object


class NCCF_array:
    '''
    data type for storing Noise Cross Correlation Function (NCCF) Experiments

    Attributes
    ----------
    NCCFs : numpy array
        2D array containing noise cross correlation functions. Array has shape
        [M,N] where M is the number of noise cross correlation functions and
        N is the number of time samples in each NCCF.
    dates : numpy array of datetime.datetime objects
        1D array containing the start times of each NCCF.
    t : numpy array
        1D array containing the time values for each time bin of a NCCF.
    stride : int
        stride of the moving average used to create each NCCF (in hours)
    avg_len : int
        length of moving average used to create each NCCF (in hours)
    peaks : dictionary
        keys are string to specify peak name. d(A or B) or sxby(A or B).
        where x is number of surface reflections, and y is the number of
        bottom reflections. A refers to lag peaks, B refers to lead peaks.
        dictionary returns numpy array of specified peak for all dates in 
        NCCFs
    peak_slices : dictionary
        keys are peak_names, returns numpy slice for tau dimension of NCCFs
    peak_window_idx : dictionary
        manual index of peak bounds
    NCCFs_c : numpy array
        NCCFs, but delay time dimension has been converted to complex signal
        using Hilbert Transform

    Methods
    -------
    '''

    def __init__(self, NCCFs, stride, avg_len, num_available, dates=None, t=None):
        self.NCCFs = NCCFs
        self.dates = dates
        self.stride = stride
        self.avg_len = avg_len
        self.num_available = num_available
        # default t (None) is Fs = 200, W = 30s
        if t == None:
            W = 30
            Fs = 200
            self.t = np.linspace(-W,W,Fs*2*W-1)
        else:
            self.t = t
        '''
        peak_windows = [
            [5535, 5635],
            [5303, 5503],
            [4969, 5169],
            [4622, 4782],
            [6367, 6467],
            [6498, 6698],
            [6833, 7033],
            [7242, 7402],
        ]

        slices = []
        for k in range(len(peak_windows)):
            slices.append(np.s_[peak_windows[k][0]:peak_windows[k][1]])

        peak_names = ['dA', 's1b0A', 's2b1A', 's3b2A', 'dB', 's1b0B', 's2b1B', 's3b2B']
        '''
        self.peak_slices = peak_slices
        
        self.peaks = {
            'dA': self.NCCFs[:,slices[0]],
            's1b0A': self.NCCFs[:,slices[1]],
            's2b1A': self.NCCFs[:,slices[2]],
            's3b2A': self.NCCFs[:,slices[3]],
            'dB': self.NCCFs[:,slices[4]],
            's1b0B': self.NCCFs[:,slices[5]],
            's2b1B': self.NCCFs[:,slices[6]],
            's3b2B': self.NCCFs[:,slices[7]]
        }

        # add NCCFs_c - complex version of NCCFs
        self.NCCFs_c = scipy.signal.hilbert(self.NCCFs)


    def phase_of_peak(self, peak_id):
        '''
        phase_of_peak - finds phase of peak identified by peak_id and returns
            the [M,] array of phase for each date

        Parameters
        ----------
        peak_id : str
            string to specify peak name. d(A or B) or sxby(A or B). where x is number
            of surface reflections, and y is the number of bottom reflections.
            A refers to lag peaks, B refers to lead peaks
        
        Returns
        -------
        phase_array : numpy array
            array of shape [M,] containing the phase (in degrees) of specified
            peak for each date m.
        '''

        tdgf_crop = self.peaks[peak_id]

        tdgf_c = self.NCCFs_c[:,self.peak_slices[peak_id]]
        tdgf_f = scipy.fft.fft(tdgf_c, axis=1)
        data_mag = np.abs(tdgf_f)
        data_ang = np.angle(tdgf_f)

        data_90 = np.real(scipy.fft.ifft(data_mag*np.exp(1j*(data_ang + np.pi/2)), axis=1))

        a1 = np.sum(tdgf_crop,axis=1)
        a2 = np.sum(data_90, axis=1)

        phase2 = np.rad2deg(np.arctan2(-a2,a1)+np.pi)
        phase2 = phase2 % 360

        return phase2


    def snr_of_peak_amp(self, peak_id):
        '''
        snr_of_peak_amp - find SNR of specified peak using amplitude method
            (peak amp / std of noise)

        Parameters
        ----------
        peak_id : str
            string to specify peak name. d(A or B) or sxby(A or B). where x is number
            of surface reflections, and y is the number of bottom reflections.
            A refers to lag peaks, B refers to lead peaks

        Returns
        -------
        snr_plot : numpy array
            array of shape [m,] containing amplitude SNR (in dB) for specified
            peak for each date instance of averaged NCCF
        '''

        noise_bounds = np.array([-1.5, 1.5])
        NCCFs_c = self.NCCFs_c

        noise_idx = np.array([
            np.argmin(np.abs(self.t - noise_bounds[0])),
            np.argmin(np.abs(self.t - noise_bounds[1]))
        ])
        
        noise_slice = np.s_[noise_idx[0]:noise_idx[1]]
        noise_std = np.std(np.abs(self.NCCFs_c[:,noise_slice]), axis=1)

        peak_amp = np.max(np.abs(NCCFs_c[:,self.peak_slices[peak_id]]),axis=1)
        peak_idx = np.argmax(np.abs(NCCFs_c[:,self.peak_slices[peak_id]]),axis=1)

        peak_time = (peak_idx + self.peak_slices[peak_id].start)/200 - 30
        SNR = 20*np.log10(peak_amp/noise_std)

        if np.max(np.abs(np.gradient(peak_idx))) > 12:
            warnings.warn('Peak index jumps more than 3 in single step')
        return SNR, peak_time, noise_std


    def SNR_plots(self):
        '''
        Calculate 24 SNR Plots (using magnitude signal) for year of data and
        save as PNGs in file directory /SNR_plots
        '''
        
        if self.avg_len != 1:
            raise Exception ('Average length must be 1 for this method. Please recalculate NCCFs_array with avg_len = 1')

        for k in range(24):
            start_time = k*356
            end_time = (k+1)*365 - 1
            print(start_time, end_time)
            NCCF_avg = []
            for n in range(start_time, end_time):
                print(start_time, n+1)
                temp = np.mean(np.abs(self.NCCFs_c[start_time:n+1,:]),axis=0)
        
                return temp


    def compute_2D_SNR_plots(self):
        '''
        takes NCCFs_array (with stride and avg_len = 1) and computes 2D SNR
        Plots. These plots consist of 2D array [dates, avg_time] where dates
        specifies the hour of the year and avg_time specifies how many hours
        are averaged. Z value is SNR in dB
        '''
        if self.stride != 1 and self.avg_len != 1:
            raise Exception('Stride and avg_len of NCCFs_array must be 1 for this method. Please recalculate NCCFs_array')


