import numpy as np
import os
import pickle
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.animation as animation
import matplotlib
import progressbar
import scipy

class NCF_analysis:

    '''
    class with tools for analyzing NCFs. NCF experiments are stored in ./NCFs.

    Attributes
    ----------
    ckpt_dir : str
        name of specific NCF directory. NCF .pkl files will be located in ckpt_dir
        where ckpt_dir is an absolute path

    Methods
    -------
    read_header()
        reads header from ckpt_dir and prints pertanent information

    NCF_plot(num_hours)
        create average NCF for num_hours (starting from beginning) and plot
        the results

    average_NCF(num_hours)
        averages NCF for num_hours (starting from beginning)

    animate(num_hours, direction, animation_len)

    available_hours(max_number)
        loops through all ckpts in ckpt_dir and reports number of valid hours.
        it also creates a plot of availability for the time frame.
    '''
    def __init__(self, ckpt_dir):
        self.ckpt_dir = ckpt_dir
        return


    def read_header(self):
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

        print(header.start_time)
        print(header.filter_cutoffs)
        print(header.htype)
        print(header.num_periods)
        print(header.node1)
        print(header.node2)
        print(header.W)
        print(header.sp_method)
        try:
            print(header.other_notes)
            return header
        except:
            return header


    def available_hours(self):
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

        fig = plt.figure(figsize=(10,0.5))
        ax = plt.imshow(NCF_available.T, aspect='auto')

        ax.axes.get_yaxis().set_visible(False)
        plt.xlabel('Hours')
        self.num_available = num_available

        print(f'Number of Available hours: {num_available}')


    def average_NCF(self, hour_start, hour_end, plot=False):
        invalid = 0
        for k in range(hour_start,hour_end):
            ckpt_name = self.ckpt_dir + '/ckpt_'+ str(k) +'.pkl'
            try:
                with open(ckpt_name, 'rb') as f:
                    xcorr_1hr = pickle.load(f)
            except:
                invalid = invalid + 1
                continue

            if np.isnan(np.sum(xcorr_1hr)):
                invalid = invalid + 1

            else:
                try:
                    xcorr = xcorr + xcorr_1hr
                except NameError:
                    xcorr = xcorr_1hr

        xcorr_avg = xcorr / self.num_available

        self.xcorr = xcorr_avg

        num_available_short = hour_end - hour_start - invalid
        if plot:
            self.NCF_plot(xcorr_avg, num_available_short)
        return xcorr_avg


    def NCF_plot(self, xcorr, num_valid, save_fig=False, file_name=None, frequency=False, xlimits=None, symetric=False):
        if 'self.xcorr' not in locals():
            Exception: 'Must Calculate average NCF first. use average_NCF(num_hours) method'
        dt = 1/self.Fs
        t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)
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
        else:
            plt.plot(t, xcorr)
            plt.title(f'{num_valid} Hours - {node1}/{node2} - Average NCF', fontsize=18, y=1.08)

            plt.xlabel('Delay τ (s)', fontsize=14)
            plt.ylabel('Average NCF Amplitude', fontsize=14)
            if xlimits is not None:
                plt.xlim(xlimits)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        if save_fig:
            f1.savefig(file_name + '.png', dpi=500)


    def animate_NCF(self, num_hours, length, file_name, xlim=None, frequency=False, interval=1, dpi=None):
        interval=int(interval)

        # set y limit
        xcorr = self.average_NCF(1,num_hours)
        ylim = np.max(np.abs(xcorr))

        # matplotlib.rcParams[figure.max_open_warning'] = 0
        if xlim == None:
            xlim = [-self.W,self.W]
        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line,time_text,

        # Hard Coded to Animate in Reverse
        def animate(i):
            i = interval*i
            dt = 1/self.Fs
            t = np.arange(-xcorr.shape[0]*dt/2,xcorr.shape[0]*dt/2,dt)
            f = np.linspace(0,self.Fs,xcorr.shape[0])
            try:
                y = self.average_NCF(1,i)
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


    def edit_header(self):
        header = self.read_header()
        header.num_periods = 1000

        header_name = os.path.dirname(os.getcwd()) + \
            'self.ckpt_dir' + '/0HEADER.pkl'

        with open(header_name, 'wb') as f:
            pickle.dump(header, f)
