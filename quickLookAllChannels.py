#!/usr/bin/env python2
"""
For sanity checking, I can look at all channels
"""

from S21analysisFunctions import *
from matplotlib import pyplot as plt

def main():
    
    #fileName = 'seaveyDataPalestine2014/S21s/140626_135748_ps_pulser_xpol_5ft_Ch1.csv'
    fileName = 'seaveyDataPalestine2014/S21s/140626_134852_ps_pulser_direct_fast_Ch1.csv'

    numChannels = 4
    channelNames = ['Ch' + str(chan) for chan in xrange(1, numChannels+1)]
    
    for chanName in channelNames:
        fileNameChan = fileName.split('_')[-1].split('.')[0]
        fileName2 = fileName[:] # Slice whole string to get copy, rather than reference.
        fileName2 = fileName2.replace(fileNameChan, chanName)
        print fileName2 # Check string is what I think it should be.
        wave, dt, t0 = getWaveform(fileName2)
        times = [t0 + i*dt for i, w in enumerate(wave)]
        plt.plot(times, wave, label = chanName)

    plt.legend()
    plt.show()




if __name__ == '__main__':
    main()
