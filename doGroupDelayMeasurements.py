#!/usr/bin/env ipython

import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import math
from glob import glob
import datetime
from matplotlib import dates
from S21analysisFunctions import *
import CableResponses as CR

def main():
    """
    Takes in csv from TDS6804B and processes waveforms.
    Interesting deconvolution stuff is in my CableResponses class.
    """

    dataDir = 'seaveyDataPalestine2014/S21s/'

    """
    vals, dt, t0 = getWaveform(dataDir + '140625_051909_rxp30_hpol_+15az+00el_Ch1.csv')
    times = [t0 + i*dt for i in xrange(len(vals))]
    
    plt.figure()
    plt.plot(times, [v*1e2 for v in vals], label = 'A waveform')

    plt.show()
    """
    waves, dts, t0s = getAllWaveformsNoiseSubtracted(ant = 30, az = 30, el = 0, pol = 'vpol')

    #print waves
    print dts
    print t0s
                               


    return 0




if __name__ == '__main__':
    main()
