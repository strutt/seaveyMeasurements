#!/usr/bin/env python2

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
    Reads in the direct data, and the farther off data.
    Then using geometry from notebook, figures out phase center location
    """

    dataDir = 'seaveyDataPalestine2014/S21s/'

    listOfChannels = ['Ch1', 'Ch2']
    listOfPols = ['hpol', 'vpol']
    listOfAnts = ['rxp25']

    fartherBools = [False, True]

    chanToPol = {'Ch1':'Aligned ', 'Ch4':'Cross-pol '}

    ant = listOfAnts[0]
    antInd = 24
    for fb in fartherBools:
        print 'Doing analysis for antenna ' + str(ant)
        # Glob is the module that string searches for files
        globString = 'seaveyDataPalestine2014/S21s/*' + ant + '*.csv'
        listOfFiles = glob(globString)

        # Now doing on an antenna by antenna basis
        for polInd, pol in enumerate(listOfPols):

            for chanInd, chan in enumerate(listOfChannels):

                fs = []
                waves = []
                dts = []
                attenFlag = -1
                normFlag = -1
                for f in listOfFiles:
                    if (chan in f or chan.lower() in f) and pol in f and 'pulser' not in f and 'az' not in f and (('farther' in f) == fb):
                        v, dt, t0 = getWaveform(f)
                        if 'atten' in f or 'off' in f:
                            attenFlag = len(fs)
                            normFlag = 1 - attenFlag
                        waves.append(v)
                        dts.append(dt)
                        fs.append(f)
                if len(dts) > 2:
                    print 'globString = ' + globString
                    print 'Must contain ' + chan + ' and ' + pol
                    print fs
                    raise Exception('Too many files matching!')
                elif len(dts) < 2:
                    print 'globString = ' + globString
                    print 'Must contain ' + chan + ' and ' + pol
                    print fs
                    raise Exception('Too few matching files!')

                # Make times, from 0 increasing in steps of dt
                times = [dts[0]*i for i in range(len(waves[0]))]

                # Subtract 60dB attenutated pulse - corrects for pulser box noise in front of antenna pulse
                newV = [v1 - v2 for v1, v2 in zip(waves[normFlag], waves[attenFlag])]
                maxInd = findPulseMaxInd(newV)

                print 'maxInd = ' + str(maxInd), chan, pol
                


    plt.show()

    return 0




if __name__ == '__main__':
    main()
