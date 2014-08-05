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
    Takes in csv from TDS6804B and processes waveforms.
    Interesting deconvolution stuff is in my CableResponses class.
    """

    dataDir = 'seaveyDataPalestine2014/S21s/'

    maxFreqMHz = 1300 #2000
    minFreqMHz = 100
    padToLength = 8192*2

    seaveyNumsVPol = [6.3, 7.7, 9.5, 9.0, 12.5]
    seaveyNumsHPol = [6.0, 8.1, 10.1, 8.0, 12.8]
    seaveyFreqs = [200, 450, 700, 950, 1200]

    savePlots = False
    printAverageVpolResponseFile = False #True #False #True
    doSqrt = False # For debugging Friis correction

    seaveySeparation = 8.89 #9.8 #8.89 #9.8 # meters

    crs = CR.CableResponses(padToLength, dataDir)

    # How much do we need?
    pulseWindow = 30 #100 #30 #100
    prePeakWindow = 1 #10 #3 #10 #ns
    postPeakWindow = pulseWindow - prePeakWindow
    
    highPass = 100
    lowPass = 1300
    
    #listOfAnts = ['rxp01', 'rxp02', 'rxp03', 'rxp04', 'rxp05', 'rxp06', 'rxp07', 'rxp08', 'rxp09', 'rxp10', 'rxp11', 'rxp12']
    #listOfAnts = ['rxp13', 'rxp14', 'rxp15', 'rxp16', 'rxp17', 'rxp18', 'rxp19', 'rxp20', 'rxp21', 'rxp22', 'rxp23', 'rxp24']
    #listOfAnts = ['rxp25', 'rxp26', 'rxp27', 'rxp28', 'rxp29', 'rxp30', 'rxp31', 'rxp32', 'rxp33', 'rxp34', 'rxp35', 'rxp36']
    #listOfAnts = ['rxp37', 'rxp38', 'rxp39', 'rxp40', 'rxp41', 'rxp42', 'rxp43', 'rxp44', 'rxp45', 'rxp46', 'rxp47', 'rxp48']
    #listOfAnts = ['rxp49', 'rxp50', 'rxp51']
    listOfAnts = ['rxp' + str(antInd+1) if antInd >= 9 else 'rxp0' + str(antInd+1) for antInd in range(51) ]
    listOfChannels = ['Ch1', 'Ch4']
    listOfPols = ['hpol', 'vpol']
    listOfAnts = ['rxp25']

    chanToPol = {'Ch1':'Aligned ', 'Ch4':'Cross-pol '}

    figA, axesA = plt.subplots(2, 2)
    plt.suptitle = 'Time domain waveforms'

    mean_vpol_gain_dB = []
    rms_vpol_gain_dB = []
    max_vpol_gain_dB = []
    min_vpol_gain_dB = []

    mean_hpol_gain_dB = []
    rms_hpol_gain_dB = []
    max_hpol_gain_dB = []
    min_hpol_gain_dB = []

    for antInd, ant in enumerate(listOfAnts):
        print 'Doing analysis for antenna ' + str(ant)
        # Glob is the module that string searches for files
        globString = 'seaveyDataPalestine2014/S21s/*' + ant + '*.csv'
        listOfFiles = glob(globString)

        # Now doing on an antenna by antenna basis
        fig, axes = plt.subplots(2)
        plt.suptitle = 'Antenna ' + str(antInd+1)

        axes[0].set_title('Aligned')
        axes[0].grid(b=True, which='major', color='black', linestyle='--')
        plt.ylabel('Gain (dBi)')
        axes[1].set_title('Cross Polarization')
        axes[1].grid(b=True, which='major', color='black', linestyle='--')
        plt.ylabel('Relative power (dB)')
        plt.xlabel('Frequency (MHz)')

        for polInd, pol in enumerate(listOfPols):

            relativeCrossPol = []
            indexOfAbsMax = 0

            sampleWindowCopol = {}

            for chanInd, chan in enumerate(listOfChannels):

                fs = []
                waves = []
                dts = []
                attenFlag = -1
                normFlag = -1
                for f in listOfFiles:
                    if (chan in f or chan.lower() in f) and pol in f and 'pulser' not in f and 'az' not in f and 'farther' not in f:
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
                if chan == 'Ch1':
                    maxInd = findPulseMaxInd(newV) 
                    sampleWindowCopol['start'] = maxInd - int(prePeakWindow/dt)
                    sampleWindowCopol['end'] = maxInd + int(postPeakWindow/dt)
                    #print sampleWindowCopol

                windowedPulse = windowPulse(newV, 
                                            startSample = sampleWindowCopol['start'], 
                                            endSample = sampleWindowCopol['end'])
                                            
                                                            

                # Limit plots...
                x0 = sampleWindowCopol['start'] - 100
                x1 = sampleWindowCopol['end'] + 100
                #axesA[chanInd, polInd].plot(times[x0:x1], newV[x0:x1], label = ant)
                axesA[chanInd, polInd].plot(times[x0:x1], windowedPulse[x0:x1], label = ant)

                axesA[chanInd, polInd].set_title(chanToPol[chan] + pol.capitalize())
                axesA[chanInd, polInd].set_xlabel('Time (ns)')
                axesA[chanInd, polInd].set_ylabel('Amplitude (mV)')
                axesA[chanInd, polInd].grid(b=True, which='major', color='black', linestyle='--')


                #windowedPulseFreqs = crs.removeCopol(windowedPulse, dt)
                antennaGain = []
                f = []
                if chan == 'Ch1': # Channel 1 on scope was copol
                    antennaGain, f = crs.removeCopolCablesAndDoFriisCorrection(wave = windowedPulse, 
                                                                               dtNs = dts[0], 
                                                                               distMeters = seaveySeparation,
                                                                               doSqrt = doSqrt)
                else:
                    antennaGain, f = crs.removeXpolCablesAndDoFriisCorrection(wave = windowedPulse, 
                                                                              dtNs = dts[0], 
                                                                              distMeters = seaveySeparation,
                                                                              doSqrt = doSqrt)

                
                #f, antennaGain, phase  = getPowerSpectrumInfo(windowedPulse, dts[0])
                antennaGain_dB  = [10*math.log10(g) if g > 0 else 0 for g in antennaGain]

                #print windowedPulse

                maxPlotInd = int(maxFreqMHz/(f[1]-f[0])) # short for maxPlotInd
                minPlotInd = int(minFreqMHz/(f[1]-f[0])) # short for maxPlotInd
                

                myLabel = pol.capitalize()
                if chan == 'Ch1': #Ch1 is direct
                    axes[chanInd].plot(f[minPlotInd:maxPlotInd], antennaGain_dB[minPlotInd:maxPlotInd], label = myLabel)
                if pol == 'vpol' and chan == 'Ch1':
                    if antInd == 0:
                        mean_vpol_gain_dB = [0 for g in antennaGain_dB]
                        rms_vpol_gain_dB = [0 for g in antennaGain_dB]
                        max_vpol_gain_dB = [-100000 for g in antennaGain_dB]
                        min_vpol_gain_dB = [1000000 for g in antennaGain_dB]
                    mean_vpol_gain_dB = [m + g for m, g in zip(mean_vpol_gain_dB, antennaGain_dB)]
                    rms_vpol_gain_dB = [r + g**2 for r, g in zip(rms_vpol_gain_dB, antennaGain_dB)]
                    max_vpol_gain_dB = [g if g > maxG else maxG for g, maxG in zip(antennaGain_dB, max_vpol_gain_dB)]
                    min_vpol_gain_dB = [g if g < minG else minG for g, minG in zip(antennaGain_dB, min_vpol_gain_dB)]
                elif pol == 'hpol' and chan == 'Ch1':
                    if antInd == 0:
                        mean_hpol_gain_dB = [0 for g in antennaGain_dB]
                        rms_hpol_gain_dB = [0 for g in antennaGain_dB]
                        max_hpol_gain_dB = [-100000 for g in antennaGain_dB]
                        min_hpol_gain_dB = [1000000 for g in antennaGain_dB]
                    mean_hpol_gain_dB = [m + g for m, g in zip(mean_hpol_gain_dB, antennaGain_dB)]
                    rms_hpol_gain_dB = [r + g**2 for r, g in zip(rms_hpol_gain_dB, antennaGain_dB)]
                    max_hpol_gain_dB = [g if g > maxG else maxG for g, maxG in zip(antennaGain_dB, max_hpol_gain_dB)]
                    min_hpol_gain_dB = [g if g < minG else minG for g, minG in zip(antennaGain_dB, min_hpol_gain_dB)]

                if chan is 'Ch1':
                    while len(relativeCrossPol) < len(antennaGain_dB):
                        relativeCrossPol.append(0)
                    relativeCrossPol = [-p for p in antennaGain_dB]
                    #print len(relativeCrossPol), f[1] - f[0], 'Ch1'
                    
                else:
                    relativeCrossPol = [rcp+p for rcp, p in zip(relativeCrossPol, antennaGain_dB) ]
                    #print len(relativeCrossPol), f[1] - f[0], 'Ch4'


            #plt.plot(f, [math.log10(rcp) for rcp in relativeCrossPol], label = 'RelativeCrossPol')
            bandPassRCP = [rcp if fVal > highPass and fVal < lowPass else -50 for rcp, fVal in zip(relativeCrossPol, f)]

            if chan == 'Ch4':
                if polInd == 1: # text selection didn't see to work here...
                    myLabel = 'Vpol to Hpol'
                else:
                    myLabel = 'Hpol to Vpol'
                axes[chanInd].plot(f[minPlotInd:maxPlotInd], relativeCrossPol[minPlotInd:maxPlotInd], label = myLabel)

        for ax in axes:
            ax.legend(loc='lower right', fancybox=True)


        if savePlots == True:
            fig.savefig('measurementSummaryDocs/' + ant + '.png',dpi=100)
    df = crs.dfMHz

    n = len(listOfAnts)

    # Finalize general calc
    mean_vpol_gain_dB = [m/n for m in mean_vpol_gain_dB]
    rms_vpol_gain_dB = [math.sqrt(r/n-m**2) for r, m in zip(rms_vpol_gain_dB, mean_vpol_gain_dB) ]
    # Finalize general calc
    mean_hpol_gain_dB = [m/n for m in mean_hpol_gain_dB]
    rms_hpol_gain_dB = [math.sqrt(r/n-m**2) for r, m in zip(rms_hpol_gain_dB, mean_hpol_gain_dB) ]

    freqs = [df*i for i in range(len(mean_vpol_gain_dB))]
    maxPlotInd = int(maxFreqMHz/(f[1]-f[0]))
    minPlotInd = int(minFreqMHz/(f[1]-f[0]))

    fig = plt.figure()
    plt.title('Vertical Polarization Antenna Gain')
    plt.plot(freqs[minPlotInd:maxPlotInd], mean_vpol_gain_dB[minPlotInd:maxPlotInd], label = 'Mean')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m+r for r, m in zip(rms_vpol_gain_dB, mean_vpol_gain_dB)][minPlotInd:maxPlotInd], label = 'Mean + RMS')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m-r for r, m in zip(rms_vpol_gain_dB, mean_vpol_gain_dB)][minPlotInd:maxPlotInd], label = 'Mean - RMS')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m for m in max_vpol_gain_dB[minPlotInd:maxPlotInd]], label = 'Bin-by-bin maximum')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m for m in min_vpol_gain_dB[minPlotInd:maxPlotInd]], label = 'Bin-by-bin minimum')
    plt.plot(seaveyFreqs, seaveyNumsVPol, 'ro', label='Seavey Measurements')
    #plt.plot(freqs[minPlotInd:maxPlotInd], [10*math.log10(m) if m > 0 else 0 for m in crs.meanVpolResponse[minPlotInd:maxPlotInd]], label = '51 antenna average')
    #plt.xticks(range(0, 2001, 200))
    #plt.yticks(range(-45, 21, 5))
    ax = plt.gca()
    ax.grid(b=True, which='major', color='black', linestyle='--')
    plt.legend(loc='lower right', fancybox=True)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Gain (dBi)')
    if savePlots == True:
        fig.savefig('measurementSummaryDocs/vpolSummary.png',dpi=100)

    fig = plt.figure()
    plt.title('Horizontal Polarization Antenna Gain')
    plt.plot(freqs[minPlotInd:maxPlotInd], mean_hpol_gain_dB[minPlotInd:maxPlotInd], label = 'Mean')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m+r for r, m in zip(rms_hpol_gain_dB, mean_hpol_gain_dB)][minPlotInd:maxPlotInd], label = 'Mean + RMS')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m-r for r, m in zip(rms_hpol_gain_dB, mean_hpol_gain_dB)][minPlotInd:maxPlotInd], label = 'Mean - RMS')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m for m in max_hpol_gain_dB[minPlotInd:maxPlotInd]], label = 'Bin-by-bin maximum')
    plt.plot(freqs[minPlotInd:maxPlotInd], [m for m in min_hpol_gain_dB[minPlotInd:maxPlotInd]], label = 'Bin-by-bin minimum')
    plt.plot(seaveyFreqs, seaveyNumsHPol, 'ro', label='Seavey Measurements')

    #plt.xticks(range(0, 2001, 200))
    #plt.yticks(range(-45, 21, 5))
    ax = plt.gca()
    ax.grid(b=True, which='major', color='black', linestyle='--')
    plt.legend(loc='lower right', fancybox=True)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Gain (dBi)')
    if savePlots == True:
        fig.savefig('measurementSummaryDocs/hpolSummary.png',dpi=100)


    if printAverageVpolResponseFile == True:
        with file('meanVpolResponse.dat', 'w') as outFile:
            outFile.write('meanVpolReponse\tFreqsMHz\n')
            for i, g in enumerate(mean_vpol_gain_dB):
                outFile.write(str(pow(10, g/10)) + '\t' + str(df*i) + '\n')


    plt.show()





    return 0




if __name__ == '__main__':
    main()
