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

    maxFreqMHz = 1300 #2000
    minFreqMHz = 100
    padToLength = 8192*2

    speedOfLight = 299792458 # m/s

    seaveyNumsVPol = [6.3, 7.7, 9.5, 9.0, 12.5]
    seaveyNumsHPol = [6.0, 8.1, 10.1, 8.0, 12.8]
    seaveyFreqs = [200, 450, 700, 950, 1200]
    
    #savePlots = True
    #printAverageVpolResponseFile = True
    #printAverageHpolResponseFile = True
    savePlots = False
    printAverageVpolResponseFile = False
    printAverageHpolResponseFile = False
    doSqrt = True #False # For debugging Friis correction

    crs = CR.CableResponses(padToLength, dataDir)
    
    # How much do we need? in nano-seconds
    pulseWindow = 30 #ns 
    prePeakWindow = 1 #ns
    postPeakWindow = pulseWindow - prePeakWindow
    
    highPass = 100 # MHz
    lowPass = 1300 # MHz
    
    listOfAnts = xrange(1, 51)
    #listOfAnts = [35]

    listOfChannels = ['Ch1', 'Ch4']
    #listOfPols = ['hpol', 'vpol']
    listOfPols = ['vpol', 'hpol']
    #listOfAnts = ['rxp25']

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
        if antInd > 0:
            continue
        #print 'Doing analysis for antenna ' + str(ant)

        # Now doing on an antenna by antenna basis
        #fig, axes = plt.subplots(3)#2)
        #plt.suptitle = 'Antenna ' + str(antInd+1)
        #
        #axes[0].set_title('Aligned')
        #axes[0].grid(b=True, which='major', color='black', linestyle='--')
        #plt.ylabel('Gain (dBi)')
        #axes[1].set_title('Cross Polarization')
        #axes[1].grid(b=True, which='major', color='black', linestyle='--')
        #plt.ylabel('Relative power (dB)')
        ##plt.xlabel('Frequency (MHz)')
        #axes[2].set_title('Phase')
        #axes[2].grid(b=True, which='major', color='black', linestyle='--')
        #plt.ylabel('Phase (rads)')
        #plt.xlabel('Frequency (MHz)')

        for polInd, pol in enumerate(listOfPols):

            waves, dts, t0s = getAllWaveformsNoiseSubtracted(ant = antInd+1, pol = pol)

            relativeCrossPol = []
            indexOfAbsMax = 0
            phaseCenterSeparation = -1 # Metres

            sampleWindowCopol = {}

            for chanInd, chan in enumerate(listOfChannels):

                newV = waves[chan]
                dt = dts[chan]
                t0 = t0s[chan]

                times = [t0 + i*dt for i, v in enumerate(newV)]

                if chan == 'Ch1':
                    maxInd = findPulseMaxInd(newV)
                    sampleWindowCopol['start'] = maxInd - int(prePeakWindow/dt)
                    sampleWindowCopol['end'] = maxInd + int(postPeakWindow/dt)
                    #print sampleWindowCopol

                windowedPulse = windowPulse(newV, 
                                            startSample = sampleWindowCopol['start'], 
                                            endSample = sampleWindowCopol['end'])
                                            #startSample = 0, 
                                            #endSample = 30000)

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
                phase = []
                if chan == 'Ch1': # Channel 1 on scope was copol
                    removedCopol = crs.removeCopol(windowedPulse,
                                                       dts[chan])

                    newV = CR.doNormalizedInvFFT(removedCopol, dtNs = dt)
                    print crs.t0s
                    print crs.dts
                    print ''

                    #plt.figure()
                    ##plt.plot([crs.t0s['Co'] + crs.dts['Co']*i for i, v in enumerate(newV)], [v*1e4 for v in newV])
                    #plt.plot([t0 + dt*i for i, v in enumerate(windowedPulse)], [v*1e2 for v in windowedPulse], label = 'windowed pulse')
                    ##plt.plot([t0 + dt*i for i, v in enumerate(newV)], newV)
                    #plt.plot([crs.t0s['Co'] + crs.dts['Co']*i for i, v in enumerate(crs.waves['Co'])], crs.waves['Co'], label='Copol')
                    #plt.plot([crs.t0s['Co5'] + crs.dts['Co5']*i for i, v in enumerate(crs.waves['Co5'])], crs.waves['Co5'], label = 'Copol+5ft')

                    
                    pulseTime = CR.doNormalizedInvFFT(crs.pulseFreqs, dtNs = 0.05) # 20Gsa
                    #plt.plot([crs.t0s['Co5'] + crs.dts['Co5']*i for i, v in enumerate(pulseTime)], pulseTime, label = 'Pure pulse')
                    #plt.plot([crs.t0s['P5'] + crs.dts['P5']*i for i, v in enumerate(crs.waves['P5'])], crs.waves['P5'], label = 'Pulse + 5ft')
                    
                    #plt.plot([crs.t0s['P5'] + crs.dts['P5']*i for i, v in enumerate(crs.pulse5ft)], crs.pulse5ft)

                    #plt.plot([t0s['Ch2'] + dts['Ch2']*i for i, v in enumerate(waves['Ch2'])], waves['Ch2'])

                    phase = CR.getPhaseFromFFT(removedCopol)
                    #dt_ab = crs.getAntToAntDelayLeadingEdge(windowedPulse, dt, t0)
                    
                    #plt.legend()
                    #phaseCenterSeparation = dt_ab*speedOfLight*1e-9 
                    faceSeparation = 8.89 # apparently
                    phaseCenterToFace = 0.20
                    phaseCenterSeparation = faceSeparation + 2*phaseCenterToFace
                    print 'Separation = ' + str(phaseCenterSeparation) + ' m'
                    #print 'Phase centre distance behind face ' + str((phaseCenterSeparation - faceSeparation)/2)
                    print 'Phase centre distance behind face ' + str((phaseCenterSeparation - faceSeparation)/2)

                    antennaGain, f = crs.removeCopolCablesAndDoFriisCorrection(wave = windowedPulse,
                                                                               dtNs = dts[chan],
                                                                               distMeters = phaseCenterSeparation,
                                                                               doSqrt = doSqrt)

                    dw = (f[1] - f[0])*2*math.pi
                    phase = [-(phase[i-1]-phase[i])/dw if i > 0 else 0 for i, p in enumerate(phase)]

                else:
                    antennaGain, f = crs.removeXpolCablesAndDoFriisCorrection(wave = windowedPulse, 
                                                                              dtNs = dts[chan], 
                                                                              distMeters = phaseCenterSeparation,
                                                                              doSqrt = doSqrt)
                antennaGain_dB  = CR.dBScale(antennaGain)
                maxPlotInd = int(maxFreqMHz/(f[1]-f[0]))
                minPlotInd = int(minFreqMHz/(f[1]-f[0]))
                

                myLabel = pol.capitalize()
                if chan == 'Ch1': #Ch1 is direct
                    #axes[chanInd].plot(f[minPlotInd:maxPlotInd], antennaGain_dB[minPlotInd:maxPlotInd], label = myLabel)
                    #axes[chanInd+2].plot(f[minPlotInd:maxPlotInd], phase[minPlotInd:maxPlotInd], label = 'phase')
                    pass
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
                    
                else:
                    relativeCrossPol = [rcp+p for rcp, p in zip(relativeCrossPol, antennaGain_dB) ]


            #plt.plot(f, [math.log10(rcp) for rcp in relativeCrossPol], label = 'RelativeCrossPol')
            bandPassRCP = [rcp if fVal > highPass and fVal < lowPass else -50 for rcp, fVal in zip(relativeCrossPol, f)]

            if chan == 'Ch4':
                if polInd == 1: # text selection didn't see to work here...
                    myLabel = 'Vpol to Hpol'
                else:
                    myLabel = 'Hpol to Vpol'
                #axes[chanInd].plot(f[minPlotInd:maxPlotInd], relativeCrossPol[minPlotInd:maxPlotInd], label = myLabel)

        #for ax in axes:
        #    ax.legend(loc='lower right', fancybox=True)


        #if savePlots == True:
        #    fig.savefig('measurementSummaryDocs/rpx' + str(ant) + '.png',dpi=100)
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
        with file('meanVPolResponse.dat', 'w') as outFile:
            outFile.write('vPolSeaveyGain_dBi\tFreqsMHz\n')
            for g, f in zip(mean_vpol_gain_dB[minPlotInd:maxPlotInd], freqs[minPlotInd:maxPlotInd]):
                outFile.write(str(g) + '\t' + str(f) + '\n')

    if printAverageHpolResponseFile == True:
        with file('meanHPolResponse.dat', 'w') as outFile:
            outFile.write('hPolSeaveyGain_dBi\tFreqsMHz\n')
            for g, f in zip(mean_hpol_gain_dB[minPlotInd:maxPlotInd], freqs[minPlotInd:maxPlotInd]):
                outFile.write(str(g) + '\t' + str(f) + '\n')

                

    plt.show()





    return 0




if __name__ == '__main__':
    main()
