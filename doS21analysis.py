import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import math
from glob import glob
import datetime
from matplotlib import dates
from S21analysisFunctions import *


def main():
    """
    Takes in csv from TDS6804B and processes waveforms
    
    
    Now have Frequencies received through entire signal chain
    Want only respone of seaveys.
    Need to remove cable responses from pulser.
    Deconvolve signal chain.
    Take FFT of pulse->scope
    Take FFT of pulser->signal chain
    Take

    """
    maxFreqMHz = 2000
    padToLength = 8192*2


    # How much do we need?
    pulseWindow = 100
    prePeakWindow = 10 #ns
    postPeakWindow = pulseWindow - prePeakWindow
    
    highPass = 100
    lowPass = 1300
    
    listOfAnts = ['rxp01', 'rxp02', 'rxp03', 'rxp04', 'rxp05', 'rxp06', 'rxp07', 'rxp08', 'rxp09', 'rxp10', 'rxp11', 'rxp12']
    #listOfAnts = ['rxp13', 'rxp14', 'rxp15', 'rxp16', 'rxp17', 'rxp18', 'rxp19', 'rxp20', 'rxp21', 'rxp22', 'rxp23', 'rxp24']
    #listOfAnts = ['rxp25', 'rxp26', 'rxp27', 'rxp28', 'rxp29', 'rxp30', 'rxp31', 'rxp32', 'rxp33', 'rxp34', 'rxp35', 'rxp36']
    listOfAnts = ['rxp37', 'rxp38', 'rxp39', 'rxp40', 'rxp41', 'rxp42', 'rxp43', 'rxp44', 'rxp45', 'rxp46', 'rxp47', 'rxp48']
    #listOfAnts = ['rxp49', 'rxp50', 'rxp51']
    #listOfAnts = ['rxp' + str(antInd+1) if antInd >= 9 else 'rxp0' + str(antInd+1) for antInd in range(48) ]
    listOfChannels = ['Ch1', 'Ch4']
    listOfPols = ['hpol', 'vpol']

    chanToPol = {'Ch1':'Aligned ', 'Ch4':'Cross-pol '}

    figA, axesA = plt.subplots(2, 2)
    plt.suptitle = 'Time domain waveforms'
    figB, axesB = plt.subplots(2, 2)
    plt.suptitle = 'Power Spectra'
    figD, axesD = plt.subplots(2, 2)
    plt.suptitle = 'Frequency Phases'

    figC, axesC = plt.subplots(2)
    plt.suptitle = 'Cross Polarizations'

    # Noise only measurement (rxp15) - ASSUMES NOISE ENVIRONMENT IS SAME IN V AND H!!!!
    nPSD = []
    noiseV, noiseDt, noiseT0 = getWaveform('seaveyDataPalestine2014/S21s/140624_065249_rxp15_hpol_pulseroff_Ch1.csv')
    midNoiseWave = len(noiseV)/2
    windowStart = midNoiseWave - int(prePeakWindow/noiseDt)
    windowEnd =  midNoiseWave + int(postPeakWindow/noiseDt)
    windowedNoise = [v1 if i > windowStart and i < windowEnd else 0 for i, v1 in enumerate(noiseV)]
    nf, npwr, nphase = getPowerSpectrumInfo(windowedNoise, noiseDt, maxFreqMHz)
    nPSD.append([10*math.log10(p*noiseDt*noiseDt) for p in npwr])
    noiseV2, noiseDt2, noiseT0 = getWaveform('seaveyDataPalestine2014/S21s/140624_065249_rxp15_hpol_pulseroff_Ch4.csv')
    windowedNoise2 = [v1 if i > windowStart and i < windowEnd else 0 for i, v1 in enumerate(noiseV2)]
    nf2, npwr2, nphase2 = getPowerSpectrumInfo(windowedNoise2, noiseDt2, maxFreqMHz)    
    nPSD.append([10*math.log10(p*noiseDt*noiseDt) for p in npwr2])



    for antInd, ant in enumerate(listOfAnts):
        print 'Doing analysis for antenna ' + str(ant)
        # Glob is the module that string searches for files
        globString = 'seaveyDataPalestine2014/S21s/*' + ant + '*.csv'
        listOfFiles = glob(globString)

        for polInd, pol in enumerate(listOfPols):

            relativeCrossPol = []
            indexOfAbsMax = 0

            for chanInd, chan in enumerate(listOfChannels):

                fs = []
                waves = []
                dts = []
                attenFlag = -1
                normFlag = -1
                for f in listOfFiles:
                    if (chan in f or chan.lower() in f) and pol in f and 'pulser' not in f and 'az' not in f and 'farther' not in f:
                        v, dt, t0 = getWaveform(f)
                        if 'atten' in f:
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
                windowedPulse = windowPulseAroundPeak(newV, int(prePeakWindow/dt), int(postPeakWindow/dt))
                

                # Limit plots...
                x0 = windowStart - 100
                x1 = windowEnd + 100
                #axesA[chanInd, polInd].plot(times[x0:x1], newV[x0:x1], label = ant)
                axesA[chanInd, polInd].plot(times[x0:x1], windowedPulse[x0:x1], label = ant)

                axesA[chanInd, polInd].set_title(chanToPol[chan] + pol.capitalize())
                axesA[chanInd, polInd].set_xlabel('Time (ns)')
                axesA[chanInd, polInd].set_ylabel('Amplitude (mV)')
                axesA[chanInd, polInd].grid(b=True, which='major', color='black', linestyle='--')


                #axesA[chanInd, polInd].xlim([390, 420])

                f, pwr, phase = getPowerSpectrumInfo(windowedPulse, dt, maxFreqMHz, normFactor = len(windowedPulse))
                pwrSpecDensity = [10*math.log10(p) for p in pwr]

                #if antInd == 0:
                    #axesB[chanInd, polInd].plot(f, nPSD[polInd], label = 'Noise only ' + pol, linestyle='dotted', color='black')
                axesB[chanInd, polInd].plot(f, pwrSpecDensity, label = ant)

                if chan is 'Ch1':
                    while len(relativeCrossPol) < len(pwr):
                        relativeCrossPol.append(0)
                    relativeCrossPol = [-p for p in pwrSpecDensity]
                    #print len(relativeCrossPol), f[1] - f[0], 'Ch1'
                    
                else:
                    relativeCrossPol = [rcp+p for rcp, p in zip(relativeCrossPol, pwrSpecDensity) ]
                    #print len(relativeCrossPol), f[1] - f[0], 'Ch4'


                axesD[chanInd, polInd].plot(f, phase, label = ant)
                #axesD[chanInd, polInd].legend()
                axesD[chanInd, polInd].set_title(pol.capitalize() + ' pulsed directly, ' + chanToPol[chan] + 'feed')
                axesD[chanInd, polInd].set_xlabel('Frequencies (MHz)')
                axesD[chanInd, polInd].set_ylabel('Phase (Degrees)')




            #plt.plot(f, [math.log10(rcp) for rcp in relativeCrossPol], label = 'RelativeCrossPol')
            bandPassRCP = [rcp if fVal > highPass and fVal < lowPass else -50 for rcp, fVal in zip(relativeCrossPol, f)]
            #bandPassRCP_dB = [10*math.log10(rcp) if rcp > 0 else -75 for rcp in bandPassRCP]
            axesC[polInd].plot(f, bandPassRCP, label='P' + str(antInd+1))
            #axesC[polInd].legend(loc='upper right')
            
            
            #plt.plot(f, relativeCrossPol, label = 'RelativeCrossPol')
            #plt.plot(f, bandPassRCP, label = 'RelativeCrossPol')
            #plt.legend()

    for polInd, pol in enumerate(listOfPols):
        for chanInd, chan in enumerate(listOfChannels):

            axesA[chanInd, polInd].legend(loc='lower right',ncol=3, fancybox=True)        
            if chanInd == 0:
                #The loc itslef can be a 2-tuple giving x,y of the lower-left corner of the legend
                #in axes coords (bbox_to_anchor is ignored).
                axesB[chanInd, polInd].legend(loc=[0.12,0.03], ncol=3, fancybox=True)
            if chanInd == 1:
                axesB[chanInd, polInd].legend(loc='upper center',ncol=3, fancybox=True)
            axesB[chanInd, polInd].set_title(pol.capitalize() + ' pulsed directly, ' + chanToPol[chan] + 'feed')
            axesB[chanInd, polInd].set_xlabel('Frequencies (MHz)')
            axesB[chanInd, polInd].set_ylabel('Power (dB)')
            y_lims = [ -40, 20 ]
            axesB[chanInd, polInd].set_ylim(y_lims)
            axesB[chanInd, polInd].grid(b=True, which='major', color='black', linestyle='--')

        axesC[polInd].set_title('Cross polarization: ' + pol.capitalize() + ' to ' + listOfPols[1-polInd].capitalize())
        x_lim = [0, 1500]
        axesC[polInd].set_xlim(x_lim)
        axesC[polInd].grid(b=True, which='major', color='black', linestyle='--')
        axesC[polInd].set_ylabel('Power fraction (dB)')
        axesC[polInd].set_xlabel('Frequency (MHz)')

    plt.show()

    return 0




if __name__ == '__main__':
    main()
