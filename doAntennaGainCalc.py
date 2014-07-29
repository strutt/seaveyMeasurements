#!/usr/bin/env python2
import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import math
from glob import glob






def main():
    """
    First pass at Antenna gain calculation...
    no deconvolution done...
    """

    maxFreqMHz = 2000
    padToLength = 8192*2
    pulseWindow = 50 #ns
    highPass = 200
    lowPass = 1200
    

    listOfAnts = ['rxp08', 'rxp11']
    listOfChannels = ['Ch1']

    #    listOfFiles = glob('S21Palestine/*_ps_pulser_xpol*Ch1.csv')
    #    listOfFiles.append(glob('S21Palestine/*_rxp11_hpol_Ch1.csv')[0])

    dataDir = 'seaveyDataPalestine2014/S21s/'
    listOfFiles = []
    listOfFiles.append(glob(dataDir + '*_ps_pulser_xpol*Ch1.csv')[0])
    listOfFiles.append(glob(dataDir + '*_rxp11_hpol_Ch1.csv')[0])

    for f in listOfFiles:
        print f
    
    waves = []
    dts = []
    maxIndices = []
    windowedPulses = []
    for f in listOfFiles:        
        w, dt = getWaveform(f, padToLength)

        # for 20dB attenuator on pulser connection
        if '_pulser_' in f:
            # FACTOR OF 10 HERE FOR 20dB ATTENUATOR ON PULSER ONLY MEASUREMENT
            w = [w1* 10 for w1 in w] 
        waves.append(w)
        dts.append(dt)

        # Take the absolute of the subtracted voltage, and find the time of the absolute maximum
        absNewV = [abs(v) for v in waves[-1]]
        maxIndices.append(absNewV.index(max(absNewV)))

        # Find pulse window in terms of sample number
        ds = pulseWindow/(2*dt) # /2 since plus minus peak

        # Zero everything not in pulse window around around the pulse peak
        windowedPulses.append([v1 if abs(i-maxIndices[-1]) < ds else 0 for i, v1 in enumerate(waves[-1])])

    if len(set(dts)) > 1:
        print 'Differing dts!!! - doing some extra zero padding...'
        lowerDtInd = dts.index(min(dts))
        higherDtInd = dts.index(max(dts))
        print int(dts.index(min(dts)))
        newLen = int(len(waves[higherDtInd])*dts[higherDtInd]/dts[lowerDtInd])
        print len(waves[lowerDtInd]), newLen
        while len(waves[lowerDtInd]) < newLen:
            waves[lowerDtInd].append(0)
            windowedPulses[lowerDtInd].append(0)
    elif len(set(dts)) > 2:
        raise Exception('More than 2 dts, you need to handle this condition')

    fig1, axes = plt.subplots(len(waves))

    for waveInd, wave in enumerate(waves):

        t0 = [dt*i for i in range(len(wave))]

        print str(waveInd) + ' plotting...'
        axes[waveInd].plot(t0, wave, label = listOfFiles[waveInd])
        axes[waveInd].plot(t0, windowedPulses[waveInd], label = listOfFiles[waveInd] + ' windowed')
        plt.xlabel('Time (ns)')
        plt.ylabel('Volts (V)')
        axes[waveInd].legend()

    fig2 = plt.figure()
    pwrs = []
    freqs = []
    for pulseInd, pulse in enumerate(windowedPulses):
        #print pulse
        f, pwr, phase = doAnalysis(pulse, dts[pulseInd], maxFreqMHz)
        print str(pulseInd) + ' df = ' + str(f[1] - f[0])
        freqs.append(f)
        pwrs.append(pwr)
        pwrSpecDensity = [10*math.log10(p*dt*dt) for p in pwr]

        plt.plot(f, pwrSpecDensity, label = 'Pow Spec')
    plt.xlabel('Frequency MHz')
    plt.ylabel('Power Spectral Density (dB)')
    plt.legend()

    print len(pwrs[0]), len(pwrs[1])

    relativePower = pwrs[1]/(pwrs[0])

    fig3 = plt.figure()
    plt.plot(freqs[0], relativePower)
    plt.xlabel('Frequency MHz')

    c = 3e8*1e-6 #m/s, *1e-6 so f goes from Hz -> MHz
    r = 8.89 #m 33.4ft correct this...
    gain = [rp*pow((4*math.pi*r*fVal)/c,2) if fVal > 0 else 0 for rp, fVal in zip(relativePower, freqs[0])]
    gain2 = [ math.sqrt(g) if fVal > 50 and fVal < 2000 else 0 for g, fVal in zip(gain, freqs[0])]
    gain2_dB = [10*math.log10(g) if g > 0 else -10 for g in gain2]

    fig4 = plt.figure()
    plt.plot(freqs[0], gain2_dB)
    plt.xlabel('Frequency MHz')
    plt.ylabel('Gain (dBi)')
    plt.grid(b=True, which='major', color='r', linestyle='--')
    
    
    plt.show()

    return 0





def getLabel(csvFile):
    """
    Extract information from file name to go on graph legends.
    """
    labelWords = csvFile.split('_')
    antNum = '?'
    antPol = '?'
    for labelWord in labelWords:
        if 'pol' in labelWord:
            antPol = labelWord[0].split('.')[0].capitalize()
        elif 'rx' in labelWord:
            antNum = labelWord[2:].split('.')[0].capitalize()
    theLabel = antNum + antPol
    return theLabel




def getWaveform(fileName, padToLength = 0):
    """
    Returns the volts and times from the waveform in the file named 'fileName'.
    """

    counter = 0
    vals = []
    times = []
    
    # Read file, skips first 5 lines, which contain header information
    for line in file(fileName):
        if counter > 5:
            words = line.split(',')
            vals.append(1*float(words[4]))
            times.append(1e9*float(words[3]))
        counter += 1

    dt = times[1]-times[0]
    trueLen = len(times)
    while len(times) < padToLength:
        times.append(times[-1] + dt)
        vals.append(0)

    return vals, dt


def doAnalysis(vals, dt, maxFreqMHz = 0):
    """
    Return the frequencies, power (dB) and phase (Deg) from the volts and times passed to the function.
    Limits the maximum frequency if 'maxFreqMHz' > 0 is passed.
    """
    #print len(times)
    #fig1 = plt.figure()
    #plt.plot(times, vals)
    #plt.ylabel('Waveform (mV)')
    #plt.xlabel('Time (ns)')
    #plt.draw()

    #fig2 = plt.figure()

    # Get frequencies from sample times
    #print 'dt = ' + str(dt) + 'ns, df = ' + str(1e3/(len(times)*dt)) + 'MHz'
    #print len(times)
    freqs = [1e3*samp/(len(vals)*dt) for samp in range(len(vals))]
    freqsBand = []

    # Limit frequency information returned by maximum freq optional arg
    if maxFreqMHz > 0:
        freqsBand = [freq for freq in freqs if freq < maxFreqMHz]
    else:
        freqsBand = freqs

    # Do FFT and convert to power spectum and phase
    theFFT = np.fft.rfft(vals)
    powSpec = getPowerSpectralDensity(theFFT, dt) 
    #powSpec_dB = [10*math.log10(samp) for samp in powSpec]
    #plt.plot(freqsBand, powSpec_dB[:len(freqsBand)])
    #plt.ylabel('Power Spectrum (dB)')
    #plt.xlabel('Frequency (MHz)')
    #plt.draw()

    #fig3 = plt.figure()
    phaseInfo = np.angle(theFFT)
    phaseInfoDeg = [samp*360/math.pi for samp in phaseInfo]
    #plt.plot(freqsBand, phaseInfoDeg[:len(freqsBand)])
    #plt.ylabel('Phase (Degrees)')
    #plt.xlabel('Frequency (MHz)')


    #plt.show()

    # Return a list of lists
    return freqsBand, powSpec[:len(freqsBand)], phaseInfoDeg[:len(freqsBand)]



def getPowerSpectralDensity(theFFT, dt):
    """
    Normalizes the power spectra according to Ryan's scheme...
    see http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf
    """
    powSpec = np.abs(theFFT)**2
    #print powSpec
    #powSpec = powSpec*dt*dt
    #print powSpec
    #print ''
    return powSpec

if __name__ == '__main__':
    main()
