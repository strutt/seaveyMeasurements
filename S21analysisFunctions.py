"""
Very general functions with I will probably use in lots of places in this analysis.
For example, get waveforms from CSV files and do FFT-like things
"""


import numpy as np
from numpy import fft
from matplotlib import pyplot as plt
import math
from glob import glob

def getWaveform(fileName, padToLength = 0):
    """
    Returns the volts and times from the waveform in the file named 'fileName'.
    """

    counter = 0
    vals = []
    times = []
    
    # Read file, skips first 5 lines, which contain header information
    for line in file(fileName):
        words = line.split(',')
        vals.append(1e3*float(words[4])) # 1e3 for V -> mV
        times.append(1e9*float(words[3]))
        counter += 1

    dt = times[1]-times[0]
    trueLen = len(times)

    meanVal = sum(vals)/len(vals)
    vals = [v - meanVal for v in vals]

    while len(times) < padToLength:
        times.append(times[-1] + dt)
        vals.append(0)


    return vals, dt, times[0]

def getAllWaveformsNoiseSubtracted(baseDir = 'seaveyDataPalestine2014/S21s/', ant = 0, pol = 'vpol', az = 0, el = 0):

    allWaves = {}
    allDts = {}
    allT0s = {}
    listOfChannels = ['Ch1', 'Ch2', 'Ch3', 'Ch4']
    antStr = 'rxp' + str(ant)
    if ant < 10:
        antStr = 'rxp0' + str(ant)

    for chan in listOfChannels:
        
        # Glob is the module that string searches for files
        globString = 'seaveyDataPalestine2014/S21s/*' + antStr + '*.csv'
        listOfFiles = glob(globString)

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
            Warning('Too few matching files! Skipping this waveform')
            continue

        # Make times, from 0 increasing in steps of dt
        times = [dts[0]*i for i in range(len(waves[0]))]

        # Subtract 60dB attenutated pulse - corrects for pulser box noise in front of antenna pulse
        newV = [v1 - v2 for v1, v2 in zip(waves[normFlag], waves[attenFlag])]

        allWaves[chan] = newV
        allDts[chan] = dt
        allT0s[chan] = t0

    return allWaves, allDts, allT0s


def getPowerSpectrumInfo(vals, dt, maxFreqMHz = 0, normFactor = 1, shortenOutput = True):
    """
    Return the frequencies, power (dB) and phase (Deg) from the volts and times passed to the function.
    Limits the maximum frequency if 'maxFreqMHz' > 0 is passed.
    """

    # Get frequencies from sample times
    freqs = makeFreqsMHz(dtNs = dt, N = len(vals))
    freqsBand = []

    # Limit frequency information returned by maximum freq optional arg
    if maxFreqMHz > 0:
        freqsBand = [freq for freq in freqs if freq < maxFreqMHz]
    else:
        freqsBand = freqs

    # Do FFT and convert to power spectum and phase
    theFFT = np.fft.fft(vals)
    powSpec = abs(theFFT)**2
    powSpec /= normFactor
    N = len(vals)
    
    if shortenOutput == True:
        powSpec = 2*powSpec[:N/2+1]
        powSpec[0]/=2
        powSpec[-1]/=2

    phaseInfo = np.angle(theFFT)
    phaseInfoDeg = [samp*360/math.pi for samp in phaseInfo]

    # Return a list of lists
    return freqsBand, powSpec[:len(freqsBand)], phaseInfoDeg[:len(freqsBand)]


def windowPulse(vals, startSample = 0, endSample = -1):
    
    vals = removeOffsetFromZero(vals)

    if endSample == -1:
        endSample = len(vals)
    vals = [v1 if i > startSample and i <= endSample else 0 for i, v1 in enumerate(vals)]
    return vals

def windowPulseAroundPeak(vals, numSampsBefore = -1, numSampsAfter = -1):
    maxInd = findPulseMaxInd(vals)
    vals = windowPulse(vals, maxInd - numSampsBefore, maxInd + numSampsAfter)
    return vals

def findPulseMaxInd(vals):
    absVals = [abs(v) for v in vals]
    maxInd = absVals.index(max(absVals))
    return maxInd

def removeOffsetFromZero(vals, numPointsAtBeginning = 100):
    offsetApprox = sum(vals[:numPointsAtBeginning])/numPointsAtBeginning
    vals = [v - offsetApprox for v in vals]
    return vals

def removeAttenuationTimeDomain(vals, atten_dB = 0):

    # e.g. 20dB in power is 10dB in voltage
    atten_v = atten_dB/2.
    mult = pow(10, atten_v/10.)
    vals = [v*mult for v in vals]
    return vals


def convolve(v1, dt1, v2, dt2):
    """
    Generic convolution function.
    Takes two waveforms in the time domain,
    checks they have the same length and dts,
    and if so, does FFTs and multiplies them,
    then inverse FFTs the result.
    (Multiplication if the Fourier domain is
    equivalent to convolution).
    """

    N1 = len(v1)
    N2 = len(v2)

    if not N1 == N2:
        print 'You must input waveforms with equal lengths'
        raise Exception('Different length waveforms')

    if not dt1 == dt2:
        print 'The dts of the waveforms must be equal'
        raise Exception('Different deltaTs')

    fft1 = np.fft.fft(v1)
    fft2 = np.fft.fft(v2)

    fft3 = [z1*z2/N1 for z1, z2 in zip(fft1, fft2)]

    v3_complex = np.fft.ifft(fft3)
    v3_imag = [float(z3.imag) for z3 in v3_complex]
    for v3i in v3_imag:
        if abs(v3i) > 0.000001:
            print "Looks like deconvolution didn't return real values"
            raise Exception('Check this convolution function!!!')
    
    v3_real = [z3.real for z3 in v3_complex]
    
    return v3_real


def makeFreqsMHz(dtNs = 0, N = 0):
    if dtNs == 0 or N == 0:
        raise Exception('You need to provide delta t in nanoseconds and a length to make the frequencies!')
    else:
        dfMHz = 1e3/(N*dtNs) # 1e3 since ns -> GHz and we want MHz
        return [dfMHz*freqInd for freqInd in xrange(N)]



if __name__ == '__main__':
    print 'This is a set of functions used by various scripts.'
    print 'You probably want doS21analysis.py.'



