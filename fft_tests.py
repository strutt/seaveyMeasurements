from S21analysisFunctions import *
import random
from scipy import signal

from matplotlib import pyplot as plt
from numpy import array
from numpy import fft as fft

import CableResponses as CR

dataDir = 'seaveyDataPalestine2014/S21s/'




def test_CableResponses2():

    pl = 8192*2 # number of samples to zeropad waveforms by
    crs = CR.CableResponses(pl, 'seaveyDataPalestine2014/S21s/')

    dataDir = 'seaveyDataPalestine2014/S21s/'
    #pulserToScope, dt_ps, t0_ps = getWaveform(dataDir + '140626_140317_ps_pulser_copol_fast_Ch1.csv')
    pulserToScope, dt_ps, t0_ps = getWaveform(dataDir + '140624_090559_rxp19_vpol_Ch1.csv')

    pulserToScope = windowPulseAroundPeak(pulserToScope, 10, 300)

    dist = 8.89
    pulseDec = crs.removeCopolAndPulse(pulserToScope, dt_ps)
    #pulseDec = [abs(x)**2 for x in pulseDec]
    df = 1e3/(len(pulseDec)*dt_ps)
    freqs = [df*i for i in range(len(pulseDec))]
    pulseDec = [x*pow(4*math.pi*dist*f/3e2, 2) for x, f in zip(pulseDec, freqs)]
    pulseDec = [math.sqrt(x) for x in pulseDec]
    pulseDec_dB = [10*math.log10(x) if x > 0 else 0 for x in pulseDec]

    plt.plot(freqs, pulseDec_dB, label = 'deconvolved pulser')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Gain (dBi)')
    plt.grid(b=True, which='major', color='r', linestyle='--')
    plt.show()


def test_CableResponses():
    pl = 8192*2 # number of samples to zeropad waveforms by
    crs = CR.CableResponses(pl, 'seaveyDataPalestine2014/S21s/')
    co_dB = [20*math.log10(abs(z)) for z in crs.cableFreqResponseCopol]
    x_dB = [20*math.log10(abs(z)) for z in crs.cableFreqResponseXpol]
    N = len(x_dB)/2+1
    plt.plot(crs.copolFreqs[:N], co_dB[:N], label='Copol cable response')
    plt.plot(crs.xpolFreqs[:N], x_dB[:N], label = 'Xpol cable response')
    plt.legend()

    plt.figure()
    pure, dt_p, t0_p = getWaveform(dataDir + '140626_140317_ps_pulser_copol_fast_Ch1.csv', padToLength = pl)
    direct, dt_d, t0_d = getWaveform(dataDir + '140626_134852_ps_pulser_direct_fast_Ch1.csv', padToLength = pl)
    times1 = [t0_p + dt_p*i for i,j in enumerate(pure)]
    times2 = [t0_d + dt_d*i for i,j in enumerate(direct)]
    plt.plot(times1, pure, label = 'w/ cable')
    plt.plot(times2, direct, label = 'w/o cable')
    plt.legend()

    w, dt, t0 = getWaveform(dataDir + '140624_042916_rxp48_vpol_Ch1.csv', padToLength = pl)
    w2, dt2, t2 = getWaveform(dataDir + '140624_042726_rxp48_vpol_60dBatten_Ch1.csv', padToLength = pl)
    
    w = [i-j for i,j in zip(w, w2)]
    w = windowPulseAroundPeak(w, 10, 300)
    gain, freqs = crs.removeCablesAndDoFriisCorrection(wave = w, dtNs = dt, distMeters = 9.8)

    plt.figure()
    plt.plot(freqs[:N], [10*math.log10(g) if g > 0 else 0 for g in gain][:N], label = 'Antenna Gain')    
    ax = plt.gca()
    plt.grid(b=True, which='major', color='r', linestyle='--')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Gain (dBi)')
    plt.show()

def test_thing():

    epsilon = 1e-10
    
    l = 8000
    tw = [random.random() for i in range(l)]
    assert l == len(tw)
    
    dt = 0.1
    f, pwr, phase = getPowerSpectrumInfo(tw, dt)
    print len(tw), len(pwr)
    print sum(v**2 for v in tw), sum(pwr)
    assert abs(sum(v**2 for v in tw) - sum(pwr)) < epsilon

        
def test_deconv():
    """
    Get two waveforms I understand and deconvolve them.
    """
    print ''

    #v0, dt0, t0 = getWaveform(dataDir + '140626_134732_ps_pulser_direct_Ch1.csv')
    #v1, dt1, t1 = getWaveform(dataDir + '140626_140219_ps_pulser_copol__Ch1.csv')

    v0, dt0, t0 = getWaveform(dataDir + '140626_134852_ps_pulser_direct_fast_Ch1.csv')
    v1, dt1, t1 = getWaveform(dataDir + '140626_135504_ps_pulser_copol_fast_5ft_Ch1.csv')

    #v0 = removeAttenuationTimeDomain(v0, atten_dB = 20)
    #v1 = removeAttenuationTimeDomain(v1, atten_dB = 20)

    v0maxInd = v0.index(max(v0))
    v1maxInd = v1.index(max(v1))
    v0 = windowPulse(v0, startSample = v0maxInd-20, endSample = v0maxInd+134)
    v1 = windowPulse(v1, startSample = v1maxInd-20, endSample = v1maxInd+134)

    v0p = sum(v**2 for v in v0)
    v1p = sum(v**2 for v in v1)

    #print v0p/v1p

    # v0 is signal, v1 is signal + filter
    quotient = CR.deconvolve(v1, v0)
    #quotient = deconvolve(v1, v0)
    #quotient, remainder = signal.deconvolve(v1, v0)

    #print 'quotient = ', len(quotient), quotient
    #print 'remainder = ', len(remainder), remainder

    #q1, r1 = signal.deconvolve([float(v+1) for v in v1], [float(v+1) for v in v0])
    #print q1
    #print r1

    plt.figure()
    plt.plot(v0, label='direct')
    plt.plot(v1, label='direct + cable')
    plt.plot(quotient, label='cable response')
    plt.legend()
    plt.figure()

    f0, pwr0, phase0 = getPowerSpectrumInfo(v0, dt0, normFactor = len(v0))
    f1, pwr1, phase1 = getPowerSpectrumInfo(v1, dt1, normFactor = len(v1))
    f2, pwr2, phase2 = getPowerSpectrumInfo(quotient, dt1, shortenOutput = False)

    #print len(v1)
    #print sum(math.sqrt(p) for p in pwr2)


    #plt.plot(f0[:len(pwr0)], [10*math.log10(p) if p > 0 else 0 for p in pwr0], label='direct')
    #plt.plot(f1[:len(pwr1)], [10*math.log10(p) if p > 0 else 0 for p in pwr1], label='direct plus cable')
    plt.plot(f0[:len(pwr0)], pwr0, label='direct')
    plt.plot(f1[:len(pwr1)], pwr1, label='direct plus cable')
    plt.legend()
    plt.figure()

    #plt.plot(f2[:len(pwr2)], [10*math.log10(p) if p > 0 else 0 for p in pwr2], label='Cable response')
    #plt.plot(f2[:len(pwr2)], pwr2, label='Cable response')
    plt.plot(f2[:len(pwr2)/2+1], [10*math.log10(p) for p in pwr2[:len(pwr2)/2+1]], label='Cable response')
    #plt.plot(f2, phase2, label='Cable response')

    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Power (dB)')
    #plt.plot(f2[:len(pwr2)], [math.sqrt(p) for p in pwr2], label='some kind of deconvolution')

    plt.legend()
    plt.show()

def test_deconv_for_dummies():

    """
    I am an idiot who does not understand convolution and deconvolution.
    This test function does things in a very simple manner so I can get to understand them.
    This is how it should, in principal, work.
    However, division by zero kills the scipy function, so I will need my own...
    """
    print ''
    # Simple signal
    sig = np.array([1., 1., 1., 1., 1., 2., 2., 2., 2.,])

    # Simple filter
    filt = np.array([2., 2., 1., 1., 1., 1., 1., 1., 1.])

    # Convolving the signal with the filter...
    res = signal.convolve(sig, filt)

    # quotient and remainder after deconvolution
    q, r = signal.deconvolve(res, filt)
    
    # Still don't have much intuition, but deconvolution removes the filter...
    # i.e. filter = quotient
    print 'signal ', sig
    print 'filter ', filt
    print 'convolution result ', res
    print 'deconvolving the filter from the result: quotient ', q
    print 'deconvolving the filter from the result: remainder ', r

    # Can we get the filter from the signal and the result of the convolution?
    q1, r1 = signal.deconvolve(res, sig)
    print 'Can we get the filter from the signal and the result of the convolution?'
    print 'q1, r1 = signal.deconvolve(res, sig)'
    print 'q1 = filter ', q1
    print 'r1 = remainder ', r1
    print 'Yes! The quotient is now the value I put in for the filter!'
