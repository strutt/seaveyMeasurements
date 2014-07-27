"""
Nice deconvolution functions and a class to hold the cable response information.


Deconvolution is just division in the fourier domain.
Deconvolution functions have undescriptive variable names like v_1 or fft_2 etc.
This is because there are two ways I like to think about the operation...

Way 1:
filterResponse = deconvolveFunction(v_1 = signalAfterFilter, v_2 = signalWithoutFilter)

Way 2:
signalWithoutFilter = deconvolveFunction(v_1 = signalAfterFilter, v_2 = filterResponse)

"""

from S21analysisFunctions import *

class CableResponses:
    """
    A simple class which will hold the cable responses and has some nicely named functions 
    """

    def __init__(self, padLength, baseDir, maxDeconvFreqMHz = 0):
        self.maxDeconvFreqMHz = maxDeconvFreqMHz
        self.dts = {}


        # From these waveforms I can get the cable response...
        # Pulser also goes through 5ft of extra cable
        # but that was also attached to the cables for this measurement
        self.pulse5ft, dt, t0 = getWaveform(baseDir + '140626_134852_ps_pulser_direct_fast_Ch1.csv', padToLength = padLength)
        self.dts['P5'] = dt
        print self.dts['P5'] 
        self.pulseThruCopol5ft, dt, t1 = getWaveform(baseDir + '140626_135504_ps_pulser_copol_fast_5ft_Ch1.csv', padToLength = padLength)
        self.dts['Co5'] = dt
        self.pulseThruXpol5ft, dt, t2 = getWaveform(baseDir + '140626_135703_ps_pulser_xpol_fast_5ft_Ch1.csv', padToLength = padLength)
        self.dts['X5'] = dt


        # Window around pulses to make the frequency domain lovely and smooth.
        # These numbers were picked by looking at the waveforms and should be fine.
        self.pulse5ft = windowPulseAroundPeak(self.pulse5ft,20, 134)
        self.pulseThruCopol5ft = windowPulseAroundPeak(self.pulseThruCopol5ft, 20, 134)
        self.pulseThruXpol5ft = windowPulse(self.pulseThruXpol5ft, 20, 134)
        
        print self.dts['X5']

        # Do the first round of deconvolution
        self.cableFreqResponseCopol = deconvolveTimeToFreq(self.pulseThruCopol5ft, self.pulse5ft)
        self.copolFreqs = makeFreqsMHz(dtNs = self.dts['Co5'], N = len(self.cableFreqResponseCopol)

        self.cableFreqResponseXpol = deconvolveTimeToFreq(self.pulseThruXpol5ft, self.dts['X5'])
        self.xpolFreqs = makeFreqsMHz(dtNs = self.dts['X5'], N = len(self.pulse5ft))


        # After getting the freq response of the cables, 
        # waveforms which don't have the additional 5ft of cable.
        self.pulseThruCopol, dt, t4 = getWaveform(baseDir + '140626_140317_ps_pulser_copol_fast_Ch1.csv', padToLength = padLength)
        self.dts['P'] = dt
        self.pulseThruCopol = windowPulseAroundPeak(self.pulseThruCopol, 20, 134)
        self.pulseThruCopolFFT = np.fft.fft(self.pulseThruCopol)

        print self.dts
        #assert len(set(self.dts)) == 1

    def removeCopol(self, wave, dtNs):
        # In order to have the same df, we need to have 
        # N_1*dt_1 == N_2*dt_2. If that's not the case, we 
        # need to take action.
        numZerosToPad = (len(self.pulseThruCopol5ft)*self.dts['Co5'])/dtNs - len(wave)
        assert numZerosToPad == int(numZerosToPad)
        for zeroInd in xrange(int(numZerosToPad)):
            wave.append(0)

        fft_wave = np.fft.fft(wave)
        withoutCables = deconvolveFreqToFreq(fft_wave, self.cableFreqResponseCopol)
        justSeaveys = deconvolveFreqToFreq(self.pulseThruCopolFFT, withoutCables)
        
        #return withoutCables 
        return justSeaveys
        



def deconvolveTimeToFreq(v_1, v_2):
    """
    Do FFT then call main deconvolver
    """
    fft_1 = np.fft.fft(v_1)
    fft_2 = np.fft.fft(v_2)
    return deconvolveFreqToFreq(fft_1, fft_2)

        
def deconvolveFreqToFreq(fft_1, fft_2):
    """
    Actually does the deconvolution.
    All other deconvolution functions end up here.
    """                
    # Deconvolution is equivalent to division in the fourier domain.
    fft_d = [fval1/fval2 for fval1, fval2 in zip(fft_1, fft_2)]
    return fft_d


def deconvolveFreqToTime(fft_1, fft_2):
    """
    Calls main deconvolver and inverse fourier transforms result
    """
    fft_d = deconvolveFreqToFreq(fft_1, fft_2)
    return self.doInvFFT(fft_d)

def deconvolveTimeToTime(v_1, v_2):
    """
    Do FFT then pass on to the next guy...
    """
    fft_1 = np.fft.fft(v_1)
    fft_2 = np.fft.fft(v_2)
    return deconvolveFreqToTime(fft_1, fft_2)


def doInvFFT(fft):
    v_deco = np.fft.ifft(fft)
    v_d = [x.real for x in v_deco]
    v_d_im = [x.imag for x in v_deco]
    if abs(sum(v_d_im)) > 0.000000001:
        print 'Something dodgy about the fft, imag values remaining after inverse!'
        raise Exception('Inverse FFT not entirely real')
                
    return v_d    
    

    




