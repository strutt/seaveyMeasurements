"""
Nice deconvolution functions and a class to hold the cable response information.


Deconvolution is just division in the fourier domain.
Deconvolution functions have undescriptive variable names like v_1 or fft_2 etc.
This is because there are two ways I like to think about the operation...

Way 1:
filterResponse = deconvolveFunction(v_1 = signalAfterFilter, v_2 = signalWithoutFilter)

Way 2:
signalWithoutFilter = deconvolveFunction(v_1 = signalAfterFilter, v_2 = filterResponse)



Friss equation for antenna transmission:

P1/P2 = G1*G2*(lambda/4*pi*R)**2
P1/P2 = G1*G2*(c/f*4*pi*R)**2

We were pulsing in thre VPOL so the gains are the same.
Lambda is the wavelength for which we are doing this calculation = c/f.
The antenna separation, R,  we actually want is the separation of the phase centers.

"""

from S21analysisFunctions import *

class CableResponses:
    """
    A simple class which will hold the cable responses and has some nicely named functions 
    """

    def __init__(self, padLength, baseDir, maxDeconvFreqMHz = 0):
        self.maxDeconvFreqMHz = maxDeconvFreqMHz
        self.dts = {}
        self.lens = {}
        self.t0s = {}

        # From these waveforms I can get the cable response...
        # Pulser also goes through 5ft of extra cable
        # but that was also attached to the cables for this measurement
        
        self.pulse5ft, dt, t0 = getWaveform(baseDir + '140626_134852_ps_pulser_direct_fast_Ch1.csv', padToLength = padLength)
        self.dts['P5'] = dt
        self.lens['P5'] = len(self.pulse5ft)
        self.t0s['P5'] = t0

        self.pulseThruCopol5ft, dt, t1 = getWaveform(baseDir + '140626_135504_ps_pulser_copol_fast_5ft_Ch1.csv', padToLength = padLength)
        self.dts['Co5'] = dt
        self.lens['Co5'] = len(self.pulseThruCopol5ft)
        self.pulseThruXpol5ft, dt, t2 = getWaveform(baseDir + '140626_135703_ps_pulser_xpol_fast_5ft_Ch1.csv', padToLength = padLength)
        self.dts['X5'] = dt
        self.lens['X5'] = len(self.pulseThruXpol5ft)
        self.t0s['X5'] = t1
        
        """
        #self.pulse5ft, dt, t0 = getWaveform(baseDir + '140626_134732_ps_pulser_direct_Ch1.csv', padToLength = padLength)
        #self.dts['P5'] = dt
        #self.lens['P5'] = len(self.pulse5ft)
        #self.pulseThruCopol5ft, dt, t1 = getWaveform(baseDir + '140626_135419_ps_pulser_copol_5ft_Ch1.csv', padToLength = padLength)
        #self.dts['Co5'] = dt
        #self.lens['Co5'] = len(self.pulseThruCopol5ft)
        #self.pulseThruXpol5ft, dt, t2 = getWaveform(baseDir + '140626_135748_ps_pulser_xpol_5ft_Ch1.csv', padToLength = padLength)
        #self.dts['X5'] = dt
        #self.lens['X5'] = len(self.pulseThruXpol5ft)
        """

        # Need to account for 20dB attenuator on pulser measurement!

        self.pulse5ft = removeAttenuationTimeDomain(self.pulse5ft, atten_dB = 20)
        self.pulseThruCopol5ft = removeAttenuationTimeDomain(self.pulseThruCopol5ft, atten_dB = 20)
        self.pulseThruXpol5ft = removeAttenuationTimeDomain(self.pulseThruXpol5ft, atten_dB = 20)

        # Window around pulses to make the frequency domain lovely and smooth.
        # These numbers were picked by looking at the waveforms and playing until samples contained just the pulse
        # I think they are fine.
        self.pulse5ft = windowPulseAroundPeak(self.pulse5ft,20, 134)
        self.pulseThruCopol5ft = windowPulseAroundPeak(self.pulseThruCopol5ft, 20, 134)
        self.pulseThruXpol5ft = windowPulseAroundPeak(self.pulseThruXpol5ft, 20, 134)

        # Do the first round of deconvolution, which gets just the cable responses and the
        # effects of the 10dB coupler on that line.
        self.cableFreqResponseCopol = deconvolveTimeToFreq(self.pulseThruCopol5ft, self.dts['Co5'],
                                                           self.pulse5ft, self.dts['P5'])
        self.lens['Co'] = len(self.cableFreqResponseCopol)
        self.copolFreqs = makeFreqsMHz(dtNs = self.dts['Co5'], N = self.lens['Co'])

        self.cableFreqResponseXpol = deconvolveTimeToFreq(self.pulseThruXpol5ft, self.dts['X5'],
                                                          self.pulse5ft, self.dts['P5'])
        self.xpolFreqs = makeFreqsMHz(dtNs = self.dts['X5'], N = len(self.cableFreqResponseXpol))

        # After getting the freq response of the cables plus 10 dB coupler I need 
        # the frequency content of the pulse.
        # To do this, I take the case where the pulse went through the signal chain, WITHOUT
        # the five feet of cable, and then remove the cable response.
        self.pulseThruCopol, dt, t4 = getWaveform(baseDir + '140626_140317_ps_pulser_copol_fast_Ch1.csv', padToLength = padLength)
        #self.pulseThruCopol, dt, t4 = getWaveform(baseDir + '140626_140219_ps_pulser_copol_Ch1.csv', padToLength = padLength)
        self.dts['Co'] = dt
        self.t0s['Co'] = t4
        self.pulseThruCopol = removeAttenuationTimeDomain(self.pulseThruCopol, atten_dB = 20)
        self.pulseThruCopol = windowPulseAroundPeak(self.pulseThruCopol, 20, 134)
        self.pulseThruCopolFFT = doNormalizedFFT(self.pulseThruCopol, self.dts['Co'])

        self.pulseFreqs = deconvolveFreqToFreq(self.pulseThruCopolFFT, self.cableFreqResponseCopol)
        
        self.pulseThruXpol, dt, t5 = getWaveform(baseDir + '140626_140010_ps_pulser_xpol_fast_Ch1.csv', padToLength = padLength)
        
        self.dts['X'] = dt
        self.t0s['X'] = t5
        self.pulseThruXpol = removeAttenuationTimeDomain(self.pulseThruXpol, atten_dB = 20)
        self.pulseThruXpol = windowPulseAroundPeak(self.pulseThruXpol, 20, 134)



        for key, dt in self.dts.items():
            # Some kind of floating point tolerance
            assert dt - self.dts['Co'] < 0.00000001

        for key, l in self.lens.items():
            # Some kind of floating point tolerance
            assert l - self.lens['Co'] < 0.00000001

        
        
        self.dfMHz = 1e3/(len(self.pulse5ft)*self.dts['P5'])
        #assert len(set(self.dts)) == 1

        # Does what is says, we need this to do the Friis correction, since we were pulsing through VPOL with p52
        self.readInMeanVpolGain()


    def getGroupDelayNs(self, wave, dtNs, t0, pol = 'copol'):
        maxInd = findPulseMaxInd(wave)
        t_max = t0 + maxInd*dtNs

        print maxInd, t_max
        
        t_max2 = 0
        if pol == 'copol':
            maxInd2 = findPulseMaxInd(self.pulseThruCopol)
            t_max2 = maxInd2*self.dts['Co'] + self.t0s['Co']
            print maxInd2, t_max2
        elif pol == 'xpol':
            maxInd2 = findPulseMaxInd(self.pulseThruXpol)
            t_max2 = maxInd2*self.dts['X'] + self.t0s['X']
        else:
            raise Exception("pol argument must be 'copol' or 'xpol'")

        return t_max - t_max2

    def removeCopol(self, wave, dtNs):
        # In order to have the same df, we need to have 
        # N_1*dt_1 == N_2*dt_2. If that's not the case, we 
        # need to take action.

        numZerosToPad = len(self.pulseThruCopol)*self.dts['Co']/dtNs - len(wave)
        assert numZerosToPad == int(numZerosToPad) # If it's not an integer, things will be tricky.
        if numZerosToPad >=0:
            for zeroInd in xrange(int(numZerosToPad)):
                wave.append(0)
        else:
            print 'Warning! Deleting things! You should probably check this is OK!'
            for zeroInd in xrange(int(abs(numZerosToPad))):
                wave.pop()

        # So now the frequencies should be the same
        assert 1e3/(len(wave)*dtNs) == 1e3/(len(self.pulseThruCopol)*self.dts['Co'])

        fft_wave = doNormalizedFFT(wave, dtNs)
        withoutCables = deconvolveFreqToFreq(fft_wave, self.cableFreqResponseCopol)
        justSeaveyToSeavey = deconvolveFreqToFreq(withoutCables, self.pulseFreqs)
        return justSeaveyToSeavey

    def removeXpol(self, wave, dtNs):
        # In order to have the same df, we need to have 
        # N_1*dt_1 == N_2*dt_2. If that's not the case, we 
        # need to take action.

        numZerosToPad = len(self.pulseThruXpol5ft)*self.dts['X5']/dtNs - len(wave)
        assert numZerosToPad == int(numZerosToPad)
        #print numZerosToPad, len(self.pulseThruCopol5ft), self.dts['Co5'], dtNs, len(wave)
        if numZerosToPad >=0:
            for zeroInd in xrange(int(numZerosToPad)):
                wave.append(0)
        else:
            print 'Warning! Deleting things!'
            for zeroInd in xrange(int(abs(numZerosToPad))):
                wave.pop()

        assert 1e3/(len(wave)*dtNs) == 1e3/(len(self.pulseThruXpol5ft)*self.dts['X5'])

        fft_wave = doNormalizedFFT(wave, dtNs)
        withoutCables = deconvolveFreqToFreq(fft_wave, self.cableFreqResponseXpol)
        justSeaveyToSeavey = deconvolveFreqToFreq(withoutCables, self.pulseFreqs)
        return justSeaveyToSeavey
        
    def doFriisCorrection(self, relativePowSpec=None, freqsMHz = None, distMeters = -1, doSqrt = False):
        if relativePowSpec == None or freqsMHz == None or distMeters < 0:
            raise Exception('Need more information!')

        # If this isn't true, someone's done somethings stupid... probably me
        assert len(relativePowSpec) == len(freqsMHz)

        # Input is P1/P2, we want G2 (assume == G1 for now)
        c = 3e2 # speed of light in m/us since input freqs are in MHz
        separationFactors = [(f*4*math.pi*distMeters/c)**2 for f in freqsMHz]

        gain = []
        if doSqrt == False:
            gain = [relPowSpec*sepFact/mvr for relPowSpec, sepFact, mvr in zip(relativePowSpec, separationFactors, self.meanVpolResponse)]
        else:
            gain = [math.sqrt(relPowSpec*sepFact) for relPowSpec, sepFact in zip(relativePowSpec, separationFactors)]
        #gain = [math.sqrt(g) for g in gain]
        return gain

    def removeCopolCablesAndDoFriisCorrection(self, wave = None, dtNs = None, distMeters = -1, doSqrt = False):
        """
        The differing deltaFs after zero padding are starting to confuse me
        so let's have one function do everything and return a set of freqs as well...
        """
        if wave == None or dtNs == None or distMeters < 0:
            raise Exception('Need more information!')
        
        justSeaveyToSeavey = self.removeCopol(wave, dtNs)
        rps = [abs(js2s)**2 for js2s in justSeaveyToSeavey]
        df = 1e3/(len(self.pulseThruCopol)*self.dts['Co5'])
        freqs = [df*i for i, r in enumerate(rps)]

        gain = self.doFriisCorrection(relativePowSpec = rps, freqsMHz = freqs, 
                                      distMeters = distMeters, doSqrt = doSqrt)

        # Everything was padded in this class so the df of the input is the same of the df of the fast measurement

        return gain, freqs

    def removeXpolCablesAndDoFriisCorrection(self, wave = None, dtNs = None, distMeters = None, doSqrt = False):
        """
        The differing deltaFs after zero padding are starting to confuse me
        so let's have one function do everything and return a set of freqs as well...
        """
        if wave == None or dtNs == None or distMeters == 0:
            raise Exception('Need more information!')
        
        justSeaveyToSeavey = self.removeXpol(wave, dtNs)
        rps = [abs(js2s)**2 for js2s in justSeaveyToSeavey]
        df = 1e3/(len(self.pulseThruCopol)*self.dts['Co5'])
        freqs = [df*i for i, r in enumerate(rps)]

        gain = self.doFriisCorrection(relativePowSpec = rps, freqsMHz = freqs, distMeters = distMeters, doSqrt = doSqrt)

        # Everything was padded in this class so the df of the input is the same of the df of the fast measurement

        return gain, freqs



    def readInMeanVpolGain(self):
        self.meanVpolResponse = []
        self.vpolReponseFreqsMHz = []
        inFileName = 'meanVpolResponse.dat'
        with file(inFileName, 'r') as inFile:
            for lineInd, line in enumerate(inFile):
                if lineInd == 0:
                    continue
                else:
                    vals = line.split('\t')
                    self.meanVpolResponse.append(float(vals[0]))
                    self.vpolReponseFreqsMHz.append(float(vals[1]))

        if len(self.meanVpolResponse) > 0:
            print 'Read in mean vpol response! ' + str(len(self.meanVpolResponse)) + ' frequencies'
        else:
            raise Exception("Couldn't find file " + inFileName)
                    

def deconvolveTimeToFreq(v_1, dt_1, v_2, dt_2):
    """
    Do FFT then call main deconvolver
    """
    fft_1 = doNormalizedFFT(v_1, dt_1) 
    fft_2 = doNormalizedFFT(v_2, dt_2)
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

def deconvolveTimeToTime(v_1, dt_1, v_2, dt_2):
    """
    Do FFT then pass on to the next guy...
    """
    fft_1 = doNormalizedFFT(v_1, dt_1)
    fft_2 = doNormalizedFFT(v_2, dt_2)
    return deconvolveFreqToTime(fft_1, fft_2)

def doNormalizedInvFFT(fft_norm, df_MHz = 0, dtNs = 0):

    # Is this right!?
    fft_unnorm = []
    if df_MHz > 0:
        N = len(fft_norm)
        fft_unnorm = [z*df_MHz*1e-3 for z in fft_norm]
    elif dtNs > 0:
        fft_unnorm = [z/dtNs for z in fft_norm]
    else:
        raise Exception('You need to provide dtNs or df_MHz to do fourier stuff!')

    v = np.fft.ifft(fft_unnorm)
    v_r = [x.real for x in v]
    v_im = [x.imag for x in v]
    if abs(sum(v_im)) > 0.000000001:
        print 'Something dodgy about the fft, imag values remaining after inverse!'
        raise Exception('Inverse FFT not entirely real')

    return v_r

def doNormalizedFFT(wave, dt):
    fft_unnorm = np.fft.fft(wave)
    fft_norm = [dt*z for z in fft_unnorm]
    return fft_norm


def makePowerSpectrum(normFFT, dtNs = 0, df_MHz = 0):
    powSpec = []
    if df_MHz > 0:
        powSpec = [(abs(z)**2)*df_MHz*1e3 for z in normFFT]
    elif dtNs > 0:
        N = len(normFFT)
        powSpec = [(abs(z)**2)/(dtNs*N) for z in normFFT]
    else:
        raise Exception('You need to provide dtNs or df_MHz to do fourier stuff!')

    return powSpec

def dBScale(powSpec, zeroVal = 0):
    return [10*math.log10(ps) if ps > 0 else zeroVal for ps in powSpec]


def getPhaseFromFFT(theFFT):
    return [math.atan2(z.imag, z.real) for z in theFFT]
    
