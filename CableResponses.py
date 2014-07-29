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


        # From these waveforms I can get the cable response...
        # Pulser also goes through 5ft of extra cable
        # but that was also attached to the cables for this measurement
        self.pulse5ft, dt, t0 = getWaveform(baseDir + '140626_134852_ps_pulser_direct_fast_Ch1.csv', padToLength = padLength)
        self.dts['P5'] = dt
        #print self.dts['P5'] 
        self.pulseThruCopol5ft, dt, t1 = getWaveform(baseDir + '140626_135504_ps_pulser_copol_fast_5ft_Ch1.csv', padToLength = padLength)
        self.dts['Co5'] = dt
        self.pulseThruXpol5ft, dt, t2 = getWaveform(baseDir + '140626_135703_ps_pulser_xpol_fast_5ft_Ch1.csv', padToLength = padLength)
        self.dts['X5'] = dt


        # Need to account for 20dB attenuator on pulser measurement!

        self.pulse5ft = removeAttenuationTimeDomain(self.pulse5ft, atten_dB = 20)
        self.pulseThruCopol5ft = removeAttenuationTimeDomain(self.pulseThruCopol5ft, atten_dB = 20)
        self.pulseThruXpol5ft = removeAttenuationTimeDomain(self.pulseThruXpol5ft, atten_dB = 20)
        

        # Window around pulses to make the frequency domain lovely and smooth.
        # These numbers were picked by looking at the waveforms and playing until samples contained just the pulse
        # I think they are fine
        self.pulse5ft = windowPulseAroundPeak(self.pulse5ft,20, 134)
        self.pulseThruCopol5ft = windowPulseAroundPeak(self.pulseThruCopol5ft, 20, 134)
        self.pulseThruXpol5ft = windowPulseAroundPeak(self.pulseThruXpol5ft, 20, 134)

        # Do the first round of deconvolution, which gets just the cable responses and the
        # effects of the 10dB coupler on that line.
        self.cableFreqResponseCopol = deconvolveTimeToFreq(self.pulseThruCopol5ft, self.pulse5ft)
        self.copolFreqs = makeFreqsMHz(dtNs = self.dts['Co5'], N = len(self.cableFreqResponseCopol))

        self.cableFreqResponseXpol = deconvolveTimeToFreq(self.pulseThruXpol5ft, self.pulse5ft)
        self.xpolFreqs = makeFreqsMHz(dtNs = self.dts['X5'], N = len(self.cableFreqResponseXpol))

        # After getting the freq response of the cables plus 10 dB coupler I need 
        # the frequency content of the pulse.
        # To do this, I take the case where the pulse went through the signal chain, WITHOUT
        # the five feet of cable, and then remove the cable response.
        self.pulseThruCopol, dt, t4 = getWaveform(baseDir + '140626_140317_ps_pulser_copol_fast_Ch1.csv', padToLength = padLength)
        self.dts['Co'] = dt

        self.pulseThruCopol = removeAttenuationTimeDomain(self.pulseThruCopol, atten_dB = 20)
        self.pulseThruCopol = windowPulseAroundPeak(self.pulseThruCopol, 20, 134)
        self.pulseThruCopolFFT = np.fft.fft(self.pulseThruCopol)

        self.pulseFreqs = deconvolveFreqToFreq(self.pulseThruCopolFFT, self.cableFreqResponseCopol)

        for key, dt in self.dts.items():
            # Some kind of floating point tolerance
            assert dt - self.dts['Co'] < 0.00000001

        self.dfMHz = 1e3/(len(self.pulse5ft)*self.dts['P5'])
        #assert len(set(self.dts)) == 1

        # Does what is says, we need this to do the Friis correction, since we were pulsing through VPOL with p52
        self.readInMeanVpolGain()

    def removeCopolAndPulse(self, wave, dtNs):
        # In order to have the same df, we need to have 
        # N_1*dt_1 == N_2*dt_2. If that's not the case, we 
        # need to take action.

        numZerosToPad = len(self.pulseThruCopol)*self.dts['Co']/dtNs - len(wave)
        assert numZerosToPad == int(numZerosToPad) # If it's not an integer, things will be tricky.
        print numZerosToPad, len(self.pulseThruCopol), self.dts['Co'], dtNs, len(wave)
        if numZerosToPad >=0:
            for zeroInd in xrange(int(numZerosToPad)):
                wave.append(0)
        else:
            print 'Warning! Deleting things! You should probably check this is OK!'
            for zeroInd in xrange(int(abs(numZerosToPad))):
                wave.pop()

        # So now the frequencies should be the same
        assert 1e3/(len(wave)*dtNs) == 1e3/(len(self.pulseThruCopol)*self.dts['Co'])

        fft_wave = np.fft.fft(wave)
        deco = deconvolveFreqToFreq(fft_wave, self.pulseThruCopolFFT)
        
        return deco


    def removeCopol(self, wave, dtNs):
        # In order to have the same df, we need to have 
        # N_1*dt_1 == N_2*dt_2. If that's not the case, we 
        # need to take action.

        numZerosToPad = len(self.pulseThruCopol)*self.dts['Co']/dtNs - len(wave)
        assert numZerosToPad == int(numZerosToPad) # If it's not an integer, things will be tricky.
        print numZerosToPad, len(self.pulseThruCopol), self.dts['Co'], dtNs, len(wave)
        if numZerosToPad >=0:
            for zeroInd in xrange(int(numZerosToPad)):
                wave.append(0)
        else:
            print 'Warning! Deleting things! You should probably check this is OK!'
            for zeroInd in xrange(int(abs(numZerosToPad))):
                wave.pop()

        # So now the frequencies should be the same
        assert 1e3/(len(wave)*dtNs) == 1e3/(len(self.pulseThruCopol)*self.dts['Co'])

        fft_wave = np.fft.fft(wave)
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

        fft_wave = np.fft.fft(wave)
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
    

    




