import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import math
from glob import glob
import datetime
from matplotlib import dates





def main():
    """
    Takes in csv from TDS6804B and processes waveform
    """

    maxFreqMHz = 2000
    padToLength = 8192*2

    pulseWindow = 100 #ns
    prePeakWindow = 10 #ns
    postPeakWindow = pulseWindow - prePeakWindow
    
    highPass = 100
    lowPass = 1300
    

    #listOfAnts = ['rxp01', 'rxp02', 'rxp03', 'rxp04', 'rxp05', 'rxp06']
    #listOfAnts = ['rxp07', 'rxp08', 'rxp09', 'rxp10', 'rxp11', 'rxp12']
    #listOfAnts = ['rxp13', 'rxp14', 'rxp15', 'rxp16', 'rxp17', 'rxp18']
    #listOfAnts = ['rxp19', 'rxp20', 'rxp21', 'rxp22', 'rxp23', 'rxp24']
    #listOfAnts = ['rxp25', 'rxp26', 'rxp27', 'rxp28', 'rxp29', 'rxp30']
    #listOfAnts = ['rxp31', 'rxp32', 'rxp33', 'rxp34', 'rxp35', 'rxp36']
    #listOfAnts = ['rxp37', 'rxp38', 'rxp39', 'rxp40', 'rxp41', 'rxp42']
    #listOfAnts = ['rxp43', 'rxp44', 'rxp45', 'rxp46', 'rxp47', 'rxp48']
    #listOfAnts = ['rxp49', 'rxp50', 'rxp51']
    listOfAnts = ['rxp30']
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

    azimuths = [-45, -30, -15, 0, 15, 30, 45]
    elevations = [-45, -30, -15, 0, 15, 30, 45]
    for elInd, el in enumerate(elevations):
        for azInd, az in enumerate(azimuths):
            print el, az


    # ax = fig.add_subplot(111)
    # image_buffer = image_buffer.split(' ')
    # image_buffer.pop()
    # image_buffer = [data_type(x) for x in image_buffer]
    # print max(image_buffer)
    # image_buffer = [image_buffer[y:y+256] for y in range(0, 262144, 256)]
    # image_array = np.array(image_buffer)
    # image_array = image_array.transpose()

    #return 0

    for antInd, ant in enumerate(listOfAnts):
        print 'Doing analysis for antenna ' + str(ant)
        # Glob is the module that string searches for files
        globString = 'S21Palestine/*' + ant + '*.csv'
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
                    if (chan in f or chan.lower() in f) and pol in f and 'pulser' not in f:
                        v, dt = getWaveform(f)
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
                #newV = [v1 for v1  in waves[normFlag]]

                # Take the absolute of the subtracted voltage, and find the time of the absolute maximum
                if chan == 'Ch1':
                    absNewV = [abs(v1) for v1 in newV]
                    indexOfAbsMax = absNewV.index(max(absNewV))

                # Find pulse window in terms of sample number
                windowStart = indexOfAbsMax - int(prePeakWindow/dt)
                windowEnd =  indexOfAbsMax + int(postPeakWindow/dt)
                #print windowStart*dt, windowEnd*dt

                #ds = pulseWindow/(2*dt) # /2 since plus minus peak

                # Zero everything not in pulse window around around the pulse peak
                #windowedPulse = [v1 if abs(i-indexOfAbsMax) < ds else 0 for i, v1 in enumerate(newV)]
                windowedPulse = [v1 if i > windowStart and i < windowEnd else 0 for i, v1 in enumerate(newV)]
                

                # Limit plots...
                x0 = windowStart - 100
                x1 = windowEnd + 100
                #axesA[chanInd, polInd].plot(times[x0:x1], newV[x0:x1], label = ant)
                axesA[chanInd, polInd].plot(times[x0:x1], windowedPulse[x0:x1], label = ant)
                axesA[chanInd, polInd].legend()
                axesA[chanInd, polInd].set_title(chanToPol[chan] + pol.capitalize())
                axesA[chanInd, polInd].set_xlabel('Time (ns)')
                axesA[chanInd, polInd].set_ylabel('Amplitude (mV)')
                #axesA[chanInd, polInd].xlim([390, 420])

                f, pwr, phase = getPowerSpectrumInfo(windowedPulse, dt, maxFreqMHz)
                pwrSpecDensity = [10*math.log10(p) for p in pwr]

                if antInd == 0:
                    axesB[chanInd, polInd].plot(f, nPSD[polInd], label = 'Noise only ' + pol, linestyle='dotted', color='black')
                axesB[chanInd, polInd].plot(f, pwrSpecDensity, label = ant)
                axesB[chanInd, polInd].legend()
                axesB[chanInd, polInd].set_title(pol.capitalize() + ' pulsed directly, ' + chanToPol[chan] + 'feed')
                axesB[chanInd, polInd].set_xlabel('Frequencies (MHz)')
                axesB[chanInd, polInd].set_ylabel('Power (dB)')

                if chan is 'Ch1':
                    while len(relativeCrossPol) < len(pwr):
                        relativeCrossPol.append(0)
                    relativeCrossPol = [-p for p in pwrSpecDensity]
                else:
                    relativeCrossPol = [rcp+p for rcp, p in zip(relativeCrossPol, pwrSpecDensity) ]


                axesD[chanInd, polInd].plot(f, phase, label = ant)
                axesD[chanInd, polInd].legend()
                axesD[chanInd, polInd].set_title(pol.capitalize() + ' pulsed directly, ' + chanToPol[chan] + 'feed')
                axesD[chanInd, polInd].set_xlabel('Frequencies (MHz)')
                axesD[chanInd, polInd].set_ylabel('Phase (Degrees)')



            #plt.plot(f, [math.log10(rcp) for rcp in relativeCrossPol], label = 'RelativeCrossPol')
            bandPassRCP = [rcp if fVal > highPass and fVal < lowPass else -50 for rcp, fVal in zip(relativeCrossPol, f)]
            #bandPassRCP_dB = [10*math.log10(rcp) if rcp > 0 else -75 for rcp in bandPassRCP]
            axesC[polInd].plot(f, bandPassRCP)
            axesC[polInd].set_title('Cross polarization: ' + pol.capitalize() + ' to ' + listOfPols[1-polInd].capitalize())
            axesC[polInd].set_ylabel('Power fraction (dB)')
            axesC[polInd].set_xlabel('Frequency (MHz)')
            
            #plt.plot(f, relativeCrossPol, label = 'RelativeCrossPol')
            #plt.plot(f, bandPassRCP, label = 'RelativeCrossPol')
            #plt.legend()

    plt.show()

    return 0




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
            vals.append(1e3*float(words[4])) # 1e3 for V -> mV
            times.append(1e9*float(words[3]))
        counter += 1

    dt = times[1]-times[0]
    trueLen = len(times)
    while len(times) < padToLength:
        times.append(times[-1] + dt)
        vals.append(0)

    return vals, dt


def getPowerSpectrumInfo(vals, dt, maxFreqMHz = 0):
    """
    Return the frequencies, power (dB) and phase (Deg) from the volts and times passed to the function.
    Limits the maximum frequency if 'maxFreqMHz' > 0 is passed.
    """

    # Get frequencies from sample times
    freqs = [1e3*samp/(len(vals)*dt) for samp in range(len(vals))]
    freqsBand = []

    # Limit frequency information returned by maximum freq optional arg
    if maxFreqMHz > 0:
        freqsBand = [freq for freq in freqs if freq < maxFreqMHz]
    else:
        freqsBand = freqs

    # Do FFT and convert to power spectum and phase
    #theFFT = np.fft.rfft(vals)
    theFFT = np.fft.fft(vals)
    powSpec = getPowerSpectrum(theFFT, dt, len(vals)) 
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



def getPowerSpectrum(theFFT, dt, N):
    """
    Normalizes the power spectra according to Ryan's scheme...
    see http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf
    """
    powSpec = np.abs(theFFT)**2
    #print powSpec
    #powSpec = powSpec*dt*dt
    #powSpec = powSpec*dt/N
    powSpec = powSpec/N
    #print powSpec
    #print ''
    return powSpec

if __name__ == '__main__':
    main()
