from doS21analysis import *
from scipy import signal
import numpy as np

def main():
    """
    Figure out phase centre from timing information with cable delays...

    I should be able to use information in the *ps* files, along with the
    measured geometry to figure out where the phase centre of the Seaveys are.

    Notes:
    Each of the cables connected to the antennas are listed here.

    The ps_pulser file 140623_052608 was just through 5ft of LMR cable on the pulse line!!!!!!!!!
    So the ones labelled copol and xpol all through entire signal chain, minus antennas, 
    plus 5ft of LMR, which needs to be subtracted out.

    Doing that subtraction gives the cable delays!

    The ps_pulse file (note the missing 'r', pulser->pulse) 
    140623_051312 should be ignored as we were triggering on the RFpulse itself, not
    the trigger.    
    """


    cableDelayCopol, cableDelayXpol = getCableDelays()
    print cableDelayCopol, cableDelayXpol

    azimuths = [-45, -30, -15, 0, 15, 30, 45]
    elevations = [-45, -30, -15, 0, 15, 30, 45]

    #    files = glob('S21Palestine/*rxp30_hpol*.csv')
    #    for f in files:
    #        print f

    # fileStr = anglesToStr(0, 0)
    # v1 = []
    # globStr = 'S21Palestine/*rxp30_hpol_' + fileStr + '*.csv'
    # print 'Found files for az = ' + str(0) + ', el = ' + str(0) + '!:'
    # files = glob(globStr)
    # for f in files:
    #     if 'Ch1' in f and 'off' not in f and '60dB' not in f :
    #         v1 = getWaveform(f)

    deltaTs = [[0 for x in xrange(len(azimuths))] for x in xrange(len(elevations))]
    sumDts=0
    counter=0

    #r = 10 #m 33.4ft correct this...
    r = 8.89 #m 33.4ft correct this...

    bl = 8.79
    br = 8.84
    tr = 8.83
    tl = 8.78
    r = (bl+br+tr+tl)/4
    r = 8.89
    print r

    timeThroughAirTheory = r/3e8*1e9 #ns
    print 'timeThroughAirTheory = ' + str(timeThroughAirTheory)

    for elInd, el in enumerate(elevations):
        for azInd, az in enumerate(azimuths):
            fileStr = anglesToStr(az, el)
            globStr = 'S21Palestine/*rxp30_hpol_' + fileStr + '*.csv'
            files = glob(globStr)

            if len(files) > 0:
                print 'Found files for az = ' + str(az) + ', el = ' + str(el) + '!:'
                for f in files:
                    #print f
                    if ('Ch1' in f or 'ch1' in f) and 'off' not in f and '60dB' not in f :
                        v, dt, t0 = getWaveform(f)
                        peakTime = findAbsPeakTime(v, dt, t0)
                        #times = [t0 + dt*i for i in range(len(v))]
                        #vEnv = getHilbertEnvelope(v)
                        #peakTime = times[vEnv.index(max(vEnv))]
                        #deltaTs[elInd][azInd] = getDeltaTFromCrossCorr(v1, v)
                        timeThroughAirMeasured = (peakTime - cableDelayCopol)
                        deltaTs[elInd][azInd] =  (timeThroughAirTheory - timeThroughAirMeasured)*3e8/1e9
                        print deltaTs[elInd][azInd]

                        if deltaTs[elInd][azInd]<0:
                            deltaTs[elInd][azInd] = 0

                        
                        sumDts += deltaTs[elInd][azInd]
                        counter += 1

                        if az == 0 and abs(el) == 45:
                            plt.plot([t0 + dt*i for i in range(len(v))],v, label = str(el))
                            #plt.plot([t0 + dt*i for i in range(len(v))],vEnv, label = str(el))
                    elif 'Ch2' in f and 'off' not in f and '60dB' not in f :
                        v, dt, t0 = getWaveform(f)
                        #print findAbsPeakTime(v, dt, t0)

    plt.legend()
    meanDeltaT = sumDts/counter
    for elInd, el in enumerate(elevations):
        for azInd, az in enumerate(azimuths):
            if deltaTs[elInd][azInd] == 0:
                deltaTs[elInd][azInd] = meanDeltaT
    
    fig2 = plt.figure()
    plt.imshow(deltaTs, origin='lower', aspect='auto')
    #self.implot = mp.imshow(self.image_array, origin='lower', aspect='auto', interpolation='nearest', extent=self.limits)
    plt.colorbar()
    plt.show()

def getDeltaTFromCrossCorr(v1, v2):
    cc = np.correlate(np.array(v1), np.array(v2), old_behavior='False')
    return cc.index(max(cc))


def getHilbertEnvelope(v):

    ht = signal.hilbert(v)
    he = [math.sqrt(vi*vi + abs(hi)**2) for vi, hi in zip(v, ht)]
    return he


def findAbsPeakTime(vals, dt, t0):
    
    times = [t0 + i*dt for i in range(len(vals))]
    maximum = max(vals)
    minimum = min(vals)
    if maximum > abs(minimum):
        return times[vals.index(maximum)]
    else:
        return times[vals.index(minimum)]

    return peakTime


def anglesToStr(az, el):
    azStr = ''
    elStr = ''
    fileStr = ''

    # azimuth part
    if az > 0:
        azStr = '+' + str(az) + 'az'
    elif az < 0:
        azStr = str(az) + 'az'
    else:
        azStr = '+00az'

    #elevation part
    if el > 0:
        elStr = '+' + str(el) + 'el'
    elif el < 0:
        elStr = str(el) + 'el'
    else:
        elStr = '+00el'


    fileStr = azStr + elStr
    if el==0 and az == 0:
        fileStr = 'Ch'

    return fileStr
    


def getCableDelays():
    """
    Open cable measurement files and get figure out cable delays
    Ch1 is the RF pulse, Ch2 is the trigger (should trigger at t=0)
    """

    files = glob('S21Palestine/*ps_pulser*.csv')

    timeJustLMR = 0
    timeCopol = 0
    timeXpol = 0

    for f in files:
        vals, dt, t0 = getWaveform(f)
        times = [t0 + i*dt for i in range(len(vals))]

        # BE AWARE... DIFFERING DTS!!
        #print len(times), dt, len(vals)

        if 'Ch1' in f:
            #print f

            peakTime = times[vals.index(max([abs(v) for v in vals]))]
            #print peakTime
            if '_051836_' in f:
                timeJustLMR = peakTime
            elif 'copol' in f:
                timeCopol = peakTime
            elif 'xpol' in f:
                timeXpol = peakTime
            #plt.plot(times, vals, label=f)

            #pass

        elif 'Ch2' in f:
            print findAbsPeakTime(vals, dt, t0)


    ################################
    # WHY!?!?!?!?!?
    ################################
    #timeCopol -= timeJustLMR
    #timeXpol -= timeJustLMR

    #plt.legend()
    #plt.show()
    

    return timeCopol, timeXpol


if __name__ == '__main__':
    main()

