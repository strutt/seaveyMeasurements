import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import math

counter = 0
vals = []
times = []
for line in file('140622_050202_picosecond_pulser_no_atten_vpol_txp52_rxp49.csv'):
    if counter > 5:
        words = line.split(',')
        vals.append(1e3*float(words[4]))
        times.append(1e9*float(words[3]))
        #print words
    counter += 1

fig1 = plt.figure()
plt.plot(times, vals)
plt.ylabel('Waveform (mV)')
plt.xlabel('Time (ns)')
plt.draw()

fig2 = plt.figure()
dt = times[1]-times[0]
freqs = [samp/(len(times)*dt) for samp in range(len(times))]
freqsBand = [freq for freq in freqs if freq < 2]

theFFT = np.fft.rfft(vals)
powSpec = np.abs(theFFT)**2
powSpec_dB = [10*math.log10(samp) for samp in powSpec]
plt.plot(freqsBand, powSpec_dB[:len(freqsBand)])
plt.ylabel('Power Spectrum (dB)')
plt.xlabel('Frequency (GHz)')
plt.draw()

fig3 = plt.figure()
phaseInfo = np.angle(theFFT)
phaseInfoDeg = [samp*360/math.pi for samp in phaseInfo]
plt.plot(freqsBand, phaseInfoDeg[:len(freqsBand)])
plt.ylabel('Phase (Degrees)')
plt.xlabel('Frequency (GHz)')


plt.show()
