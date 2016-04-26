from matplotlib import pyplot as plt
import math


"""
From Pk-Pk numbers read off screen by Abby during data taking for P30.
"""

angles = [-45, -30, -15, 0, 15, 30, 45] # degrees
PkPk = [16, 28, 48, 60, 43, 27, 15] #mV
PkPkSqr = [x**2 for x in PkPk]
PkPkSqr_dB = [10*math.log10(x) for x in PkPkSqr]
PkPkSqr_dB2 = [x - max(PkPkSqr_dB) for x in PkPkSqr_dB]

plt.plot(angles, PkPkSqr_dB2, marker='o')
plt.xlabel('Angle off boresight (Degrees)')
plt.ylabel('Power (dB)')
plt.grid(b=True, which='major', color='r', linestyle='--')
plt.show()
