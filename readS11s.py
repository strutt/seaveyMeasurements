#!/usr/bin/env python2
from matplotlib import pyplot as plt
import math
import numpy as np
from glob import glob

"""
Simple script to read and draw the S11 data.
"""

def main():

    filesV = glob('seaveyDataPalestine2014/S11s/*V.csv')
    filesH = glob('seaveyDataPalestine2014/S11s/*H.csv')
    filesV = orderGlobResults(filesV, 'V')
    filesH = orderGlobResults(filesH, 'H')

    nAnts = len(filesV)
    nAntsPerAxis = 1

    mean_vpol_s11_dB = []
    rms_vpol_s11_dB = []
    max_vpol_s11_dB = []
    min_vpol_s11_dB = []

    mean_vpol_phase = []
    rms_vpol_phase = []
    max_vpol_phase = []
    min_vpol_phase = []

    mean_hpol_s11_dB = []
    rms_hpol_s11_dB = []
    max_hpol_s11_dB = []
    min_hpol_s11_dB = []

    mean_hpol_phase = []
    rms_hpol_phase = []
    max_hpol_phase = []
    min_hpol_phase = []

    savePlots = False

    counter = 0
    for antInd, (fV, fH) in enumerate(zip(filesV, filesH)):
        if counter % nAntsPerAxis == 0:
            fig, axes = plt.subplots(2)
            plt.suptitle('Antenna P' + str(antInd+1) + ' S11')

            for ax in axes:
                ax.grid(b=True, which='major', color='black', linestyle='--')


        mags, phases, freqs = readS11(fV)
        mags_dB = [10*math.log10(m) for m in mags]
        #phases = [p*180./math.pi for p in phases]

        axes[0].plot(freqs, mags_dB, label='VPol')
        axes[1].plot(freqs, phases, label='VPol')

        if antInd == 0:
            mean_vpol_s11_dB = [0 for g in mags_dB]
            rms_vpol_s11_dB = [0 for g in mags_dB]
            max_vpol_s11_dB = [-100000 for g in mags_dB]
            min_vpol_s11_dB = [1000000 for g in mags_dB]
            mean_vpol_phase = [0 for g in phases]
            rms_vpol_phase = [0 for g in phases]
            max_vpol_phase = [-100000 for g in phases]
            min_vpol_phase = [1000000 for g in phases]
        mean_vpol_s11_dB = [m + g for m, g in zip(mean_vpol_s11_dB, mags_dB)]
        rms_vpol_s11_dB = [r + g**2 for r, g in zip(rms_vpol_s11_dB, mags_dB)]
        max_vpol_s11_dB = [g if g > maxG else maxG for g, maxG in zip(mags_dB, max_vpol_s11_dB)]
        min_vpol_s11_dB = [g if g < minG else minG for g, minG in zip(mags_dB, min_vpol_s11_dB)]
        mean_vpol_phase = [m + g for m, g in zip(mean_vpol_phase, phases)]
        rms_vpol_phase = [r + g**2 for r, g in zip(rms_vpol_phase, phases)]
        max_vpol_phase = [g if g > maxG else maxG for g, maxG in zip(phases, max_vpol_phase)]
        min_vpol_phase = [g if g < minG else minG for g, minG in zip(phases, min_vpol_phase)]


        mags, phases, freqs = readS11(fH)
        mags_dB = [10*math.log10(m) for m in mags]
        #phases = [p*180./math.pi for p in phases]

        axes[0].plot(freqs, mags_dB, label='HPol')
        axes[1].plot(freqs, phases, label='HPol')

        if antInd == 0:
            mean_hpol_s11_dB = [0 for g in mags_dB]
            rms_hpol_s11_dB = [0 for g in mags_dB]
            max_hpol_s11_dB = [-100000 for g in mags_dB]
            min_hpol_s11_dB = [1000000 for g in mags_dB]
            mean_hpol_phase = [0 for g in phases]
            rms_hpol_phase = [0 for g in phases]
            max_hpol_phase = [-100000 for g in phases]
            min_hpol_phase = [1000000 for g in phases]
        mean_hpol_s11_dB = [m + g for m, g in zip(mean_hpol_s11_dB, mags_dB)]
        rms_hpol_s11_dB = [r + g**2 for r, g in zip(rms_hpol_s11_dB, mags_dB)]
        max_hpol_s11_dB = [g if g > maxG else maxG for g, maxG in zip(mags_dB, max_hpol_s11_dB)]
        min_hpol_s11_dB = [g if g < minG else minG for g, minG in zip(mags_dB, min_hpol_s11_dB)]
        mean_hpol_phase = [m + g for m, g in zip(mean_hpol_phase, phases)]
        rms_hpol_phase = [r + g**2 for r, g in zip(rms_hpol_phase, phases)]
        max_hpol_phase = [g if g > maxG else maxG for g, maxG in zip(phases, max_hpol_phase)]
        min_hpol_phase = [g if g < minG else minG for g, minG in zip(phases, min_hpol_phase)]

        counter += 1
        if counter % nAntsPerAxis == 0 or counter == nAnts - 1:
            axes[0].legend(loc='lower right', fancybox = True)
            axes[1].legend(loc='lower right', fancybox = True)
            for ax in axes:
                ax.set_xlabel('Frequency (MHz)')
                ax.set_xlim([0, 1500])
            axes[0].set_ylabel('Power (dB)')
            axes[1].set_ylabel('Group delay (ns)')
            #axes[1].set_ylim([-180, 180])
            if savePlots == True:
                antNum = antInd + 1
                fileName = ''
                if antNum <= 9:
                    fileName = 's11p0' + str(antNum)
                else:
                    fileName = 's11p' + str(antNum)
                fig.savefig('measurementSummaryDocs/'+fileName+'.png',dpi=100)




    n = counter


    # Finalize general calc
    mean_vpol_s11_dB = [m/n for m in mean_vpol_s11_dB]
    rms_vpol_s11_dB = [math.sqrt(r/n-m**2) for r, m in zip(rms_vpol_s11_dB, mean_vpol_s11_dB) ]
    mean_hpol_s11_dB = [m/n for m in mean_hpol_s11_dB]
    rms_hpol_s11_dB = [math.sqrt(r/n-m**2) for r, m in zip(rms_hpol_s11_dB, mean_hpol_s11_dB) ]
    mean_vpol_phase = [m/n for m in mean_vpol_phase]
    rms_vpol_phase = [math.sqrt(r/n-m**2) for r, m in zip(rms_vpol_phase, mean_vpol_phase) ]
    mean_hpol_phase = [m/n for m in mean_hpol_phase]
    rms_hpol_phase = [math.sqrt(r/n-m**2) for r, m in zip(rms_hpol_phase, mean_hpol_phase) ]

    fig = plt.figure()
    plt.title('Vertical Polarization Antenna S11')
    plt.plot(freqs, mean_vpol_s11_dB, label = 'Mean')
    plt.plot(freqs, [m+r for r, m in zip(rms_vpol_s11_dB, mean_vpol_s11_dB)], label = 'Mean + RMS')
    plt.plot(freqs, [m-r for r, m in zip(rms_vpol_s11_dB, mean_vpol_s11_dB)], label = 'Mean - RMS')
    plt.plot(freqs, [m for m in max_vpol_s11_dB], label = 'Bin-by-bin maximum')
    plt.plot(freqs, [m for m in min_vpol_s11_dB], label = 'Bin-by-bin minimum')
    plt.grid(b=True, which='major', color='black', linestyle='--')
    plt.legend(loc='lower right', fancybox=True)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('S11 (dB)')
    if savePlots == True:
        fig.savefig('measurementSummaryDocs/s11VPolSummary.png',dpi=100)


    plt.figure()
    plt.title('Horizontal Polarization Antenna S11')
    plt.plot(freqs, mean_hpol_s11_dB, label = 'Mean')
    plt.plot(freqs, [m+r for r, m in zip(rms_hpol_s11_dB, mean_hpol_s11_dB)], label = 'Mean + RMS')
    plt.plot(freqs, [m-r for r, m in zip(rms_hpol_s11_dB, mean_hpol_s11_dB)], label = 'Mean - RMS')
    plt.plot(freqs, [m for m in max_hpol_s11_dB], label = 'Bin-by-bin maximum')
    plt.plot(freqs, [m for m in min_hpol_s11_dB], label = 'Bin-by-bin minimum')
    plt.grid(b=True, which='major', color='black', linestyle='--')
    plt.legend(loc='lower right', fancybox=True)
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('S11 (dB)')
    if savePlots == True:
        fig.savefig('measurementSummaryDocs/s11HPolSummary.png',dpi=100)

    plt.show()


def orderGlobResults(fileList, polChar):

    antNumber = []
    for f in fileList:
        a = f.split('P')[-1].split(polChar)[0]
        antNumber.append(int(a))
        

    fileListSorted = fileList[:]
    for i, a in enumerate(antNumber):
        fileListSorted[a-1] = fileList[i] 
    return fileListSorted
    
        

def readS11(fileName):

    freqs = []
    reals = []
    imags = []

    for lineInd, line in enumerate(file(fileName)):
    #print lineInd, line
        if 'END' in line:
            break
        if lineInd < 18:
            continue
        else:
            vals = line.split(',')
            freqs.append(1e-6*float(vals[0]))
            reals.append(float(vals[1]))
            imags.append(float(vals[2]))

    dw = 2*math.pi*1e6*(freqs[1] - freqs[0])
    mags = [re*re+im*im for re, im in zip(reals, imags)]
    #phases = [math.asin(im/math.sqrt(mag)) for re, im, mag in zip(reals, imags, mags)]
    phases = [math.atan2(im, re) for im, re in zip(imags, reals)] #math.sqrt(mag)) for re, im, mag in zip(reals, imags, mags)]
    phases = np.unwrap(phases)
    for i, (p, f) in enumerate(zip(phases, freqs)):
        if i > 0:
            if p > phases[i-1]:
                print f
    #phases = [-p/dw for p in phases]

    return mags, phases, freqs

if __name__ == '__main__':
    main()



