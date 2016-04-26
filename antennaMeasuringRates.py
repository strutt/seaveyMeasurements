import matplotlib.pyplot as plt
import numpy as np
from numpy import fft
import math
from glob import glob
import datetime
from matplotlib import dates





def main():
    """
    Takes in csv from TDS6804B and plots time stamps in the file name.
    Just a bit of fun ;)
    """
    
    dstamps = []
    tstamps = []

    listOfAnts = ['rxp' + str(antInd+1) if antInd+1 > 9 else 'rxp0' + str(antInd+1) for antInd in range(48)]

    for antInd, ant in enumerate(listOfAnts):
        # Glob is the module that string searches for files
        globString = 'S21Palestine/*' + ant + '*.csv'
        listOfFiles = glob(globString)

        if len(listOfFiles) > 0:

            dstamp = listOfFiles[-1].split('_')[0]
            tstamp = listOfFiles[-1].split('_')[1]
            secs = int(tstamp[-1]) + 10*int(tstamp[-2])
            mins = int(tstamp[-3]) + 10*int(tstamp[-4])
            # +4 since scope clock is 4 hours behind Palestine time
            hours = int(tstamp[-5]) + 10*int(tstamp[-6]) + 4
            day = int(dstamp[-1]) + 10*int(dstamp[-2])
            month = int(dstamp[-3]) + 10*int(dstamp[-4])
            year = 2000 + int(dstamp[-5]) + 10*int(dstamp[-6])

            tstamps.append(datetime.datetime(year, month, day, hours, mins, secs))

    figT, axT = plt.subplots(1)
    axT.plot(sorted(tstamps), range(len(tstamps)))
    plt.suptitle('Work rate')
    axT.set_xlabel('Time')
    axT.set_ylabel('Antenna S21 measurements completed')
    hfmt = dates.DateFormatter('%m/%d %H:%M')
    axT.xaxis.set_major_locator(dates.DayLocator())
    axT.xaxis.set_major_formatter(hfmt)
    #plt.gcf().autofmt_xdate()
    plt.xticks(rotation='vertical')
    plt.subplots_adjust(bottom=.3)

    plt.show()

    return 0

if __name__ == '__main__':
    main()
