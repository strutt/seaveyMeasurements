from glob import glob
import copy
from os import rename
from matplotlib import pyplot as plt
from doS21analysis import *

def main():
    """
    Quick script to look for common typos in the naming of S21Palestine data
    """
    # csvFiles = glob('S21Palestine/*.csv')
    
    # for fileName in csvFiles:

    #     if 'rpx' in fileName:
    #         print fileName
    #         fileName2 = fileName.replace('rpx', 'rxp')
    #         rename(fileName, fileName2)
    #         print 'Renamed ' + fileName + ' as ' +  fileName2


    # csvFiles = glob('S21Palestine/*.csv')
    # for fileName in csvFiles:

    #     if 'rdp' in fileName:
    #         print fileName
    #         fileName2 = fileName.replace('rdp', 'rxp')
    #         rename(fileName, fileName2)
    #         print 'Renamed ' + fileName + ' as ' +  fileName2

    # csvFiles = glob('S21Palestine/*rxp52*.csv')
    # for fileName in csvFiles:

    #     if 'rxp52' in fileName:
    #         print fileName
    #         fileName2 = fileName.replace('rxp52', 'rxp51')
    #         rename(fileName, fileName2)
    #         print 'Renamed ' + fileName + ' as ' +  fileName2


    # csvFiles = glob('S21Palestine/*rxp51*.csv')
    # for fileName in csvFiles:

    #     if 'pol' not in fileName:
    #         print fileName
    #         fileName2 = fileName.replace('rxp51', 'rxp51_vpol')
    #         rename(fileName, fileName2)
    #         print 'Renamed ' + fileName + ' as ' +  fileName2

    #csvFiles = glob('S21Palestine/*_ch*.*')
    
    # rxp29_vpol_-30az+00el_off_Ch
    #for fileName in csvFiles:

#        if '_ch' in fileName:
#            print fileName
#            fileName2 = fileName.replace('_ch', '_Ch')
#            rename(fileName, fileName2)
#            print 'Renamed ' + fileName + ' as ' +  fileName2
            
    checkOffAxisMeasurementFiles(1)

def checkOffAxisMeasurementFiles(moveFlag = 0):

    ants = ['rxp25']
    pols = ['_vpol_', '_hpol_']
    #pols = ['_hpol_']
    puls = ['_','_off_']
    puls2 = ['','60dBatten_']
    chs = ['Ch1', 'Ch2', 'Ch3', 'Ch4']
    azs = ['-45az', '-30az', '-15az', '+00az', '+15az', '+30az', '+45az']
    els = ['-45el', '-30el', '-15el', '+00el', '+15el', '+30el', '+45el']

    for ant in ants:
        for pol in pols:
            for azInd, az in enumerate(azs):
                for elInd, el in enumerate(els):
                    specialCase = 0
                    for pulInd, pul in enumerate(puls):
                        for ch in chs:
                            globString = '*' + ant + pol + az + el + pul + ch + '.csv'

                            # Special boresight case
                            if azInd == 3 and azInd ==3:
                                globString = '*' + ant + pol + puls2[pulInd] + ch + '.csv'

                            # Allow to carry on if along 00 in either axis
                            elif azInd == 3:
                                pass
                            elif elInd == 3:
                                pass

                            # Allows to carry on if one of the off axis measurements
                            elif elInd > 2 and elInd < 4 and azInd > 2 and azInd < 4:
                                pass

                            # Otherwise skip this combo
                            else:
                                continue

                            l = glob('S21Palestine/' + globString)
                            if len(l) != 1:
                                print globString + ' finds:' 
                                for m in l:
                                    print m                        
                                print ''

                            if moveFlag > 0:
                                fileName = l[0]
                                fileName2 = l[0].split('/')[0] + '/offAxisMeasurements/' + l[0].split('/')[-1]
                                #print fileName, fileName2
                                #rename(fileName, fileName2)

    return 0




if __name__ == '__main__':
    main()

