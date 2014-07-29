#!/usr/bin/env python2
"""
Because copying and pasting code is for chumps
"""

for antInd in range(51):
    antNum = antInd + 1

    for plotInd in range(2):
        baseFileName = ''
        if plotInd == 0:
            baseFileName = 's11p'
        else:
            baseFileName = 'rxp'

        fileName = ''
        if antNum <= 9:
            fileName = baseFileName + '0' + str(antNum)
        else:
            fileName = baseFileName + str(antNum)
        print '\\begin{figure}[H]'
        print '\t\\begin{center}'
        print '\t\t\\includegraphics[width=11cm]{' + fileName + '.png}'
        if plotInd == 0:
            print '\t\t\\caption{S11 plot of power and phase for antenna P' + str(antNum) + ', 100-1300 MHz.}'
        else:
            print '\t\t\\caption{S21 gain and fractional cross polarization for antenna P' + str(antNum) + ', 100-1300 MHz.}'
        print '\t\t\\label{'+fileName+'plot}'
        print '\t\\end{center}'
        print '\\end{figure}'
        print ''
        print ''
