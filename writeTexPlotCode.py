#!/usr/bin/env python2
"""
Because copying and pasting code is for chumps
"""

for antInd in range(51):

    antNum = antInd + 1
    fileName = ''
    if antNum <= 9:
        fileName = 's11p0' + str(antNum)
        #fileName = 'rxp0' + str(antNum)
    else:
        fileName = 's11p' + str(antNum)
        #fileName = 'rxp' + str(antNum)
    print '\\begin{figure}[H]'
    print '\t\\begin{center}'
    print '\t\t\\includegraphics[width=11cm]{' + fileName + '.png}'
    print '\t\t\\caption{S21 plot of power and phase for antenna P' + str(antNum) + ', 100-1300 MHz.}'
    print '\t\t\\label{'+fileName+'plot}'
    print '\t\\end{center}'
    print '\\end{figure}'
    print ''
    print ''
