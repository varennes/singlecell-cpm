# functions related to energy and work
import math

# output cell 'COM' at time-step 't' to file 'filename'
def outputCOM( com, t, filename):
    strCOM =  ''
    strCOM += '{: >10.5f}'.format( com[0])
    strCOM += ' '
    strCOM += '{: >10.5f}'.format( com[1])
    strCOM += ' '
    strCOM += '{: >4d}'.format( t)
    strCOM += '\n'

    with open( filename, 'a') as f:
        f.write(strCOM)


def outputRCELL( rCell, t, filename):
    for lattice in rCell:
        s  = '{: >6d}'.format( lattice[0])
        s += '{: >6d}'.format( lattice[1])
        s += '{: >5d}'.format( t)
        s += '\n'
        with open( filename, 'a') as f:
            f.write(s)


# setup the file for the COM output, return filename
def fileSetupCOM( iRun):
    if iRun < 10:
        fnCOM = 'com00' + str(iRun) + '.dat'
    elif iRun < 100:
        fnCOM = 'com0'  + str(iRun) + '.dat'
    else:
        fnCOM = 'com'   + str(iRun) + '.dat'
    with open( fnCOM, 'w') as f:
        pass
    return fnCOM


# setup the file for the COM output, return filename
def fileSetupRCELL():
    fnRCELL = 'rcell.dat'
    with open( fnRCELL, 'w') as f:
        pass
    return fnRCELL
