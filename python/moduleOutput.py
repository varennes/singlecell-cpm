# functions related to energy and work
import math

# output cell 'vector' at time-step 't' to file 'filename'
def outputVector( vector, t, filename):
    if len(vector) != 2:
        print 'ERROR with outpurVector: len != 2'
        return

    s =  ''
    s += '{: >10.5f}'.format( vector[0])
    s += ' '
    s += '{: >10.5f}'.format( vector[1])
    s += ' '
    s += '{: >4d}'.format( t)
    s += '\n'

    with open( filename, 'a') as f:
        f.write(s)


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


# setup the file for the polarization vector output, return filename
def fileSetupPLR():
    fnPLR = 'polar.dat'
    with open( fnPLR, 'w') as f:
        pass
    return fnPLR


# setup the file for the rCell lattice output, return filename
def fileSetupRCELL():
    fnRCELL = 'rcell.dat'
    with open( fnRCELL, 'w') as f:
        pass
    return fnRCELL
