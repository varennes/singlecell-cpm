# functions related to energy and work
import math

# output cell 'COM' at time-step 't' to file 'filename'
def outputCOM( com, t, filename):
    strCOM  = str(com[0]) +  ' ' + str(com[0]) + ' '
    strCOM += str(t) + '\n'
    with open( filename, 'a') as f:
        f.write(strCOM)


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
