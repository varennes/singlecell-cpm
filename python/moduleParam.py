# contains all parameters

def getParameters():
    param = dict()
    # energy related parameters
    param['alpha']  = 0.8
    param['lArea']  = 0.5
    param['lPerim'] = 0.01
    param['w0']     = 1.0
    # cell related parameters
    param['aCell'] = 400.0
    param['rVec']  = 0.5
    param['eVec']  = 0.1
    # concentration profile parameters
    param['c0'] = 0.0
    param['g']  = 0.5
    # simulation parameters
    param['pxReal']   = 5.0
    param['lFinish']  = 9
    param['runTotal'] = 1
    param['tMCmax']   = 2

    return param
