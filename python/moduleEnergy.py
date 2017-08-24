# functions related to energy, work, polarization, and chemical concentration
import math
import numpy as np

# mean chemical concentration, linear profile
def cLinear( x, param):
    return param['c0'] + (param['g'] * float(x) / param['pxReal'])


# calculate probability of adding/removing a lattice site
def calcProb( ui, uf, w):
    uDelta = uf - ui
    return math.exp( min( 0.0, w-uDelta ) )


def calcEnergy( rCell, param):
    energy = 0.0

    area0 = param['aCell'] / param['pxReal']**2
    area  = float( len(rCell))
    energy += param['lArea'] * ( area - area0)**2

    perim0 = param['pCell'] / param['pxReal']
    perimList, perim = calcPerimeter( rCell)
    energy += param['lPerim'] * ( float(perim) - perim0)**2

    return energy


def calcWork( comi, comf, plr):
    work = 0.0
    deltaCOM = [ comf[0] - comi[0], comf[1] - comi[1]]
    norm = sum( x*x for x in deltaCOM)
    if norm > 1e-10:
        deltaCOM = [ x/norm for x in deltaCOM]
        work = deltaCOM[0]*plr[0] + deltaCOM[1]*plr[1]

    return work


def evolvePlr( plr, rCell, com, deltaCOM, param):
    # get cell perimeter
    perimList, perim = calcPerimeter( rCell)
    # calculate local and mean cell chemical concentration signals
    localSignal, meanSignal = calcSignal( rCell, perimList, param)
    # make polarization vector update
    q = [ 0.0, 0.0]
    for lattice in rCell:
        if lattice in perimList:
            i = perimList.index(lattice)
            signal = ( localSignal[i] - meanSignal) / meanSignal
            qtemp  = [ float(i)-j for i,j in zip( lattice, com)]
            norm   = ( qtemp[0]**2 + qtemp[1]**2)**0.5
            q[0]  += signal * qtemp[0] / norm
            q[1]  += signal * qtemp[1] / norm
    q = [ x / float(len(perimList)) for x in q]
    for i in range(2):
        plr[i] = (1.0 - param['rVec'])*plr[i] + param['eVec']*q[i] + deltaCOM[i]
    return plr


def calcSignal( rCell, perimList, param):
    local = []
    mean  = 0.0
    for lattice in rCell:
        c = cLinear( lattice[0], param)
        if c > 0.0:
            local.append( np.random.poisson(c))
        else:
            local.append( 0.0)
    mean = np.mean( local)
    return local, mean


def calcPerimeter( rCell):
    nnMove = [[1,0], [0,1], [-1,0], [0,-1]]
    perim = 0
    perimList = []
    for lattice in rCell:
        perimCheck = 0
        for nn in nnMove:
            latticeNN = []
            latticeNN.append( lattice[0] + nn[0])
            latticeNN.append( lattice[1] + nn[1])
            if latticeNN in rCell:
                pass
            else:
                perimCheck += 1
        if perimCheck > 0:
            perim += perimCheck
            perimList.append( lattice)
    return perimList, perim
