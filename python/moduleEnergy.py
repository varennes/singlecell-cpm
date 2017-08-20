# functions related to energy and work
import math

def calcEnergy( rCell, param):
    energy = 0.0

    area0 = param['aCell'] / param['pxReal']**2
    area  = float( len(rCell))
    energy += param['lArea'] * ( area - area0)**2

    perim0 = param['pCell'] / param['pxReal']
    perimList, perim = calcPerimeter( rCell)
    energy += param['lPerim'] * ( float(perim) - perim0)**2

    return energy


def calcProb( ui, uf, w):
    uDelta = uf - ui
    return math.exp( min( 0.0, w-uDelta ) )


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
