# functions related to energy and work
import math

def calcEnergy( rCell, param):
    energy = 0.0

    area0 = param['aCell'] / param['pxReal']**2
    area  = float( len(rCell))
    energy += param['lArea'] * ( area - area0)**2

    return energy


def calcProb( ui, uf, w):
    uDelta = uf - ui
    return math.exp( min( 0.0, w-uDelta ) )
