# functions for simulation initialization
import random as rd
import moduleEnergy as mE


def evolveCell( rCell, com, plr, lCell, param):
    elemTime = int(9*lCell**2)

    lmin = int(lCell*1.5)
    lmax = 3*lCell - lmin - 1
    xlist = [ int(com[0]) - lmin, int(com[0]) + lmax]
    ylist = [ int(com[1]) - lmin, int(com[1]) + lmax]

    for iElem in xrange( elemTime):
        a, b = pickLatticePair( xlist, ylist)
        a.append( (a in rCell) )
        b.append( (b in rCell) )
        # print 'after in check', len(rCell)
        if a[2] != b[2]:
            if a[2]:
                # cell is adding a lattice site
                rTemp = rCell[:]
                rTemp.append( b[:2])
                # rTemp is simply connected
                connectedCheck = True
            else:
                # cell is removing a lattice site
                rTemp = rCell[:]
                rTemp.remove( b[:2])
                # check if rTemp is simply connected
                floodList = []
                floodList = floodFill( rCell[0], floodList, rCell)
                if len(floodList) == len(rCell):
                    connectedCheck = True
                else:
                    connectedCheck = False

            if connectedCheck:
                # calculate change in energy
                ui = mE.calcEnergy( rCell, param)
                uf = mE.calcEnergy( rTemp, param)
                # calculate work
                comTemp = calcCOM( rTemp)
                w  = mE.calcWork( com, comTemp, plr)
                # get probability of accepting change
                prob = mE.calcProb( ui, uf, w)
                # print 'prob = %s, ui = %s, uf = %s' %(prob,ui,uf)
                if ( rd.random() < prob):
                    # print 'lattice change accepted'
                    rCell = rTemp

    return rCell


def relaxCell( rCell, param):
    lCell = len(rCell)**0.5
    elemTime = int(9*lCell**2)
    print 'Number of initialization time steps: %i' % elemTime

    xlist = [-lCell, 2*lCell-1]
    ylist = [-lCell, 2*lCell-1]

    for iRelax in xrange( 1*elemTime):
        a, b = pickLatticePair( xlist, ylist)
        a.append( (a in rCell) )
        b.append( (b in rCell) )
        # print 'after in check', len(rCell)
        if a[2] != b[2]:
            if a[2]:
                # cell is adding a lattice site
                rTemp = rCell[:]
                rTemp.append( b[:2])
                # rTemp is simply connected
                connectedCheck = True
            else:
                # cell is removing a lattice site
                rTemp = rCell[:]
                rTemp.remove( b[:2])
                # check if rTemp is simply connected
                floodList = []
                floodList = floodFill( rCell[0], floodList, rCell)
                if len(floodList) == len(rCell):
                    connectedCheck = True
                else:
                    connectedCheck = False

            if connectedCheck:
                # calculate change in energy
                ui = mE.calcEnergy( rCell, param)
                uf = mE.calcEnergy( rTemp, param)
                # get probability of accepting change
                prob = mE.calcProb( ui, uf, 0.0)
                # print 'prob = %s, ui = %s, uf = %s' %(prob,ui,uf)
                if ( rd.random() < prob):
                    # print 'lattice change accepted'
                    rCell = rTemp

    return rCell


# initialize list of cell lattice site locations
def initSystem( param):
    # calculate cell length in terms of simulation lattices
    lCell = int( (param['aCell'] / param['pxReal']**2)**0.5 )
    if lCell == 0:
        lCell = 1

    # set initial cell location
    rCell = []
    for ix in xrange(lCell):
        for iy in xrange(lCell):
            rCell.append( [ix, iy])

    return rCell


# calculate center of mass (COM) of cell
def calcCOM( rCell):
    com = [ 0.0, 0.0]
    for lattice in rCell:
        com[0] += float(lattice[0])
        com[1] += float(lattice[1])
    com = [ x/float(len(rCell)) for x in com]

    return com


def floodFill( lattice, floodList, rCell):
    # exit if lattice is already in floodList
    if lattice in floodList:
        return floodList
    # exit if lattice is not in rCell
    if lattice in rCell:
        pass
    else:
        return floodList
    # exit if floodList has more lattices than rCell
    if len(floodList) > len(rCell):
        return floodList

    floodList.append( lattice)
    nnMove = [[1,0], [0,1], [-1,0], [0,-1]]
    for nn in nnMove:
        newLattice = [ i+j for i,j in zip(lattice, nn)]

    return floodList


def pickLatticePair( x, y):
    a = [0,0]
    a[0] = rd.randint( x[0], x[1])
    a[1] = rd.randint( y[0], y[1])

    nnMove = [[1,0], [0,1], [-1,0], [0,-1]]
    inn = rd.randint(0,3)
    b = [ i+j for i,j in zip(a,nnMove[inn])]

    return a, b
