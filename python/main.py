import random as rd
import moduleParam as m0
import moduleInit  as m1

if __name__ == '__main__':
    # print rd.random()

    param = m0.getParameters()

    for iRun in xrange( param['runTotal']):
        comCell = []
        rCell = m1.initSystem( param)
        comCell.append( m1.calcCOM( rCell))
        print len(rCell), comCell

        # relax cell shape and its polarization vector before starting movement tracking
        m1.relaxCell( rCell, param)
