import random as rd
import moduleParam  as m0
import moduleInit   as m1
import moduleOutput as m2
import moduleEnergy as m3

if __name__ == '__main__':
    # get parameters
    param = m0.getParameters()

    for iRun in xrange( param['runTotal']):
        # setup output files
        fnCOM = m2.fileSetupCOM( iRun)
        fnRCL = m2.fileSetupRCELL()

        comCell = []
        rCell = m1.initSystem( param)
        lCell = len(rCell)**0.5
        comCell.append( m1.calcCOM( rCell))

        # # perimeter test
        # perimList, perim = m3.calcPerimeter( rCell)
        # print len(perimList), perim

        # # test output
        # m2.outputCOM( comCell[-1], 0, fnCOM)

        # relax cell shape and its polarization vector before starting movement tracking
        rCell = m1.relaxCell( rCell, param)
        comCell.append( m1.calcCOM( rCell))
        m2.outputCOM( comCell[-1], 0, fnCOM)
        m2.outputRCELL( rCell, 0, fnRCL)

        # time evoluation of cell
        plr = [ 0.0, 0.0]
        for tMC in xrange( param['tMCmax']):
            # update cell lattices
            rCell = m1.evolveCell( rCell, comCell[-1], lCell, param)
            comCell.append( m1.calcCOM( rCell))
            # update polarization vector
            deltaCOM = [ i-j for i,j in zip( comCell[-1], comCell[-2])]
            plr = m3.evolvePlr( plr, rCell, com, deltaCOM, perimList, param)
            # output data
            m2.outputCOM( comCell[-1], tMC+1, fnCOM)
            m2.outputRCELL( rCell, tMC+1, fnRCL)
