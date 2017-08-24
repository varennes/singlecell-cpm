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
        fnPLR = m2.fileSetupPLR()
        fnRCL = m2.fileSetupRCELL()

        comCell = []
        rCell = m1.initSystem( param)
        lCell = len(rCell)**0.5
        comCell.append( m1.calcCOM( rCell))

        # relax cell shape
        rCell = m1.relaxCell( rCell, param)
        comCell.append( m1.calcCOM( rCell))
        # initialize polarization vector
        plr = [ 0.0, 0.0]
        for i in xrange( int( 5.0 / param['rVec'])):
            plr = m3.evolvePlr( plr, rCell, comCell[-1], [ 0.0, 0.0], param)
            # print plr

        # output initial cell state
        m2.outputRCELL( rCell, 0, fnRCL)
        m2.outputVector( comCell[-1], 0, fnCOM)
        m2.outputVector( plr, 0, fnPLR)

        # time evoluation of cell
        for tMC in xrange( param['tMCmax']):
            # update cell lattices
            rCell = m1.evolveCell( rCell, comCell[-1], plr, lCell, param)
            comCell.append( m1.calcCOM( rCell))
            # update polarization vector
            deltaCOM = [ i-j for i,j in zip( comCell[-1], comCell[-2])]
            plr = m3.evolvePlr( plr, rCell, comCell[-1], deltaCOM, param)
            # output data
            m2.outputRCELL( rCell, tMC+1, fnRCL)
            m2.outputVector( comCell[-1], tMC+1, fnCOM)
            m2.outputVector( plr, tMC+1, fnPLR)
