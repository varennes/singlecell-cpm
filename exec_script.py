from tempfile import mkstemp
from shutil import move
from os import remove
import sys
import os

# change a parameter value in file 'mod0param.f90'
# function adapted from https://gist.github.com/kirang89/6478017
def replaceParam( param, newValue):
    paramList = [
    "alpha", "lArea", "lPerim",
    "c0", "g",
    "rVec", "eVec",
    "aCell", "pxReal",
    "lfinish", "runTotal", "tMCmax"
    ]
    if paramList.index(param) > 8:
        sOld = "integer,  parameter :: " + param
        sNew = sOld
        sNew += " =  " + str(int(newValue)) + "\n"
    else:
        sOld = "real(b8), parameter :: " + param
        sNew = sOld
        sNew += " =  " + str(float(newValue)) + "_b8" + "\n"

    source_file_path = "mod0param.f90"
    fh, target_file_path = mkstemp()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
                if line.startswith(sOld):
                    target_file.write(sNew)
                else:
                    target_file.write(line)
    remove(source_file_path)
    move(target_file_path, source_file_path)

# compiles as fortran files using gfortran
def compile():
    os.system("gfortran -c mod0* mod1* mod2* mod3*")
    os.system("rm *.o")
    os.system("gfortran main.f90 mod*")
    os.system("rm *.mod")
    pass

'''
List of parameters that can be changed:

"alpha", "lArea", "lPerim"
"c0", "g"
"rVec", "eVec"
"aCell", "pxReal"
"lfinish", "runTotal", "tMCmax"
'''

if __name__ == '__main__':
    # all possible parameters and their default values
    paramList = [
    "alpha", "lArea", "lPerim",
    "c0", "g",
    "rVec", "eVec",
    "aCell", "pxReal",
    "lfinish", "runTotal", "tMCmax"
    ]
    paramDefault = [
    0.8, 0.5, 0.01,
    00.0, 0.05,
    0.5, 0.01,
    400, 2.0,
    9, 50, 200
    ]
    # iterate over psList of parameters with values contained in pvList
    psList = ["runTotal"]
    pvList = [ [1]]

    # psList = [ "rVec", "eVec"]
    # pvList = [ [0.5, 0.1, 0.05, 0.01], [0.0, 0.01, 0.1]]

    # set all paramters to their default values
    for param in paramList:
        i = paramList.index(param)
        replaceParam( param, paramDefault[i])

    # check if directory for output data exists
    data_dir_location = 'data/test'
    if not os.path.exists(data_dir_location):
        os.makedirs(data_dir_location)

    # compile and run program for all parameters and values
    for ps in psList:

        ps_subdir_location = 'data/test/' + ps
        if not os.path.exists(ps_subdir_location):
            os.makedirs(ps_subdir_location)

        i = psList.index(ps)
        # reset previous parameters to their default value
        for j in range(i):
            param = psList[j]
            jDefault = paramList.index(param)
            replaceParam( param, paramDefault[jDefault])

        # run program over all values of current parameter
        for pv in pvList[i]:
            replaceParam( ps, pv)
            compile()
            print "EXECUTING w/ " + ps + " = " + str(pv) + "\n"
            os.system("./a.out")
            print "\n" + "\n"

            # move output files to desired directory
            data_subdir = ps_subdir_location + '/' + str(pv)
            if not os.path.exists(data_subdir):
                os.makedirs(data_subdir)
            os.system("mv *.dat " + data_subdir)
