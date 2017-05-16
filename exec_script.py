from tempfile import mkstemp
from shutil import move
from os import remove
import sys
import os

# change a parameter value in file 'mod0param.f90'
# function adapted from https://gist.github.com/kirang89/6478017
def replaceParam( param, newValue):
    paramList = [
    "alpha", "lArea", "lPerim", "w0", "aCell", "c0", "g",
    "pxReal", "lfinish", "runTotal", "tMCmax"
    ]
    if paramList.index(param) > 7:
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

### List of parameters that can be changed ###
#
#  "alpha", "lArea", "lPerim", "w0", "aCell"
#  "c0", "g"
#  "pxReal", "lfinish", "runTotal", "tMCmax"
#
###

if __name__ == '__main__':
    # list of paramater names and values to iterate over
    psList = ["tMCmax", "alpha"]
    pvList = [ [10, 20, 30], [0.1, 0.5]]

    # check if directory for output data exists
    data_dir_location = 'data/test'
    if not os.path.exists(data_dir_location):
        os.makedirs(data_dir_location)

    # compile and run program for all parameters and values
    for ps in psList:

        data_dir_location = 'data/test/' + ps
        if not os.path.exists(data_dir_location):
            os.makedirs(data_dir_location)

        i = psList.index(ps)
        for pv in pvList[i]:
            replaceParam( ps, pv)
            compile()
            print "EXECUTING w/ " + ps + " = " + str(pv) + "\n"
            os.system("./a.out")
            print "\n" + "\n"

            # move output files to desired directory
            data_subdir = data_dir_location + '/' + str(pv)
            if not os.path.exists(data_subdir):
                os.makedirs(data_subdir)
            os.system("mv *.dat " + data_subdir)
