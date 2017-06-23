# singlecell-cpm
Cellular Potts Model (CPM) Simulations of single cell chemotaxis.

<img src="data/cpm_mov1.gif" style="width: 300px;"/>



## Running singlecell-cpm

The python script `exec_script.py` can compile and run the program for you. You can also compile and run the program manually by compiling module files and main program file in the following manner:
```
gfortran -c mod0* mod1* mod2* mod3*
rm *.o
gfortran main.f90 mod*
rm *.mod
./a.out
```
In order to use `exec_script.py` the [GNU Fortran compiler](https://gcc.gnu.org/fortran/) must be installed on your machine. A different compiler can be used if the `compile()` function in `exec_script.py` is appropriately modified.

`exec_script.py` allows you to run the program for various changing parameters. For each respective parameter and value the output from the program is placed in its designated directory.

## Output Files

The module file `mod3wrt.f90` contains all functions and subroutines which create the output files listed below.

### `ci.dat`

Outputs the Chemotactic Index (CI) of the cell.

- First column is the CI value.
- Second column is the run of corresponding to that row's CI value.

### `cr.dat`

Outputs the Chemotactic Ratio (CR) of the cell.

- First column is the CR value.
- Second column is the run of corresponding to that row's CR value.

### `v_mean.dat`

Outputs mean instantaneous speed over each run.

- First column is the mean speed in units of pixels per time-steps.
- Third column is the run of the recorded mean speed.

### `v_inst.dat`

Outputs instantaneous speed of cell as a function of time for each run.

- First column is the speed in units of pixels per time-steps.
- Second column is the time-steps of the recorded speed.
- Third column is the run of the recorded speed.
