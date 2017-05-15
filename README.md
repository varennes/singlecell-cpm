# singlecell-cpm
CPM Simulations of Single Cell Chemotaxis

Compilation:
`gfortran -c mod0* mod1* mod2* mod3*; rm *.o; gfortran main.f90 mod*; rm *.mod`

## Output Files

### `v_inst.dat`

Outputs instantaneous speed of cell as a function of time for each run.

- First column is the speed in units of pixels per time-steps.
- Second column is the time-steps of the recorded speed.
- Third column is the run of the recorded speed.

### `ci.dat`

Outputs the Chemotactic Index (CI) of the cell.

- First column is the CI value.
- Second column is the run of corresponding to that row's CI value.

### `cr.dat`

Outputs the Chemotactic Ratio (CR) of the cell.

- First column is the CR value.
- Second column is the run of corresponding to that row's CR value.
