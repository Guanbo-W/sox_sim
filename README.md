# sox_sim

## Simulation code for "Structured learning in Cox model with time-dependent covariates".
The paper shows how to structurally select variables in time-dependent Cox models.

- File "Derive_SD_TDCox.R" is the R code for deriving the selection dictionary in the simulation study (low-dim case).

* File "Low_Dim.R" is example code for the low-dim case-Senario 2 when the sample size is 500, which calls "Low_Dim_Function.R" containing necessary functions. The code does not call the package ```sox```, for the sake of showing how to perform selecting time-dependent variables step by step. Of course, the package ```sox``` can do exactly the same job.

+ File "High_Dim.R" is the imulation code for the high-dimensional setting, calls the package sox.

- The corresponding .sh file is the files for submitting jobs to high performance computing.

Please see [here](https://cran.r-project.org/web/packages/sox/index.html) for the details about the ```R``` package ```sox```.
