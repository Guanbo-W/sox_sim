# sox_sim

The repository contains the simulation code for "Structural learning in Cox model with time-dependent covariates".

File "Derive_SD_TDCox.R" is the R code for deriving the selection dictionary in the simulation study (low-dim case).

File "Low_Dim.R" is an example code to run the low-dim case-Senario 2 when the sample size is 500, which calls "Function.R" containing necessary functions to generate covariates, and to perform proximal gradient descent, and cross-validation. The code does not call the package sox, for the sake of showing how to perform selecting time-dependent variables step by step. Of course, the package soc can do exactly the same job.

The corresponding .sh file is the files for submitting jobs to high performance computing.
