# StructuralTDCox

The repository contains the example simulation code for "Structural learning in Cox model with time-dependent covariates".

File "Derive_SD_TDCox.R" is the R code for deriving the selection dictionary in the simulation study.

File "500_ABAB.R" is an example code to run Senario 2 when the sample size is 500.

The above two file calls "Function.R", which contains necessary functions to generate covariates, and to perform proximal gradient descent, and cross-validation.

The corresponding .sh file is the files for submitting jobs to high performance computing.
