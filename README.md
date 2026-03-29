# KIQR: Knowledge Integration Quantile Regression

**Knowledge Integration Quantile Regression for Ultra-High-Dimensional Data**

[![R](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)

## Overview

This repository provides the R implementation and simulation scripts for the **Knowledge Integration Quantile Regression (KIQR)** method.

KIQR is a penalized quantile regression approach designed for ultra-high-dimensional data. It enables simultaneous variable selection and estimation while flexibly incorporating prior knowledge. The method uses a Huber loss approximation for the quantile check loss and combines it with the local linear approximation (LLA) of the SCAD penalty, achieving computational efficiency and desirable theoretical properties.

KIQR is particularly useful when interest lies in modeling specific conditional quantiles (e.g., upper quantiles) rather than the conditional mean.

## Key Features

- Simultaneous variable selection and parameter estimation
- Flexible integration of prior knowledge (prior set)
- Huber loss approximation for quantile regression
- LLA-SCAD penalty with cyclic coordinate descent algorithm
- Support for modeling at different quantile levels (τ)
- Extensive numerical simulations for performance evaluation

## Repository Structure

- `R/`: Core implementation files (API, loss and penalty functions, solvers, tuning criteria, and utilities)
- `paper_simulations/`: Simulation scripts for the numerical studies
  - `simulation_exp1_n200p1500_norm.R`
  - `simulation_exp1_n200p1500_t3.R`
  - `simulation_exp2_n200p1500_heterX.R`
- `compile_scripts/`: Scripts to run simulations and compile results
  - `compile_res_exp1.R`
  - `generate_table_for_exp2.R`
- `Rdata/`: Intermediate results and simulation outputs (for result replication).

## How to Use

### Install Dependencies

Install the required R packages (for simulations):

```r
install.packages(c("quantreg", "glmnet", "Matrix", "foreach", "doParallel", "rqPen",
                   ,"mvtnorm", "MASS", "ggplot2", "dplyr"))
```

## License

This project is licensed under the MIT License.
