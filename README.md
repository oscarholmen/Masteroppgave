This repository contains the R scripts used in the development, testing, and evaluation of Bayesian inference methods for a superposed renewal process. The code supports the Gibbs sampling procedure, dataset generation, and an experimental Kalman filter–based approach.

## **Files**

- **Gibbs Sampler.R**  
  Core implementation of the Gibbs sampler used for posterior inference on the two renewal processes.

- **KalmanFilterMethod.R**  
  Implementation of the filter Method, implemented with a MCMCmetrop1r.

- **Testing Blocks.R**  
  A bad name, but this contains all the code for testing the Gibbs Sampling methods discussed. 

- **TestingFilterMet.R**  
  Supplementary testing code for evaluating the Kalman filter method.

- **UpdatedGeneratingDataset.R**  
  Script for generating synthetic datasets consistent with the model structure used in the thesis.

## **How to run**

All scripts can be sourced directly in R. Dependencies are limited to standard statistical and plotting packages (e.g., `ggplot2`, `dplyr`, `coda`, etc.).
