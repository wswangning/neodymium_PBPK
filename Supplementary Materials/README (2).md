# Neodymium PBPK Model - Supplementary Materials

This repository contains the supplementary materials, model code, and data for the manuscript:

**"Development and Validation of a Physiologically Based Pharmacokinetic Model for Neodymium in Mice"**

## Contents
- **Model code**: Berkeley Madonna and R (mrgsolve) implementations.
- **Model parameters**: physiological, partition coefficients, and pharmacokinetic parameters.
- **Experimental data**: time-course concentrations (processed).
- **Sensitivity analysis**: global Sobol indices.
- **Supplementary tables**: experimental design, NCA parameters, half-life comparison, model parameters, and sensitivity results.

## How to Use
- The `Code/` directory contains the PBPK model source files.  
  - For Berkeley Madonna: open `Nd_PBPK_Model_MAIN.mmd`.  
  - For R: source `Nd_PBPK_mrgsolve.R` and run the analysis script.

- The `Data/` directory provides the individual animal data used for model calibration.

- Supplementary tables are provided as CSV files and referenced in the manuscript.

## Citation
If you use this model, please cite the original manuscript (under review) and this repository.

## License
MIT License. See LICENSE file.