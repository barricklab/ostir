# OSTIR Calibration

This directory contains the input file and R script used to calibrate the energy model parameters used by OSTIR with ViennaRNA 2.4.17 and the default Turner2004 energy parameters.

The data is from: Salis HM, Mirsky EA, Voigt CA. 2009. Automated design of synthetic ribosome binding sites to control protein expression. _Nat. Biotechnol._ **27**:946–950.

The fitting procedure used by OSTIR also follows the procedure from Salis et al. 2009, except with some modifications described in the R script.

To run the calibration:
```bash
Rscript ostir_calibration.R
```

## Fitting Δ_G_(spacing)

![Spacing Model Plot](output.spacing_dG.png)

## Total Δ_G_ versus measured translation initiation rate
![Spacing Model Plot](output.total_dG_versus_measured_rate.png)

## Predicted versus measured translation initiation rates
![Predicted Versus Measured Rate](output.predicted_rate_versus_measured_rate.png)

## Fold error predicted versus measured
![Log2 Fold Error](output.log2_fold_error.png)
