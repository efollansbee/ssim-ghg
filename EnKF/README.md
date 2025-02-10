# Kalman filtering and ensemble methods

Notebooks demonstrating various methods for analytical Kalman
filtering, demonstrating the balance between agreeing with
observational constraints and priors, and solving the toy problem
using three different schemes.

## Available notebooks

* EnKF/welch/welch.ipynb - analytical Kalman filtering example
* EnKF/L-curve/lcurve.ipynb - graphical demonstration of finding optimal balance between agreeing with observational constraints and priors
* EnKF/base/common_problem.ipynb - solving toy problem (no time propagation)
* EnKF/time-varying/persist.ipynb - solving toy problem using persistence time propagation
* EnKF/time-varying/modpersist.ipynb - solving toy problem using modified persistence time propagation

## Dependencies

### Support codes located in EnKF/tools:

* enkf.r
* find.indir.r
* inversion_032024_nobgd.R
* load.ncdf4.r
* normality.test.r
* progress.bar.r
* time.r

### Required third party R packages:

* svd
* plotrix
* ncdf4
* stats
* e1071
* yaml
* txtplot
* MASS
* gplots
* repr
* Matrix
* EnvStats

## Setup

* You'll need to make sure that `input_folder` and `output_folder` are correctly set in `site_settings.yml`.

## Executing notebooks

* In a shell:

```
jupyter lab L-curve/lcurve.ipynb
```

* It is good practice to restart the kernel before executing new
  notebooks. Persistent data can cause unexpected results.
