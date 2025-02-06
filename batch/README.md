# Toy Retrieval

Simple Jupyter Notebook toy retrieval where we generate synthetic radiances and then try to retrieve them.

## Dependencies

* build_sensitivity_matrices_032024.R: Code to generate Jacobians from offline sensitivities (from GEOS-Chem pulses)
* collect_HEMCO_output_write_netcdf_prior_032024.R: Code to collect the GeosCHEM HEMCO flux output to create a "prior" which can be convoluted across inversion solution to create posterior flux est
* generate_transcom_flux_ensemble_from_inversion.R: code to pull MC sample from posterior flux distribution (could be deprecated since write_inversion*R does same thing)
* inversion_032024.R:  core inversion code
* plot_concentrations.R: code to plot concentration data
* util_code_032024.R: a variety of utility functions
* write_inversion_2_netcdf_032024.R:  code to write inversion output to a netcdf file

## Setup

* You'll need to make sure that `input_folder` in `site_settings.ini` points to the top level of the data folder acquired from Zenodo, wherever that resides on your file system. For example, if one unzips ssim-ghg-data.tgz to ssim-ghg-data folder than input should point to that folder


## Executing program

There are currently 4 different Jupyter notebooks to play with: batch_demo.ipynb, chi_square.ipynb, impose_bias_case.ipynb, and kalman_gain.ipynb

* To launch, in a shell:

```
jupyter notebook batch_demo.ipynb
```

Or just,


```
jupyter notebook
```
And select notebook to open after the fact.

* You can also run the Jupyter Notebook in a GUI such as Visual Studio

