# Variation inversion

There are different Jupyter notebooks here for exploring different aspects of the variational inversion technique.

* `bits_and_pieces.ipynb` for understanding preconditioning
* `explore_biases.ipynb` for understanding what happens to estimated fluxes if there are persistent biases in a subset of obs
* `explore_errors.ipynb` for understanding posterior covariances estimated by running an ensemble of variational inversions
* `var4d_demo.ipynb` is a sampling of the variational technique with different input datasets

You probably want to start with `var4d_demo.ipynb` if this is your first rodeo.

## Setup

* You'll need to make sure that `input_folder` in `site_settings.ini` points to the top level of the data folder acquired from Zenodo, wherever that resides on your file system. For example, if one unzips ssim-ghg-data.tgz to ssim-ghg-data folder than input should point to that folder

## Executing example program

In whatever system you have to execute Jupyter notebooks, load `var4d_demo.ipynb`. The first couple of commands in the notebook,

```python
from var4d_components import Var4D_Components
from visualize_results import Visualize_Obs, Visualize_Fluxes, Diagnostic_Plots
```

is for loading the necessary python classes, which are in turn coded in `var4d_components.py` and `visualize_results.py`. After that, the world is your oyster. The first example code block,

```python
var4d = Var4D_Components('only_noaa_observatories', verbose=True, store_intermediate=True)
flux_corr_structure = {'temp_corr': 2.0}
obs_assim_dict = {'sites': ['mlo', 'spo', 'brw', 'smo']}
prior_flux_unc_dict = {'prior_unc_source': 'nee', 'prior_unc_scale': {'land': 0.25, 'ocean': 0.5}}
var4d.var4d_setup(obs_to_assim=obs_assim_dict, corr_structure=flux_corr_structure, **prior_flux_unc_dict)
```

sets up a basic flux inversion with the following specifications:

* It is given a name `only_noaa_observatories` that helps you remember what that inversion does. The code doesn't care. You could name it `mindolluin` or `employee_benefits` as long as you can remember what it is. The inversion output is placed in an output folder with this name.
* `verbose=True` prints the cost function and gradient for every single iteration, while `verbose=False` prints a set of progress bars that tracks the number of function evaluations, gradient evaluations, etc.
* `store_intermediate=True` stores all the intermediate states as the optimized iterates. Setting it to `False` will only store the initial (prior) and final (posterior) emissions.
* `obs_assim_dict` is a dictionary that specifies which observations are assimilated. There are too many options to list here, check `setup_obs_errors` of the `Vard4D_Components` class in `var4d_components.py`. The specification above assimilates only *flask* samples from the four insitu sites specified by the three-letter codes. There are options to assimilate OCO2 soundings, TCCON data, some combination of data, etc.
* The error specification in `prior_flux_unc_dict` says that the prior flux error is 0.25 x abs(NEE) over land and 0.5 x abs(flux) over ocean. The land flux specification can be changed to key off GPP or respiration, see code to see how.
* `flux_corr_structure` basically says that the prior flux error is correlated spatially with an e-folding time of 2 months. There is no spatial correlation imposed because we're solving for TRANSCOM regions.
* The call to `var4d_setup` sets up the prior and obs covariances, and creates the pseudo-obs to assimilate by running transport with the "true" fluxes.

The next call should be to `var4d.var4d_chain(gradnorm=1.0E-5)` which says "iterate until the L2 norm of the gradient goes down by $10^5$". This will print a bunch of diagnostic messages until the optimization converges. The following lines of code are keyed off the name of the inversion you provided, `only_noaa_observatories`, and plots the convergence trajectory, fits to obs, fluxes, etc.

```python
diag = Diagnostic_Plots('only_noaa_observatories')
diag.plot_convergence()
po1 = Visualize_Obs('only_noaa_observatories')
po1.plot_site(['spo','mlo','smo', 'brw'])
po1.plot_site(['kum','amt','lef','crz'])
vf1 = Visualize_Fluxes('only_noaa_observatories')
vf1.plot_region(['North American Boreal', 'North American Temperate', 'South American Tropical', 'South American Temperate'])
```