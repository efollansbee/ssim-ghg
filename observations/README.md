### Observation Analysis

This notebook demonstrates some analyses of different observation datasets at the LEF tower in Park Falls, Wisconsin. There are three in situ analyzer altitudes as well as aircraft data and TCCON total column CO2. The basic concepts covered are:
1. Reading in netCDF files with xarray
2. Creating Pandas DataFrames and using them for time series analysis
3. Seasonal cycles and trends in atmospheric data
4. Diurnal cycles in surface in situ data
5. Eddy covariance CO2 flux measurements

To run the notebook, first configure the variables in `site_settings.yml` to point at the folders and files you want to read and analyze. 

From the command line, run `aws s3 sync s3://ghg-ssim/ssim-ghg-data /tmp/ssim-ghg-data` to download the needed data to your local temp directory. If you are doing this work outside of the US GHGHub, you will need to place the data in a different place and modify the `site_settings.yml` file accordingly.

Then run Jupyter-lab and open the `lef_observations.ipynb` notebook. You can select "Run all cells" under the Run menu or step through each cell one by one. 

Libraries needed:
- pandas
- numpy
- xarray
- pytz
- scipy
- pyyaml
- matplotlib