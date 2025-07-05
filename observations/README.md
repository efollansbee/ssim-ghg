### Observation Analysis

This notebook demonstrates some analyses of different observation datasets at the LEF tower in Park Falls, Wisconsin. There are three in situ analyzer altitudes as well as aircraft data and TCCON total column CO2. The basic concepts covered are:
1. Reading in netCDF files with xarray
2. Creating Pandas DataFrames and using them for time series analysis
3. Seasonal cycles and trends in atmospheric data
4. Diurnal cycles in surface in situ data
5. Eddy covariance CO2 flux measurements

FFrom the command line, run `aws s3 sync s3://ghg-ssim/ssim-ghg-data /tmp/ssim-ghg-data` to download the needed data to your local temp directory. Navigate to `~/ssim-ghg/` and copy `site_settings.yml.tmp` to `site_settings.yml`. Open `site_settings.yml` in an editor like `vi` or Jupyter itself and add `/tmp/ssim-ghg-data/` after the colon on the `input_folder` line under `global_paths`. 

 *If you are doing this work outside of the US GHGHub, you will need to place the data in a different place and modify the `site_settings.yml` file accordingly.*

To run the notebook, open the `lef_observations.ipynb` notebook within Jupyter. You can select "Run all cells" under the Run menu or step through each cell one by one. 

Libraries needed (should all be present in the GHGHub SSIM image):
- pandas
- numpy
- xarray
- pytz
- scipy
- pyyaml
- matplotlib
- seaborn
