# Band Generation for RRTMGP

As part of generating the bands used in RRTMGP computations, we need to validate the calculated fluxes and heating rates with those from bands in a reference model (LBLRTM). The scripts in this repository perform this comparison via profile (line) and statistics (scatter) plots, which are generated for each band (including broadband) and profile. In the profile plots, each profile has a separate plotting page with 6 panels: raw downwelling flux, upwelling flux, and heating rates along with their respective test-reference differences. Every band corresponds to a separate PDF file. In the stats plots, the test-reference errors at certain layers in the atmosphere are plotted as a function of the LBLRTM value.

The plotting was always subject to trial and error, so many iterations have been done and conventions have changed, leading to some pretty ad-hoc code. Some cleanup was done with regard to this -- obsolete functions and code snippets were removed, Python 3 was invoked, and naming conventions were changed (caveat: `band_lims_wvn` was changed to `band_wavenumbers`, and the array dimensions were flipped, but this was based on old input netCDFs and we will likely have to revert back the convention when working with newer netCDFs).

Plotting also takes a noticeable amount of time (profiles: ~10s per band). The bottleneck is creating a multi-page PDF. No solution is known currently.

To run out of the box, the following netCDF files are expected to exist:

```
lblrtm-lw-flux-inputs-outputs-garand-all.nc
new-rrtm-lw-flux-inputs-outputs-garand-all.nc
```

Additionally, the default [configuration file](https://docs.python.org/3/library/configparser.html#quick-start) is in version control (`LBLRTM_RRTMGP_config.ini`) and assumes these two netCDFs are in the working directory. It does not need to be provided as an argument, but other configuration files with different names can be specified with the `--config_file` keyword. Then:

```
% ./LBLRTM_RRTMGP_compare.py -p # profile plots
% ./LBLRTM_RRTMGP_compare.py -p -b 2 # plot a single band
% ./LBLRTM_RRTMGP_compare.py -p --broad # plot profiles for broadband only
% ./LBLRTM_RRTMGP_compare.py -s -s1 # without `-s1`, each band is plotted to its own page
```

Other options:

```
usage: LBLRTM_RRTMGP_compare.py [-h] [--config_file CONFIG_FILE] [--plot_profiles] [--plot_stats] [--single_stat]
                                [--tropopause_pressure TROPOPAUSE_PRESSURE] [-csv STATS_CSV] [--band BAND [BAND ...]] [--log_y] [--broad_only]

Generate comparative plots for calculations from two radiative transfer models.

optional arguments:
  -h, --help            show this help message and exit
  --config_file CONFIG_FILE, -i CONFIG_FILE
                        Path to configuration file that contains the values to keyword arguments (so this method can be used instead of providing
                        the keywords; config_file data supercede any keyword argument input).
  --plot_profiles, -p   Generates an N-band page PDF, each with a 3x2 array of plots that compare the upward and downwward fluxes from the
                        reference and test models. Also plots heating rate for each model.
  --plot_stats, -s      Generates a single page with RRTMGP-LBLRTM vs. LBLRTM statistics for TOA upwelling, SFC downwelling, and Max Net flux
                        cases as well as Max Tropospheric and Stratospheric HR.
  --single_stat, -s1    Provide a single PDF for all bands in the stats plots (one band per page).
  --tropopause_pressure TROPOPAUSE_PRESSURE, -t TROPOPAUSE_PRESSURE
                        Pressure threshold that separates troposphere and stratosphere.
  -csv STATS_CSV, --stats_csv STATS_CSV
                        Filename for CSV file that contains band/column stats.
  --band BAND [BAND ...], -b BAND [BAND ...]
                        Number of band to plot. Default (band=None) is all bands.
  --log_y               Generate a semilog-y plot.
  --broad_only          Only generate a broadband plot.
```

# Parallelizing Jobs

To speed up the very time consuming job of plotting the profile errors, a module that uses multithreading has been written. The number of available cores, _N_, that can be used can be deduced with:

```
./pool_LBLRTM_RRTMGP_compare.py --config_file LBLRTM_RRTMGP_config.ini -ac
```

Then, something less than _N_ can be used in farming out all of the jobs:

```
./pool_LBLRTM_RRTMGP_compare.py --config_file LBLRTM_RRTMGP_config.ini -p -c 10
```

**Note**: jobs are only parallelized over bands (including broadband), not profiles.

**Note**: stats plots are not parallelized.

# Forcing

To be continued...
